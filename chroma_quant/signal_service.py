"""可编辑色谱信号服务（基线扣除 + 区间面积）。

职责边界：本服务只认识采样点与基线锚点，**不接受任何组分方法**；
组分定量是 :class:`~chroma_quant.repository.QuantRepository` 的职责。
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from fractions import Fraction
from typing import Sequence

from . import limits
from .errors import DomainError, LimitError, NotFoundError, ValidationError
from .quant import SegmentProduct, build_product
from .rational import frac_str, stable_id, to_fraction
from .segments import PLF, Point, make_plf
from .store import connect, get_frac, init_db


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class SignalService:
    def __init__(self, db_path: str):
        self.conn: sqlite3.Connection = connect(db_path)
        init_db(self.conn)

    def close(self) -> None:
        self.conn.close()

    # ------------------------------------------------------------------
    # 导入
    # ------------------------------------------------------------------
    def create_run(self, run_id: str,
                   samples: Sequence[tuple], note: str | None = None) -> None:
        """导入一条信号：t 严格递增的有理 (t, y) 序列，线性插值。"""
        if self._run_exists(run_id):
            raise ValidationError(f"run_id 已存在: {run_id}")
        pts = [Point(to_fraction(t, "t"), to_fraction(y, "y"))
               for t, y in samples]
        if len(pts) < 2:
            raise ValidationError("每条信号至少需要两个采样点")
        if len(pts) > limits.MAX_SAMPLES_PER_RUN:
            raise LimitError(
                f"每条信号最多 {limits.MAX_SAMPLES_PER_RUN} 个采样点，"
                f"收到 {len(pts)}")
        for p, q in zip(pts, pts[1:]):
            if q.t <= p.t:
                raise ValidationError("采样时间必须严格递增")
        self.conn.execute(
            "INSERT INTO runs(run_id, created_at, note) VALUES (?,?,?)",
            (run_id, _now(), note))
        self.conn.executemany(
            "INSERT INTO samples(run_id, idx, t_num, t_den, y_num, y_den)"
            " VALUES (?,?,?,?,?,?)",
            [(run_id, i, p.t.numerator, p.t.denominator,
              p.y.numerator, p.y.denominator) for i, p in enumerate(pts)])
        self.conn.commit()

    def _run_exists(self, run_id: str) -> bool:
        row = self.conn.execute(
            "SELECT 1 FROM runs WHERE run_id=?", (run_id,)).fetchone()
        return row is not None

    def get_signal(self, run_id: str) -> PLF:
        rows = self.conn.execute(
            "SELECT * FROM samples WHERE run_id=? ORDER BY idx",
            (run_id,)).fetchall()
        if not rows:
            raise NotFoundError(f"未知 run_id: {run_id}")
        return make_plf(
            (get_frac(r, "t_num", "t_den"),
             get_frac(r, "y_num", "y_den")) for r in rows)

    # ------------------------------------------------------------------
    # 基线：至少两个互异锚点、锚点位于采样范围内；更正只追加新版本
    # ------------------------------------------------------------------
    def add_baseline(self, run_id: str, anchors: Sequence[tuple]) -> int:
        signal = self.get_signal(run_id)  # 同时校验 run 存在
        pts = [Point(to_fraction(t, "t"), to_fraction(y, "y"))
               for t, y in anchors]
        self._validate_anchors(pts, signal)
        ver_row = self.conn.execute(
            "SELECT COALESCE(MAX(version), 0) AS v FROM baseline_versions"
            " WHERE run_id=?", (run_id,)).fetchone()
        version = int(ver_row["v"]) + 1
        self.conn.execute(
            "INSERT INTO baseline_versions(run_id, version, created_at,"
            " supersedes) VALUES (?,?,?,?)",
            (run_id, version, _now(), version - 1 if version > 1 else None))
        self.conn.executemany(
            "INSERT INTO baseline_anchors(run_id, version, idx,"
            " t_num, t_den, y_num, y_den) VALUES (?,?,?,?,?,?,?)",
            [(run_id, version, i, p.t.numerator, p.t.denominator,
              p.y.numerator, p.y.denominator) for i, p in enumerate(pts)])
        self.conn.commit()
        return version

    def correct_baseline(self, run_id: str,
                         anchors: Sequence[tuple]) -> int:
        """更正一个基线锚点（或整组重设）：返回新版本号，旧版保留。"""
        return self.add_baseline(run_id, anchors)

    @staticmethod
    def _validate_anchors(pts: list[Point], signal: PLF) -> None:
        if len(pts) < 2:
            raise ValidationError("基线至少需要两个锚点")
        if len(pts) > limits.MAX_ANCHORS_PER_BASELINE:
            raise LimitError(
                f"基线锚点最多 {limits.MAX_ANCHORS_PER_BASELINE} 个，"
                f"收到 {len(pts)}")
        for p, q in zip(pts, pts[1:]):
            if q.t <= p.t:
                raise ValidationError("基线锚点时间必须互异且递增")
        lo, hi = signal.domain
        for p in pts:
            if p.t < lo or p.t > hi:
                raise DomainError(
                    f"基线锚点 t={p.t} 超出采样范围 [{lo}, {hi}]，"
                    "锚点必须位于采样范围内")

    def list_baseline_versions(self, run_id: str) -> list[int]:
        rows = self.conn.execute(
            "SELECT version FROM baseline_versions WHERE run_id=?"
            " ORDER BY version", (run_id,)).fetchall()
        return [int(r["version"]) for r in rows]

    def get_baseline(self, run_id: str, version: int | None = None) -> PLF:
        version = self._resolve_version(run_id, version)
        rows = self.conn.execute(
            "SELECT * FROM baseline_anchors WHERE run_id=? AND version=?"
            " ORDER BY idx", (run_id, version)).fetchall()
        if not rows:
            raise NotFoundError(
                f"基线版本不存在: run={run_id} version={version}")
        return make_plf(
            (get_frac(r, "t_num", "t_den"),
             get_frac(r, "y_num", "y_den")) for r in rows)

    def _resolve_version(self, run_id: str,
                         version: int | None) -> int:
        versions = self.list_baseline_versions(run_id)
        if not versions:
            raise NotFoundError(f"run {run_id} 尚无基线")
        if version is None:
            return versions[-1]
        if version not in versions:
            raise NotFoundError(f"基线版本 {version} 不存在")
        return version

    def coverage(self, run_id: str, version: int | None = None
                 ) -> tuple[Fraction, Fraction]:
        """锚点覆盖的闭范围即有效域；域外拒绝而非外推。"""
        return self.get_baseline(run_id, version).domain

    # ------------------------------------------------------------------
    # 独立区间查询：原信号 / 基线 / 净响应 / 有效覆盖
    # ------------------------------------------------------------------
    def query(self, run_id: str, a, b, version: int | None = None) -> dict:
        a, b = to_fraction(a, "a"), to_fraction(b, "b")
        if a >= b:
            raise ValidationError("查询区间要求 a < b")
        signal = self.get_signal(run_id)
        baseline = self.get_baseline(run_id, version)
        dom_lo, dom_hi = baseline.domain
        if a < dom_lo or b > dom_hi:
            raise DomainError(
                f"查询区间 [{a}, {b}] 超出有效域 [{dom_lo}, {dom_hi}]，"
                "拒绝外推")
        from .segments import net_response  # 避免循环导入
        net = net_response(signal, baseline)
        return {
            "run_id": run_id,
            "version": self._resolve_version(run_id, version),
            "interval": (a, b),
            "coverage": (dom_lo, dom_hi),
            "raw_signal": signal.clip(a, b),
            "baseline": baseline.clip(a, b),
            "net_response": net.clip(a, b),
        }

    def interval_area(self, run_id: str, a, b,
                      version: int | None = None) -> Fraction:
        """信号服务独立交付：任意区间的净响应面积（精确有理）。"""
        q = self.query(run_id, a, b, version)
        return q["net_response"].integral()

    # ------------------------------------------------------------------
    # 产品发布与复读
    # ------------------------------------------------------------------
    def publish_product(self, run_id: str,
                        version: int | None = None) -> SegmentProduct:
        signal = self.get_signal(run_id)
        version = self._resolve_version(run_id, version)
        baseline = self.get_baseline(run_id, version)
        product_id = stable_id(
            "product", run_id, version,
            [(frac_str(k.t), frac_str(k.y)) for k in signal.knots],
            [(frac_str(k.t), frac_str(k.y)) for k in baseline.knots],
        )
        existing = self.conn.execute(
            "SELECT product_id FROM products WHERE product_id=?",
            (product_id,)).fetchone()
        if existing:
            return self.get_product(product_id)
        product = build_product(product_id, run_id, version, signal, baseline)
        self._persist_product(product)
        return product

    def _persist_product(self, product: SegmentProduct) -> None:
        p = product
        self.conn.execute(
            "INSERT INTO products(product_id, run_id, baseline_version,"
            " fingerprint, created_at) VALUES (?,?,?,?,?)",
            (p.product_id, p.run_id, p.baseline_version,
             p.fingerprint(), _now()))
        for kind, plf in (("signal", p.signal),
                          ("baseline", p.baseline),
                          ("net", p.net)):
            self.conn.executemany(
                "INSERT INTO product_segments(product_id, kind, idx,"
                " t_num, t_den, y_num, y_den) VALUES (?,?,?,?,?,?,?)",
                [(p.product_id, kind, i, k.t.numerator, k.t.denominator,
                  k.y.numerator, k.y.denominator)
                 for i, k in enumerate(plf.knots)])
        self.conn.commit()

    def get_product(self, product_id: str) -> SegmentProduct:
        row = self.conn.execute(
            "SELECT * FROM products WHERE product_id=?",
            (product_id,)).fetchone()
        if row is None:
            raise NotFoundError(f"未知产品: {product_id}")
        loaded = self._load_product(product_id)
        # 复读校验：用原采样与锚点重新计算净响应，必须与存储的切段一致
        rebuilt = build_product(
            product_id, loaded.run_id, loaded.baseline_version,
            loaded.signal, loaded.baseline)
        if rebuilt.net.fingerprint() != loaded.net.fingerprint():
            raise ValidationError("产品复读校验失败：净响应无法由快照复算")
        return loaded

    def _load_product(self, product_id: str) -> SegmentProduct:
        rows = self.conn.execute(
            "SELECT * FROM product_segments WHERE product_id=?"
            " ORDER BY kind, idx", (product_id,)).fetchall()
        groups: dict[str, list[Point]] = {"signal": [], "baseline": [],
                                          "net": []}
        for r in rows:
            groups[r["kind"]].append(Point(
                get_frac(r, "t_num", "t_den"),
                get_frac(r, "y_num", "y_den")))
        meta = self.conn.execute(
            "SELECT * FROM products WHERE product_id=?",
            (product_id,)).fetchone()
        return SegmentProduct(
            product_id=product_id,
            run_id=meta["run_id"],
            baseline_version=int(meta["baseline_version"]),
            signal=make_plf(groups["signal"]),
            baseline=make_plf(groups["baseline"]),
            net=make_plf(groups["net"]),
        )

    def list_products(self, run_id: str | None = None) -> list[str]:
        if run_id is None:
            rows = self.conn.execute(
                "SELECT product_id FROM products ORDER BY created_at"
            ).fetchall()
        else:
            rows = self.conn.execute(
                "SELECT product_id FROM products WHERE run_id=?"
                " ORDER BY created_at", (run_id,)).fetchall()
        return [r["product_id"] for r in rows]
