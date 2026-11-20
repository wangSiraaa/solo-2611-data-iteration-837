"""组分定量仓：消费公开分段产品，管理方法、分析缓存与复读。

与 :class:`~chroma_quant.signal_service.SignalService` 物理上共用一个
SQLite 文件、逻辑上互不调用：仓库只认识分段产品，因此也能分析
``product_from_dict`` 读入的外部等价分段产品。
"""
from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from fractions import Fraction
from typing import Sequence

from . import limits
from .errors import LimitError, NotFoundError, ValidationError
from .quant import (ABS_TOL, AnalysisResult, MethodSpec, SegmentProduct,
                    WeightFunction, analyze, make_weight_function,
                    product_from_dict, product_to_dict)
from .rational import frac_str, stable_id
from .segments import Point, make_plf
from .store import connect, get_frac, init_db


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class QuantRepository:
    def __init__(self, db_path: str):
        self.conn: sqlite3.Connection = connect(db_path)
        init_db(self.conn)

    def close(self) -> None:
        self.conn.close()

    # ------------------------------------------------------------------
    # 方法
    # ------------------------------------------------------------------
    def create_method(self, method_id: str,
                      components: Sequence[dict],
                      note: str | None = None) -> str:
        """components: [{"name","factor","nodes":[(t,w),...]}]，顺序即身份。"""
        if self.conn.execute("SELECT 1 FROM methods WHERE method_id=?",
                             (method_id,)).fetchone():
            raise ValidationError(f"method_id 已存在: {method_id}")
        if len(components) > limits.MAX_COMPONENTS_PER_METHOD:
            raise LimitError(
                f"每方法最多 {limits.MAX_COMPONENTS_PER_METHOD} 个组分，"
                f"收到 {len(components)}")

        wfs: list[WeightFunction] = []
        for comp in components:
            wf = make_weight_function(
                comp["name"], comp["factor"], comp["nodes"])
            if len(wf.nodes) > limits.MAX_NODES_PER_COMPONENT:
                raise LimitError(
                    f"组分 {wf.name} 权重节点最多 "
                    f"{limits.MAX_NODES_PER_COMPONENT} 个")
            wfs.append(wf)
        spec = MethodSpec(method_id=method_id, components=tuple(wfs))

        self.conn.execute(
            "INSERT INTO methods(method_id, created_at, note)"
            " VALUES (?,?,?)", (method_id, _now(), note))
        for i, wf in enumerate(spec.components):
            self.conn.execute(
                "INSERT INTO method_components(method_id, idx, name,"
                " factor_num, factor_den) VALUES (?,?,?,?,?)",
                (method_id, i, wf.name,
                 wf.factor.numerator, wf.factor.denominator))
            self.conn.executemany(
                "INSERT INTO weight_nodes(method_id, comp_idx, idx,"
                " t_num, t_den, w_num, w_den) VALUES (?,?,?,?,?,?,?)",
                [(method_id, i, j, p.t.numerator, p.t.denominator,
                  p.y.numerator, p.y.denominator)
                 for j, p in enumerate(wf.nodes)])
        self.conn.commit()
        return method_id

    def get_method(self, method_id: str) -> MethodSpec:
        row = self.conn.execute(
            "SELECT 1 FROM methods WHERE method_id=?",
            (method_id,)).fetchone()
        if row is None:
            raise NotFoundError(f"未知方法: {method_id}")
        comp_rows = self.conn.execute(
            "SELECT * FROM method_components WHERE method_id=?"
            " ORDER BY idx", (method_id,)).fetchall()
        comps: list[WeightFunction] = []
        for cr in comp_rows:
            node_rows = self.conn.execute(
                "SELECT * FROM weight_nodes WHERE method_id=? AND comp_idx=?"
                " ORDER BY idx", (method_id, cr["idx"])).fetchall()
            nodes = tuple(
                Point(get_frac(n, "t_num", "t_den"),
                      get_frac(n, "w_num", "w_den"))
                for n in node_rows)
            factor = get_frac(cr, "factor_num", "factor_den")
            comps.append(WeightFunction(
                name=cr["name"], factor=factor, nodes=nodes))
        return MethodSpec(method_id=method_id, components=tuple(comps))

    def list_methods(self) -> list[str]:
        return [r["method_id"] for r in self.conn.execute(
            "SELECT method_id FROM methods ORDER BY created_at")]

    # ------------------------------------------------------------------
    # 外部等价分段产品：只存公开分段，可直接来自 JSON
    # ------------------------------------------------------------------
    def import_product(self, product_dict: dict) -> str:
        product = product_from_dict(product_dict)
        # 复算校验：外部净响应必须与其原信号、基线分段一致
        from .quant import build_product
        rebuilt = build_product(
            product.product_id, product.run_id,
            product.baseline_version, product.signal, product.baseline)
        if rebuilt.net.fingerprint() != product.net.fingerprint():
            raise ValidationError(
                "外部产品校验失败：净响应分段与信号/基线不一致")
        existing = self.conn.execute(
            "SELECT product_id FROM products WHERE fingerprint=?",
            (product.fingerprint(),)).fetchall()
        if existing:
            # 同函数产品等价：复用第一条，返回其 id
            return existing[0]["product_id"]
        self._persist(product)
        return product.product_id

    def _persist(self, product: SegmentProduct) -> None:
        p = product
        self.conn.execute(
            "INSERT OR IGNORE INTO runs(run_id, created_at, note)"
            " VALUES (?,?,?)", (p.run_id, _now(), "外部产品引用"))
        self.conn.execute(
            "INSERT OR IGNORE INTO baseline_versions(run_id, version,"
            " created_at) VALUES (?,?,?)",
            (p.run_id, p.baseline_version, _now()))
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

    def load_product(self, product_id: str) -> SegmentProduct:
        rows = self.conn.execute(
            "SELECT * FROM product_segments WHERE product_id=?"
            " ORDER BY kind, idx", (product_id,)).fetchall()
        if not rows:
            raise NotFoundError(f"未知产品: {product_id}")
        groups: dict[str, list] = {"signal": [], "baseline": [], "net": []}
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

    # ------------------------------------------------------------------
    # 分析：缓存键 = 方法 × 净响应函数指纹（同函数不同分段命中同一缓存）
    # ------------------------------------------------------------------
    def run_analysis(self, method_id: str, product_id: str,
                     *, use_cache: bool = True) -> AnalysisResult:
        method = self.get_method(method_id)
        product = self.load_product(product_id)
        fp = product.fingerprint()

        if use_cache:
            hit = self.conn.execute(
                "SELECT analysis_id FROM analyses WHERE method_id=?"
                " AND fingerprint=? ORDER BY created_at LIMIT 1",
                (method_id, fp)).fetchone()
            if hit is not None:
                return self.get_analysis(hit["analysis_id"])

        n = self.conn.execute(
            "SELECT COUNT(*) AS c FROM analyses").fetchone()["c"]
        if n >= limits.MAX_ANALYSES:
            raise LimitError(
                f"分析次数上限 {limits.MAX_ANALYSES} 已达；"
                "缓存命中不计入新分析")

        result = analyze(product, method)
        analysis_id = stable_id("analysis", method_id, fp)
        self._persist_analysis(analysis_id, method_id, product, fp, result)
        return result

    def _persist_analysis(self, analysis_id: str, method_id: str,
                          product: SegmentProduct, fp: str,
                          result: AnalysisResult) -> None:
        a, b = result.domain
        for name in result.component_area:
            self.conn.execute(
                "INSERT INTO analyses(analysis_id, method_id, product_id,"
                " created_at, fingerprint, component_name, area,"
                " concentration, unallocated_area, total_area,"
                " conservation_residual, tolerance, domain_t0_num,"
                " domain_t0_den, domain_t1_num, domain_t1_den)"
                " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (analysis_id, method_id, product.product_id, _now(), fp,
                 name, result.component_area[name],
                 result.concentration[name], result.unallocated_area,
                 result.total_area, result.conservation_residual,
                 result.tolerance, a.numerator, a.denominator,
                 b.numerator, b.denominator))
        self.conn.commit()

    def get_analysis(self, analysis_id: str) -> AnalysisResult:
        rows = self.conn.execute(
            "SELECT * FROM analyses WHERE analysis_id=? ORDER BY component_name",
            (analysis_id,)).fetchall()
        if not rows:
            raise NotFoundError(f"未知分析: {analysis_id}")
        r0 = rows[0]
        area = {r["component_name"]: r["area"] for r in rows}
        conc = {r["component_name"]: r["concentration"] for r in rows}
        dom = (Fraction(r0["domain_t0_num"], r0["domain_t0_den"]),
               Fraction(r0["domain_t1_num"], r0["domain_t1_den"]))
        return AnalysisResult(
            method_id=r0["method_id"],
            product_id=r0["product_id"],
            domain=dom,
            component_area=area,
            concentration=conc,
            unallocated_area=r0["unallocated_area"],
            total_area=r0["total_area"],
            conservation_residual=r0["conservation_residual"],
            tolerance=r0["tolerance"],
            trace=(),
        )

    def list_analyses(self) -> list[str]:
        return [r["analysis_id"] for r in self.conn.execute(
            "SELECT DISTINCT analysis_id FROM analyses"
            " ORDER BY created_at")]

    def find_analysis(self, method_id: str,
                      product_id: str) -> AnalysisResult:
        """按方法与产品取回分析（基于净响应函数指纹的缓存键）。"""
        product = self.load_product(product_id)
        row = self.conn.execute(
            "SELECT analysis_id FROM analyses WHERE method_id=?"
            " AND fingerprint=? ORDER BY created_at LIMIT 1",
            (method_id, product.fingerprint())).fetchone()
        if row is None:
            raise NotFoundError(
                f"无分析记录: method={method_id} product={product_id}")
        return self.get_analysis(row["analysis_id"])

    # ------------------------------------------------------------------
    # 解释入口：重新计算（不落库）以返回逐格轨迹
    # ------------------------------------------------------------------
    def explain(self, method_id: str, product_id: str) -> dict:
        method = self.get_method(method_id)
        product = self.load_product(product_id)
        result = analyze(product, method)
        cells = []
        for c in result.trace:
            cells.append({
                "t_lo": frac_str(c.t_lo),
                "t_hi": frac_str(c.t_hi),
                "net_integral": str(c.net_integral),
                "allocated": str(c.allocated_integral),
                "unallocated": str(c.unallocated_integral),
                "W_lo": frac_str(c.W_lo),
                "W_hi": frac_str(c.W_hi),
                "component_areas": {
                    k: repr(v) for k, v in c.component_areas.items()},
            })
        return {
            "method_id": method_id,
            "product_id": product_id,
            "domain": [frac_str(product.domain[0]),
                       frac_str(product.domain[1])],
            "component_order": method.names,
            "cells": cells,
            "summary": {
                "total_area": result.total_area,
                "unallocated_area": result.unallocated_area,
                "component_area": result.component_area,
                "concentration": result.concentration,
                "conservation_residual": result.conservation_residual,
            },
        }

    def export_product_json(self, product_id: str) -> str:
        return json.dumps(product_to_dict(self.load_product(product_id)),
                          ensure_ascii=False, indent=2)
