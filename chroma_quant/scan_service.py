"""往返扫描与印样拼接服务（SQLite 持久化，追加式不可变历史）。

职责
----
* 管理电机行程（u(t)、g、起始台面位置）与台面轨迹；
* 管理刻度原点订正版本——订正只新增偏移版本，**不重算电机运动**；
* 管理扫行说明（哪条行程的哪次单调行程、哪版原点、哪个净响应产品）；
* 管理印样单（重叠列与禁用位置）与已定位图带——**换印样单只是换约束，
  既不触碰行程也不触碰图带**；
* 编排 :mod:`motion` 的逐列重建与 :mod:`stitch` 的接缝动态规划。

净响应产品复用既有 ``products`` / ``product_segments`` 表（kind='net'），
本服务只读取公开净响应分段，不改变旧定量命令的任何行为。
"""
from __future__ import annotations

import sqlite3
from datetime import datetime, timezone
from fractions import Fraction
from typing import Sequence

from . import limits
from .errors import LimitError, NotFoundError, ValidationError
from .motion import (OriginVersion, ReconstructedRow, RowSpec, Stroke,
                     build_stage_trajectory, check_motor_knots,
                     check_scan_limits, reconstruct_rows,
                     split_monotonic_runs)
from .rational import to_fraction
from .segments import PLF, make_plf
from .stitch import (LocatedStrip, PrintSheet, SeamResult, find_seams)
from .store import connect, get_frac, init_db


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class ScanService:
    def __init__(self, db_path: str):
        self.conn: sqlite3.Connection = connect(db_path)
        init_db(self.conn)
        self._traj_cache: dict[str, object] = {}

    def close(self) -> None:
        self.conn.close()

    # ------------------------------------------------------------------
    # 电机行程
    # ------------------------------------------------------------------
    def create_stroke(self, stroke_id: str,
                      motor_points: Sequence[tuple], g,
                      initial_stage, note: str | None = None) -> Stroke:
        if self._stroke_exists(stroke_id):
            raise ValidationError(f"stroke_id 已存在: {stroke_id}")
        motor = make_plf(motor_points)
        check_motor_knots(len(motor.knots))
        stroke = Stroke(stroke_id=stroke_id, motor=motor,
                        backlash=to_fraction(g, "g"),
                        initial_stage=to_fraction(initial_stage,
                                                  "initial_stage"))
        self.conn.execute(
            "INSERT INTO scan_strokes(stroke_id, g_num, g_den, s0_num,"
            " s0_den, created_at, note) VALUES (?,?,?,?,?,?,?)",
            (stroke_id, stroke.backlash.numerator,
             stroke.backlash.denominator,
             stroke.initial_stage.numerator,
             stroke.initial_stage.denominator, _now(), note))
        self.conn.executemany(
            "INSERT INTO scan_motor_knots(stroke_id, idx, t_num, t_den,"
            " u_num, u_den) VALUES (?,?,?,?,?,?)",
            [(stroke_id, i, k.t.numerator, k.t.denominator,
              k.y.numerator, k.y.denominator)
             for i, k in enumerate(motor.knots)])
        self.conn.commit()
        # 首版原点偏移为 0（实物坐标 = 台面坐标）
        self._insert_origin(stroke_id, Fraction(0), supersedes=None)
        return stroke

    def _stroke_exists(self, stroke_id: str) -> bool:
        return self.conn.execute(
            "SELECT 1 FROM scan_strokes WHERE stroke_id=?",
            (stroke_id,)).fetchone() is not None

    def get_stroke(self, stroke_id: str) -> Stroke:
        row = self.conn.execute(
            "SELECT * FROM scan_strokes WHERE stroke_id=?",
            (stroke_id,)).fetchone()
        if row is None:
            raise NotFoundError(f"未知行程: {stroke_id}")
        knots = self.conn.execute(
            "SELECT * FROM scan_motor_knots WHERE stroke_id=?"
            " ORDER BY idx", (stroke_id,)).fetchall()
        motor = make_plf(
            (get_frac(r, "t_num", "t_den"),
             get_frac(r, "u_num", "u_den")) for r in knots)
        return Stroke(
            stroke_id=stroke_id, motor=motor,
            backlash=Fraction(row["g_num"], row["g_den"]),
            initial_stage=Fraction(row["s0_num"], row["s0_den"]))

    def list_strokes(self) -> list[str]:
        return [r["stroke_id"] for r in self.conn.execute(
            "SELECT stroke_id FROM scan_strokes ORDER BY created_at")]

    def list_monotonic_runs(self, stroke_id: str):
        return split_monotonic_runs(self.get_stroke(stroke_id))

    def get_stage_trajectory(self, stroke_id: str):
        if stroke_id not in self._traj_cache:
            stroke = self.get_stroke(stroke_id)
            self._traj_cache[stroke_id] = build_stage_trajectory(stroke)
        return self._traj_cache[stroke_id]

    # ------------------------------------------------------------------
    # 刻度原点订正（只追加；台面轨迹与偏移无关，立即返回且不重算运动）
    # ------------------------------------------------------------------
    def _insert_origin(self, stroke_id: str, offset: Fraction,
                       supersedes: int | None) -> int:
        ver_row = self.conn.execute(
            "SELECT COALESCE(MAX(version),0) AS v FROM scan_origin_versions"
            " WHERE stroke_id=?", (stroke_id,)).fetchone()
        version = int(ver_row["v"]) + 1
        if version > limits.MAX_ORIGIN_VERSIONS:
            raise LimitError(
                f"原点版本上限 {limits.MAX_ORIGIN_VERSIONS} 已达")
        self.conn.execute(
            "INSERT INTO scan_origin_versions(stroke_id, version, off_num,"
            " off_den, created_at, supersedes) VALUES (?,?,?,?,?,?)",
            (stroke_id, version, offset.numerator, offset.denominator,
             _now(), supersedes))
        self.conn.commit()
        return version

    def correct_origin(self, stroke_id: str, offset) -> int:
        """追加一个原点偏移版本（实物坐标 = 台面位置 + offset）。"""
        self.get_stroke(stroke_id)  # 校验存在
        versions = self.list_origin_versions(stroke_id)
        return self._insert_origin(stroke_id,
                                   to_fraction(offset, "offset"),
                                   supersedes=versions[-1])

    def list_origin_versions(self, stroke_id: str) -> list[int]:
        self.get_stroke(stroke_id)
        return [int(r["version"]) for r in self.conn.execute(
            "SELECT version FROM scan_origin_versions WHERE stroke_id=?"
            " ORDER BY version", (stroke_id,))]

    def get_origin_version(self, stroke_id: str,
                           version: int) -> OriginVersion:
        row = self.conn.execute(
            "SELECT * FROM scan_origin_versions WHERE stroke_id=?"
            " AND version=?", (stroke_id, version)).fetchone()
        if row is None:
            raise NotFoundError(
                f"原点版本不存在: stroke={stroke_id} version={version}")
        return OriginVersion(
            stroke_id=stroke_id, version=version,
            offset=Fraction(row["off_num"], row["off_den"]))

    # ------------------------------------------------------------------
    # 扫行说明
    # ------------------------------------------------------------------
    def add_scan_row(self, row_id: str, stroke_id: str, run_index: int,
                     origin_version: int,
                     product_id: str | None = None,
                     baseline_version: int | None = None) -> None:
        if self.conn.execute("SELECT 1 FROM scan_rows WHERE row_id=?",
                             (row_id,)).fetchone():
            raise ValidationError(f"row_id 已存在: {row_id}")
        stroke = self.get_stroke(stroke_id)
        self.get_origin_version(stroke_id, origin_version)  # 校验版本
        runs = split_monotonic_runs(stroke)
        if not (0 <= run_index < len(runs)):
            raise ValidationError(
                f"行程 {stroke_id} 单调行程编号 {run_index} 越界"
                f"（共 {len(runs)} 次）")
        if product_id is not None and not self._product_exists(product_id):
            raise NotFoundError(f"净响应产品不存在: {product_id}")
        self.conn.execute(
            "INSERT INTO scan_rows(row_id, stroke_id, run_index,"
            " origin_version, product_id, baseline_version, created_at)"
            " VALUES (?,?,?,?,?,?,?)",
            (row_id, stroke_id, run_index, origin_version,
             product_id, baseline_version, _now()))
        self.conn.commit()

    def _product_exists(self, product_id: str) -> bool:
        return self.conn.execute(
            "SELECT 1 FROM products WHERE product_id=?",
            (product_id,)).fetchone() is not None

    def list_scan_rows(self) -> list[RowSpec]:
        rows = self.conn.execute(
            "SELECT * FROM scan_rows ORDER BY created_at").fetchall()
        return [RowSpec(
            row_id=r["row_id"], stroke_id=r["stroke_id"],
            run_index=int(r["run_index"]),
            origin_version=int(r["origin_version"]),
            product_id=r["product_id"],
            baseline_version=(int(r["baseline_version"])
                              if r["baseline_version"] is not None else None))
            for r in rows]

    # ------------------------------------------------------------------
    # 重建：逐整数列采样，附完整出处
    # ------------------------------------------------------------------
    def reconstruct(self, columns_by_row: dict[str, Sequence[int]]
                    ) -> tuple[ReconstructedRow, ...]:
        """按列计划重建其中点名的扫行（不必覆盖库里全部扫行）。"""
        all_specs = self.list_scan_rows()
        if not all_specs:
            raise ValidationError("尚无扫行说明，无法重建")
        known_ids = {s.row_id for s in all_specs}
        missing = [rid for rid in columns_by_row if rid not in known_ids]
        if missing:
            raise ValidationError(f"列计划引用了未知扫行: {missing}")
        specs = [s for s in all_specs
                 if s.row_id in columns_by_row]
        n_overlap = max((len(v) for v in columns_by_row.values()),
                        default=0)
        check_scan_limits(len(specs), n_overlap)
        strokes: list[Stroke] = []
        trajectories = {}
        for sid in {s.stroke_id for s in specs}:
            stroke = self.get_stroke(sid)
            strokes.append(stroke)
            trajectories[sid] = self.get_stage_trajectory(sid)
        origins = {}
        for sid in {s.stroke_id for s in specs}:
            for v in self.list_origin_versions(sid):
                ov = self.get_origin_version(sid, v)
                origins[(sid, v)] = ov
        nets: dict[str, PLF] = {}
        for pid in {s.product_id for s in specs
                    if s.product_id is not None}:
            nets[pid] = self._load_net(pid)
        return reconstruct_rows(strokes, trajectories, origins, nets,
                                specs, columns_by_row)

    def _load_net(self, product_id: str) -> PLF:
        rows = self.conn.execute(
            "SELECT * FROM product_segments WHERE product_id=?"
            " AND kind='net' ORDER BY idx", (product_id,)).fetchall()
        if not rows:
            raise NotFoundError(f"产品 {product_id} 缺少净响应分段")
        return make_plf(
            (get_frac(r, "t_num", "t_den"),
             get_frac(r, "y_num", "y_den")) for r in rows)

    # ------------------------------------------------------------------
    # 印样单：重叠列范围 + 禁用位置（换单不触发运动重算）
    # ------------------------------------------------------------------
    def create_print_sheet(self, sheet_id: str, n_rows: int,
                           col_lo: int, col_hi: int,
                           disabled: Sequence[tuple[int, int]] = ()
                           ) -> PrintSheet:
        if self.conn.execute(
                "SELECT 1 FROM print_sheets WHERE sheet_id=?",
                (sheet_id,)).fetchone():
            raise ValidationError(f"sheet_id 已存在: {sheet_id}")
        cols = tuple(range(int(col_lo), int(col_hi) + 1))
        sheet = PrintSheet(
            sheet_id=sheet_id, n_rows=int(n_rows),
            overlap_columns=cols, disabled=frozenset(
                (int(r), int(c)) for r, c in disabled))
        self.conn.execute(
            "INSERT INTO print_sheets(sheet_id, n_rows, col_lo, col_hi,"
            " created_at) VALUES (?,?,?,?,?)",
            (sheet_id, sheet.n_rows, cols[0], cols[-1], _now()))
        self.conn.executemany(
            "INSERT INTO print_sheet_disabled(sheet_id, row_idx, col)"
            " VALUES (?,?,?)",
            [(sheet_id, r, c) for r, c in sorted(sheet.disabled)])
        self.conn.commit()
        return sheet

    def get_print_sheet(self, sheet_id: str) -> PrintSheet:
        row = self.conn.execute(
            "SELECT * FROM print_sheets WHERE sheet_id=?",
            (sheet_id,)).fetchone()
        if row is None:
            raise NotFoundError(f"未知印样单: {sheet_id}")
        disabled = frozenset(
            (int(r["row_idx"]), int(r["col"]))
            for r in self.conn.execute(
                "SELECT row_idx, col FROM print_sheet_disabled"
                " WHERE sheet_id=?", (sheet_id,)))
        cols = tuple(range(int(row["col_lo"]), int(row["col_hi"]) + 1))
        return PrintSheet(sheet_id=sheet_id, n_rows=int(row["n_rows"]),
                          overlap_columns=cols, disabled=disabled)

    # ------------------------------------------------------------------
    # 已定位图带：可独立导入，不依赖行程；缺测列直接缺行
    # ------------------------------------------------------------------
    def import_located_strip(self, strip_id: str,
                             rows: Sequence[dict[int, Fraction]]) -> LocatedStrip:
        if self.conn.execute(
                "SELECT 1 FROM located_strips WHERE strip_id=?",
                (strip_id,)).fetchone():
            raise ValidationError(f"strip_id 已存在: {strip_id}")
        check_scan_limits(len(rows), max((len(m) for m in rows), default=0))
        normalized: list[dict[int, Fraction]] = []
        payload = []
        for r, mapping in enumerate(rows):
            clean = {}
            for c, y in mapping.items():
                col = int(c)
                clean[col] = to_fraction(y, "y")
                payload.append((strip_id, r, col,
                                clean[col].numerator, clean[col].denominator))
            normalized.append(clean)
        self.conn.executemany(
            "INSERT INTO located_strips(strip_id, row_idx, col, y_num,"
            " y_den) VALUES (?,?,?,?,?)", payload)
        self.conn.commit()
        return LocatedStrip(strip_id=strip_id, rows=tuple(normalized))

    def get_located_strip(self, strip_id: str) -> LocatedStrip:
        if self.conn.execute(
                "SELECT 1 FROM located_strips WHERE strip_id=?",
                (strip_id,)).fetchone() is None:
            raise NotFoundError(f"未知图带: {strip_id}")
        n_rows = 1 + int(self.conn.execute(
            "SELECT COALESCE(MAX(row_idx),0) AS m FROM located_strips"
            " WHERE strip_id=?", (strip_id,)).fetchone()["m"])
        rows: list[dict[int, Fraction]] = [dict() for _ in range(n_rows)]
        for r in self.conn.execute(
                "SELECT * FROM located_strips WHERE strip_id=?"
                " ORDER BY row_idx, col", (strip_id,)):
            rows[int(r["row_idx"])][int(r["col"])] = get_frac(
                r, "y_num", "y_den")
        return LocatedStrip(strip_id=strip_id, rows=tuple(rows))

    # ------------------------------------------------------------------
    # 接缝拼接（纯内存 DP；换印样单只重跑这里，电机运动完全不参与）
    # ------------------------------------------------------------------
    def stitch(self, sheet_id: str, left_strip_id: str,
               right_strip_id: str,
               full_columns: Sequence[int] | None = None) -> SeamResult:
        sheet = self.get_print_sheet(sheet_id)
        left = self.get_located_strip(left_strip_id)
        right = self.get_located_strip(right_strip_id)
        return find_seams(left, right, sheet, full_columns)
