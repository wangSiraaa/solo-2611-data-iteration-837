"""双图带拼接缝：动态规划选缝（纯函数、精确有理代价）。

印样端独立导入两条**已定位**图带（左图、右图），只在连续重叠列上拼接：

* 每行在重叠列中选一个两侧均有值、且未被印样单禁用的位置作接缝；
  缝列归左图，其后的列用右图；
* 相邻两行的缝列之差至多 1（连通接缝）；
* 最小化全部缝列上两图读数之差绝对值之和；
* 平局取从首行起字典序最小的列号序列；
* 某行找不到可与上一行接续的合法接缝时，指出首个无法接续的行；
* 任一图负责的区域若有缺测则保持缺测，不用另一图或邻列填补。

规模：至多 80 行、每行 80 个重叠列，DP 状态数 ≤ 6400，直接逐行递推。
所有读数为 ``Fraction``，代价比较无浮点；字典序平局用完整缝列路径
元组比较，精确满足「从首行起字典序最小」。
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from .errors import ValidationError
from .motion import check_scan_limits

# 拼图单元格状态
CELL_LEFT = "left"     # 取自左图（含缝列）
CELL_RIGHT = "right"   # 取自右图
CELL_MISSING = "missing"


@dataclass(frozen=True)
class LocatedStrip:
    """一条已定位图带：每行给出「实物列号 → 读数」，缺测列不出现在映射中。"""

    strip_id: str
    rows: tuple[dict[int, Fraction], ...]

    def __post_init__(self) -> None:
        if not self.strip_id:
            raise ValidationError("strip_id 不能为空")


@dataclass(frozen=True)
class PrintSheet:
    """印样单：禁用的 (行号, 重叠列) 位置；换单不影响电机运动与图带本身。"""

    sheet_id: str
    n_rows: int
    overlap_columns: tuple[int, ...]          # 连续整数列，严格递增
    disabled: frozenset[tuple[int, int]]

    def __post_init__(self) -> None:
        cols = self.overlap_columns
        if len(cols) < 1:
            raise ValidationError("至少需要一个重叠列")
        for a, b in zip(cols, cols[1:]):
            if b != a + 1:
                raise ValidationError("重叠列必须是连续整数列")
        check_scan_limits(self.n_rows, len(cols))
        for r, c in self.disabled:
            if not (0 <= r < self.n_rows):
                raise ValidationError(f"禁用位置行号越界: {r}")
            if c not in cols:
                raise ValidationError(f"禁用位置列 {c} 不在重叠列内")


@dataclass(frozen=True)
class CompositeCell:
    column: int
    status: str
    value: Fraction | None
    source_strip: str | None     # "left" / "right"；缺测为 None


@dataclass(frozen=True)
class RowSeam:
    row_index: int
    seam_column: int
    cost: Fraction               # 该行缝列上 |左-右|


@dataclass(frozen=True)
class SeamResult:
    feasible: bool
    seams: tuple[RowSeam, ...]
    total_cost: Fraction
    composite: tuple[tuple[CompositeCell, ...], ...]
    first_blocked_row: int | None     # 首个无法接续的行（0 基）；可行时为 None


def legal_columns(row_index: int, sheet: PrintSheet,
                  left_row: dict[int, Fraction],
                  right_row: dict[int, Fraction]) -> tuple[int, ...]:
    """该行重叠列中两侧均有值且未被印样单禁用的缝候选列。"""
    out: list[int] = []
    for c in sheet.overlap_columns:
        if (row_index, c) in sheet.disabled:
            continue
        if c in left_row and c in right_row:
            out.append(c)
    return tuple(out)


def find_seams(left: LocatedStrip, right: LocatedStrip,
               sheet: PrintSheet,
               full_columns: Sequence[int] | None = None) -> SeamResult:
    """求最优连通接缝并输出实际拼图。

    ``full_columns`` 给出拼图输出的实物列范围（默认取重叠列）。
    缝列及其左侧用左图、右侧用右图；被选来源缺测的单元格保持缺测。
    """
    n = sheet.n_rows
    if len(left.rows) != n or len(right.rows) != n:
        raise ValidationError("两条图带的行数必须与印样单一致")
    check_scan_limits(n, len(sheet.overlap_columns))

    legals = [legal_columns(r, sheet, left.rows[r], right.rows[r])
              for r in range(n)]
    if not legals[0]:
        return _blocked(0, sheet, full_columns)

    # DP：状态 = 缝列 → (最小累计代价, 字典序最小的完整缝列路径)
    # 元组 (cost, path) 直接比较：先比代价，平局逐列字典序。
    dp: dict[int, tuple[Fraction, tuple[int, ...]]] = {}
    for c in legals[0]:
        cost = abs(left.rows[0][c] - right.rows[0][c])
        dp[c] = (cost, (c,))

    last_feasible_row = 0
    for r in range(1, n):
        nxt: dict[int, tuple[Fraction, tuple[int, ...]]] = {}
        for c in legals[r]:
            step = abs(left.rows[r][c] - right.rows[r][c])
            best: tuple[Fraction, tuple[int, ...]] | None = None
            for prev, (prev_cost, prev_path) in dp.items():
                if abs(prev - c) > 1:
                    continue
                cand = (prev_cost + step, prev_path + (c,))
                if best is None or cand < best:
                    best = cand
            if best is not None:
                nxt[c] = best
        if not nxt:
            return _blocked(r, sheet, full_columns,
                            seams_so_far=_best_path(dp),
                            left=left, right=right)
        dp = nxt
        last_feasible_row = r

    best_cost, best_path = _best_state(dp)
    seams = tuple(
        RowSeam(row_index=r, seam_column=best_path[r],
                cost=abs(left.rows[r][best_path[r]]
                         - right.rows[r][best_path[r]]))
        for r in range(n))
    composite = _compose(left, right, best_path, sheet, full_columns)
    return SeamResult(
        feasible=True, seams=seams, total_cost=best_cost,
        composite=composite, first_blocked_row=None)


def _best_state(dp: dict[int, tuple[Fraction, tuple[int, ...]]]
                ) -> tuple[Fraction, tuple[int, ...]]:
    return min(dp.values())


def _best_path(dp: dict[int, tuple[Fraction, tuple[int, ...]]]
               ) -> tuple[int, ...]:
    return _best_state(dp)[1]


def _blocked(row_index: int, sheet: PrintSheet,
             full_columns: Sequence[int] | None,
             seams_so_far: tuple[int, ...] = (),
             left: LocatedStrip | None = None,
             right: LocatedStrip | None = None) -> SeamResult:
    """构造不可行结果：指出首个无法接续的行，已确定行仍给出缝与代价。"""
    seams = tuple(
        RowSeam(row_index=r, seam_column=seams_so_far[r],
                cost=abs(left.rows[r][seams_so_far[r]]
                         - right.rows[r][seams_so_far[r]]))
        for r in range(len(seams_so_far)))
    total = sum((s.cost for s in seams), Fraction(0))
    return SeamResult(
        feasible=False, seams=seams, total_cost=total,
        composite=(), first_blocked_row=row_index)


def _compose(left: LocatedStrip, right: LocatedStrip,
             path: tuple[int, ...], sheet: PrintSheet,
             full_columns: Sequence[int] | None
             ) -> tuple[tuple[CompositeCell, ...], ...]:
    cols = tuple(full_columns) if full_columns is not None \
        else sheet.overlap_columns
    out_rows: list[tuple[CompositeCell, ...]] = []
    for r, seam in enumerate(path):
        cells: list[CompositeCell] = []
        for c in cols:
            if c <= seam:
                # 缝列归左图
                if c in left.rows[r]:
                    cells.append(CompositeCell(
                        c, CELL_LEFT, left.rows[r][c], left.strip_id))
                else:
                    cells.append(CompositeCell(c, CELL_MISSING, None, None))
            else:
                if c in right.rows[r]:
                    cells.append(CompositeCell(
                        c, CELL_RIGHT, right.rows[r][c], right.strip_id))
                else:
                    cells.append(CompositeCell(c, CELL_MISSING, None, None))
        out_rows.append(tuple(cells))
    return tuple(out_rows)
