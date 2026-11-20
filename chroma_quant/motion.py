"""往返扫描的回差坐标重建（纯函数、全程 ``Fraction``，不碰数据库）。

问题背景
--------
直接把编码器列号（电机位置 u(t)）当作扫描位置，会在电机换向后的一段时间
内重复画同一点：机械台面通过非负回差半宽 g 的间隙与电机耦合，只有当
``[u-g, u+g]`` 的带边越过台面当前位置时台面才贴着该边移动，否则台面不动。

台面运动规则（对每个单调电机段统一成立）
----------------------------------------
设台面在前一时刻的位置为 s，电机位置为 u：

* 若 ``s < u - g``：台面落在带外（左侧），贴右边移动到 ``u - g``；
* 若 ``s > u + g``：台面落在带外（右侧），贴左边移动到 ``u + g``；
* 否则台面保持不动。

即一步更新 ``s ← min(max(s, u-g), u+g)``（对持续到达的电机位置逐点施加）。
验证：g=1、电机从 3 折返到 0、折返时台面在 2 —— 电机先退到 1（带左边
u+g=2 才碰到台面 2），电机继续退向 0 时台面才跟着回到 1、0。

实物坐标 = 台面位置 + 刻度原点偏移。一条行程的台面轨迹只由
（电机轨迹, g, 起始台面位置）决定，**原点偏移订正不需要重算电机运动**。

本模块同时负责：把电机轨迹分解为若干次「单调行程」（平台/换向处分界），
在每次行程的连续运动部分按实物坐标反查时刻，再从该扫行的净响应 PLF
插值读取目标整数列；纯停顿、缺覆盖、净响应无效一律给缺测，绝不造数。
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from . import limits
from .errors import DomainError, LimitError, ValidationError
from .rational import to_fraction
from .segments import PLF, Point, make_plf

# 列采样状态（写入 ColumnSample.status）
STATUS_OK = "ok"                 # 连续运动中反查、净响应有效
STATUS_DWELL = "dwell"           # 目标实物坐标只在纯停顿（台面不动）期间出现
STATUS_NO_COVERAGE = "no_coverage"        # 该行程连续运动部分不覆盖该列
STATUS_INVALID_NET = "invalid_net"        # 反查时刻落在净响应有效域之外

MISSING_STATUSES = frozenset(
    {STATUS_DWELL, STATUS_NO_COVERAGE, STATUS_INVALID_NET})


# ----------------------------------------------------------------------
# 数据结构
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Stroke:
    """一条扫描行程：分段线性电机位置 u(t)、回差半宽 g、起始台面位置。

    ``initial_stage`` 必须落在起始容许带 ``[u(t0)-g, u(t0)+g]`` 内。
    原点偏移不在此处：它属于可订正的刻度版本（见 OriginVersion），
    使原点订正与电机运动解耦。
    """

    stroke_id: str
    motor: PLF
    backlash: Fraction
    initial_stage: Fraction

    def __post_init__(self) -> None:
        if not self.stroke_id:
            raise ValidationError("stroke_id 不能为空")
        g = to_fraction(self.backlash, "backlash")
        if g < 0:
            raise ValidationError("回差半宽 g 必须非负")
        s0 = to_fraction(self.initial_stage, "initial_stage")
        u0 = self.motor.knots[0].y
        if not (u0 - g <= s0 <= u0 + g):
            raise ValidationError(
                f"起始台面位置 {s0} 不在起始容许带 "
                f"[{u0 - g}, {u0 + g}] 内")
        object.__setattr__(self, "backlash", g)
        object.__setattr__(self, "initial_stage", s0)


@dataclass(frozen=True)
class MonotonicRun:
    """电机轨迹上一次极大严格单调行程（平台与换向处分界）。

    ``direction`` 为 +1（电机递增）或 -1（电机递减）；
    ``motor_segment`` 是该行程的电机子 PLF（闭区间，端点与相邻行程共享）。
    """

    index: int
    direction: int
    motor_segment: PLF


@dataclass(frozen=True)
class StageTrajectory:
    """台面位置 s(t) 的分段线性轨迹（平台段显式保留）。"""

    stage_plf: PLF

    def __call__(self, t: Fraction) -> Fraction:
        return self.stage_plf(t)


@dataclass(frozen=True)
class OriginVersion:
    """刻度原点偏移的一个不可变版本：实物坐标 x = s(t) + offset。"""

    stroke_id: str
    version: int
    offset: Fraction

    def __post_init__(self) -> None:
        object.__setattr__(self, "offset",
                           to_fraction(self.offset, "offset"))


@dataclass(frozen=True)
class RowSpec:
    """一个扫行：取哪条行程的哪次单调行程、哪一版原点与净响应产品。"""

    row_id: str
    stroke_id: str
    run_index: int
    origin_version: int
    product_id: str
    baseline_version: int   # 净响应来源版本（出处记录用）


@dataclass(frozen=True)
class ColumnSample:
    """一个目标整数列的采样结果与完整出处。"""

    column: int
    status: str
    value: Fraction | None          # 缺测时为 None
    stroke_id: str
    run_index: int
    time: Fraction | None           # 反查时刻；缺测时为 None
    origin_version: int
    product_id: str | None
    baseline_version: int | None    # 净响应源版本


# ----------------------------------------------------------------------
# 台面轨迹：clamp 规则
# ----------------------------------------------------------------------
def build_stage_trajectory(stroke: Stroke) -> StageTrajectory:
    """按回差规则逐电机段推进台面，返回台面位置 s(t) 的 PLF。

    电机在相邻折点间线性移动，带边 ``u±g`` 也是线性的。台面在一段内
    只有两种状态：

    * 仍在带内：保持不动（平台）；
    * 已贴上某条带边：随该边线性移动。

    两种状态的切换时刻是「带边越过静止台面」的真实交点，该交点可能落在
    电机段内部，必须作为台面折点显式插入（与净响应真实切段同一思想）。
    逐段更新：

    * 上升段（u 递增）：若当前台面 s 低于终点右边 u1-g，则在右边追上 s
      的时刻前保持 s，之后贴右边移动到 u1-g；否则整段不动；
    * 下降段（u 递减）：对称地处理左边 u+g。
    """
    g = stroke.backlash
    knots = stroke.motor.knots
    s = stroke.initial_stage
    out: list[Point] = [Point(knots[0].t, s)]
    for p, q in zip(knots, knots[1:]):
        dt = q.t - p.t
        du = q.y - p.y
        if du == 0:
            # 电机平台：带不动，台面不可能在段内新贴上带边
            out.append(Point(q.t, s))
            continue
        if du > 0:
            edge_lo = p.y - g          # 段起点右边
            edge_hi = q.y - g          # 段终点右边
            if s < edge_hi:
                # 右边将越过 s：求交点 r∈(0,1]，交点前静止、之后贴边
                if s <= edge_lo:
                    r = Fraction(0)    # 段起点已经（或仍）贴边
                else:
                    r = (s - edge_lo) / (edge_hi - edge_lo)
                if 0 < r < 1:
                    out.append(Point(p.t + dt * r, s))
                s = edge_hi
            # s >= edge_hi：整段台面都在带内，不动
        else:
            edge_lo = p.y + g          # 段起点左边
            edge_hi = q.y + g          # 段终点左边（递减，故 edge_hi < edge_lo）
            if s > edge_hi:
                if s >= edge_lo:
                    r = Fraction(0)
                else:
                    # 左边从 edge_lo 降到 edge_hi，求其等于 s 的位置
                    r = (edge_lo - s) / (edge_lo - edge_hi)
                if 0 < r < 1:
                    out.append(Point(p.t + dt * r, s))
                s = edge_hi
        out.append(Point(q.t, s))
    return StageTrajectory(make_plf(out))


def split_monotonic_runs(stroke: Stroke) -> tuple[MonotonicRun, ...]:
    """把电机轨迹分解为极大严格单调行程。

    相邻电机折点间斜率为正/负的归入当前行程，斜率为零（平台）或变号
    （换向）时结束当前行程；平台结束后以平台末端为起点开启新行程。
    首尾的纯平台不产生行程（没有单调运动）。
    """
    knots = stroke.motor.knots
    runs: list[MonotonicRun] = []
    start = 0
    direction = 0
    for i in range(1, len(knots)):
        dy = knots[i].y - knots[i - 1].y
        if dy == 0:
            sign = 0
        else:
            sign = 1 if dy > 0 else -1
        if sign == 0:
            # 平台：若已有进行中的行程，在平台起点闭合它
            if direction != 0:
                runs.append(_make_run(len(runs), direction,
                                      knots[start:i]))
                direction = 0
            start = i  # 平台末端作为下一行程的候选起点
        elif direction == 0:
            direction = sign
            start = i - 1
        elif sign != direction:
            # 换向：在上一折点闭合当前行程，开启新行程
            runs.append(_make_run(len(runs), direction,
                                  knots[start:i]))
            direction = sign
            start = i - 1
    if direction != 0:
        runs.append(_make_run(len(runs), direction, knots[start:]))
    return tuple(runs)


def _make_run(index: int, direction: int,
              knots: Sequence[Point]) -> MonotonicRun:
    return MonotonicRun(index=index, direction=direction,
                        motor_segment=make_plf(knots))


# ----------------------------------------------------------------------
# 单调行程上的连续运动区间与实物坐标反查
# ----------------------------------------------------------------------
def run_motion_interval(stroke: Stroke, trajectory: StageTrajectory,
                        run: MonotonicRun
                        ) -> tuple[Fraction, Fraction, int] | None:
    """返回该单调行程中台面「连续运动部分」的 (t_lo, t_hi, direction)。

    台面在该行程内沿运动方向单调；行程开头若存在回差带动的纯停顿
    （台面平台），其末端才是连续运动的起点。整段台面都不动（纯停顿
    行程）返回 None —— 这种行程不能为任何整数列造数。
    """
    t_lo, t_hi = run.motor_segment.domain
    s_lo = trajectory(t_lo)
    s_hi = trajectory(t_hi)
    if s_lo == s_hi:
        return None
    # 台面 s(t) 是 PLF；找到最后一个取值等于 s_lo 的折点时刻，
    # 即开头平台段（若有）的末端、连续运动的起点。
    stage_knots = trajectory.stage_plf.knots
    t_motion_start = t_lo
    for k in stage_knots:
        if t_lo <= k.t <= t_hi and k.y == s_lo:
            t_motion_start = k.t
        if k.t > t_hi:
            break
    return t_motion_start, t_hi, run.direction


def invert_monotone(plf: PLF, value: Fraction,
                    lo: Fraction, hi: Fraction) -> Fraction:
    """在严格单调（按当前折点值）的 PLF 子区间 [lo,hi] 上反查 value 的时刻。

    value 必须落在 ``[plf(lo), plf(hi)]``（含端点）；闭区间端点按
    「运动一侧」取值：递增行程取下界时刻、递减行程取上界时刻，由调用
    方保证 lo/hi 已经裁到连续运动部分。
    """
    value = to_fraction(value, "value")
    y_lo, y_hi = plf(lo), plf(hi)
    if not (min(y_lo, y_hi) <= value <= max(y_lo, y_hi)):
        raise DomainError(
            f"反查值 {value} 不在运动区间值域 [{min(y_lo, y_hi)},"
            f" {max(y_lo, y_hi)}] 内")
    knots = plf.knots
    # 定位包含 value 的线性子段（严格单调，折点值互异）
    for p, q in zip(knots, knots[1:]):
        if q.t < lo or p.t > hi:
            continue
        a = max(p.t, lo)
        b = min(q.t, hi)
        ya, yb = plf(a), plf(b)
        if min(ya, yb) <= value <= max(ya, yb):
            if yb == ya:  # 平台段：只可能在与运动段的衔接端点命中
                continue
            r = (value - ya) / (yb - ya)
            return a + (b - a) * r
    # 防御：理论上不可达
    raise DomainError(f"反查 {value} 失败")


# ----------------------------------------------------------------------
# 扫行重建：整数列采样 + 出处
# ----------------------------------------------------------------------
def sample_row(stroke: Stroke, run: MonotonicRun,
               trajectory: StageTrajectory, offset: Fraction,
               net: PLF | None, product_id: str | None,
               baseline_version: int | None,
               columns: Sequence[int]) -> tuple[ColumnSample, ...]:
    """对一个扫行的目标整数列逐一采样。

    判定顺序（严格不造数）：

    1. 该单调行程无台面连续运动（纯停顿）→ 全部 ``no_coverage``；
    2. 实物列坐标落在连续运动值域之外 → ``no_coverage``；
    3. 反查时刻落在净响应有效域之外（或未提供净响应）→ ``invalid_net``；
    4. 否则在净响应 PLF 上精确有理插值，状态 ``ok``。

    纯停顿期间出现的实物坐标天然被排除在连续运动值域之外，因而只会得到
    缺测而不会读到停顿期间仍在变化的信号——这正是「纯停顿不能造数」。
    """
    offset = to_fraction(offset, "offset")
    interval = run_motion_interval(stroke, trajectory, run)
    samples: list[ColumnSample] = []
    for col in columns:
        col = int(col)
        base = _missing(col, stroke, run, None)
        if interval is None:
            samples.append(_with_source(base, STATUS_NO_COVERAGE,
                                        product_id, baseline_version))
            continue
        t_lo, t_hi, _direction = interval
        # 实物坐标 x = s + offset；运动区间上的实物值域
        x_lo = trajectory(t_lo) + offset
        x_hi = trajectory(t_hi) + offset
        target = Fraction(col)
        if not (min(x_lo, x_hi) <= target <= max(x_lo, x_hi)):
            samples.append(_with_source(base, STATUS_NO_COVERAGE,
                                        product_id, baseline_version))
            continue
        t = invert_monotone(trajectory.stage_plf, target - offset,
                            t_lo, t_hi)
        base = ColumnSample(
            column=col, status=STATUS_OK, value=None,
            stroke_id=stroke.stroke_id, run_index=run.index, time=t,
            origin_version=0, product_id=None, baseline_version=None)
        if net is None or t < net.t0 or t > net.t1:
            samples.append(_with_source(base, STATUS_INVALID_NET,
                                        product_id, baseline_version))
            continue
        ok = ColumnSample(
            column=col, status=STATUS_OK, value=net(t),
            stroke_id=stroke.stroke_id, run_index=run.index, time=t,
            origin_version=0, product_id=product_id,
            baseline_version=baseline_version)
        samples.append(ok)
    return tuple(samples)


def _missing(column: int, stroke: Stroke, run: MonotonicRun,
             t: Fraction | None) -> ColumnSample:
    return ColumnSample(
        column=column, status=STATUS_NO_COVERAGE, value=None,
        stroke_id=stroke.stroke_id, run_index=run.index, time=t,
        origin_version=0, product_id=None, baseline_version=None)


def _with_source(sample: ColumnSample, status: str,
                 product_id: str | None,
                 baseline_version: int | None) -> ColumnSample:
    return ColumnSample(
        column=sample.column, status=status, value=None,
        stroke_id=sample.stroke_id, run_index=sample.run_index,
        time=sample.time, origin_version=0,
        product_id=product_id, baseline_version=baseline_version)


# ----------------------------------------------------------------------
# 重建计划：一次给出整组扫行（出处中的原点版本在此回填）
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class ReconstructedRow:
    """一个扫行的重建结果。"""

    row_id: str
    samples: tuple[ColumnSample, ...]

    def value_map(self) -> dict[int, Fraction]:
        """有效列 → 净响应读数（缺测列不在字典中）。"""
        return {s.column: s.value for s in self.samples
                if s.status == STATUS_OK and s.value is not None}


def reconstruct_rows(strokes: Sequence[Stroke],
                     trajectories: dict[str, StageTrajectory],
                     origin_versions: dict[tuple[str, int], OriginVersion],
                     nets: dict[str, PLF],
                     rows: Sequence[RowSpec],
                     columns_by_row: dict[str, Sequence[int]]
                     ) -> tuple[ReconstructedRow, ...]:
    """按扫行说明批量重建；每条行程的台面轨迹只算一次并复用。

    ``origin_versions`` 键为 ``(stroke_id, version)``；
    ``nets`` 键为 ``product_id``。
    """
    stroke_by_id = {s.stroke_id: s for s in strokes}
    out: list[ReconstructedRow] = []
    for spec in rows:
        stroke = stroke_by_id.get(spec.stroke_id)
        if stroke is None:
            raise ValidationError(f"未知行程: {spec.stroke_id}")
        ov = origin_versions.get((spec.stroke_id, spec.origin_version))
        if ov is None:
            raise ValidationError(
                f"行程 {spec.stroke_id} 原点版本 "
                f"{spec.origin_version} 不存在")
        all_runs = split_monotonic_runs(stroke)
        if not (0 <= spec.run_index < len(all_runs)):
            raise ValidationError(
                f"行程 {spec.stroke_id} 不存在单调行程 "
                f"#{spec.run_index}（共 {len(all_runs)} 次）")
        run = all_runs[spec.run_index]
        traj = trajectories[spec.stroke_id]
        net = nets.get(spec.product_id)
        raw = sample_row(stroke, run, traj, ov.offset, net,
                         spec.product_id, spec.baseline_version,
                         columns_by_row[spec.row_id])
        # 回填真实原点版本到出处
        filled = tuple(
            ColumnSample(
                column=s.column, status=s.status, value=s.value,
                stroke_id=s.stroke_id, run_index=s.run_index, time=s.time,
                origin_version=ov.version, product_id=s.product_id,
                baseline_version=s.baseline_version)
            for s in raw)
        out.append(ReconstructedRow(row_id=spec.row_id, samples=filled))
    return tuple(out)


def check_scan_limits(n_rows: int, n_overlap: int) -> None:
    """题面规模上限：至多 80 行、每行 80 个重叠列。"""
    if n_rows > limits.MAX_SCAN_ROWS:
        raise LimitError(
            f"扫行上限 {limits.MAX_SCAN_ROWS}，收到 {n_rows}")
    if n_overlap > limits.MAX_OVERLAP_COLUMNS:
        raise LimitError(
            f"每行重叠列上限 {limits.MAX_OVERLAP_COLUMNS}，"
            f"收到 {n_overlap}")


def check_motor_knots(n_knots: int) -> None:
    if n_knots > limits.MAX_MOTOR_KNOTS:
        raise LimitError(
            f"单条行程电机折点上限 {limits.MAX_MOTOR_KNOTS}，"
            f"收到 {n_knots}")
