"""组分定量引擎（纯函数，不碰数据库）。

归属规则
--------
每组分有响应因子 f_j（正有理数）和归属权重节点 ``(t,w)``（w≥0），
节点间线性插值、节点范围外权重恒为 0。同一时刻：

* W(t) = Σ_j n_j(t)；W>0 时组分 j 份额 = n_j(t)/W(t)；
* W=0 时该时刻全部净响应记为「未分配」。

事件网格 = 净响应折点 ∪ 各组分权重节点 ∪ 各组分支撑边界 ∩ 有效域。
每个格子内 u、n_j、W 全部线性，故：

* 总面积、未分配面积用精确有理梯形积分；
* 各组分面积用 :mod:`integration` 的闭式解析式；
* 守恒用精确值核对，绝不为凑数反扣某组分。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Sequence

from .errors import ValidationError
from .integration import component_interval_area
from .rational import to_fraction
from .segments import PLF, Point, net_response, plf_from_dict

ABS_TOL = Fraction(1, 10 ** 8)  # 题面统一绝对误差 1e-8


@dataclass(frozen=True)
class WeightFunction:
    """组分归属权重：节点间线性，支撑范围外为零。"""

    name: str
    factor: Fraction          # 固定正响应因子
    nodes: tuple[Point, ...]  # 时间严格递增，权重非负

    def __post_init__(self) -> None:
        if not self.name:
            raise ValidationError("组分名称不能为空")
        if self.factor <= 0:
            raise ValidationError(f"组分 {self.name} 响应因子必须为正有理数")
        if len(self.nodes) < 2:
            raise ValidationError(f"组分 {self.name} 至少需要两个权重节点")
        for i, p in enumerate(self.nodes):
            if i > 0 and p.t <= self.nodes[i - 1].t:
                raise ValidationError(f"组分 {self.name} 权重节点时间必须严格递增")
            if p.y < 0:
                raise ValidationError(f"组分 {self.name} 权重必须非负")

    @property
    def t0(self) -> Fraction:
        return self.nodes[0].t

    @property
    def t1(self) -> Fraction:
        return self.nodes[-1].t

    def weight(self, t: Fraction) -> Fraction:
        """范围外恒为 0（不外推）。"""
        return _eval_nodes(self.nodes, t)


# 节点序列不大（≤30），每次求值直接线性扫描，简单且零外部状态
def _eval_nodes(nodes: tuple[Point, ...], t: Fraction) -> Fraction:
    if t < nodes[0].t or t > nodes[-1].t:
        return Fraction(0)
    for p, q in zip(nodes, nodes[1:]):
        if p.t <= t <= q.t:
            if q.t == p.t:  # 结构上不可能，防御
                return p.y
            r = (t - p.t) / (q.t - p.t)
            return p.y + (q.y - p.y) * r
    return nodes[-1].y


def _trapezoid(plf: PLF, a: Fraction, b: Fraction) -> Fraction:
    """PLF 在 [a,b] 上的精确梯形积分（a、b 可落在子区间内部）。"""
    if a == b:
        return Fraction(0)
    return (plf(a) + plf(b)) * (b - a) / 2


@dataclass(frozen=True)
class MethodSpec:
    """一个定量方法：若干组分。组分顺序即身份，不做任何字典序重排。"""

    method_id: str
    components: tuple[WeightFunction, ...]

    def __post_init__(self) -> None:
        if not self.components:
            raise ValidationError("方法至少包含一个组分")
        names = [c.name for c in self.components]
        if len(set(names)) != len(names):
            raise ValidationError("组分名称必须唯一（名称即身份）")

    @property
    def names(self) -> list[str]:
        return [c.name for c in self.components]


@dataclass(frozen=True)
class SegmentProduct:
    """对外发布的分段产品：原采样、基线锚点快照 + 可复算净响应。"""

    product_id: str
    run_id: str
    baseline_version: int
    signal: PLF
    baseline: PLF
    net: PLF

    @property
    def domain(self) -> tuple[Fraction, Fraction]:
        return self.net.domain

    def fingerprint(self) -> str:
        """净响应函数内容指纹（去共线），等价产品同指纹。"""
        return self.net.fingerprint()


@dataclass(frozen=True)
class CellRecord:
    """解释用：一个事件格子上的归属明细。"""

    t_lo: Fraction
    t_hi: Fraction
    net_integral: Fraction          # 格子内净响应精确面积
    allocated_integral: Fraction    # W>0 部分（梯形精确值）
    unallocated_integral: Fraction  # W=0 部分（梯形精确值）
    component_areas: dict[str, float]
    W_lo: Fraction
    W_hi: Fraction


@dataclass(frozen=True)
class AnalysisResult:
    method_id: str
    product_id: str
    domain: tuple[Fraction, Fraction]
    component_area: dict[str, float]      # 净响应×份额 的积分（float，闭式）
    concentration: dict[str, float]       # area / 响应因子
    unallocated_area: float               # 未分配（精确有理转 float）
    total_area: float                     # 有效域内净响应总面积（精确有理）
    conservation_residual: float          # 守恒核对残差（解析和 vs 精确值）
    tolerance: float
    trace: tuple[CellRecord, ...] = field(repr=False, default_factory=tuple)


def build_product(product_id: str, run_id: str, baseline_version: int,
                  signal: PLF, baseline: PLF) -> SegmentProduct:
    """构造产品并自校验：重新计算净响应，确认与定义一致。"""
    net = net_response(signal, baseline)
    return SegmentProduct(product_id, run_id, baseline_version,
                          signal, baseline, net)


def product_from_dict(d: dict) -> SegmentProduct:
    """消费外部等价分段产品：只需分段本身，不依赖本库数据库。"""
    return SegmentProduct(
        product_id=d["product_id"],
        run_id=d["run_id"],
        baseline_version=int(d["baseline_version"]),
        signal=plf_from_dict(d["signal"]),
        baseline=plf_from_dict(d["baseline"]),
        net=plf_from_dict(d["net"]),
    )


def product_to_dict(p: SegmentProduct) -> dict:
    return {
        "product_id": p.product_id,
        "run_id": p.run_id,
        "baseline_version": p.baseline_version,
        "signal": p.signal.to_dict(),
        "baseline": p.baseline.to_dict(),
        "net": p.net.to_dict(),
    }


def _event_grid(product: SegmentProduct,
                components: Sequence[WeightFunction]) -> list[Fraction]:
    """所有真实折点的并集，裁剪到有效域；这是不依赖固定采样步长的关键。"""
    a, b = product.domain
    grid = {a, b}
    for k in product.net.knots:
        if a <= k.t <= b:
            grid.add(k.t)
    for c in components:
        for p in c.nodes:
            if a <= p.t <= b:
                grid.add(p.t)
        # 支撑边界必须进网格：边界两侧权重从 0 变正（或反之），
        # 即使节点本身在域外，域边界也已在网格中完成切段。
        if a <= c.t0 <= b:
            grid.add(c.t0)
        if a <= c.t1 <= b:
            grid.add(c.t1)
    return sorted(grid)


def analyze(product: SegmentProduct, method: MethodSpec,
            *, abs_tol: Fraction = ABS_TOL) -> AnalysisResult:
    """对一个分段产品执行定量分析。"""
    a, b = product.domain
    comps = method.components
    grid = _event_grid(product, comps)
    u = product.net

    comp_area: dict[str, float] = {c.name: 0.0 for c in comps}
    trace: list[CellRecord] = []
    total_exact = Fraction(0)
    unallocated_exact = Fraction(0)
    allocated_exact = Fraction(0)

    for lo, hi in zip(grid, grid[1:]):
        dt = hi - lo
        u0, u1 = u(lo), u(hi)
        cell_net = (u0 + u1) * dt / 2
        total_exact += cell_net

        # 网格已在每个组分的支撑边界切段。对活跃（整格位于支撑内）的组分
        # 取其线性端权重；仅在边界点接触的组分在格内恒 0（半开间断，
        # 边界点零测度）。
        n_lo: list[Fraction] = []
        n_hi: list[Fraction] = []
        for c in comps:
            if c.t0 <= lo and hi <= c.t1:
                n_lo.append(_eval_nodes(c.nodes, lo))
                n_hi.append(_eval_nodes(c.nodes, hi))
            else:
                n_lo.append(Fraction(0))
                n_hi.append(Fraction(0))
        W0, W1 = sum(n_lo), sum(n_hi)
        if W0 < 0 or W1 < 0:  # 防御：权重非负
            raise ValidationError("权重和出现负值")

        # 求 W 在格内的零点（W 线性，至多一个），切成 W 恒正/恒零子格。
        # 统一处理：支撑间断、内部零权重间隙、W 从零变正/归零。
        cuts = [lo, hi]
        if W0 > 0 and W1 > 0 and W0 != W1:
            root = lo + dt * (-W0) / (W1 - W0)
            if lo < root < hi:
                cuts = [lo, root, hi]

        cell_comp: dict[str, float] = {c.name: 0.0 for c in comps}
        cell_alloc = Fraction(0)
        cell_unalloc = Fraction(0)

        for s_lo, s_hi in zip(cuts, cuts[1:]):
            r_lo = (s_lo - lo) / dt
            r_hi = (s_hi - lo) / dt
            sW0 = W0 + (W1 - W0) * r_lo
            sW1 = W0 + (W1 - W0) * r_hi
            sub_net = _trapezoid(u, s_lo, s_hi)
            if sW0 == 0 and sW1 == 0:
                # 子格 W 恒零：未分配
                cell_unalloc += sub_net
                continue
            # 子格 W 恒正（端点零测度不影响）：组分权重线性，闭式积分
            cell_alloc += sub_net
            sdt = s_hi - s_lo
            su0, su1 = u(s_lo), u(s_hi)
            for j, c in enumerate(comps):
                v0 = n_lo[j] + (n_hi[j] - n_lo[j]) * r_lo
                v1 = n_lo[j] + (n_hi[j] - n_lo[j]) * r_hi
                if v0 == 0 and v1 == 0:
                    continue
                area = component_interval_area(
                    su0, su1, v0, v1, sW0, sW1, sdt)
                if area < -1e-12:
                    raise ValidationError(
                        f"组分 {c.name} 出现负面积 {area}")
                cell_comp[c.name] += area
                comp_area[c.name] += area

        allocated_exact += cell_alloc
        unallocated_exact += cell_unalloc
        trace.append(CellRecord(
            t_lo=lo, t_hi=hi, net_integral=cell_net,
            allocated_integral=cell_alloc,
            unallocated_integral=cell_unalloc,
            component_areas=cell_comp, W_lo=W0, W_hi=W1,
        ))

    # ---- 守恒核对（精确有理真值，禁止反扣组分）----
    total_f = float(total_exact)
    sum_comp = sum(comp_area.values())
    residual = abs(sum_comp + float(unallocated_exact) - total_f)
    if residual > float(abs_tol):
        raise ValidationError(
            f"守恒核对失败: 残差 {residual:.3e} 超过容差 {float(abs_tol):.0e}；"
            "不得通过反扣组分凑守恒，请检查积分实现")
    if abs(allocated_exact + unallocated_exact - total_exact) > Fraction(1, 10 ** 18):
        raise ValidationError("精确分配/未分配之和不等于总面积")

    concentration = {
        c.name: comp_area[c.name] / float(c.factor) for c in comps
    }
    return AnalysisResult(
        method_id=method.method_id,
        product_id=product.product_id,
        domain=(a, b),
        component_area=dict(comp_area),
        concentration=concentration,
        unallocated_area=float(unallocated_exact),
        total_area=total_f,
        conservation_residual=residual,
        tolerance=float(abs_tol),
        trace=tuple(trace),
    )


def make_weight_function(name: str, factor,
                         nodes: Sequence[tuple]) -> WeightFunction:
    pts = tuple(Point(to_fraction(t, "t"), to_fraction(w, "weight"))
                for t, w in nodes)
    return WeightFunction(name=name, factor=to_fraction(factor, "factor"),
                          nodes=pts)
