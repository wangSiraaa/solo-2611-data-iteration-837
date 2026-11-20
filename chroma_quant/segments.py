"""分段线性函数（Piecewise Linear Function, PLF）。

系统内一切连续曲线——原信号、基线、净响应、归属权重——都表示为
``(Fraction x, Fraction y)`` 折点序列上的分段线性函数，折点严格递增。

关键约定
--------
* 折点即真实切段点：净响应 ``max(y-b,0)`` 的交点会被真实地插入切段，
  而不是在采样点处截零后再连线。
* 净响应在整个有效域内显式保留零段（值为 0 的折点不丢弃），
  因此任何区间查询与积分都不需要外推或猜值。
"""
from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

from .errors import DomainError, ValidationError
from .rational import frac_str, stable_id, to_fraction


@dataclass(frozen=True)
class Point:
    t: Fraction
    y: Fraction


@dataclass(frozen=True)
class PLF:
    """不可变分段线性函数，``knots`` 至少两个且时间严格递增。"""

    knots: tuple[Point, ...]

    def __post_init__(self) -> None:
        if len(self.knots) < 2:
            raise ValidationError("PLF 至少需要两个折点")
        for i, k in enumerate(self.knots):
            if not isinstance(k.t, Fraction) or not isinstance(k.y, Fraction):
                raise ValidationError("PLF 折点必须是 Fraction")
            if i > 0 and k.t <= self.knots[i - 1].t:
                raise ValidationError("PLF 折点时间必须严格递增")

    # ---- 基本性质 -------------------------------------------------
    @property
    def t0(self) -> Fraction:
        return self.knots[0].t

    @property
    def t1(self) -> Fraction:
        return self.knots[-1].t

    @property
    def domain(self) -> tuple[Fraction, Fraction]:
        return self.t0, self.t1

    def in_domain(self, t: Fraction, *, closed: bool = True) -> bool:
        if closed:
            return self.t0 <= t <= self.t1
        return self.t0 < t < self.t1

    # ---- 求值 -----------------------------------------------------
    def __call__(self, t: Fraction) -> Fraction:
        """闭域内精确有理求值；域外抛 DomainError，绝不外推。"""
        t = to_fraction(t, "t")
        if t < self.t0 or t > self.t1:
            raise DomainError(f"t={t} 落在有效域 [{self.t0}, {self.t1}] 之外")
        xs = [k.t for k in self.knots]
        i = bisect_left(xs, t)
        if i < len(xs) and xs[i] == t:
            return self.knots[i].y
        lo, hi = self.knots[i - 1], self.knots[i]
        r = (t - lo.t) / (hi.t - lo.t)
        return lo.y + (hi.y - lo.y) * r

    # ---- 变换 -----------------------------------------------------
    def clip(self, a: Fraction, b: Fraction) -> "PLF":
        """裁到闭区间 [a,b]（含端点插值）；域外拒绝。"""
        a, b = to_fraction(a, "a"), to_fraction(b, "b")
        if a >= b:
            raise ValidationError("clip 要求 a < b")
        if a < self.t0 or b > self.t1:
            raise DomainError("clip 区间超出有效域，拒绝外推")
        out: list[Point] = [Point(a, self(a))]
        for k in self.knots:
            if k.t > a and k.t < b:
                out.append(k)
        out.append(Point(b, self(b)))
        return PLF(tuple(out))

    def canonical(self) -> "PLF":
        """去掉斜率不变的冗余共线折点（端点保留）。

        同一条直线无论用多少中间点表示，规范化后完全一致，
        这是“不同输入分段但同函数”等价判定的基础。
        """
        kept: list[Point] = [self.knots[0]]
        kept.extend(self.knots[1:-1])
        kept.append(self.knots[-1])
        # 三点共线则中间点冗余：
        # (q.y-p.y)/(q.t-p.t) == (r.y-q.y)/(r.t-q.t)，交叉相乘避免除法
        result: list[Point] = [kept[0]]
        for i in range(1, len(kept) - 1):
            p, q, r = kept[i - 1], kept[i], kept[i + 1]
            if (q.y - p.y) * (r.t - q.t) != (r.y - q.y) * (q.t - p.t):
                result.append(q)
        result.append(kept[-1])
        return PLF(tuple(result))

    def integral(self, a: Fraction | None = None, b: Fraction | None = None) -> Fraction:
        """精确有理定积分（梯形），默认整条函数。"""
        a = self.t0 if a is None else to_fraction(a, "a")
        b = self.t1 if b is None else to_fraction(b, "b")
        if a < self.t0 or b > self.t1:
            raise DomainError("积分区间超出有效域")
        if a > b:
            raise ValidationError("积分要求 a <= b")
        if a == b:
            return Fraction(0)
        piece = self.clip(a, b)
        total = Fraction(0)
        for p, q in zip(piece.knots, piece.knots[1:]):
            total += (p.y + q.y) * (q.t - p.t) / 2
        return total

    def fingerprint(self) -> str:
        """规范化（去共线）后的内容哈希；同函数必同哈希。"""
        c = self.canonical()
        return stable_id(
            "PLF",
            [(frac_str(k.t), frac_str(k.y)) for k in c.knots],
        )

    def to_dict(self) -> dict:
        return {
            "knots": [
                {"t": frac_str(k.t), "y": frac_str(k.y)} for k in self.knots
            ]
        }


def make_plf(points: Iterable[tuple[Any, Any] | Point]) -> PLF:
    """由 (t, y) 序列构造 PLF，要求时间严格递增、至少两点。"""
    pts: list[Point] = []
    for item in points:
        if isinstance(item, Point):
            pts.append(item)
        else:
            t, y = item
            pts.append(Point(to_fraction(t, "t"), to_fraction(y, "y")))
    return PLF(tuple(pts))


def plf_from_dict(d: dict) -> PLF:
    return make_plf((p["t"], p["y"]) for p in d["knots"])


# ----------------------------------------------------------------------
# 净响应：真实切段
# ----------------------------------------------------------------------
def net_response(signal: PLF, baseline: PLF) -> PLF:
    """返回 ``max(signal - baseline, 0)``，交点处真实切段。

    有效域取两者域的交集（基线锚点覆盖的闭范围即有效域）。
    交点用一次线性求根精确（有理）定位，净函数在整个域上显式保留零段。
    """
    a = max(signal.t0, baseline.t0)
    b = min(signal.t1, baseline.t1)
    if a >= b:
        raise ValidationError("信号与基线没有正向长度的公共有效域")

    # 合并两者所有折点（裁剪到公共域），保证每个子区间上信号、基线均线性
    grid = {a, b}
    for k in signal.knots:
        if a < k.t < b:
            grid.add(k.t)
    for k in baseline.knots:
        if a < k.t < b:
            grid.add(k.t)
    xs = sorted(grid)

    def d(t: Fraction) -> Fraction:
        return signal(t) - baseline(t)

    knots: list[Point] = [Point(a, max(d(a), Fraction(0)))]
    for t_prev, x in zip(xs, xs[1:]):
        d0, d1 = d(t_prev), d(x)
        # 在 (t_prev, x) 内找 d 的唯一根（线性至多一个）
        if d0 == 0 or (d0 > 0) != (d1 > 0):
            if d0 != d1:  # 恒零段无根可插
                root = t_prev + (x - t_prev) * (-d0) / (d1 - d0)
                if t_prev < root < x:
                    knots.append(Point(root, Fraction(0)))
        knots.append(Point(x, max(d1, Fraction(0))))
    return PLF(tuple(knots))
