"""积分策略。

主路径
------
份额 ``weight/W`` 相除后一般不是多项式（W 线性、可为非常数），
在事件网格的每个格子内：

* 净响应 u、组分权重 n_j、权重和 W 都是线性的；
* 被积函数 ``u*n_j/W`` 是「二次多项式 / 一次多项式」，存在闭式原函数
  （含一个对数项）。:func:`component_interval_area` 返回该闭式值（float64）。

独立核对路径（不共享任何推导代码）
----------------------------------
* :func:`adaptive_simpson_float`：经典递归自适应 Simpson，float 运算；
* :func:`composite_simpson_fraction`：纯 ``Fraction`` 复合 Simpson，
  逐级加密直到两次精确结果之差小于有理容差，全程无浮点。

题面允许统一绝对误差 1e-8；验证脚本要求三条路径两两一致。
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Callable

from .errors import ChromaError


def component_interval_area(
    u0: Fraction,
    u1: Fraction,
    n0: Fraction,
    n1: Fraction,
    W0: Fraction,
    W1: Fraction,
    dt: Fraction,
) -> float:
    """一个事件格子内 ``∫ u(t)·n(t)/W(t) dt`` 的闭式值。

    归一坐标 s∈[0,1]（物理长度 dt）：
      u(s) = u0 + du·s，  n(s) = n0 + dn·s，  W(s) = W0 + B·s
    分子 u·n 是 **二次** 多项式（两个线性之积），分母 W 线性。
    W 在格子内部严格为正（调用方保证），端点允许恰好为 0（零测度）。
    """
    if dt <= 0:
        raise ChromaError("事件格子长度必须为正")
    if W0 < 0 or W1 < 0 or (W0 == 0 and W1 <= 0) or (W1 == 0 and W0 <= 0):
        raise ChromaError(
            f"事件格子 W 结构非法: W0={W0}, W1={W1}")

    du, dn, B = u1 - u0, n1 - n0, W1 - W0
    # u·n = c0 + c1·s + c2·s²
    c0, c1, c2 = u0 * n0, u0 * dn + du * n0, du * dn

    if B == 0:
        # W 恒为正常数：二次多项式精确积分
        # ∫_0^1 (c0+c1 s+c2 s²) ds = c0 + c1/2 + c2/3
        return float(dt * (c0 + c1 / 2 + c2 / 3) / W0)

    # W(s)=W0+B·s。做多项式除法，把二次分子除以一次分母：
    #   (c0+c1 s+c2 s²)/W = α·s + β + γ/W，解出 α,β,γ（均为精确有理数）
    alpha = c2 / B
    beta = (c1 - alpha * W0) / B
    gamma = c0 - beta * W0
    # 积分（s∈[0,1]，再乘物理长度 dt）：
    #   ∫ α s ds = α/2，∫ β ds = β，∫ γ/W ds = γ/B · ln(W1/W0)
    if W0 == 0 or W1 == 0:
        # 端点 W=0：γ/B·W·ln(W) 项按极限 → 0，对数贡献消失（零测度）
        log_term = 0.0
    else:
        log_term = float(gamma / B) * math.log(float(W1) / float(W0))
    return float(dt) * (float(alpha / 2 + beta) + log_term)


# ----------------------------------------------------------------------
# 独立核对路径 1：float 递归自适应 Simpson
# ----------------------------------------------------------------------
def _simpson(f: Callable[[float], float], a: float, b: float,
             fa: float, fm: float, fb: float) -> float:
    return (b - a) * (fa + 4.0 * fm + fb) / 6.0


def adaptive_simpson_float(f: Callable[[float], float], a: float, b: float,
                           tol: float = 1e-12, max_depth: int = 40) -> float:
    fa, fb = f(a), f(b)
    m = (a + b) / 2.0
    fm = f(m)
    whole = _simpson(f, a, b, fa, fm, fb)

    def rec(lo: float, hi: float, flo: float, fmid: float,
            fhi: float, whole_v: float, eps: float, depth: int) -> float:
        mid = (lo + hi) / 2.0
        lm = (lo + mid) / 2.0
        rm = (mid + hi) / 2.0
        flm, frm = f(lm), f(rm)
        left = _simpson(f, lo, mid, flo, flm, fmid)
        right = _simpson(f, mid, hi, fmid, frm, fhi)
        delta = left + right - whole_v
        if depth <= 0 or abs(delta) <= 15.0 * eps:
            return left + right + delta / 15.0
        return (rec(lo, mid, flo, flm, fmid, left, eps / 2, depth - 1)
                + rec(mid, hi, fmid, frm, fhi, right, eps / 2, depth - 1))

    return rec(a, b, fa, fm, fb, whole, tol, max_depth)


# ----------------------------------------------------------------------
# 独立核对路径 2：纯 Fraction 复合 Simpson（无任何浮点）
# ----------------------------------------------------------------------
def composite_simpson_fraction(
    f: Callable[[Fraction], Fraction],
    a: Fraction,
    b: Fraction,
    tol: Fraction = Fraction(1, 10 ** 10),
    max_n: int = 8192,
) -> Fraction:
    """偶数等分复合 Simpson，逐级加倍直到相邻两次精确值之差 < tol。"""
    prev: Fraction | None = None
    n = 4
    while n <= max_n:
        h = (b - a) / n
        s = f(a) + f(b)
        for i in range(1, n):
            s += (4 if i % 2 else 2) * f(a + i * h)
        cur = s * h / 3
        if prev is not None and abs(cur - prev) < tol:
            return cur
        prev = cur
        n *= 2
    if prev is None:  # pragma: no cover
        raise ChromaError("Fraction Simpson 未能启动")
    return prev
