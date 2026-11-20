"""有理数工具：所有时间、响应、权重在系统内部均为 ``Fraction``。

只接受本地合成的有理输入；float 也会按其二进制精确值转成 Fraction
（建议调用方使用 int 或 "p/q" 字符串以避免十进制歧义）。
"""
from __future__ import annotations

import hashlib
from fractions import Fraction
from typing import Any

from .errors import ValidationError


def to_fraction(x: Any, name: str = "value") -> Fraction:
    """把 int/str/Fraction/float 严格转为 Fraction。"""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, bool):  # bool 是 int 的子类，语义上拒绝
        raise ValidationError(f"{name} 不接受布尔值")
    if isinstance(x, int):
        return Fraction(x)
    if isinstance(x, float):
        return Fraction(x)
    if isinstance(x, str):
        s = x.strip()
        try:
            return Fraction(s)
        except (ZeroDivisionError, ValueError) as exc:
            raise ValidationError(f"{name} 不是合法有理数: {x!r}") from exc
    raise ValidationError(f"{name} 类型无法转为有理数: {type(x).__name__}")


def frac_str(x: Fraction) -> str:
    """稳定的有理数字符串形式，用于序列化与哈希。"""
    return f"{x.numerator}/{x.denominator}"


def parse_fraction(s: str) -> Fraction:
    return to_fraction(s)


def stable_id(*parts: Any) -> str:
    """由若干可读部件生成确定性 sha256 标识（缓存与复读的基础）。"""
    h = hashlib.sha256()
    for p in parts:
        h.update(repr(p).encode("utf-8"))
        h.update(b"\x1f")
    return h.hexdigest()
