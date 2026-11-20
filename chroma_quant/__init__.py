"""色谱基线产品与重叠峰组分定量仓（本地合成有理时间序列，无真实仪器）。"""
from __future__ import annotations

from .errors import (ChromaError, DomainError, LimitError, NotFoundError,
                     ValidationError)
from .integration import (adaptive_simpson_float,
                          composite_simpson_fraction,
                          component_interval_area)
from .quant import (ABS_TOL, AnalysisResult, MethodSpec, SegmentProduct,
                    WeightFunction, analyze, build_product,
                    make_weight_function, product_from_dict, product_to_dict)
from .repository import QuantRepository
from .segments import PLF, Point, make_plf, net_response
from .signal_service import SignalService

__all__ = [
    "ChromaError", "ValidationError", "LimitError", "DomainError",
    "NotFoundError",
    "Point", "PLF", "make_plf", "net_response",
    "WeightFunction", "MethodSpec", "SegmentProduct", "AnalysisResult",
    "make_weight_function", "analyze", "build_product",
    "product_from_dict", "product_to_dict", "ABS_TOL",
    "component_interval_area", "adaptive_simpson_float",
    "composite_simpson_fraction",
    "SignalService", "QuantRepository",
]
