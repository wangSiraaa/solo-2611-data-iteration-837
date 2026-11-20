"""色谱基线产品与重叠峰组分定量仓（本地合成有理时间序列，无真实仪器）。"""
from __future__ import annotations

from .errors import (ChromaError, DomainError, LimitError, NotFoundError,
                     ValidationError)
from .integration import (adaptive_simpson_float,
                          composite_simpson_fraction,
                          component_interval_area)
from .motion import (MISSING_STATUSES, STATUS_DWELL, STATUS_INVALID_NET,
                     STATUS_NO_COVERAGE, STATUS_OK, ColumnSample,
                     MonotonicRun, OriginVersion, ReconstructedRow, RowSpec,
                     StageTrajectory, Stroke, build_stage_trajectory,
                     invert_monotone, reconstruct_rows, run_motion_interval,
                     sample_row, split_monotonic_runs)
from .quant import (ABS_TOL, AnalysisResult, MethodSpec, SegmentProduct,
                    WeightFunction, analyze, build_product,
                    make_weight_function, product_from_dict, product_to_dict)
from .repository import QuantRepository
from .scan_service import ScanService
from .segments import PLF, Point, make_plf, net_response
from .signal_service import SignalService
from .stitch import (CELL_LEFT, CELL_MISSING, CELL_RIGHT, CompositeCell,
                     LocatedStrip, PrintSheet, RowSeam, SeamResult,
                     find_seams)

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
    # 往返扫描：回差坐标重建
    "Stroke", "MonotonicRun", "StageTrajectory", "OriginVersion",
    "RowSpec", "ColumnSample", "ReconstructedRow",
    "build_stage_trajectory", "split_monotonic_runs",
    "run_motion_interval", "invert_monotone", "sample_row",
    "reconstruct_rows",
    "STATUS_OK", "STATUS_DWELL", "STATUS_NO_COVERAGE",
    "STATUS_INVALID_NET", "MISSING_STATUSES",
    # 双图带拼接缝
    "LocatedStrip", "PrintSheet", "CompositeCell", "RowSeam",
    "SeamResult", "find_seams",
    "CELL_LEFT", "CELL_RIGHT", "CELL_MISSING",
    "ScanService",
]
