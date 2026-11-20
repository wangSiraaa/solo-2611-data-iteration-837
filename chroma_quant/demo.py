"""可复跑端到端演示：合成信号 → 基线产品 → 重叠峰定量 → 基线更正 → 复读。

运行：``python -m chroma_quant.demo``
演示库写入临时目录，结束后打印路径；所有输入均为本地合成有理数，
不涉及任何真实仪器或化学物质识别。
"""
from __future__ import annotations

import json
import os
import tempfile
from fractions import Fraction as F

from .quant import product_to_dict
from .repository import QuantRepository
from .signal_service import SignalService


def main() -> None:
    workdir = tempfile.mkdtemp(prefix="chroma_demo_")
    db = os.path.join(workdir, "demo.sqlite")
    print(f"演示数据库: {db}\n")

    svc = SignalService(db)
    # 一条含两个重叠峰的合成信号（整数有理采样，线性插值）
    samples = [(t, _signal_y(t)) for t in range(0, 9)]
    svc.create_run("DEMO", samples, note="本地合成重叠峰信号")

    # 基线 v1：略微上倾
    v1 = svc.add_baseline("DEMO", [(0, 1), (8, 3)])
    p1 = svc.publish_product("DEMO", v1)
    print(f"基线 v{v1} 产品: {p1.product_id[:16]}…")
    print(f"  有效域: {p1.net.domain}")
    print(f"  净响应总面积(精确): {p1.net.integral()} "
          f"≈ {float(p1.net.integral()):.6f}\n")

    repo = QuantRepository(db)
    # 两个重叠组分 + 一个时间很窄的归属窗口，名称顺序即身份（非字典序）
    repo.create_method("OVL", [
        {"name": "zeta", "factor": F(3, 2),
         "nodes": [(0, 2), (3, 2), (5, 0), (8, 0)]},
        {"name": "alpha", "factor": F(7, 3),
         "nodes": [(0, 0), (2, 0), (4, 3), (6, 3), (8, 1)]},
        {"name": "tiny", "factor": 1,
         "nodes": [(F(7, 3), 1), (F(8, 3), 1)]},  # 宽度 1/3 的窄窗
    ])
    r1 = repo.run_analysis("OVL", p1.product_id)
    _print_result("基线 v1 定量结果", r1)

    # 独立交付：信号服务自己的区间净面积（不经过组分方法）
    area = svc.interval_area("DEMO", F(1), F(7), v1)
    print(f"信号服务独立区间面积 [1,7] = {float(area):.6f}\n")

    # 基线更正：锚点整体下移，新版本；旧产品/结果保留
    v2 = svc.correct_baseline("DEMO", [(0, 0), (8, 1)])
    p2 = svc.publish_product("DEMO", v2)
    r2 = repo.run_analysis("OVL", p2.product_id)
    _print_result(f"基线 v{v2} 更正后定量结果", r2)

    # 复读旧结果
    old = repo.find_analysis("OVL", p1.product_id)
    print("旧结果复读: zeta 面积 = "
          f"{old.component_area['zeta']:.6f}（仍为更正前值）")

    # 外部等价分段产品：导出 JSON 到新库，按净响应函数指纹去重复用
    ext_db = os.path.join(workdir, "external.sqlite")
    ext_repo = QuantRepository(ext_db)
    pid = ext_repo.import_product(product_to_dict(p2))
    ext_repo.create_method("OVL", [
        {"name": "zeta", "factor": F(3, 2),
         "nodes": [(0, 2), (3, 2), (5, 0), (8, 0)]},
        {"name": "alpha", "factor": F(7, 3),
         "nodes": [(0, 0), (2, 0), (4, 3), (6, 3), (8, 1)]},
        {"name": "tiny", "factor": 1,
         "nodes": [(F(7, 3), 1), (F(8, 3), 1)]},
    ])
    r_ext = ext_repo.run_analysis("OVL", pid)
    print("\n外部库分析与原库一致: "
          f"{abs(r_ext.total_area - r2.total_area) < 1e-12}")
    print(f"演示文件保留在: {workdir}")
    svc.close()
    repo.close()
    ext_repo.close()


def _signal_y(t: int) -> int:
    # 两个三角峰叠加的本地合成整数序列
    peak1 = max(0, 4 - abs(t - 3))
    peak2 = max(0, 3 - abs(t - 6))
    return 2 + peak1 + peak2


def _print_result(title: str, res) -> None:
    print(f"== {title} ==")
    print(json.dumps({
        "component_area": {k: round(v, 6)
                           for k, v in res.component_area.items()},
        "concentration": {k: round(v, 6)
                          for k, v in res.concentration.items()},
        "unallocated_area": round(res.unallocated_area, 6),
        "total_area": round(res.total_area, 6),
        "conservation_residual": f"{res.conservation_residual:.2e}",
    }, ensure_ascii=False, indent=2))
    assert abs(sum(res.component_area.values())
               + res.unallocated_area - res.total_area) <= res.tolerance
    print()


if __name__ == "__main__":
    main()
