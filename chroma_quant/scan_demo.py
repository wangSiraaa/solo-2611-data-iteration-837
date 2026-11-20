"""往返扫描 / 双图带拼接的可复跑端到端演示。

运行：``python -m chroma_quant.scan_demo``
演示库写入临时目录。全部为本地合成有理数行程与净响应，不涉及真实仪器。

演示主线
--------
1. 建立带换向停顿的真实电机行程，展示回差带动与「停顿期间信号仍变化、
   但停顿不造数」；
2. 正扫、反扫分别按实物列反查时刻、从净响应产品插值，附完整出处；
3. 用真实重建的两条图带重放题面三个接缝结果；
4. 原点订正只追加偏移版本，电机运动不重算；换印样单不触碰运动。
"""
from __future__ import annotations

import os
import tempfile
from fractions import Fraction as F

from .scan_service import ScanService
from .signal_service import SignalService


def _publish(svc: SignalService, run_id: str, signal: list,
             baseline: list | None = None) -> str:
    svc.create_run(run_id, [(t, signal[t]) for t in range(len(signal))])
    anchors = baseline if baseline is not None else [
        (0, 0), (len(signal) - 1, 0)]
    version = svc.add_baseline(run_id, anchors)
    return svc.publish_product(run_id, version).product_id


def main() -> None:
    workdir = tempfile.mkdtemp(prefix="chroma_scan_")
    db = os.path.join(workdir, "scan.sqlite")
    print(f"演示数据库: {db}\n")
    svc = SignalService(db)
    scan = ScanService(db)

    # ------------------------------------------------------------------
    # 1. 换向停顿行程：电机 0→3（t=3）端点停顿到 t=5 →0（t=9），g=1
    # ------------------------------------------------------------------
    scan.create_stroke(
        "ST-DWELL",
        [(0, 0), (3, 3), (5, 3), (9, 0)],
        g=F(1), initial_stage=F(0),
        note="端点换向停顿；停顿期间净响应仍变化")
    traj = scan.get_stage_trajectory("ST-DWELL")
    print("台面轨迹折点 (t, s):")
    print("  " + ", ".join(f"({k.t},{k.y})"
                            for k in traj.stage_plf.knots))
    runs = scan.list_monotonic_runs("ST-DWELL")
    print("单调行程: " + ", ".join(f"#{r.index}(方向{r.direction:+d})"
                                    for r in runs))

    # 净响应在停顿段 t∈[3,5] 从 30 变到 60 —— 信号在停顿期间确实变化
    pid_dwell = _publish(svc, "SIG-DWELL", [0, 10, 20, 30, 60, 70, 80, 90,
                                            90, 90])
    scan.add_scan_row("dwell-fwd", "ST-DWELL", 0, 1, pid_dwell, 1)
    scan.add_scan_row("dwell-rev", "ST-DWELL", 1, 1, pid_dwell, 1)
    rebuilt = {r.row_id: r for r in scan.reconstruct({
        "dwell-fwd": [0, 1, 2, 3],
        "dwell-rev": [0, 1, 2, 3],
    })}
    print("\n换向停顿例（停顿段 t=3..5，净值 30→60）：")
    for rid in ("dwell-fwd", "dwell-rev"):
        for s in rebuilt[rid].samples:
            t = f"t={s.time}" if s.time is not None else "—"
            v = str(s.value) if s.value is not None else f"<{s.status}>"
            print(f"  {rid:10s} 列{s.column}: {v:>6}  ({t},"
                  f" 行程{s.stroke_id}#{s.run_index},"
                  f" 原点v{s.origin_version})")
    rev2 = next(s for s in rebuilt["dwell-rev"].samples if s.column == 2)
    print(f"  => 反扫列2取 t={rev2.time}（在停顿末 t=5 之后），"
          f"绝不读停顿值 30..60\n")

    # ------------------------------------------------------------------
    # 2. 真实图带：左图两行恒零；右图行1=(0,8,1)、行2=(8,8,0,0)
    # ------------------------------------------------------------------
    pid_l1 = _publish(svc, "SIG-L1", [0] * 7)
    pid_l2 = _publish(svc, "SIG-L2", [0] * 7)
    pid_r1 = _publish(svc, "SIG-R1", [0, 0, 0, 8, 1, 0, 0])
    pid_r2 = _publish(svc, "SIG-R2", [0, 0, 8, 8, 0, 0, 0])
    for sid in ("K-L1", "K-L2", "K-R1", "K-R2"):
        scan.create_stroke(sid, [(t, t) for t in range(7)],
                           g=F(1), initial_stage=F(0))
    scan.add_scan_row("row-L1", "K-L1", 0, 1, pid_l1, 1)
    scan.add_scan_row("row-L2", "K-L2", 0, 1, pid_l2, 1)
    scan.add_scan_row("row-R1", "K-R1", 0, 1, pid_r1, 1)
    scan.add_scan_row("row-R2", "K-R2", 0, 1, pid_r2, 1)
    plan = {rid: [1, 2, 3]
            for rid in ("row-L1", "row-L2", "row-R1", "row-R2")}
    rb = {r.row_id: r.value_map() for r in scan.reconstruct(plan)}
    scan.import_located_strip("STRIP-L", (rb["row-L1"], rb["row-L2"]))
    scan.import_located_strip("STRIP-R", (rb["row-R1"], rb["row-R2"]))
    scan.create_print_sheet("SHEET-A", 2, 1, 3)
    scan.create_print_sheet("SHEET-B", 2, 1, 3, disabled=((0, 3),))

    def show(title: str, sheet: str, right: str) -> None:
        res = scan.stitch(sheet, "STRIP-L", right,
                          full_columns=(0, 1, 2, 3, 4))
        if not res.feasible:
            print(f"{title}: 首个无法接续的行 = "
                  f"{res.first_blocked_row}")
            return
        seams = [s.seam_column for s in res.seams]
        print(f"{title}: 缝={seams} 总代价={res.total_cost} "
              f"各行代价={[str(s.cost) for s in res.seams]}")

    print("题面三个接缝重放（真实行程+净响应产品重建）：")
    show("  重放1 全通行        ", "SHEET-A", "STRIP-R")
    show("  重放2 禁首行第3列    ", "SHEET-B", "STRIP-R")

    # 原点订正：K-R2 向左移一列（offset=-1，新版本 v2，运动不重算）
    fp_before = scan.get_stage_trajectory("K-R2").stage_plf.fingerprint()
    v2 = scan.correct_origin("K-R2", F(-1))
    scan.add_scan_row("row-R2b", "K-R2", 0, v2, pid_r2, 1)
    rb2 = {r.row_id: r.value_map() for r in
           scan.reconstruct({"row-R2b": [1, 2, 3, 4]})}
    fp_after = scan.get_stage_trajectory("K-R2").stage_plf.fingerprint()
    scan.import_located_strip(
        "STRIP-R2", (rb["row-R1"], rb2["row-R2b"]))
    show("  重放3 原点左移一列后 ", "SHEET-A", "STRIP-R2")
    print(f"\n原点订正后台面轨迹指纹不变（电机运动未重算）: "
          f"{fp_before == fp_after}")
    print(f"演示文件保留在: {workdir}")
    svc.close()
    scan.close()


if __name__ == "__main__":
    main()
