"""往返扫描 / 双图带拼接的验证入口。

运行：``python -m chroma_quant.scan_validation``
与旧 :mod:`chroma_quant.validation` 完全独立，旧 54 项断言不受影响。

覆盖：
* S1 回差力学（含题面 g=1 折返原话例、段内真实交点切段）；
* S2 换向停顿期间信号变化但停顿不造数、净响应无效不造数；
* S3 正扫/反扫出处（行程、时刻、源版本）与原点偏移；
* S4-S6 题面三个接缝重放（全通行 / 禁用 / 原点订正），全部用真实
  合成行程与净响应产品经完整重建管线产出图带；
* S7 不可接续行定位、缺测保持缺测；
* S8 持久化、原点订正不重算电机运动、换印样单不触碰运动。
"""
from __future__ import annotations

import os
import tempfile
from fractions import Fraction as F

from chroma_quant import (LocatedStrip, PrintSheet, Stroke,
                          build_stage_trajectory, sample_row,
                          split_monotonic_runs)
from .errors import ValidationError
from .motion import (STATUS_INVALID_NET, STATUS_NO_COVERAGE, STATUS_OK)
from .scan_service import ScanService
from .segments import make_plf
from .signal_service import SignalService
from .stitch import CELL_LEFT, CELL_MISSING, CELL_RIGHT, find_seams

results: list[tuple[str, bool, str]] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    results.append((name, bool(cond), detail))
    flag = "PASS" if cond else "FAIL"
    print(f"[{flag}] {name}" + (f"  {detail}" if detail else ""))
    if not cond:
        raise AssertionError(name + " " + detail)


# ----------------------------------------------------------------------
# S1：回差力学（题面原话例 + 段内交点）
# ----------------------------------------------------------------------
def case_backlash_mechanics() -> None:
    # 题面：g=1、电机从 3 折返到 0、折返时台面在 2；电机先退到 1，
    # 台面才会继续往回运动。
    st_rev = make_plf([(0, 3), (3, 0)])
    stroke = Stroke("rev", st_rev, backlash=F(1), initial_stage=F(2))
    traj = build_stage_trajectory(stroke)
    # 左边 u+1 在 t=2（u=1）才碰到台面 2：(0,2)-(2,2) 平台，(2,2)-(3,1) 贴边
    check("S1 折返原话例折点",
          [(k.t, k.y) for k in traj.stage_plf.knots]
          == [(F(0), F(2)), (F(2), F(2)), (F(3), F(1))],
          str([(str(k.t), str(k.y)) for k in traj.stage_plf.knots]))
    check("S1 电机退到1台面仍不动", traj(F(2)) == F(2))
    check("S1 电机到0台面回到1", traj(F(3)) == F(1))

    # 完整往返：0→3（停顿）→0，g=1，台面起 0
    full = Stroke(
        "full", make_plf([(0, 0), (3, 3), (6, 3), (9, 0)]),
        backlash=F(1), initial_stage=F(0))
    t_full = build_stage_trajectory(full)
    expect = [(F(0), 0), (F(1), 0), (F(3), 2), (F(6), 2),
              (F(8), 2), (F(9), 1)]
    check("S1 往返台面轨迹(段内交点)",
          [(k.t, k.y) for k in t_full.stage_plf.knots] == expect)
    runs = split_monotonic_runs(full)
    check("S1 平台/换向拆成两次单调行程",
          [(r.index, r.direction) for r in runs] == [(0, 1), (1, -1)])

    # g=0：台面严格跟随电机（无回差）
    tight = Stroke(
        "tight", make_plf([(0, 0), (2, 5)]),
        backlash=F(0), initial_stage=F(0))
    check("S1 g=0 台面即电机",
          build_stage_trajectory(tight)(F(1)) == F(5, 2))

    # 起始台面不在容许带内 => 拒绝
    try:
        Stroke(
            "bad", make_plf([(0, 0), (2, 1)]),
            backlash=F(1), initial_stage=F(5))
        check("S1 起始容许带校验", False)
    except ValidationError:
        check("S1 起始容许带校验", True)


# ----------------------------------------------------------------------
# S2：换向停顿期间信号变化；纯停顿不造数、净响应无效不造数
# ----------------------------------------------------------------------
def case_reversal_dwell() -> None:
    # 电机 0→3（t=3）端点停顿到 t=5 →0（t=9）
    stroke = Stroke(
        "dwell", make_plf([(0, 0), (3, 3), (5, 3), (9, 0)]),
        backlash=F(1), initial_stage=F(0))
    traj = build_stage_trajectory(stroke)
    fwd, rev = split_monotonic_runs(stroke)
    # 净响应在停顿段 t∈[3,5] 仍从 30 变到 60
    net = make_plf([(0, 0), (3, 30), (5, 60), (9, 90)])
    rf = sample_row(stroke, fwd, traj, F(0), net, "P-DW", 1, [0, 1, 2, 3])
    rr = sample_row(stroke, rev, traj, F(0), net, "P-DW", 1, [0, 1, 2, 3])
    df = {s.column: (s.status, s.value) for s in rf}
    dr = {s.column: (s.status, s.value) for s in rr}
    check("S2 正扫列0..2有效、列3缺覆盖",
          df == {0: (STATUS_OK, F(10)), 1: (STATUS_OK, F(20)),
                 2: (STATUS_OK, F(30)),
                 3: (STATUS_NO_COVERAGE, None)})
    check("S2 反扫只覆盖台面回退值域1..2",
          dr[0][0] == STATUS_NO_COVERAGE and dr[3][0]
          == STATUS_NO_COVERAGE and dr[1][0] == STATUS_OK
          and dr[2][0] == STATUS_OK)
    t_rev_2 = next(s.time for s in rr if s.column == 2)
    check("S2 停顿不造数：反扫列2时刻在停顿之后",
          t_rev_2 > F(5), f"t={t_rev_2}")
    check("S2 反扫列2读数取运动一侧而非停顿值",
          dr[2][1] == F(80))  # 停顿末值 60 被排除
    # 净响应有效域只到 t=7：反扫时刻 23/3≈7.67、9 均在域外
    net_short = make_plf([(0, 0), (3, 30), (5, 60), (7, 80)])
    rr_bad = sample_row(stroke, rev, traj, F(0), net_short,
                        "P-DW2", 2, [1, 2])
    check("S2 净响应无效=>invalid_net 不造数",
          all(s.status == STATUS_INVALID_NET and s.value is None
              for s in rr_bad))
    # 无净响应产品（product_id=None）同样不得造数
    rr_none = sample_row(stroke, rev, traj, F(0), None, None,
                         None, [1, 2])
    check("S2 缺净响应产品=>invalid_net",
          all(s.status == STATUS_INVALID_NET for s in rr_none))


# ----------------------------------------------------------------------
# S3：出处（行程/时刻/源版本）与原点偏移；正扫反扫
# ----------------------------------------------------------------------
def case_provenance_and_offset() -> None:
    stroke = Stroke(
        "prov", make_plf([(0, 0), (4, 4)]),
        backlash=F(0), initial_stage=F(0))
    traj = build_stage_trajectory(stroke)
    run = split_monotonic_runs(stroke)[0]
    net = make_plf([(0, 100), (4, 500)])  # net(t)=100+100t
    # 偏移 +10：实物列 = 台面 +10；列 10..14 ↔ t=0..4
    res = sample_row(stroke, run, traj, F(10), net, "P-PV", 3,
                     range(10, 15))
    check("S3 偏移后列覆盖10..14",
          all(s.status == STATUS_OK for s in res))
    check("S3 出处完整(行程/产品/基线版本/时刻)",
          all(s.stroke_id == "prov" and s.product_id == "P-PV"
              and s.baseline_version == 3 and s.time == F(s.column - 10)
              for s in res))
    check("S3 插值读数精确",
          {s.column: s.value for s in res}
          == {10: F(100), 11: F(200), 12: F(300), 13: F(400), 14: F(500)})
    # 偏移左移一列（-1）：原列 c 的值现在出现在 c-1
    res2 = sample_row(stroke, run, traj, F(-1), net, "P-PV", 3,
                      range(-1, 4))
    check("S3 原点左移一列后读数迁移",
          {s.column: s.value for s in res2 if s.status == STATUS_OK}
          == {-1: F(100), 0: F(200), 1: F(300), 2: F(400), 3: F(500)})
    # 反扫出处方向标记
    rev_stroke = Stroke(
        "provR", make_plf([(0, 4), (4, 0)]),
        backlash=F(0), initial_stage=F(4))
    rev_traj = build_stage_trajectory(rev_stroke)
    rev_run = split_monotonic_runs(rev_stroke)[0]
    check("S3 反扫方向=-1", rev_run.direction == -1)
    rrev = sample_row(rev_stroke, rev_run, rev_traj, F(0), net,
                      "P-PV", 3, [0, 1, 2, 3, 4])
    # 台面 4→0：列4 ↔ t=0，列0 ↔ t=4
    check("S3 反扫时刻反向对应",
          {s.column: s.time for s in rrev}
          == {4: F(0), 3: F(1), 2: F(2), 1: F(3), 0: F(4)})


# ----------------------------------------------------------------------
# 真实合成行程 + 净响应产品，经完整重建管线重放题面三个接缝结果
# ----------------------------------------------------------------------
def _publish_net(svc: SignalService, run_id: str, signal: list,
                 baseline: list) -> str:
    svc.create_run(run_id, signal)
    version = svc.add_baseline(run_id, baseline)
    return svc.publish_product(run_id, version).product_id


def _forward_stroke(scan: ScanService, sid: str, g=F(1)) -> None:
    # 电机 t=0..6 从 0 线性到 6，台面起 0；g=1 => 台面 t=1 起贴右边，
    # 连续运动覆盖台面 [0,5]（t∈[1,6]），足以覆盖实物列 1..4。
    scan.create_stroke(
        sid, [(t, t) for t in range(7)], g=g, initial_stage=F(0))


def case_three_replays(db: str) -> None:
    svc = SignalService(db)
    scan = ScanService(db)

    # ---- 真实净响应产品（基线恒 0，净响应=信号）----
    # 左图两行恒 0
    pid_l1 = _publish_net(svc, "SIG-L1", [(t, 0) for t in range(7)],
                          [(0, 0), (6, 0)])
    pid_l2 = _publish_net(svc, "SIG-L2", [(t, 0) for t in range(7)],
                          [(0, 0), (6, 0)])
    # 右图行1：实物列1,2,3 → 0,8,1（台面 s=t-1、偏移0 => 列c ↔ t=c+1，
    # 故 net(2)=0、net(3)=8、net(4)=1）
    sig_r1 = [0, 0, 0, 8, 1, 0, 0]  # 索引即 t
    pid_r1 = _publish_net(svc, "SIG-R1",
                          [(t, sig_r1[t]) for t in range(7)],
                          [(0, 0), (6, 0)])
    # 右图行2：实物列1..4 → 8,8,0,0（列c ↔ t=c+1）
    sig_r2 = [0, 0, 8, 8, 0, 0, 0]
    pid_r2 = _publish_net(svc, "SIG-R2",
                          [(t, sig_r2[t]) for t in range(7)],
                          [(0, 0), (6, 0)])

    # 四条正扫行程（每行一条，取第 0 次单调行程，原点 v1=偏移0）
    for sid in ("K-L1", "K-L2", "K-R1", "K-R2"):
        _forward_stroke(scan, sid)
    scan.add_scan_row("row-L1", "K-L1", 0, 1, pid_l1, 1)
    scan.add_scan_row("row-L2", "K-L2", 0, 1, pid_l2, 1)
    scan.add_scan_row("row-R1", "K-R1", 0, 1, pid_r1, 1)
    scan.add_scan_row("row-R2", "K-R2", 0, 1, pid_r2, 1)
    cols = [1, 2, 3]
    plan = {rid: cols for rid in ("row-L1", "row-L2", "row-R1", "row-R2")}
    rebuilt = {r.row_id: r for r in scan.reconstruct(plan)}

    # 先校验重建管线本身给出的读数（台面 s=t-1：列 c ↔ t=c+1）
    def vmap(rid):
        return rebuilt[rid].value_map()
    check("S4 重建左图两行均零",
          vmap("row-L1") == {1: F(0), 2: F(0), 3: F(0)}
          and vmap("row-L2") == {1: F(0), 2: F(0), 3: F(0)})
    check("S4 重建右图行1=(0,8,1)",
          vmap("row-R1") == {1: F(0), 2: F(8), 3: F(1)})
    check("S4 重建右图行2重叠=(8,8,0)",
          vmap("row-R2") == {1: F(8), 2: F(8), 3: F(0)})
    # 出处抽查：列 2 的时刻与源版本
    s_r1_c2 = next(s for s in rebuilt["row-R1"].samples if s.column == 2)
    check("S4 出处时刻/行程/版本",
          s_r1_c2.time == F(3) and s_r1_c2.stroke_id == "K-R1"
          and s_r1_c2.origin_version == 1
          and s_r1_c2.baseline_version == 1
          and s_r1_c2.product_id == pid_r1)

    # ---- 导入为两条已定位图带（印样端独立导入）----
    left = scan.import_located_strip(
        "STRIP-L", (vmap("row-L1"), vmap("row-L2")))
    right = scan.import_located_strip(
        "STRIP-R", (vmap("row-R1"), vmap("row-R2")))
    sheet = scan.create_print_sheet("SHEET-A", 2, 1, 3)

    # 重放 1：全通行 => 缝 (3,3)，代价 1
    r1 = scan.stitch("SHEET-A", "STRIP-L", "STRIP-R",
                     full_columns=(0, 1, 2, 3, 4))
    check("S4 重放1 缝=(3,3) 代价1",
          [s.seam_column for s in r1.seams] == [3, 3]
          and r1.total_cost == F(1)
          and [s.cost for s in r1.seams] == [F(1), F(0)],
          f"seams={[s.seam_column for s in r1.seams]} "
          f"cost={r1.total_cost}")
    check("S4 重放1 可行且无阻断行",
          r1.feasible and r1.first_blocked_row is None)

    # 重放 2：固定图带只禁首行第 3 列 => 缝 (1,1)，代价 8
    scan.create_print_sheet("SHEET-B", 2, 1, 3, disabled=((0, 3),))
    r2 = scan.stitch("SHEET-B", "STRIP-L", "STRIP-R")
    check("S5 重放2 禁(0,3)=>缝=(1,1) 代价8",
          [s.seam_column for s in r2.seams] == [1, 1]
          and r2.total_cost == F(8),
          f"seams={[s.seam_column for s in r2.seams]} "
          f"cost={r2.total_cost}")

    # 重放 3：第二行行程原点订正为向左移一列（新偏移 -1，新版本 v2）
    v2 = scan.correct_origin("K-R2", F(-1))
    check("S6 原点订正追加新版本v2", v2 == 2
          and scan.list_origin_versions("K-R2") == [1, 2])
    # 台面轨迹对象在订正前后不变（电机运动未重算）
    traj_before = scan.get_stage_trajectory("K-R2")
    scan.add_scan_row("row-R2b", "K-R2", 0, 2, pid_r2, 1)
    rebuilt2 = {r.row_id: r for r in
                scan.reconstruct({"row-R2b": [1, 2, 3, 4]})}
    traj_after = scan.get_stage_trajectory("K-R2")
    check("S6 原点订正不重算电机运动(同一轨迹)",
          traj_before.stage_plf.fingerprint()
          == traj_after.stage_plf.fingerprint())
    # 订正后（偏移 -1）：实物列 c 取原台面列 c+1，即原实物列 c+1 的值；
    # 原实物列1..4=(8,8,0,0)，故新列1..4 = 原列2..5 = (8,0,0,0)，
    # 重叠列 1..3 上即题面的 (8,0,0)。
    check("S6 订正后行2重叠=(8,0,0)",
          rebuilt2["row-R2b"].value_map()
          == {1: F(8), 2: F(0), 3: F(0), 4: F(0)})
    right2 = scan.import_located_strip(
        "STRIP-R2", (vmap("row-R1"), rebuilt2["row-R2b"].value_map()))
    r3 = scan.stitch("SHEET-A", "STRIP-L", "STRIP-R2")
    check("S6 重放3 缝=(1,2) 代价0",
          [s.seam_column for s in r3.seams] == [1, 2]
          and r3.total_cost == F(0),
          f"seams={[s.seam_column for s in r3.seams]} "
          f"cost={r3.total_cost}")
    svc.close()
    scan.close()


# ----------------------------------------------------------------------
# S7：无法接续的首行定位 + 缺测保持缺测
# ----------------------------------------------------------------------
def case_blocked_and_missing() -> None:
    cols = (1, 2, 3)
    L = LocatedStrip(
        "L", ({1: F(0), 2: F(0), 3: F(0)},
              {1: F(0), 2: F(0), 3: F(0)}))
    R = LocatedStrip(
        "R", ({1: F(1), 2: F(1), 3: F(1)},
              {3: F(5)}))  # 第二行只有列3有值
    sheet = PrintSheet("OK", 2, cols, frozenset())
    # 第一行缝只能在 1/2/3；第二行合法列只有 3，要求首行缝≥2 —— 可接续
    r_ok = find_seams(L, R, sheet)
    check("S7 可接续时给出缝",
          r_ok.feasible and r_ok.seams[1].seam_column == 3
          and abs(r_ok.seams[0].seam_column - 3) <= 1)
    # 禁用首行全部候选 => 首行即阻断
    sheet_block0 = PrintSheet(
        "B0", 2, cols, frozenset({(0, 1), (0, 2), (0, 3)}))
    r_b0 = find_seams(L, R, sheet_block0)
    check("S7 首行无合法列=>阻断行0",
          not r_b0.feasible and r_b0.first_blocked_row == 0)
    # 第二行合法列与首行无法接续（首行只留列1，第二行只有列3，差2）
    sheet_block1 = PrintSheet(
        "B1", 2, cols, frozenset({(0, 2), (0, 3)}))
    r_b1 = find_seams(L, R, sheet_block1)
    check("S7 第二行无法接续=>阻断行1",
          not r_b1.feasible and r_b1.first_blocked_row == 1
          and [s.seam_column for s in r_b1.seams] == [1])

    # 缺测保持：左图行1缺列2，缝选列1时列2应由右图，但右图也缺 => missing
    L_gap = LocatedStrip("L", ({1: F(0)}, {1: F(0), 2: F(0), 3: F(0)}))
    R_gap = LocatedStrip(
        "R", ({1: F(0), 3: F(0)}, {1: F(0), 2: F(0), 3: F(0)}))
    sheet_g = PrintSheet("G", 2, cols, frozenset())
    r_g = find_seams(L_gap, R_gap, sheet_g, full_columns=(1, 2, 3))
    row0 = r_g.composite[0]
    by_col = {c.column: c for c in row0}
    check("S7 来源缺测保持缺测(不跨图填补)",
          by_col[2].status == CELL_MISSING and by_col[2].value is None)
    check("S7 缝列归左图、其后归右图",
          by_col[1].status == CELL_LEFT
          and by_col[3].status == CELL_RIGHT)


# ----------------------------------------------------------------------
# S8：换印样单不重算电机运动（持久化层面）
# ----------------------------------------------------------------------
def case_sheet_switch_decoupling(db: str) -> None:
    scan = ScanService(db)
    scan.create_stroke(
        "K", [(t, t) for t in range(5)], g=F(1), initial_stage=F(0))
    traj_fp = scan.get_stage_trajectory("K").stage_plf.fingerprint()
    # 同两条图带，配两张不同印样单；切换印样单只改约束。
    # 左图恒 0，右图 (1,2,9)：无禁用时最优缝取最小代价列 1；
    # 禁用列 1 后最优缝被迫移到列 2。
    scan.import_located_strip(
        "A", ({1: F(0), 2: F(0), 3: F(0)},))
    scan.import_located_strip(
        "B", ({1: F(1), 2: F(2), 3: F(9)},))
    scan.create_print_sheet("SH1", 1, 1, 3)
    scan.create_print_sheet("SH2", 1, 1, 3, disabled=((0, 1),))
    r_a = scan.stitch("SH1", "A", "B")
    r_b = scan.stitch("SH2", "A", "B")
    check("S8 两张印样单给出不同缝",
          r_a.seams[0].seam_column == 1
          and r_b.seams[0].seam_column == 2
          and r_a.total_cost == F(1) and r_b.total_cost == F(2),
          f"{r_a.seams[0].seam_column}/{r_a.total_cost} vs "
          f"{r_b.seams[0].seam_column}/{r_b.total_cost}")
    check("S8 换单后台面轨迹指纹不变",
          scan.get_stage_trajectory("K").stage_plf.fingerprint() == traj_fp)
    # 旧印样单仍可复读（不可变历史）
    check("S8 旧印样单可复读",
          scan.get_print_sheet("SH1").disabled == frozenset())
    scan.close()


def main() -> int:
    cases = [
        case_backlash_mechanics,
        case_reversal_dwell,
        case_provenance_and_offset,
        case_blocked_and_missing,
    ]
    for fn in cases:
        fn()
    with tempfile.TemporaryDirectory() as tmp:
        case_three_replays(os.path.join(tmp, "replay.sqlite"))
        case_sheet_switch_decoupling(os.path.join(tmp, "sheet.sqlite"))
    passed = sum(1 for _, ok, _ in results if ok)
    print(f"\n{passed}/{len(results)} 项断言通过")
    return 0 if passed == len(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
