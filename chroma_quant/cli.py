"""命令行入口。

示例
----
::

  python -m chroma_quant.cli init --db demo.sqlite
  python -m chroma_quant.cli import-run --db demo.sqlite --run R \
      --samples '0:1 1:5 2:2 3:6'
  python -m chroma_quant.cli add-baseline --db demo.sqlite --run R \
      --anchors '0:3 3:1'
  python -m chroma_quant.cli publish --db demo.sqlite --run R
  python -m chroma_quant.cli query --db demo.sqlite --run R --a 0 --b 3
  python -m chroma_quant.cli area --db demo.sqlite --run R --a 0 --b 3
  python -m chroma_quant.cli add-method --db demo.sqlite --method M \
      --components 'A:2:(0,1;3,1) B:5:(2,1;3,2)'
  python -m chroma_quant.cli analyze --db demo.sqlite --method M --product PID
  python -m chroma_quant.cli explain --db demo.sqlite --method M --product PID
  python -m chroma_quant.cli validate

  # 往返扫描 / 双图带拼接
  python -m chroma_quant.cli create-stroke --db demo.sqlite --stroke K \
      --motor '0:0 1:1 2:2 3:3' --g 1 --stage0 0
  python -m chroma_quant.cli correct-origin --db demo.sqlite --stroke K --offset -1
  python -m chroma_quant.cli add-scan-row --db demo.sqlite --row r1 --stroke K \
      --run-index 0 --origin-version 1 --product PID --baseline-version 1
  python -m chroma_quant.cli reconstruct --db demo.sqlite --row-cols 'r1:1,2,3'
  python -m chroma_quant.cli create-sheet --db demo.sqlite --sheet S \
      --rows 2 --col-lo 1 --col-hi 3 --disabled '0,3'
  python -m chroma_quant.cli import-strip --db demo.sqlite --strip L \
      --rows-data 'r1:1=0,2=0,3=0' 'r2:1=0,2=0,3=0'
  python -m chroma_quant.cli stitch --db demo.sqlite --sheet S \
      --left L --right R --full-cols '0,1,2,3,4'
  python -m chroma_quant.cli validate-scan
"""
from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction

from .errors import ChromaError
from .rational import frac_str
from .repository import QuantRepository
from .scan_service import ScanService
from .signal_service import SignalService
from . import validation as validation_mod
from . import scan_validation as scan_validation_mod


def _parse_pairs(text: str) -> list[tuple[Fraction, Fraction]]:
    """解析 "t:y t:y" 形式的有理点序列。"""
    out = []
    for tok in text.split():
        if not tok:
            continue
        t, y = tok.split(":", 1)
        out.append((Fraction(t), Fraction(y)))
    return out


def _parse_components(text: str) -> list[dict]:
    """解析 "名称:因子:(t,w;t,w) 名称:因子:(...)"。"""
    comps = []
    for chunk in text.split():
        name, factor, nodes_raw = chunk.split(":", 2)
        inner = nodes_raw.strip()
        if inner.startswith("(") and inner.endswith(")"):
            inner = inner[1:-1]
        nodes = []
        for pair in inner.split(";"):
            if pair:
                t, w = pair.split(",")
                nodes.append((Fraction(t), Fraction(w)))
        comps.append({"name": name, "factor": Fraction(factor),
                      "nodes": nodes})
    return comps


def _plf_to_list(plf) -> list[list[str]]:
    return [[frac_str(k.t), frac_str(k.y)] for k in plf.knots]


def cmd_init(args) -> None:
    svc = SignalService(args.db)
    svc.close()
    print(f"已初始化数据库 {args.db}")


def cmd_import_run(args) -> None:
    svc = SignalService(args.db)
    svc.create_run(args.run, _parse_pairs(args.samples), note=args.note)
    print(f"已导入 run={args.run}，{len(args.samples.split())} 个采样点")
    svc.close()


def cmd_add_baseline(args) -> None:
    svc = SignalService(args.db)
    if args.correct:
        v = svc.correct_baseline(args.run, _parse_pairs(args.anchors))
        print(f"已更正基线，新版本号 v{v}（旧版本保留）")
    else:
        v = svc.add_baseline(args.run, _parse_pairs(args.anchors))
        print(f"已写入基线版本 v{v}")
    svc.close()


def cmd_versions(args) -> None:
    svc = SignalService(args.db)
    print(json.dumps(svc.list_baseline_versions(args.run)))
    svc.close()


def cmd_publish(args) -> None:
    svc = SignalService(args.db)
    product = svc.publish_product(args.run, args.version)
    print(product.product_id)
    svc.close()


def cmd_query(args) -> None:
    svc = SignalService(args.db)
    q = svc.query(args.run, Fraction(args.a), Fraction(args.b),
                  args.version)
    payload = {
        "run_id": q["run_id"],
        "baseline_version": q["version"],
        "interval": [str(q["interval"][0]), str(q["interval"][1])],
        "coverage": [str(q["coverage"][0]), str(q["coverage"][1])],
        "raw_signal": _plf_to_list(q["raw_signal"]),
        "baseline": _plf_to_list(q["baseline"]),
        "net_response": _plf_to_list(q["net_response"]),
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    svc.close()


def cmd_area(args) -> None:
    svc = SignalService(args.db)
    area = svc.interval_area(args.run, Fraction(args.a),
                             Fraction(args.b), args.version)
    print(frac_str(area), f"≈ {float(area):.12g}")
    svc.close()


def cmd_add_method(args) -> None:
    repo = QuantRepository(args.db)
    repo.create_method(args.method, _parse_components(args.components),
                       note=args.note)
    print(f"已创建方法 {args.method}")
    repo.close()


def cmd_import_product(args) -> None:
    repo = QuantRepository(args.db)
    with open(args.file, encoding="utf-8") as fh:
        pid = repo.import_product(json.load(fh))
    print(pid)
    repo.close()


def cmd_export_product(args) -> None:
    repo = QuantRepository(args.db)
    print(repo.export_product_json(args.product))
    repo.close()


def cmd_analyze(args) -> None:
    repo = QuantRepository(args.db)
    res = repo.run_analysis(args.method, args.product,
                            use_cache=not args.no_cache)
    payload = {
        "method_id": res.method_id,
        "product_id": res.product_id,
        "domain": [str(res.domain[0]), str(res.domain[1])],
        "component_area": res.component_area,
        "concentration": res.concentration,
        "unallocated_area": res.unallocated_area,
        "total_area": res.total_area,
        "conservation_residual": res.conservation_residual,
        "tolerance": res.tolerance,
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    repo.close()


def cmd_explain(args) -> None:
    repo = QuantRepository(args.db)
    print(json.dumps(repo.explain(args.method, args.product),
                     ensure_ascii=False, indent=2))
    repo.close()


# ----------------------------------------------------------------------
# 往返扫描 / 双图带拼接
# ----------------------------------------------------------------------
def _parse_motor(text: str) -> list[tuple[Fraction, Fraction]]:
    """解析 "t:u t:u ..." 电机折点序列。"""
    return _parse_pairs(text)


def cmd_create_stroke(args) -> None:
    scan = ScanService(args.db)
    scan.create_stroke(args.stroke, _parse_motor(args.motor),
                       g=Fraction(args.g),
                       initial_stage=Fraction(args.stage0),
                       note=args.note)
    runs = scan.list_monotonic_runs(args.stroke)
    print(f"已写入行程 {args.stroke}，分解出 {len(runs)} 次单调行程: "
          + ", ".join(f"#{r.index}(方向{r.direction:+d})" for r in runs))
    scan.close()


def cmd_correct_origin(args) -> None:
    scan = ScanService(args.db)
    v = scan.correct_origin(args.stroke, Fraction(args.offset))
    print(f"已追加原点版本 v{v}（台面/电机运动未重算）")
    scan.close()


def cmd_add_scan_row(args) -> None:
    scan = ScanService(args.db)
    scan.add_scan_row(args.row, args.stroke, args.run_index,
                      args.origin_version, args.product,
                      args.baseline_version)
    print(f"已写入扫行 {args.row}")
    scan.close()


def cmd_reconstruct(args) -> None:
    scan = ScanService(args.db)
    # 每行的列计划：'行号:c1,c2,...'；所有行共用同一组列时可用 --cols
    plan: dict[str, list[int]] = {}
    if args.cols is not None:
        cols = [int(x) for x in args.cols.split(",") if x != ""]
        for r in scan.list_scan_rows():
            plan[r.row_id] = cols
    for chunk in args.row_cols or []:
        rid, cols_raw = chunk.split(":", 1)
        plan[rid] = [int(x) for x in cols_raw.split(",") if x != ""]
    rows = scan.reconstruct(plan)
    payload = []
    for r in rows:
        payload.append({
            "row_id": r.row_id,
            "samples": [
                {"column": s.column, "status": s.status,
                 "value": (frac_str(s.value) if s.value is not None
                           else None),
                 "stroke_id": s.stroke_id, "run_index": s.run_index,
                 "time": (frac_str(s.time) if s.time is not None
                          else None),
                 "origin_version": s.origin_version,
                 "product_id": s.product_id,
                 "baseline_version": s.baseline_version}
                for s in r.samples],
        })
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    scan.close()


def cmd_create_sheet(args) -> None:
    scan = ScanService(args.db)
    disabled = []
    for tok in args.disabled or []:
        r, c = tok.split(",", 1)
        disabled.append((int(r), int(c)))
    scan.create_print_sheet(args.sheet, args.rows, args.col_lo,
                            args.col_hi, disabled)
    print(f"已写入印样单 {args.sheet}（{args.rows} 行，"
          f"重叠列 {args.col_lo}..{args.col_hi}，禁用 {len(disabled)} 处）")
    scan.close()


def cmd_import_strip(args) -> None:
    scan = ScanService(args.db)
    # 每行：'行号:列=值,列=值,...'；缺测列直接不写
    rows: list[dict[int, Fraction]] = []
    for chunk in args.rows_data:
        _rid, body = chunk.split(":", 1)
        mapping: dict[int, Fraction] = {}
        for pair in body.split(","):
            if not pair:
                continue
            c, y = pair.split("=", 1)
            mapping[int(c)] = Fraction(y)
        rows.append(mapping)
    scan.import_located_strip(args.strip, rows)
    print(f"已导入图带 {args.strip}（{len(rows)} 行）")
    scan.close()


def cmd_stitch(args) -> None:
    scan = ScanService(args.db)
    full = ([int(x) for x in args.full_cols.split(",")]
            if args.full_cols else None)
    res = scan.stitch(args.sheet, args.left, args.right, full)
    if not res.feasible:
        print(json.dumps({
            "feasible": False,
            "first_blocked_row": res.first_blocked_row,
            "seams_so_far": [
                {"row_index": s.row_index,
                 "seam_column": s.seam_column,
                 "cost": frac_str(s.cost)} for s in res.seams],
        }, ensure_ascii=False, indent=2))
        scan.close()
        return
    payload = {
        "feasible": True,
        "total_cost": frac_str(res.total_cost),
        "seams": [
            {"row_index": s.row_index, "seam_column": s.seam_column,
             "row_cost": frac_str(s.cost)} for s in res.seams],
        "composite": [
            [{"column": c.column, "status": c.status,
              "value": (frac_str(c.value) if c.value is not None else None),
              "source_strip": c.source_strip} for c in row]
            for row in res.composite],
    }
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    scan.close()


def cmd_validate_scan(args) -> None:
    raise SystemExit(scan_validation_mod.main())


def cmd_validate(args) -> None:
    rc = validation_mod.main()
    raise SystemExit(rc)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="chroma_quant",
        description="色谱基线产品与重叠峰组分定量仓")
    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("init"); sp.set_defaults(func=cmd_init)
    sp.add_argument("--db", required=True)

    sp = sub.add_parser("import-run"); sp.set_defaults(func=cmd_import_run)
    sp.add_argument("--db", required=True)
    sp.add_argument("--run", required=True)
    sp.add_argument("--samples", required=True, help='"t:y t:y ..."')
    sp.add_argument("--note")

    sp = sub.add_parser("add-baseline"); sp.set_defaults(func=cmd_add_baseline)
    sp.add_argument("--db", required=True)
    sp.add_argument("--run", required=True)
    sp.add_argument("--anchors", required=True, help='"t:y t:y ..."')
    sp.add_argument("--version", type=int)
    sp.add_argument("--correct", action="store_true",
                    help="更正基线：追加新版本，旧版本保留")

    sp = sub.add_parser("versions"); sp.set_defaults(func=cmd_versions)
    sp.add_argument("--db", required=True); sp.add_argument("--run", required=True)

    sp = sub.add_parser("publish"); sp.set_defaults(func=cmd_publish)
    sp.add_argument("--db", required=True)
    sp.add_argument("--run", required=True)
    sp.add_argument("--version", type=int)

    sp = sub.add_parser("query"); sp.set_defaults(func=cmd_query)
    sp.add_argument("--db", required=True)
    sp.add_argument("--run", required=True)
    sp.add_argument("--a", required=True); sp.add_argument("--b", required=True)
    sp.add_argument("--version", type=int)

    sp = sub.add_parser("area"); sp.set_defaults(func=cmd_area)
    sp.add_argument("--db", required=True)
    sp.add_argument("--run", required=True)
    sp.add_argument("--a", required=True); sp.add_argument("--b", required=True)
    sp.add_argument("--version", type=int)

    sp = sub.add_parser("add-method"); sp.set_defaults(func=cmd_add_method)
    sp.add_argument("--db", required=True)
    sp.add_argument("--method", required=True)
    sp.add_argument("--components", required=True)
    sp.add_argument("--note")

    sp = sub.add_parser("import-product")
    sp.set_defaults(func=cmd_import_product)
    sp.add_argument("--db", required=True)
    sp.add_argument("--file", required=True)

    sp = sub.add_parser("export-product")
    sp.set_defaults(func=cmd_export_product)
    sp.add_argument("--db", required=True)
    sp.add_argument("--product", required=True)

    sp = sub.add_parser("analyze"); sp.set_defaults(func=cmd_analyze)
    sp.add_argument("--db", required=True)
    sp.add_argument("--method", required=True)
    sp.add_argument("--product", required=True)
    sp.add_argument("--no-cache", action="store_true")

    sp = sub.add_parser("explain"); sp.set_defaults(func=cmd_explain)
    sp.add_argument("--db", required=True)
    sp.add_argument("--method", required=True)
    sp.add_argument("--product", required=True)

    # ---- 往返扫描 / 双图带拼接 ----
    sp = sub.add_parser("create-stroke",
                        help="写入电机行程（u(t)、回差半宽、起始台面）")
    sp.set_defaults(func=cmd_create_stroke)
    sp.add_argument("--db", required=True)
    sp.add_argument("--stroke", required=True)
    sp.add_argument("--motor", required=True, help='"t:u t:u ..."')
    sp.add_argument("--g", required=True, help="非负回差半宽")
    sp.add_argument("--stage0", required=True, help="起始台面位置")
    sp.add_argument("--note")

    sp = sub.add_parser("correct-origin",
                        help="追加刻度原点偏移版本（不重算电机运动）")
    sp.set_defaults(func=cmd_correct_origin)
    sp.add_argument("--db", required=True)
    sp.add_argument("--stroke", required=True)
    sp.add_argument("--offset", required=True)

    sp = sub.add_parser("add-scan-row", help="写入扫行说明")
    sp.set_defaults(func=cmd_add_scan_row)
    sp.add_argument("--db", required=True)
    sp.add_argument("--row", required=True)
    sp.add_argument("--stroke", required=True)
    sp.add_argument("--run-index", type=int, required=True)
    sp.add_argument("--origin-version", type=int, required=True)
    sp.add_argument("--product")
    sp.add_argument("--baseline-version", type=int)

    sp = sub.add_parser("reconstruct",
                        help="按整数列重建扫行，输出读数与完整出处")
    sp.set_defaults(func=cmd_reconstruct)
    sp.add_argument("--db", required=True)
    sp.add_argument("--cols", help="对全部扫行共用的列，如 '1,2,3'")
    sp.add_argument("--row-cols", nargs="*",
                    help="逐行列计划，如 'row1:1,2,3 row2:1,2'")

    sp = sub.add_parser("create-sheet", help="写入印样单（重叠列+禁用位）")
    sp.set_defaults(func=cmd_create_sheet)
    sp.add_argument("--db", required=True)
    sp.add_argument("--sheet", required=True)
    sp.add_argument("--rows", type=int, required=True)
    sp.add_argument("--col-lo", type=int, required=True)
    sp.add_argument("--col-hi", type=int, required=True)
    sp.add_argument("--disabled", nargs="*",
                    help="禁用位置 '行,列'，如 '0,3'")

    sp = sub.add_parser("import-strip", help="独立导入一条已定位图带")
    sp.set_defaults(func=cmd_import_strip)
    sp.add_argument("--db", required=True)
    sp.add_argument("--strip", required=True)
    sp.add_argument("--rows-data", nargs="+", required=True,
                    help="每行 '行号:列=值,列=值'，缺测列不写")

    sp = sub.add_parser("stitch", help="双图带接缝动态规划与拼图输出")
    sp.set_defaults(func=cmd_stitch)
    sp.add_argument("--db", required=True)
    sp.add_argument("--sheet", required=True)
    sp.add_argument("--left", required=True)
    sp.add_argument("--right", required=True)
    sp.add_argument("--full-cols", help="拼图输出列范围，如 '0,1,2,3,4'")

    sp = sub.add_parser("validate-scan",
                        help="往返扫描/拼接验证（38 项断言）")
    sp.set_defaults(func=cmd_validate_scan)

    sp = sub.add_parser("validate"); sp.set_defaults(func=cmd_validate)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        args.func(args)
    except ChromaError as exc:
        print(f"错误: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
