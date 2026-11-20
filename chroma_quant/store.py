"""SQLite 持久层：建表与分数（Fraction ↔ 整数分子/分母）读写。

所有历史对象不可变：更正基线只会新增版本与新产品，旧版本、旧产品、
旧分析永远保留并可复读。
"""
from __future__ import annotations

import sqlite3
from fractions import Fraction
from typing import Any, Iterable

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS runs (
    run_id      TEXT PRIMARY KEY,
    created_at  TEXT NOT NULL,
    note        TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    run_id   TEXT NOT NULL REFERENCES runs(run_id),
    idx      INTEGER NOT NULL,
    t_num    INTEGER NOT NULL,
    t_den    INTEGER NOT NULL,
    y_num    INTEGER NOT NULL,
    y_den    INTEGER NOT NULL,
    PRIMARY KEY (run_id, idx)
);

CREATE TABLE IF NOT EXISTS baseline_versions (
    run_id      TEXT NOT NULL REFERENCES runs(run_id),
    version     INTEGER NOT NULL,
    created_at  TEXT NOT NULL,
    supersedes  INTEGER,
    PRIMARY KEY (run_id, version)
);

CREATE TABLE IF NOT EXISTS baseline_anchors (
    run_id   TEXT NOT NULL,
    version  INTEGER NOT NULL,
    idx      INTEGER NOT NULL,
    t_num    INTEGER NOT NULL,
    t_den    INTEGER NOT NULL,
    y_num    INTEGER NOT NULL,
    y_den    INTEGER NOT NULL,
    PRIMARY KEY (run_id, version, idx),
    FOREIGN KEY (run_id, version) REFERENCES baseline_versions(run_id, version)
);

CREATE TABLE IF NOT EXISTS products (
    product_id        TEXT PRIMARY KEY,
    run_id            TEXT NOT NULL,
    baseline_version  INTEGER NOT NULL,
    fingerprint       TEXT NOT NULL,
    created_at        TEXT NOT NULL,
    FOREIGN KEY (run_id, baseline_version)
        REFERENCES baseline_versions(run_id, version)
);

CREATE TABLE IF NOT EXISTS product_segments (
    product_id  TEXT NOT NULL REFERENCES products(product_id),
    kind        TEXT NOT NULL CHECK (kind IN ('signal','baseline','net')),
    idx         INTEGER NOT NULL,
    t_num       INTEGER NOT NULL,
    t_den       INTEGER NOT NULL,
    y_num       INTEGER NOT NULL,
    y_den       INTEGER NOT NULL,
    PRIMARY KEY (product_id, kind, idx)
);

CREATE TABLE IF NOT EXISTS methods (
    method_id   TEXT PRIMARY KEY,
    created_at  TEXT NOT NULL,
    note        TEXT
);

CREATE TABLE IF NOT EXISTS method_components (
    method_id    TEXT NOT NULL REFERENCES methods(method_id),
    idx          INTEGER NOT NULL,
    name         TEXT NOT NULL,
    factor_num   INTEGER NOT NULL,
    factor_den   INTEGER NOT NULL,
    PRIMARY KEY (method_id, idx),
    UNIQUE (method_id, name)
);

CREATE TABLE IF NOT EXISTS weight_nodes (
    method_id  TEXT NOT NULL,
    comp_idx   INTEGER NOT NULL,
    idx        INTEGER NOT NULL,
    t_num      INTEGER NOT NULL,
    t_den      INTEGER NOT NULL,
    w_num      INTEGER NOT NULL,
    w_den      INTEGER NOT NULL,
    PRIMARY KEY (method_id, comp_idx, idx),
    FOREIGN KEY (method_id, comp_idx)
        REFERENCES method_components(method_id, idx)
);

CREATE TABLE IF NOT EXISTS analyses (
    analysis_id           TEXT NOT NULL,
    method_id             TEXT NOT NULL REFERENCES methods(method_id),
    product_id            TEXT NOT NULL REFERENCES products(product_id),
    created_at            TEXT NOT NULL,
    fingerprint           TEXT NOT NULL,
    component_name        TEXT NOT NULL,
    area                  REAL NOT NULL,
    concentration         REAL NOT NULL,
    unallocated_area      REAL NOT NULL,
    total_area            REAL NOT NULL,
    conservation_residual REAL NOT NULL,
    tolerance             REAL NOT NULL,
    domain_t0_num         INTEGER NOT NULL,
    domain_t0_den         INTEGER NOT NULL,
    domain_t1_num         INTEGER NOT NULL,
    domain_t1_den         INTEGER NOT NULL,
    PRIMARY KEY (analysis_id, component_name)
);

CREATE INDEX IF NOT EXISTS idx_products_fp ON products(fingerprint);
CREATE INDEX IF NOT EXISTS idx_analyses_lookup
    ON analyses(method_id, product_id);

-- ======================================================================
-- 往返扫描 / 双图带拼接（全部为追加式不可变历史，与旧表无外键耦合）
-- ======================================================================

-- 电机行程：u(t) 分段线性折点 + 回差半宽 + 起始台面位置
CREATE TABLE IF NOT EXISTS scan_strokes (
    stroke_id        TEXT PRIMARY KEY,
    g_num            INTEGER NOT NULL,
    g_den            INTEGER NOT NULL,
    s0_num           INTEGER NOT NULL,
    s0_den           INTEGER NOT NULL,
    created_at       TEXT NOT NULL,
    note             TEXT
);

CREATE TABLE IF NOT EXISTS scan_motor_knots (
    stroke_id  TEXT NOT NULL REFERENCES scan_strokes(stroke_id),
    idx        INTEGER NOT NULL,
    t_num      INTEGER NOT NULL,
    t_den      INTEGER NOT NULL,
    u_num      INTEGER NOT NULL,
    u_den      INTEGER NOT NULL,
    PRIMARY KEY (stroke_id, idx)
);

-- 刻度原点订正：只追加新版本；台面轨迹不依赖偏移，订正不重算电机运动
CREATE TABLE IF NOT EXISTS scan_origin_versions (
    stroke_id    TEXT NOT NULL REFERENCES scan_strokes(stroke_id),
    version      INTEGER NOT NULL,
    off_num      INTEGER NOT NULL,
    off_den      INTEGER NOT NULL,
    created_at   TEXT NOT NULL,
    supersedes   INTEGER,
    PRIMARY KEY (stroke_id, version)
);

-- 扫行说明：取哪条行程的哪次单调行程、哪版原点、哪个净响应产品
CREATE TABLE IF NOT EXISTS scan_rows (
    row_id             TEXT PRIMARY KEY,
    stroke_id          TEXT NOT NULL REFERENCES scan_strokes(stroke_id),
    run_index          INTEGER NOT NULL,
    origin_version     INTEGER NOT NULL,
    product_id         TEXT,
    baseline_version   INTEGER,
    created_at         TEXT NOT NULL
);

-- 印样单：重叠列范围 + 禁用 (行,列)；换单只是换一组约束，不触发任何运动重算
CREATE TABLE IF NOT EXISTS print_sheets (
    sheet_id     TEXT PRIMARY KEY,
    n_rows       INTEGER NOT NULL,
    col_lo       INTEGER NOT NULL,
    col_hi       INTEGER NOT NULL,
    created_at   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS print_sheet_disabled (
    sheet_id  TEXT NOT NULL REFERENCES print_sheets(sheet_id),
    row_idx   INTEGER NOT NULL,
    col       INTEGER NOT NULL,
    PRIMARY KEY (sheet_id, row_idx, col)
);

-- 已定位图带：每行每实物列一个读数；缺测列直接没有行（不造数）
CREATE TABLE IF NOT EXISTS located_strips (
    strip_id   TEXT NOT NULL,
    row_idx    INTEGER NOT NULL,
    col        INTEGER NOT NULL,
    y_num      INTEGER NOT NULL,
    y_den      INTEGER NOT NULL,
    PRIMARY KEY (strip_id, row_idx, col)
);

CREATE INDEX IF NOT EXISTS idx_strips_lookup
    ON located_strips(strip_id, row_idx);
"""


def connect(path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


def put_frac(num_col: str, den_col: str, f: Fraction) -> dict[str, int]:
    return {num_col: f.numerator, den_col: f.denominator}


def get_frac(row: sqlite3.Row | dict[str, Any],
             num_col: str, den_col: str) -> Fraction:
    return Fraction(int(row[num_col]), int(row[den_col]))


def insert_many(conn: sqlite3.Connection, sql: str,
                rows: Iterable[dict[str, Any]]) -> None:
    conn.executemany(sql, list(rows))
