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
