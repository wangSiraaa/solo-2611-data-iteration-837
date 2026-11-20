"""题面规定的规模上限。"""

MAX_SAMPLES_PER_RUN = 500
MAX_ANCHORS_PER_BASELINE = 20
MAX_COMPONENTS_PER_METHOD = 10
MAX_NODES_PER_COMPONENT = 30
MAX_ANALYSES = 50

# ---- 往返扫描 / 双图带拼接（题面给定上限）----
MAX_SCAN_ROWS = 80          # 一次拼图至多 80 行
MAX_OVERLAP_COLUMNS = 80    # 每行至多 80 个连续重叠列
MAX_MOTOR_KNOTS = 500       # 单条行程电机位置折点数（与采样上限同量级）
MAX_ORIGIN_VERSIONS = 50    # 每条行程原点订正版本数
