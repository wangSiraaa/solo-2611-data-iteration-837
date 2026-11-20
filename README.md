# 色谱基线产品与重叠峰组分定量仓

纯 Python 3.11 标准库实现（SQLite + `fractions.Fraction`，无第三方依赖），
只处理本地合成的有理数时间序列，**不涉及真实仪器或化学物质识别**。

## 解决的问题

同一条色谱中两个组分的峰会重叠。若对每个组分的归属区间分别积分，共同部分
会被重复取走，使总量超过真实观测。本系统把「归属」建模为**同一时刻的份额
归一**：W(t)=Σ 权重，份额 = 权重/W，重叠处按份额瓜分同一份净响应，从机制上
杜绝重复计数。

## 核心设计

### 1. 全程有理、分段线性（`segments.py`）

* 采样点 `(t,y)`、基线锚点、权重节点全部以 `Fraction` 精确存储；相邻点线性插值。
* 基线由 ≥2 个位于采样范围内、时间互异的锚点给出，锚点覆盖的闭范围即有效域；
  域外查询/锚点一律拒绝，**不外推**。
* 净响应 `max(signal-baseline, 0)`：在信号与基线折点并集上逐段线性求根，
  交点被**真实插入为折点**（不是采样点截零后连线），零段显式保留。
* `PLF.canonical()` 去除共线冗余折点；`fingerprint()` 对规范化后的函数内容
  取哈希——同一条函数无论用多少分段表示，指纹相同。

### 2. 定量引擎（`quant.py`）

* 每组分：固定正有理响应因子 + 一组非负权重节点；节点间线性、范围外恒 0。
* 事件网格 = 净响应折点 ∪ 权重节点 ∪ 各组分支撑边界（裁剪到有效域），
  因此**不依赖固定采样步长恰好命中窄峰**。
* 每个事件格内再求权重和 W 的内部零点并切成子格，统一处理：
  支撑边界的半开间断、内部零权重间隙、W 从零变正/归零。
* W=0 子格的净响应计入「未分配」；W>0 子格按份额积分。

### 3. 积分策略（`integration.py`）

份额相除后 `u·n/W` 是「二次多项式 / 一次多项式」，非多项式。主路径用闭式
解析解：多项式除法分解为 `αs + β + γ/W(s)`，积分含一个对数项，端点 W=0 时
`W·ln W → 0` 取极限。两条**代码完全独立**的核对路径：

* `adaptive_simpson_float`：float 递归自适应 Simpson；
* `composite_simpson_fraction`：纯 `Fraction` 复合 Simpson，全程无浮点。

总面积、未分配面积用精确有理梯形；守恒以精确值核对，残差超 1e-8 直接报错，
**禁止反扣任何组分凑守恒**。

### 4. 持久化与不可变历史（`store.py` / `signal_service.py` / `repository.py`）

* 分数以整数分子/分母存入 SQLite。
* 基线「更正」只追加新版本；产品与分析结果不可变，旧方法/旧产品/旧结果
  永远可复读。
* 分析缓存键 = 方法 × 净响应函数指纹：不同输入分段但同函数命中同一缓存。
* 信号服务只做基线扣除与区间面积，**不接受组分方法**；定量仓只消费公开分段
  产品，也能分析 `product_from_dict` 读入的外部等价产品（导入时复算校验）。
* 组分名称只是身份：结果按方法定义顺序返回，不做字典序重排。

## 规模上限（`limits.py`）

每条信号 500 采样点、20 基线锚点；每方法 10 组分、每组分 30 权重节点；
50 次分析（缓存命中不计入）。

## 入口

```bash
# 完整自动化验证（54 项断言，含手算例与三路径积分核对）
python -m chroma_quant.validation

# 可复跑端到端演示
python -m chroma_quant.demo

# CLI（在仓库根目录执行，或设置 PYTHONPATH）
python -m chroma_quant.cli validate
python -m chroma_quant.cli init --db demo.sqlite
python -m chroma_quant.cli import-run --db demo.sqlite --run R --samples '0:1 1:5 2:2'
python -m chroma_quant.cli add-baseline --db demo.sqlite --run R --anchors '0:3 2:1'
python -m chroma_quant.cli publish --db demo.sqlite --run R
python -m chroma_quant.cli query  --db demo.sqlite --run R --a 0 --b 2
python -m chroma_quant.cli area   --db demo.sqlite --run R --a 0 --b 2
python -m chroma_quant.cli add-method --db demo.sqlite --method M \
    --components 'alpha:2:(0,1;2,0) beta:5:(1,1;2,2)'
python -m chroma_quant.cli analyze --db demo.sqlite --method M --product <PID>
python -m chroma_quant.cli explain --db demo.sqlite --method M --product <PID>
python -m chroma_quant.cli export-product --db demo.sqlite --product <PID> > p.json
python -m chroma_quant.cli import-product --db other.sqlite --file p.json
```

## 验证案例（`validation.py`）

| 案例 | 覆盖点 |
|---|---|
| C1 | 可手算常数份额（纯有理、零残差）、名称身份顺序 |
| C2 | 基线跨零、交点不落采样点的真实切段、子采样点插值 |
| C3 | 交点恰好落在采样点 |
| C4 | 宽度 1/3 的窄归属窗，不依赖步长命中；两条独立积分器核对 |
| C5 | 重叠峰不重复取走、零权重间隙→未分配、精确守恒 |
| C6 | 不同输入分段同函数→同指纹、导入去重、缓存命中、坏产品拒绝 |
| C7 | 更正基线锚点后结果改变；旧方法/产品/结果仍可复读 |
| C8 | 有效域外拒绝外推、超范围/单锚点拒绝 |
| C9 | 500/20/10/30/50 全部规模上限 |
| C10 | 严格递增、权重非负、正因子、名称唯一等输入校验 |
| C11 | 净响应与权重折点交错的通用案例，三条积分路径两两一致 |

---

# 往返扫描的回差坐标与双图带拼接缝

全部为**新增模块与新增表**，旧定量命令（C1–C11、CLI、产品格式）完全不变。
仍全程 `Fraction`、分段线性、追加式不可变历史，无第三方依赖。

## 解决的问题

直接把编码器列号（电机位置 u(t)）当作扫描位置，会在电机换向后的一段
时间内重复画同一点：机械台面经回差间隙（非负半宽 g）与电机耦合。

## 回差力学（`motion.py`）

* 台面更新规则：前一台面位置 s 仍落在容许带 `[u-g, u+g]` 内则**不动**，
  带边越过才贴该边移动，即连续意义下 `s ← clamp(s, u-g, u+g)`。
* 台面在电机段内部就可能追上带边，交点被**真实插入为台面折点**
  （与净响应真实切段同思想），不是只在电机折点处 clamp。
* 电机轨迹在平台/换向处分解为极大严格**单调行程**；每个行程剔除开头
  回差带动的纯停顿，得到台面「连续运动部分」。
* 实物坐标 = 台面位置 + 刻度原点偏移。**偏移是可订正的独立版本**，
  订正只追加新版本，台面轨迹与电机运动完全不重算。
* 每个目标整数列在指定行程的连续运动部分按实物坐标反查时刻（运动端点
  按运动一侧取值），再从该扫行的净响应 PLF 精确插值；纯停顿、缺覆盖、
  净响应无效（域外或缺产品）一律缺测，**绝不造数**。
* 每列附完整出处：行程、单调行程号、反查时刻、原点版本、净响应产品、
  基线版本；正扫（+1）反扫（−1）统一处理。

题面回差例（g=1、电机 3→0、折返时台面在 2）：电机先退到 1（带左边
u+g=2 才碰到台面 2）台面才开始回退，轨迹折点为 `(0,2)-(2,2)-(3,1)`。

## 拼接缝（`stitch.py`）

印样端**独立导入**两条已定位图带，只在连续重叠列上拼接：

* 合法缝列：两侧均有值且未被印样单禁用；缝列归左图、其后归右图；
* 相邻行缝列之差至多 1（连通），违例即定位首个无法接续的行；
* 动态规划最小化全部缝列 `|左−右|` 之和（精确有理，无浮点），
  平局以完整缝列路径的字典序最小者为准（从首行起逐列比较）；
* 被选来源缺测的单元格保持缺测，不跨图、不跨列填补。
* 规模上限 80 行 × 每行 80 个重叠列（`limits.py`）。

换印样单只换一组约束，重新跑内存 DP 即可，**不触碰电机运动与图带**。

## 持久化（`scan_service.py` + `store.py` 新表）

新表 `scan_strokes / scan_motor_knots / scan_origin_versions /
scan_rows / print_sheets / print_sheet_disabled / located_strips`
与旧表无外键耦合；净响应复用既有 `products` 公开分段。

## 入口

```bash
# 往返扫描/拼接验证（38 项断言，与旧 54 项独立）
python -m chroma_quant.scan_validation
python -m chroma_quant.cli validate-scan
# 可复跑端到端演示（含换向停顿信号变化例与三个接缝重放）
python -m chroma_quant.scan_demo

# 电机行程（u(t)、g、起始台面），自动分解单调行程
python -m chroma_quant.cli create-stroke --db d.sqlite --stroke K \
    --motor '0:0 1:1 2:2 3:3' --g 1 --stage0 0
# 原点订正（只追加版本，不重算运动）
python -m chroma_quant.cli correct-origin --db d.sqlite --stroke K --offset -1
# 扫行说明：哪条行程的哪次单调行程、哪版原点、哪个净响应产品
python -m chroma_quant.cli add-scan-row --db d.sqlite --row r1 --stroke K \
    --run-index 0 --origin-version 1 --product <PID> --baseline-version 1
# 按整数列重建，输出读数与完整出处（含缺测状态）
python -m chroma_quant.cli reconstruct --db d.sqlite --row-cols 'r1:1,2,3'
# 印样单（重叠列范围 + 禁用位）与两条已定位图带
python -m chroma_quant.cli create-sheet --db d.sqlite --sheet S \
    --rows 2 --col-lo 1 --col-hi 3 --disabled '0,3'
python -m chroma_quant.cli import-strip --db d.sqlite --strip L \
    --rows-data 'r1:1=0,2=0,3=0' 'r2:1=0,2=0,3=0'
# 接缝与实际拼图（不可行时给出首个无法接续的行）
python -m chroma_quant.cli stitch --db d.sqlite --sheet S \
    --left L --right R --full-cols '0,1,2,3,4'
```

## 扫描验证案例（`scan_validation.py`）

| 案例 | 覆盖点 |
|---|---|
| S1 | 题面 g=1 折返原话例、段内真实交点、平台/换向拆行程、g=0、起始容许带校验 |
| S2 | 换向停顿期间信号变化但停顿不造数、净响应无效/缺产品不造数 |
| S3 | 正扫反扫出处（行程/时刻/源版本）、原点偏移与左移一列、反扫反向时刻 |
| S4 | 真实行程+净响应产品重放：全通行缝 (3,3) 代价 1 |
| S5 | 禁用首行第 3 列：缝 (1,1) 代价 8 |
| S6 | 原点订正左移一列（不重算运动）：缝 (1,2) 代价 0 |
| S7 | 首个无法接续行定位、来源缺测保持缺测、缝列归属 |
| S8 | 换印样单不重算电机运动、旧印样单可复读 |
