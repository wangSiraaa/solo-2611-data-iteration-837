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
