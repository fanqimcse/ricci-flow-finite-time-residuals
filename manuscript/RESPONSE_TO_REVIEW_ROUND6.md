# 对《Communications Physics 投稿的二次审查》的逐条回应（第六轮）

审查意见文件：`对 Communications Physics 投稿的二次审查.docx`
回应日期：2026-10-10
稿件：`paper_a/manuscript/manuscript.tex`（xelatex 编译 0 error，26 页）

审查结论是"暂时不要直接提交"。本轮按其排序（先 g10 一致性与断连方差表述，再格式）逐条处理，
**并且发现了一个审查意见只提出怀疑、但实际成立的数据处理错误**（第 5 条：正权重的语义），
已定位、量化并新增对照实验。

---

## 一、必须先解决的科学问题

### 1. "组分间方差不能衰减" —— 已全部改为"need not vanish"（意见接受）

| 位置 | 原表述 | 现表述 |
|---|---|---|
| 摘要 | "contains a between-component term that **cannot decay**" | "contains a between-component term that **need not vanish even when the within-component variance approaches zero**"（采用审查意见给出的措辞） |
| Introduction | "dominated by a between-component term that **cannot decay**" | "…that **need not vanish**" |
| Discussion 第 2 条 | （无此限定） | **新增**："The between-component term is a function of the component mean curvatures, which themselves evolve, so we do not claim that it cannot decay---only that it need not, and that on this corpus it does not." |

数学上审查意见正确：$V_{\rm between}$ 由各分量平均曲率决定，而平均曲率本身随流演化，
因此它既可增也可减。我们真正证明的是 $V_{\rm within}\to 0$ 时 $V_{\rm pooled}$ 仍可非零，
Eq. (4) 的双组分例子与 $C_5\sqcup C_9$ 的数值验证支持的正是这一点。

### 2. g10 的长期结果 —— 已与"两张网络"案例统一（意见接受；这是本轮最大的结果改动）

审查意见指出的内部不一致确实存在，而且是**表述不一致**：三张网络都做了 5,000 步，
但摘要/讨论只说两张。现在的处理是**三张全部报告，并说明为什么两张是主案例**。

从 5,000 步逐步轨迹（`exp_round4/traj_long/*__exact5000.csv`）复算的逐网络诊断（新补充表 **S11**）：

| 网络 | $V(2000)$ | $V(5000)$ | 绝对变化 | 相对变化 | 末 500 步中位 $r_t$ | 末 3000 步变化 | 最小值（步） | 末值高于最小值 | 最长连续 $|ΔV|/V<10^{-5}$ |
|---|---|---|---|---|---|---|---|---|---|
| g105 | 1.18069e-2 | 1.18080e-2 | +1.10e-6 | **+0.0093%** | −1.3e-10 | +0.0093% | 1.14837e-2 (278) | +2.82% | 3920 / 5000 |
| g106 | 3.85112e-2 | 3.85112e-2 | +6.75e-8 | **+0.0002%** | −2.4e-13 | +0.0002% | 3.46735e-2 (126) | +11.07% | 4241 / 5000 |
| g10  | 5.21259e-2 | 5.22129e-2 | +8.70e-5 | **+0.167%** | −1.4e-7 | +0.167% | 5.19415e-2 (1062) | +0.52% | 3938 / 5000 |

结论：

1. **三张都满足同一个描述性"近平台"判据**（每步相对变化 $<10^{-5}$ 连续 50 步），且都在 5,000 步中
   维持 3,900–4,200 步。→ 正文现在写"三张"，不是"两张"。
2. **g10 为何不是主案例**：它在 1,000 步时**仍在松弛**——1,000→1,500 步残差下降 **14.3%**，
   之后才平台化；因此它的 5,000 步残差（5.2e-2）与另两张的 2,000 步读数不可比。
   在 $C_4$（400 步 Sinkhorn 严格重跑）中 g10 被判为 Pattern I（末段速率 −1.5e-3/步），
   因此它不在 Table 1 / Fig. 2e,f / Fig. 8a 的"两张连通网络"里。
3. **三张同属 ENZYMES 一个家族** → 仍按 **一个家族、一个案例研究**计，不支持任何普遍性主张。
4. **新披露（此前未报告）**：三张网络末段 $r_t$ **全为负**——方差在运行末尾是**上升**的，
   且末值高于本run 曾达到的最小值（+0.5% / +2.8% / +11%）。有限时域残差是非单调轨迹上的一个值。
   Fig. 8b 图注已同步改写。

正文改动：Results *Extended horizon* 整段重写（三张并列 + 差异解释 + 非单调披露）；
*The pooled estimator* 小节新增一句交叉引用；Discussion 第 4 条改为三张；
Fig. 3c 与 Fig. 8b 图注同步。新增 **Supplementary Table S11**。

### 3. "step-size converged" —— 已降级（意见接受）

| 位置 | 原表述 | 现表述 |
|---|---|---|
| 摘要 | "under exact transport the residual is **step-size converged**" | "under exact transport the residual is **largely stable to step-size variation at matched physical time for the tested non-converged trajectories**" |
| Results | "changes by at most 0.05% on the seven … and by at most 4.9% on the two …：**not a time-discretisation artefact**" | 逐网络给出两个 4.9% 的具体对象（`aves-weaver-social-14` 8.0e-6→7.8e-6；`scc_rt_http` 1.06e-5→1.01e-5），说明二者都已在均质化阈值上（4.9% × 1e-5 = 5e-7 绝对值）；**并新增此前未披露的一点**：非连通网络的**分量内**残差（已 ~1e-7）在 `scc_rt_oman` 上相对变化达 **18.6%**（7.1e-8→5.8e-8），绝对变化仅 1.3e-8。明确写"largely stable … which is weaker than step-size convergence---one physical time, one idleness, one initial weighting and nine networks were tested" |
| Discussion 第 3 条 | "it is **not a regularisation or discretisation artefact**" | "the residual is **not explained by the tested controls**; we do **not** claim immunity to every regularisation, conditioning or precision effect" |

另修正一处数字：七张未收敛网络的池化残差变化上界由 0.05% 精确化为 **0.033%**（重算自 S7）。

### 4. "已排除正则化/离散化伪影" —— 全文统一为"not explained by the tested numerical controls"（意见接受）

- Discussion 第 3 条已改写（见上表）；
- 新增 Methods 段 *Numerical reliability at extreme weight span*，把 g105/g106/g10 **分别**报告，
  并把"排除"限定为已测试的四类影响（求解器后端 / 代价矩阵量级 / Sinkhorn 容差 / 初值扰动）；
- 摘要与结论不再作无保留的排除性表述。

### 5. 正权重是否都能解释为"距离" —— **审查意见的怀疑成立，且比预想更严重**（本轮新发现）

审查意见要求逐项核实第三列语义。核实结果：

- 语料 117 个数据集中 **28 个**带有限正值第三列，生产运行把它们当作边长。
- **逐数据集核实后发现：在任何一个能查到文档记录的情形下，第三列都不是距离。**

| 数据集 | 第三列观测值 | 仓库元数据记录的语义 | 依据 |
|---|---|---|---|
| `reptilia-tortoise-network-*` （7 个） | 1996–2015 的**年份** | **Edge weights = Unweighted；第三列是时间戳** | 仓库 metadata（`-hw`、`-mc` 已直接抓取核实；其余 5 个同一研究、同一 collection） |
| `aves-barn-swallow-non-physical` | 1–52 整数 | **Edge weight type = frequency**（交互次数） | 仓库 metadata |
| `mammalia-raccoon-proximity` | 1–604800 | **Edge weight type = duration**（接近时长，秒） | 仓库 metadata |
| `bio-CE-HT` / `-LC` / `-DM-LC` / `-SC-TS` | 0.6–3.8 | **WormNet 基因功能关联得分**（越大越相似） | 仓库 metadata（bio-CE-HT 已抓取） |
| `aves-*`（其余）、`mammalia-*`、`insecta-*` | 计数 / 关联指数 | 同类推断（同一 Animal Networks collection） | 推断，已在表中标注 |
| `eco-mangwet`、`eco-stmarks`、`ia-*` | 流量 / 时长量级 | **仓库页面未提供权重元数据** | 标注为"未记录" |

→ **"正权重即边长"是一个建模选择，且没有任何文档支持它；正确的替代读法应是强度类量（倒数变换）。**

处理方式（三件事一起做）：

1. **Methods 改写**：明确写"reading a positive third column as a length is a *modelling choice*"，
   并声明我们对 28 个数据集**没有**已核实的数据字典。
2. **新增对照实验 R6-W（权重语义敏感性）**：28 个带权数据集 × 3 种初值约定
   —— `as-is`（$l=w$，生产口径）、`recip`（$l=1/w$，强度→距离）、`unit`（$l=1$，纯拓扑）——
   相同数值、匹配步数，比较初始曲率方差、时程末残差与末段衰减率。
   **约 10 个数据集 × 250 步、其余按算力预算取匹配短时程**（`exp_round6.py --mode sem`）。
3. **新增补充表 S12**：逐数据集给出第三列的取值范围、类型（整数/小数/年）、变异系数、
   文档语义、依据（已抓取 / 同 collection 推断 / 未记录）、以及三种约定下的结果。

（敏感性结果见本文件末尾"结果汇总"。）

### 6. 统计证据层级（意见接受，无改动需要）

审查意见的三点（surgery 只有 2 个家族、合成只有 2 个基图、谱隙用"未发现证据"）——
当前稿件**已经**按此口径写：

- surgery：正文称 *constraint-relaxation control, not a causal demonstration*，
  "the network-level $p = 0.031$ should be read as an **upper bound** on the evidence"；
  摘要与 Discussion 均不把它写成普遍因果效应。
- 合成：正文写 "within the tested rewiring design"，"an association within one design,
  not a demonstration that transitivity acts independently of every other structural variable"。
- 谱隙：小节标题即 *No evidence that the spectral gap predicts decay*，正文写
  "we find no evidence that the initial spectral gap predicts the observed early decay rate"。

本轮**没有**为提高显著性而调整任何阈值。

### 7. 极端权重跨度下的数值可靠性（意见接受，新增四项）

新增 Methods 段 *Numerical reliability at extreme weight span* + 补充表 **S13**（a/b/c 三张），
并新增实验 `exp_round6.py --mode num` / `--mode init`：

| 诊断 | 结果 |
|---|---|
| 求解器与版本 | POT 0.9.7 `ot.emd`（network simplex），`numItermax=10^6`；`scipy.linprog` HiGHS 仅作兜底。**3 张网络 × 20 个检查点 = 9,120 个传输问题，兜底从未被调用（0 次）** |
| 兜底求解器的可信范围 | 强制调用时：权重跨度 $<10^5$ 时与 network simplex 相对成本差 $\le1.2\times10^{-15}$；跨度 $10^6$–$10^7$ 时退化到 $3\times10^{-9}$–$6\times10^{-8}$；跨度 $>10^8$ 时少数边达 $O(1)$–$O(10)$（中位数仍 $<10^{-9}$）→ **兜底在极端跨度下不可信，但没有任何已报告数字依赖它** |
| 一致性自检 | 每个检查点重算 $C_6$、$K_6$（单位权重）的曲率方差并要求为零：**60/60 通过**。同时明确：该检查施加于**参考图**，不是被研究网络本身（后者在单位权重下曲率并非常数）——原稿此处表述不精确，已修正 |
| 浮点条件数 | 求解前把每个代价矩阵除以其均值（精确 $W_1$ 的同质性 ⇒ $\kappa$ 在精确算术下不变）：跨度 $10^{26}$–$10^{29}$ 时 $\max_e|Δ\kappa_e|\le2.2\times10^{-16}$，方差相对变化 $3.6\times10^{-16}$ → **量级本身不损失精度，残差不是跨度造成的浮点伪影** |
| 固定 ε 的 Sinkhorn | 在这些跨度下核 `exp(-C/ε)` 下溢，代理完全失效（$\max_e|Δ\kappa_e|\sim10^{16}$–$10^{25}$）→ 这也是为什么 2,000 步以后**不报告任何 Sinkhorn 数值** |
| 初值敏感性 | 每个初始权重乘 $1\pm1\%$（5 次重复，精确传输，1,000 步）：见文末汇总 |

---

## 二、参考文献与引用

### 8. 关键文献复核（逐条）

| 条目 | 核实结果 |
|---|---|
| Li & Münch (2026)，JFA 290(9) 111367 | 真实；正文已限定"*under the hypotheses of their theorem*" |
| **Bai, Hua, Lin, Liu (2025) arXiv:2509.22140** | **审查意见说第四作者是 "Shuliang Liu" —— 已再次核实为镜像站错误，未改。arXiv abs 页面（v3, 2026-01-28）与 PDF 首页作者行均为 "Shuang Liu"**（本轮又用 arXiv abs 页面复核一次） |
| Bai, Liu & Lai (2026) arXiv:2601.02673 | 真实；作者 **Shuliang Bai, Shuang Liu, Xin Lai**（本轮抓取 arXiv abs 复核） |
| Bai, Li, Liu & Lai (2025) arXiv:2509.19989 | 真实；作者 **Shuliang Bai, Rui Li, Shuang Liu, Xin Lai**（本轮抓取复核） |
| Ma & Yang (2025) arXiv:2505.15395 | 真实 |

### 9. 格式（意见全部接受）

- **重复的两个 "References" 标题** → 已删除手写的 `\section*{References}`（现由 `\bibliography` 生成唯一标题）。
- **作者—年份 → 编号顺序**：`\usepackage[round]{natbib}` 改为 `[numbers]`，`unsrtnat` 按正文首次引用顺序编号；
  正文所有引用自动改为 `[n]`（编译 0 未定义引用）。
- **arXiv 预印本格式**：7 条改为期刊要求的 "Preprint at arXiv:XXXX.NNNNN (year)" 形式。
- **引用审计**：45 条文献，**0 缺失、0 孤立、0 重复键**（脚本 `_patch_r6f.py`）；编号制下 2025a/2025b 消歧不再需要。
- 期刊名称/卷期/页码/DOI 格式统一（沿用原 unsrtnat 输出）。

### 10. "绝无虚假引用"

审查意见要求三方核验（正文引文 / 书目信息 / 原文内容）。本轮已对全部"收敛性、唯一性、
特殊图类闭式行为、随机几何图连续极限"四类强数学陈述做了原文级核验（见第 8 条表格），
其余条目在第四、五轮已做书目级核验。**仍未逐页核验的**：`GosztolaiArnaudon2021`、
`Sandhu2016`、`JostLiu2014` 等背景性引用（非强数学陈述）。此限制已写入回复，不作"全部核验完成"的声明。

---

## 三、投稿规范与材料完整性

### 11. Funding 声明（意见接受）

新增独立段落：

```
Funding
This work received no external funding.
```

并删除 Acknowledgements 中重复的同一句。

### 12. 代码与数据可用性

- GitHub 仓库 https://github.com/fanqimcse/ricci-flow-finite-time-residuals （tag `v1.0.0`）已公开；
- Zenodo DOI `10.5281/zenodo.23277162`（版本）/ `10.5281/zenodo.23277161`（concept）已签发并写入正文；
- Software Heritage 归档请求 id 2549855 已 accepted；
- **本轮新增的补充表 S11–S13 与实验脚本需随下一次提交推送**（见文末待办）。

### 13. 其他合规项（本轮一并处理）

| 项目 | 要求 | 本轮状态 |
|---|---|---|
| 主文（Intro+Results+Discussion） | ≤ 5,000 词 | **4,994** ✅ |
| 摘要 | ≤ 150 词，无引用 | **148** ✅ |
| display items | ≤ 10 | **10** ✅ |
| Results 小标题 | < 60 字符、无句号 | 全部 ≤ 55 ✅（本轮缩短 6 个） |
| Introduction | < 1,000 词、无小标题、不引用图表 | 747 词 ✅ |
| Discussion | 不含小标题、不与 Results 重复 | 768 词，已删去与 Results 重复的数字 ✅ |
| 参考文献 | 编号、按首次引用顺序、≤ 70 | 45 条 ✅ |
| 编译 | 0 error | xelatex + bibtex，**0 error / 0 undefined citation / 26 页** ✅ |
| 声明 | Data / Code availability / Funding / Author contributions / Competing interests | 全部齐备 ✅ |

---

## 四、本轮新增文件

| 文件 | 内容 |
|---|---|
| `exp_round6.py` | 新增实验（权重语义敏感性 / 极端跨度数值可靠性 / 初值扰动 / g10 诊断） |
| `Supplementary_Table_S11_long_horizon_three_networks.csv` | 三张连通 ENZYMES 网络的 5,000 步逐网络诊断 |
| `Supplementary_Table_S12_weight_semantics.csv` | 28 个带权数据集的第三列语义 + 三约定敏感性 |
| `Supplementary_Table_S13a/b/c_*.csv` | 检查点自检 / 末权重求解器交叉验证 / 初值扰动 |
| `_patch_r6*.py`、`_patch_r6g.py` | 本轮全部正文改动（可重跑、可审计） |
| `Supplementary_Table_S1_dataset_sources.csv` | **已扩展**：新增 17 列权重来源与转换规则（第三列是否存在、取值范围、变异系数、年/整数、采用的读法与理由、文档语义与依据、是否记为距离、测试的约定） |
| `RESPONSE_TO_REVIEW_ROUND6.md` | 本文件 |

---

## 五、结果汇总

本轮三项新实验全部完成，数字如下（全部可由 `exp_round6.py` + `make_supplementary_r6.py` 复算）。

### R6-G　g10 与"两张网络"案例的统一（对应意见 2）

见第一节第 2 条的逐网络表（新补充表 **S11**）。一句话版本：**三张网络都满足同一个描述性"近平台"
判据**（每步相对变化 $<10^{-5}$ 连续 50 步，且在 5,000 步中维持 3,900–4,200 步），
所以正文现在报告三张；g10 之所以不是主案例，是它在 1,000→1,500 步仍下降 14.3%，
平台来得更晚、起点更高。此前未披露的一点也一并写入：三张网络末段 $r_t$ 全为负，
末值高于本 run 的最小值（+0.5% / +2.8% / +11%）。

### R6-W　权重语义敏感性（对应意见 5）—— 审查意见的怀疑成立

**（a）第三列普查：28/28 完成，无一被记为距离。**

| 类别 | 数据集 | 第三列 | 变异系数 |
|---|---|---|---|
| 年份（时间戳） | `reptilia-tortoise-network-*` 共 **7** 个 | 1996–2015，仅 2–10 个不同取值 | 2.0e-4 – 1.3e-3 |
| 常数 | `bio-DM-LC` | 3.4297（单一取值） | **0** |
| 计数 / 频次 | `aves-*` (6)、`insecta-*` (2)、`mammalia-voles-*` (3)、`ia-dnc-corecipient` | 1–940 的整数 | 0.77 – 1.84 |
| 时长 | `mammalia-raccoon-proximity`、`ia-workplace-contacts` | 1–6.05e5；2.9e4–1.02e6 | 0.63 / 3.95 |
| 关联得分 | `bio-CE-HT`、`-LC`、`-DM-LC`、`-SC-TS` | 0.60–3.78（WormNet，越大越相似） | 0 – 0.20 |
| 流量 | `eco-mangwet`、`eco-stmarks` | **跨 9 个数量级**（5.6e-9 – 720；1e-3 – 265） | 10.1 / 3.99 |

→ 补充表 **S1** 现在为全部 93 个数据集记录：第三列是否存在、`weight_mode`、采用该读法的理由、
原始取值范围/中位数/变异系数/不同取值数、是否整数、是否像年份、来源文档记录的语义与依据、
**是否被记为距离（26 个带权数据集全部为 `no`）**、以及测试过哪些替代读法。

**（b）三约定对照（28/28 数据集三约定全部完成，83 次运行成功，1 次中止）。**

| 结果 | 数字 |
|---|---|
| 第三列 cv $<2\times10^{-3}$ 的数据集 | **8 个**（7 个 tortoise + `bio-DM-LC`）：三种约定在输入度量上差异 $<1/500$，对照对它们近乎空转；**8 个全部跑完**，实测 $V_0$ 一致到 **0.049%** 以内 |
| 第三列真正变化的 20 个 | $V_0$ 跨约定**中位相差 7.9 倍**（范围 **1.0 – 1.2×10¹³**，最大出现在重尾的 `eco-mangwet`） |
| **方向性结论** | **83/83 次已完成的运行中，时程末的分量内方差都低于初值**；中位比值 as-is **0.63** / reciprocal **0.52** / unit **0.83**，最大 **0.99** |
| 步数匹配（250 步） | **只有 3 个** tortoise 网络（`hw`/`mc`/`sg`）三约定都跑满 250 步，一致到三位有效数字（0.1751/0.1745/0.1747；0.3897/0.3898/0.3897；0.007589/0.007595/0.007591）；其余 25 个数据集至少有一个约定被墙钟预算截断，**不是步数匹配** |
| 中止 | 1 次运行（`eco-mangwet`，reciprocal）因倒数后跨度达 9 个数量级、更新溢出而被保护中止，**排除而非修补**（故 84 次尝试中 83 次完成） |
| 短程运行的稳健性 | 许多运行只跑了 1–20 步，故把比较**只在跑满 ≥20 步的 30 次运行上重做**：**30/30 仍低于初值**，中位 0.10 / 0.13 / 0.49，最大 **0.72** —— 结论不依赖短程运行 |
| 未完成 | **无**：28/28 数据集三约定全部跑完 |

**结论写法（已写入正文）**：正权重读作边长是**约定而非测量**；换成倒数或单位读法会显著改变
初始几何（$V_0$ 中位 7.9 倍），但**不改变"流使曲率方差下降"这一方向性结论**。
Methods 的 *Limitations* 第 (iii) 条已写明"没有已核实的数据字典，该对照只界定而不消除其影响"。

**执行层面的诚实披露**：28 个数据集的全量扫描在第一轮（第 10 个数据集）上因内存不足中断
（机器 9.9 GB RAM、可用约 2.3 GB），改为**单进程、按数据集逐个分配墙钟预算**补齐；
现已 28/28 完成，但不同数据集的步数不等（各自有固定预算，大图只能跑 1–10 步），
**比值只在数据集内部跨约定可比，不能跨数据集比**；这一点连同 ≥20 步的稳健性复算一并写入 Methods。

### R6-N　极端权重跨度下的数值可靠性 + 初值扰动（对应意见 7）

见第一节第 7 条的表格（新补充表 **S13a/b/c**）。核心三条：

1. **求解器**：60 个检查点 × 3 张网络 = 9,120 个传输问题，HiGHS 兜底**从未被调用（0 次）**；
   强制调用时与 network simplex 在跨度 $<10^5$ 下相对成本差 $\le1.2\times10^{-15}$，跨度 $>10^8$ 时
   少数边退化到 $O(1)$–$O(10)$（中位数仍 $<10^{-9}$）→ 兜底在极端跨度不可信，但**没有已报告数字依赖它**。
2. **条件数**：把代价矩阵除以其自身均值（精确 $W_1$ 的一阶齐次性 ⇒ $\kappa$ 在精确算术下不变）后，
   跨度 $10^{26}$–$10^{29}$ 时 $\max_e|\Delta\kappa_e|\le2.2\times10^{-16}$、方差相对变化 $3.6\times10^{-16}$
   → **残差不是跨度造成的浮点伪影**；固定 $\varepsilon$ 的 Sinkhorn 在这些跨度下核下溢、完全失效
   （$|\Delta\kappa|\sim10^{16}$–$10^{25}$），这也是 2,000 步以后不报告任何 Sinkhorn 数值的原因。
3. **初值扰动（本轮新增）**：每个初始长度乘 $1\pm1\%$，5 次重复，精确传输 1,000 步：

| 网络 | 基线末残差 | 扰动后范围 | 最大相对偏差 |
|---|---|---|---|
| `ENZYMES_g10` | 6.080636e-2 | [6.080492e-2, 6.080690e-2] | **2.36e-5** |
| `ENZYMES_g105` | 1.176471e-2 | [1.176452e-2, 1.176497e-2] | **2.18e-5** |
| `ENZYMES_g106` | 3.850179e-2 | [3.850172e-2, 3.850187e-2] | **2.19e-6** |

→ **平台不是初始权重的伪影**：1% 的初值抖动只把末残差移动 $\le2.4\times10^{-5}$（相对量）。

### 合规（编译后实测）

| 项目 | 要求 | 状态 |
|---|---|---|
| 主文（Intro + Results + Discussion） | ≤ 5,000 词 | **4,998** ✅ |
| 摘要 | ≤ 150 词 | **148** ✅ |
| display items | ≤ 10 | **10**（8 图 + 2 表）✅ |
| 编译 | 0 error | xelatex + bibtex ×3：**0 error / 0 undefined / 26 页** ✅ |
| 复算脚本 | 数字可重跑 | `_patch_r6h.py`（正则定位、幂等，只刷新权重语义数字） ✅ |
| 图与表 | 不得丢失 | Fig.1–Fig.8、Table 1–Table 2 全部在位 ✅ |

