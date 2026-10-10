# 对《Communications Physics 投稿的 3 次审查》的逐条回应（第七轮）

审查日期：2026-10-11 ｜ 审查对象：`manuscript(5).pdf`（26 页）
本轮结论：**阻断项（第 1 条）已从原始轨迹逐位重算并解决**；其余 8 条全部落地。
编译状态：**0 error / 0 undefined reference / 27 页 / 主文 5,000 词 / 摘要 150 词 / 10 个 display items / 45 条参考文献**。

---

## 第 1 条（阻断）：停止规则与 16 个严格重跑的终止窗口统计互相矛盾

**审查的意见**：Methods 定义 $|V_t-V_{t-1}|/V_t<10^{-5}$ 连续 50 步；但正文又报告最后 50 步的局部率中位数 $4.7\times10^{-3}$、14/16 超过 $10^{-4}$、窗口内方差变化 21%。按 $e^{r_t}-1$ 的换算，这两组数字不可能同时成立；"至少有一项统计定义、数据汇总、终止状态标记或文字报告需要纠正"。

**处理：直接读取逐步轨迹重算，不做任何推断。** 新脚本 `_r7_stop.py` 从 `ricci_flow_v2/per_dataset/*.csv`（即 Supplementary Table S10 的同一批逐步记录；配置 `tol=1e-5`、`tol_mode=current`、`patience=50`、`max_steps=400`、`budget=700s`）逐步重算 $|V_t-V_{t-1}|/V_t$ 与 $r_t=-\log(V_t/V_{t-1})$，生成 **Supplementary Table S15**（新增）。

**重算结果（这是本轮最重要的一个数字）**：

| 量 | 值 |
|---|---|
| 判据真正触发（连续 ≥50 步）的运行数 | **1 / 16** |
| 触发的网络 | `scc2034_kilifi...`（2 个分量，池化方差 **100% 为分量间项**） |
| 该运行的细节 | 345 步中**恰好 50 步**低于容差，且**正好是最后 50 步**；运行在此结束 |
| 其余 15 个运行的终止方式 | 步数地平线（245–400 步）或 700 s 墙钟预算 |
| 其中**一步都没进入容差**的 | **11 / 15** |
| 其中间歇满足但从未连续 50 步的 | 4 / 15（最长连续段 1–33 步） |
| 让判据在这 15 个运行的终止窗口内触发所需的最小容差（中位） | **$5.1\times10^{-3}$/步**（范围 $4.2\times10^{-5}$–$1.7\times10^{-2}$），约为所配置 $10^{-5}$ 的 **500 倍** |
| 3 个 5,000 步精确传输运行的同一量 | $1.3\times10^{-7}$、$8.9\times10^{-11}$、$1.1\times10^{-13}$（**新表 S16**） |

**结论：审查者指出的矛盾是真实存在的，性质是"文字报告"错了，而不是统计量错了。** 原稿的写法会让读者以为这 16 个运行是被判据终止的；实际上它们几乎全部是被地平线/预算终止的，而正文报告的 $4.7\times10^{-3}$、14/16、21% 描述的正是**判据没有触发的那些运行**——两组数字并不冲突，此前是标签把它们放在了一起。

**稿件改动**：

1. Methods「Convergence」段整段重写：明确写出判据在 $C_4$ 上**只触发 1 次**、其余 15 个的终止方式、11/15 从未进入容差、所需容差的中位数与范围，以及与 3 个 5,000 步运行的对照；并明确"这些终止窗口统计不是与规则冲突，而是规则没有触发的那些运行的记录"。
2. **术语三分、不再混用**（审查要求的"统一定义"）：
   - *rule-terminated*：判据触发（1 个）；
   - *horizon-terminated / budget-terminated*：地平线或墙钟终止（15 个）；
   - *still falling at the horizon*：终止窗口中位 $r_t>10^{-4}$/步（14 个）。
   - 字符串 *near-plateau* **只**用于 Supplementary Table S16 中三个已验证满足判据的 5,000 步精确传输运行；不再从任何一个"跑到地平线就结束"的运行上推断出来。
3. 全稿删除 "operationally stopped" 作为 $C_4$ 的运行标签（保留 "operational stopping criterion" 作为规则名称）。
4. Results 段同步：写明判据"触发 1 次"、触发者是那个分量间项占 100% 的断连网络，并**明确标注这是"与该论证一致的一个观察，不是对该论证的演示"**（n=1）。
5. 图 2b 图注改为"re-run to the longer horizon"并加注"the corrected rule fired on one of them"；图 8a 的红/蓝改为按终止窗口中位率（Pattern II / Pattern I）定义，不再用"near-plateau"。
6. 摘要加入 "the corrected rule fired once in 16 re-runs"。

**未改动**：$4.7\times10^{-3}$、14/16、21%、58%、两个负率——这些数字全部按逐步数据重算复核，与 S10 一致。

---

## 第 2 条：修正后的停止规则不能普遍区分正极限与零极限

**处理：接受，结论已被限定。** 审查指出 $\lambda$ 可以小到 $1-e^{-\lambda}<10^{-5}$，此时即使真实极限为零，缓慢衰减也会满足判据，连续 50 步无法消除这一逻辑漏洞。

- 原句 "The corrected rule fires on runs settling onto a positive value and not on runs still decaying towards zero" **删除**，改为：
  > "The direction is the one we want, but the inference is not valid in general: whenever $\lambda$ is small enough that $1-e^{-\lambda}<10^{-5}$, a trajectory decaying all the way to zero also satisfies the test, and demanding 50 consecutive steps does not remove that. The corrected rule is a local diagnostic…"
- 并按审查建议，把规则定性为"局部诊断"，不声称能判定极限是否为正、也不声称轨迹已收敛。

**式（2）的下标偏差已修正**：分子 $A(1-e^{-\lambda})e^{-\lambda t}\to A(1-e^{-\lambda})e^{-\lambda(t-1)}$，与 $V_t-V_{t-1}=Ae^{-\lambda(t-1)}(1-e^{-\lambda})$ 一致。渐近论证不变。

---

## 第 3 条："三个网络达到平台"需与停止规则明确区分

**处理：接受，两件事在文中被显式拆开。**

- 摘要 "barely decays over 5,000 steps" → **"moves <0.2% from 2,000 to 5,000 steps"**（按审查建议改为对窗口的描述，而不是对"衰减"的描述）；同时保留"terminal value above the run minimum"这一事实。
- Results「Extended horizon」段改写为：三个运行满足描述性判据 $3{,}920$–$4{,}241$ 个连续步（占全部 5,000 步的 79–85%，新表 S16），**并明确"由此得到的只是关于我们观测到的 2,000–5,000 步窗口的陈述，与 $t\to\infty$ 无关，后者我们不讨论"**。
- 明写"它们是本文唯一满足该判据的运行"，杜绝与 $C_4$ 的 15 个地平线终止运行混读。

---

## 第 4 条：exact-transport 的 0.1 小时耗时需用原始日志复核

**处理：已复核，0.1 h 属实，并在文中补上计算方式。** 逐次运行的 `seconds` 字段：

| 阶段 | 步数 | 实测 |
|---|---|---|
| $C_5$ 九个网络 × 2,000 步 | 18,000 | **173.4 s** |
| 其中三个网络续跑至 5,000 步 | 15,000 | **213.5 s** |
| 合计 | 33,000 | **386.9 s = 0.107 h** |

**为什么会这么小**（原稿未解释，这是审查疑虑的真正来源）：$C_5$ 的九个网络只有 **5–33 个节点、6–69 条边**，单步成本 **5–21 ms**；而 3.3 s 的单步中位数是全 110 个数据集语料库的中位数，量级完全不同。这一解释已写入 Methods。

**同时纠正了同段两个复核后发现不一致的数字**（属同一句话组，一并修正以免留下新的可检验矛盾）：

| 项目 | 原稿 | 实测（逐步计时器求和） |
|---|---|---|
| 首轮 110 网络 × 10 步 | 13 h | **34 h** |
| 严格重跑（本文报告的 16 网络集） | 0.6 h | **1.5 h**（另有 0.5 h 的扩展地平线辅助集） |
| 步长扫描 | 0.1 h | **0.08 h** ✅ 一致 |
| ε 轨迹扫描 | 1 h | **1.05 h** ✅ 一致 |
| 34 个断连网络的分量内重跑 | 15 h | 未按运行计时，**已声明不报告** |

并在段首写明"每个数字都是该运行自己写下的逐步计时器之和"，消除"wall clock / CPU time / 子任务时间"的歧义。

---

## 第 5 条：极端权重跨度下的数值稳定性，结论继续保持克制

**处理：接受，并显式列出这些控制"没有"建立什么。** 在 Methods「Numerical reliability at extreme weight span」段末新增：

> 这些控制只界定了三种具体失效模式——求解器一致性、双精度下的条件数、对初始权重的敏感性。它们**没有**证明该残差在这些跨度下不含数值误差。极端权重比与传输子问题的条件数仍是核心限制；常曲率检查施加在两个参考图上，而非被研究的网络上；两个独立优化器没有找到更低方差的解，不等于全局最优已被证明。我们所主张的是：**已测试的数值控制无法解释该残差**。

原有的 "not an artefact of the span / of the initial weighting" 均为限定性表述（"of X"），保留。

---

## 第 6 条：权重的物理语义不能再扩展解释

**处理：接受，逐处对齐。**

- 摘要中的 "geometric infimum" 已与正文定义严格对齐：摘要写作 "the infimum $V_G^*$ (not computed)"，正文 Eq.(1) 定义后紧跟 "Nothing below is a lower bound on $V_G^*$"，Methods 前沿同样声明。**不写入任何"已计算"或"被实验下界约束"的读法。**
- Methods 中"把正值读成长度是约定而非测量"的表述保留（28/28 数据集无一被文档记为距离）；Results 与摘要保持一致口径。
- 摘要末句由 "structure does not predict the residual across families" 改为 **"the tested descriptors do not predict it across families"**（同时回应第 9 条）。

---

## 第 7 条：摘要第一句对离散 Ricci flow 的概括过宽

**处理：接受，按审查建议改写。**

> 原：*Discrete Ricci flow on a connected weighted graph converges to constant curvature; finite computations differ.*
> 现：**Certain discrete Ollivier–Ricci flows converge to constant curvature under assumptions; finite computations need a separate reading.**

与正文中"现有收敛结果只针对特定构造"以及 Table 2（Li & Münch 假设对照表）的边界完全一致。

---

## 第 8 条："Homogenisation itself is not in doubt" 过于绝对

**处理：接受，按审查建议改成关联表述。** 讨论段现为：

> The variance reduction we observe is **associated with** the curvature–edge correspondence **under the intervention we tested**: destroying it in the sham control is associated with suppressed homogenisation in 8/8 families. One intervention in one finite window is **evidence of association, not a proof that the correspondence is necessary**, nor that every network or flow scheme must homogenise.

引言中并列的 "Homogenisation itself survives every control and depends on the curvature–edge correspondence" 同步改为 "The variance reduction we observe survives every control and is **associated with** the curvature–edge correspondence **under the intervention we tested**"。

---

## 第 9 条：合成网络与跨家族预测的结论不要超过样本设计

**处理：接受，两处收窄。**

- 摘要："structure does not predict the residual across families" → **"the tested descriptors do not predict it across families"**（不写成"任何结构指标都无法预测"）。
- 讨论："transitivity moves the relaxation time and not the endpoint" → 加限定 **"in the synthetic design we tested"**。
- 合成实验的"两个基础图/探索性"标注、手术控制的"两个家族"标注、跨家族 LOFO 的负 $R^2$ 与"这些模型不迁移 ≠ 残差本质上家族特异"的区分，均保留未动。

---

## 投稿规范

| 项目 | 状态 |
|---|---|
| 主文词数（不含 Methods） | **5,000**（≤5,000） |
| 摘要 | **150** 词 |
| Display items | **10**（≤10） |
| 参考文献 | 45 条，编号按出现顺序，0 漏引 / 0 未被引 / 0 重复 |
| 独立章节 | Data availability ✅ / Code availability ✅ / Funding ✅ / Author contributions ✅ / Competing interests ✅ |
| 章节顺序 | Introduction–Results–Discussion–Methods ✅（该刊允许） |
| 补充表 | S1–S16（本轮新增 **S15 停止判据诊断**、**S16 描述性平台三网络**） |
| 编译 | xelatex + bibtex，**0 error / 0 undefined**，27 页 |

---

## 仍待作者完成（需 token，非本轮可完成）

1. 把 `github_release` 推送到远端并发 v1.0.1（本地已同步 S15/S16、`_r7_stop.py` 与最新稿件）。
2. 向 Zenodo 申请新版本 DOI（稿件中的 DOI 仍指向 v1.0.0，S15–S16 尚未进归档）。
3. Software Heritage 归档完成后补 SWHID。
