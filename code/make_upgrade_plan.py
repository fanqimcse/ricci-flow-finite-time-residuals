# -*- coding: utf-8 -*-
"""Build the 'how do we raise the tier of paper A' strategy document."""

import os

BASE = r"D:\2026codes\2026821100data"
OUT = r"D:\2026codes\outputs\论文A_提级方案.html"

CSS = open(os.path.join(BASE, "_report_css.txt"), encoding="utf-8").read()

H = []
A = H.append

A("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
A("<meta name='viewport' content='width=device-width,initial-scale=1'>")
A("<title>论文 A 提级方案</title><style>{}</style></head><body>".format(CSS))
A("<div class='wrap'>")

A("<h1>论文 A 提级方案：从「观察」到「机制」还差什么</h1>")
A("<div class='sub'>离散 Ricci 流下的跨域曲率均质化 · 证据账本与行动排序 · "
  "2026-10-03</div>")

A("<div class='note'><b>一句话结论。</b>论文目前最硬的两块证据（sham 对照、"
  "合成图因果）已经足够支撑 <b>PRE / Communications Physics / EPJ Data "
  "Science / Network Science</b> 这一档；但本轮审计暴露出<b>两条会直接被"
  "审稿人抓住的软肋</b>——真实网络上的同配性效应是<b>家族层面的生态学相关</b>"
  "，样本外预测能力≈0；以及<b>92× 的量级差距没有解释</b>。"
  "想再往上一档，靠加数据没用，必须补「机制性解释」或「理论」。</div>")

A("<div class='kpis'>")
A("<div class='kpi ok'><div class='v'>22/22</div><div class='l'>sham 对照："
  "真实臂残余显著低于置换臂（p=2.4e-7）</div></div>")
A("<div class='kpi ok'><div class='v'>+0.34</div><div class='l'>合成图同配性"
  "固定效应 β（p=2.7e-10），度序列钉死</div></div>")
A("<div class='kpi warn'><div class='v'>+0.04</div><div class='l'>留一家族"
  "样本外 R²（同配性）——几乎为零</div></div>")
A("<div class='kpi bad'><div class='v'>92×</div><div class='l'>合成图与真实"
  "网络残余量级差距，未解释</div></div>")
A("<div class='kpi bad'><div class='v'>6.1×</div><div class='l'>surgery 后残余"
  "中位下降倍数（删边≥6 的 3/3 全降）</div></div>")
A("</div>")

A("<div class='warnbox'><b>本次最关键的新结果：surgery 把残余打掉了。</b><br>"
  "这是审稿人的第一问（Ma–Yang 2025 证明带 surgery 的流收敛到常数曲率），"
  "现在有了答案，而且是「是」："
  "<table style='margin-top:8px'>"
  "<tr><th>数据集</th><th>无 surgery</th><th>有 surgery</th><th>删边</th>"
  "<th>下降倍数</th></tr>"
  "<tr><td>ENZYMES_g1</td><td>0.1740</td><td>0.0284</td><td>6</td>"
  "<td>6.1×</td></tr>"
  "<tr><td>ENZYMES_g10</td><td>0.3603</td><td>0.0027</td><td>9</td>"
  "<td>131.8×</td></tr>"
  "<tr><td>ENZYMES_g105</td><td>0.2393</td><td>0.0361</td><td>14</td>"
  "<td>6.6×</td></tr>"
  "</table>"
  "删边 ≥6 的 3 个数据集<b>全部</b>下降，中位 6.4×。"
  "（只删了 1 条边的 reptilia-sg 反而上升 0.3×，属未完成的瞬时扰动；"
  "预算截断的 g104/g107 不计入。）<br><br>"
  "<b>这意味着主主张必须改写</b>：不能再把平台说成「Ricci 流的结构性障碍」，"
  "而要说成「<i>拓扑保持</i>的 Ricci 流的可均质化上限——"
  "即一个图在不改变边集的前提下离常数曲率有多远」。"
  "好处是这个量 surgery 类方法会摧毁，因而是互补而非重复；"
  "坏处是 novelty 从「发现新障碍」降为「量化已知收敛理论的代价」。</div>")

# ---------- 1. evidence ledger
A("<h2>1. 证据账本：哪些能写，哪些不能</h2>")
A("<table>")
A("<tr><th>结论</th><th>证据强度</th><th>n</th><th>能否写进论文</th>"
  "<th>必须附的限制</th></tr>")
rows = [
    ("均质化是<b>部分的</b>：跑到饱和后残余≈初始的 20%",
     "<span class='tag ok'>强</span>", "93 网络",
     "可，作为主结果", "必须说明原 10 步只到瞬态（0.781→0.204）"),
    ("残余<b>需要真实的曲率—边映射</b>：置换后方差不降反升",
     "<span class='tag ok'>强</span>", "22 配对 / 8 家族",
     "可，作为核心 novelty", "两臂步数严格配对；剔除发散与惰性臂"),
    ("同配性<b>因果性</b>地抬高残余（度序列固定）",
     "<span class='tag ok'>强</span>", "60 合成图",
     "可，但限定在合成总体", "量级不可外推（差 92×）"),
    ("同配性在<b>真实网络家族层面</b>与残余相关",
     "<span class='tag warn'>中</span>", "7 家族",
     "可，但只能说「家族层面倾向」",
     "ρ=+0.857 p=0.014，但仅 7 点；4 个家族 n≤2"),
    ("同配性可<b>预测</b>单个网络的残余",
     "<span class='tag bad'>否</span>", "24 网络",
     "<b>不可</b>——样本外 R²=+0.04",
     "必须把这条零结果写进论文"),
    ("家族<b>内部</b>同配性越高残余越高",
     "<span class='tag bad'>否</span>", "ENZYMES n=11",
     "<b>不可</b>——ρ=−0.345（ns）", "唯一可估家族内为负"),
    ("社区结构（模块度）解释量级差距",
     "<span class='tag warn'>未定</span>", "—", "待实验 C",
     "家族层面 ρ=+0.786 p=0.036，家族内中位 −0.007"),
    ("残余是「没做 surgery」造成的",
     "<span class='tag bad'>是</span>", "6 配对 / 4 完成",
     "<b>必须改写主主张</b>",
     "删边 ≥6 的 3/3 全部下降 6.1–131.8×"),
]
for r in rows:
    A("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>"
      .format(*r))
A("</table>")

# ---------- 2. the three new negatives
A("<h2>2. 本轮审计新发现的三条坏消息（必须自己先写出来）</h2>")

A("<div class='warnbox'><b>（1）同配性效应是生态学相关，不是预测因子。</b><br>"
  "在统一的已收敛样本（24 网络 / 7 家族）上："
  "<ul>"
  "<li>家族层面：Spearman(平均同配性, 平均残余) = <b>+0.857, p=0.014</b>，"
  "留一家族 jackknife 稳定在 +0.771…+0.943 —— 这条是真的。</li>"
  "<li>家族内部：唯一可估的家族 ENZYMES（n=11）ρ = <b>−0.345（p=0.298）</b>，方向相反。</li>"
  "<li>留一家族外推：同配性模型样本外 <b>R² = +0.037（p=0.11）</b>，"
  "模块度 −0.144，四变量模型 −1.725。没有任何模型跑赢「用留出家族自己的均值」。</li>"
  "</ul>"
  "结论：可以写「<i>在生成源家族的层面上</i>，同配性越高的网络家族残余越大」；"
  "<b>不能</b>写「同配性预测残余」。把这条零结果写进去反而加分——它显示你知道"
  "生态学相关和网络级预测的区别。</div>")

A("<div class='warnbox'><b>（2）样本口径一换，效应就变。</b><br>"
  "同一批数据，三种口径给出三个数：已收敛子集 ρ=+0.857（n=7 家族）、"
  "全部有比值网络 ρ=+0.291（n=13 家族）、LOFO 加状态过滤 ρ 不显著（n=25）。"
  "这正是弱效应被误当强效应的典型路径。<b>论文必须锁定一个口径并全文不换</b>，"
  "其余口径放附录做敏感性分析。</div>")

A("<div class='warnbox'><b>（3）92× 量级差距仍然没有解释。</b><br>"
  "合成 configuration-model 图中位残余 0.00214，真实网络 0.197。"
  "这说明 configuration model 抓不到的结构（社区、层级、几何嵌入）才是量级的"
  "主控因素，同配性只决定<b>序</b>。好消息：家族层面模块度 ρ=+0.786（p=0.036）"
  "比同配性还强一点；坏消息：家族内部模块度中位 ρ = −0.007，同样是生态学相关。</div>")

# ---------- 3. what caps the tier
A("<h2>3. 决定天花板的是什么</h2>")
A("<table>")
A("<tr><th>档次</th><th>还需要什么</th><th>当前状态</th></tr>")
A("<tr><td>J. Complex Networks /<br>Applied Network Science /<br>"
  "JPhys Complexity / Chaos</td>"
  "<td>已完成（sham 对照 + 饱和长程）</td>"
  "<td><span class='tag ok'>已够</span></td></tr>")
A("<tr><td>PRE / EPJ Data Science /<br>Communications Physics /<br>"
  "Network Science</td>"
  "<td>P0 对照 + 合成图机制验证</td>"
  "<td><span class='tag ok'>已够</span>（本轮完成）</td></tr>")
A("<tr><td>Nature Communications /<br>PNAS / PRL</td>"
  "<td>解析结果（下界 / 收敛率 / 相变）+ 一个应用端验证</td>"
  "<td><span class='tag bad'>差距大</span>：目前无任何定理，"
  "且预测性已被证伪</td></tr>")
A("</table>")
A("<p class='sub'>说明：PRL 与 NatComm 的分野在于「是否改变了人们对这类系统的"
  "理解」。目前这项工作的形状是「一个干净的实验 + 一个受控的因果验证」，"
  "这正好是 PRE/Comm Phys 的标准形态，不是 NatComm 的形态。</p>")

# ---------- 4. ranked actions
A("<h2>4. 行动清单（按「能不能改变论文的档次」排序）</h2>")

A("<h3>P0 — 不做就无法定稿</h3>")
A("<table>")
A("<tr><th>#</th><th>行动</th><th>为什么</th><th>成本</th></tr>")
A("<tr><td>1</td><td><b>改写主主张</b>：残余 = 拓扑保持流的可均质化上限</td>"
  "<td>surgery 实验已答（删边≥6 的 3/3 下降 6.1–131.8×）。"
  "继续声称「结构性障碍」会被 Ma–Yang 2025 直接驳回。</td>"
  "<td>1 d（写作）</td></tr>")
A("<tr><td>2</td><td><b>锁定单一分析口径</b>，其余进附录</td>"
  "<td>见第 2 节（2）。口径漂移是当前最大的可信度风险。</td><td>1 h</td></tr>")
A("<tr><td>3</td><td><b>把零结果写进论文</b>（LOFO R²≈0、家族内 ρ&lt;0）</td>"
  "<td>主动报告比被审稿人挖出来好；也是「我们做了真正的样本外检验」的凭证。"
  "</td><td>2 h</td></tr>")
A("<tr><td>4</td><td>把 surgery 配对数从 3 提到 8–10</td>"
  "<td>目前 3 个完成配对，符号检验 p=0.062 不足以独立支撑；"
  "需要更多「删边≥3 且跑完」的配对。</td><td>~3 h 算力</td></tr>")
A("</table>")

A("<h3>P1 — 决定能否再上一档</h3>")
A("<table>")
A("<tr><th>#</th><th>行动</th><th>为什么</th><th>成本</th></tr>")
A("<tr><td>4</td><td><b>实验 C：在合成图上加社区结构</b><br>"
  "固定度序列，用度保持双边交换推向目标混合参数 μ，"
  "测残余是否从 0.002 升到 ~0.2</td>"
  "<td>唯一能解释 92× 差距的实验。若成立，论文就有一条"
  "「残余 = f(同配性, 社区强度)」的双因子定律，"
  "从「相关」升级为「机制性定量规律」——这是提级的真正杠杆。</td>"
  "<td>~3 h 算力</td></tr>")
A("<tr><td>5</td><td><b>解析性下限</b>：对一类图（如完全二部、"
  "规则树、双社区）推导残余方差的闭式或下界</td>"
  "<td>纯数学，但一旦有一个可检验的闭式，档次立刻不同。"
  "建议从「双社区图」切入，因为实验 C 会给它数值校准。</td>"
  "<td>不确定（人月级）</td></tr>")
A("<tr><td>6</td><td><b>把残余做成可用量</b>：证明它能区分生成源家族"
  "（ICC=0.51）优于现成指标</td>"
  "<td>应用端。若成立可加一节「Ricci 残余作为网络指纹」。</td>"
  "<td>1 d</td></tr>")
A("</table>")

A("<h3>P2 — 加分但不改档次</h3>")
A("<table>")
A("<tr><th>#</th><th>行动</th><th>说明</th></tr>")
A("<tr><td>7</td><td>补第二枚 sham 种子</td>"
  "<td>目前只跑了 seed 1；第二枚种子把 p 从 2.4e-7 变成可复现的稳健结论。"
  "</td></tr>")
A("<tr><td>8</td><td>把 45 个未收敛数据集的截断写清楚</td>"
  "<td>600 s 预算截断主要打在大图和多分量图上，必须声明。</td></tr>")
A("<tr><td>9</td><td>统一 attack/collapse 定义（为 Track B 预留）</td>"
  "<td>暂时不影响论文 A。</td></tr>")
A("</table>")

# ---------- 5. do not
A("<h2>5. 不要做的事</h2>")
A("<ul>")
A("<li><b>不要把 κ₀ 当预测因子。</b>它与分母纠缠（ρ(κ₀, 初始方差)=−0.477），"
  "绝对值尺度上消失（−0.140, p=0.38），OLS 中不显著（p=0.394）。已撤回。</li>")
A("<li><b>不要把 110 个网络当 110 个独立样本。</b>其实只有约 16 个生成源家族，"
  "ENZYMES 一家就占 12–16 个。</li>")
A("<li><b>不要外推合成图的效应量（92×）。</b>只外推方向与因果归属。</li>")
A("<li><b>不要声称「方差普遍下降」。</b>主主张必须是「部分均质化 + 结构决定的"
  "残余平台」。</li>")
A("<li><b>不要在论文里出现原始隐喻</b>，只在致谢/引言作为灵感来源提一句。</li>")
A("</ul>")

A("<h2>6. 执行顺序建议</h2>")
A("<pre>现在  → surgery 跑完 → 定稿主主张（是否改写为「拓扑保持流的障碍」）\n"
  "  ↓\n"
  "锁定口径 → 把零结果写进论文 → 投稿 J. Complex Networks 档（保底）\n"
  "  ↓\n"
  "实验 C（社区强度）\n"
  "  ├── 成功：双因子定律 → 改投 PRE / Comm Phys / Network Science\n"
  "  └── 失败：差距归因为「未知结构」→ 维持上一档，明确列为 open problem\n"
  "  ↓\n"
  "解析下界（若有人力）→ 才谈 NatComm / PRL</pre>")

A("<p class='sub'>本文件由 make_upgrade_plan.py 生成；所有数字来自 "
  "<code>ricci_flow_longrun_results/analysis/FINAL_PREDICTIVE_AUDIT.txt</code>、"
  "<code>sham_results/analysis/ROBUST_REPORT.txt</code>、"
  "<code>synthetic_results/analysis/SYNTH_REPORT.txt</code>。</p>")

A("</div></body></html>")

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(H))
print("written:", OUT, "| bytes:", os.path.getsize(OUT))
