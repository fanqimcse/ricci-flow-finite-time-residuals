# -*- coding: utf-8 -*-
"""Report the completed L1 -> L2 -> L3 sequence."""

import os

BASE = r"D:\2026codes\2026821100data"
OUT = r"D:\2026codes\outputs\论文A_L1L2L3执行报告.html"
CSS = open(os.path.join(BASE, "_report_css.txt"), encoding="utf-8").read()

H = []
A = H.append
A("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
A("<meta name='viewport' content='width=device-width,initial-scale=1'>")
A("<title>论文 A · L1–L3 执行报告</title><style>{}</style></head><body>"
  .format(CSS))
A("<div class='wrap'>")
A("<h1>论文 A · L1 → L2 → L3 执行报告</h1>")
A("<div class='sub'>2026-10-03 · 全部数字可追溯到 "
  "<code>surgery_results/</code>、<code>synthetic_results/</code>、"
  "<code>ricci_flow_longrun_results/analysis/</code></div>")

A("<div class='note'><b>三条路线的净结果。</b>"
  "L1（surgery）把主主张逼着改写；"
  "L2（社区结构）给出了本次唯一的<b>新机制</b>，并把 92× 量级差距压到 5.7×；"
  "L3（定量定律 + 样本外）给出一条受控总体内成立的标度关系，同时再次确认"
  "<b>在真实网络上什么都预测不了</b>。"
  "论文现在可以声称：一个受控的因果机制 + 一个有对照的主结果；"
  "不能声称：预测。</div>")

A("<div class='kpis'>")
A("<div class='kpi bad'><div class='v'>5/5</div><div class='l'>L1：surgery "
  "后残余全部下降，符号检验 p=0.031，中位 6.6×</div></div>")
A("<div class='kpi ok'><div class='v'>3/3</div><div class='l'>L2：三个基图"
  "残余随社区强度单调上升，P=0.0046</div></div>")
A("<div class='kpi ok'><div class='v'>5.7×</div><div class='l'>L2：量级差距"
  "从 92× 收窄到 5.7×（下界）</div></div>")
A("<div class='kpi bad'><div class='v'>+0.04</div><div class='l'>L3：真实网络"
  "留一家族样本外 R²，仍≈0</div></div>")
A("</div>")

# ---- L1
A("<h2>L1 · surgery 检查 —— 答案是「是」，主主张必须改写</h2>")
A("<p>审稿人的第一问：Ma–Yang 2025 证明带 surgery 的流收敛到常数曲率，"
  "你的残余是不是只是没删边？现在有答案了。</p>")
A("<table><tr><th>数据集</th><th>无 surgery</th><th>有 surgery</th>"
  "<th>删边</th><th>下降</th></tr>")
for name, a_, b_, rm in (("ENZYMES_g10", "0.3603", "0.0027", 9),
                         ("ENZYMES_g111", "0.2236", "0.0023", 11),
                         ("ENZYMES_g105", "0.2393", "0.0361", 14),
                         ("ENZYMES_g1", "0.1740", "0.0284", 6),
                         ("scc_rt_damascus", "0.2510", "0.2382", 4)):
    A("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td><b>{:.1f}×</b>"
      "</td></tr>".format(name, a_, b_, rm, float(a_) / float(b_)))
A("</table>")
A("<p class='sub'>删边 ≥3 且手术臂跑完的 <b>5/5</b> 全部下降，"
  "单侧精确符号检验 <b>p = 0.0312</b>，中位 6.6×（范围 1.1×–131.8×）。"
  "被排除的：从未触发删边的（与无 surgery 臂完全相同）、只删 1–2 条边的"
  "（属扰动，reptilia-sg 反而上升 0.3×）、以及被墙钟预算截断的 6 个"
  "（g104/g107/g109/g112/g113/mammalia）。</p>")
A("<div class='warnbox'><b>主主张改写。</b>不能再说「Ricci 流的结构性障碍」，"
  "要说「<b>拓扑保持</b>流的均质化上限——一个图在不改变边集的前提下离常数"
  "曲率还有多远」。这个量 surgery 类方法会摧毁，所以与 Ma–Yang 2025 "
  "<i>互补</i>而非冲突；代价是 novelty 从「发现新障碍」降为"
  "「量化已知收敛理论的代价」。</div>")

# ---- L2
A("<h2>L2 · 实验 C——社区结构才是量级的主控因素</h2>")
A("<p>固定度序列的 configuration model 基图，只把连线推向目标跨社区边比例 μ"
  "（4 个等大组，μ=0.75 即无结构）。度序列一字不变。</p>")
A("<table><tr><th>基图</th><th>μ=0.75（无社区）</th><th>μ=0.45</th>"
  "<th>μ=0.15（强社区）</th><th>趋势</th></tr>")
A("<tr><td>base 0</td><td>0.00705</td><td>0.01308</td><td>0.01670</td>"
  "<td>递增</td></tr>")
A("<tr><td>base 1</td><td>0.00860</td><td>0.01376</td><td>—</td>"
  "<td>递增</td></tr>")
A("<tr><td>base 2</td><td>0.00645</td><td>0.00918</td><td>0.05277</td>"
  "<td>递增</td></tr>")
A("</table>")
A("<ul>")
A("<li>Spearman(log 残余, Q) = <b>+0.905, p=0.002</b>；"
  "Spearman(log 残余, μ) = −0.945, p=0.0004</li>")
A("<li><b>3/3 基图单调上升</b>，随机下概率 (1/6)³ = <b>0.0046</b></li>")
A("<li><b>同配性 ρ = +0.000 (p=1.0)</b> —— 社区效应与同配性效应是"
  "两个完全不同的机制，不是同一个东西的两面</li>")
A("<li>量级：最弱 0.00705 → 最强 0.03474（4.9×）；对真实网络中位 0.19668 的"
  "差距从 <b>92× 收窄到 5.7×</b></li>")
A("<li>边数无关（ρ=−0.056, p=0.90）</li>")
A("</ul>")
A("<div class='warnbox'><b>必须随每个数字一起出现的限制。</b>"
  "8 张图全部停在<b>匹配时程 95 步</b>（status = horizon_reached，"
  "<b>未饱和</b>）。所以严格成立的陈述是「在相同流时间下，社区结构让残余"
  "更高且下降更慢」，而<b>不是</b>「渐近平台更高」。合成值若继续跑还会降，"
  "因此 5.7× 只是真实差距的<b>下界</b>。</div>")

# ---- L3
A("<h2>L3 · 定量定律 + 样本外验证</h2>")
A("<h3>Part 1：受控总体内的标度关系</h3>")
A("<pre>log10(残余) = -2.164 + 1.135 * Q      R2 = 0.778,  p(Q) = 0.0038\n"
  "  加基图固定效应：beta(Q) = +1.177 (p = 0.016), R2 = 0.808\n"
  "  即：植入模块度每增加 0.1，残余乘以 1.30×</pre>")
A("<h3>Part 2：同一对因子放到真实网络上</h3>")
A("<table><tr><th>层面</th><th>同配性</th><th>模块度</th><th>两者联合</th></tr>")
A("<tr><td>家族层面（n=7 家族）</td><td>+0.857 (p=0.014)</td>"
  "<td>+0.786 (p=0.036)</td><td>—</td></tr>")
A("<tr><td>留一家族样本外 R²</td><td><b>+0.037</b></td><td><b>−0.144</b></td>"
  "<td><b>−0.159</b></td></tr>")
A("</table>")
A("<div class='warnbox'><b>L3 的净结论是分裂的，而且两边都要写。</b><br>"
  "在受控总体里，社区强度→残余是一条<b>可定量的因果律</b>（固定度序列、"
  "基图固定效应下 β=+1.18）；<br>"
  "在真实网络语料里，<b>没有任何模型</b>的样本外 R² 跑赢「用留出家族自己的"
  "均值」。所以论文可以声称受控因果机制与家族层面倾向，"
  "<b>永远不能声称预测</b>。</div>")

# ---- what the paper can say
A("<h2>现在论文能写什么 / 不能写什么</h2>")
A("<table><tr><th>可写</th><th>不可写</th></tr>")
A("<tr><td>拓扑保持 Ricci 流产生<b>部分</b>均质化并停在平台"
  "（长程跑到饱和，93 网络）</td>"
  "<td>「Ricci 流的结构性障碍」（会被 surgery 结果驳回）</td></tr>")
A("<tr><td>均质化<b>需要真实的曲率—边映射</b>（sham 22/22，p=2.4e-7）</td>"
  "<td>「方差普遍下降」（瞬态与平台必须区分）</td></tr>")
A("<tr><td>同配性在<b>固定度序列</b>下因果地抬高残余（β=+0.341, p=2.7e-10）</td>"
  "<td>合成图的效应<b>量级</b>外推到真实网络</td></tr>")
A("<tr><td>社区强度独立地抬高残余（ρ(Q)=+0.905；同配性 ρ=+0.000）</td>"
  "<td>「同配性/模块度预测残余」（LOFO R²≈0）</td></tr>")
A("<tr><td>家族层面的同配性倾向（ρ=+0.857, p=0.014, jackknife 稳定）</td>"
  "<td>把 110 个网络当 110 个独立样本</td></tr>")
A("</table>")

A("<h2>还差什么才能再上一档</h2>")
A("<ol>")
A("<li><b>把实验 C 跑到饱和</b>——目前只到匹配时程，是最大的方法学缺口。"
  "需要 ~8 张图各跑到饱和（约 3–4 h 算力，本机只能前台分批跑）。</li>")
A("<li><b>解析结果</b>：残余方差的下界或收敛率。这是从 PRE/Comm Phys 档"
  "迈向 NatComm/PRL 的唯一通道，目前完全没有。</li>")
A("<li><b>surgery 配对从 5 提到 8–10</b>，并把「删边后残余如何随删边比例变化」"
  "做成一条曲线而不是两点比较。</li>")
A("<li>锁定单一分析口径（已收敛子集），其余口径进附录。</li>")
A("</ol>")

A("<p class='sub'>生成：make_l123_report.py · "
  "数据源：surgery_results/SURGERY_REPORT_v2.txt、"
  "synthetic_results/COMMUNITY_REPORT.txt、COMMUNITY_MONOTONICITY.txt、"
  "L3_REPORT.txt、ricci_flow_longrun_results/analysis/"
  "FINAL_PREDICTIVE_AUDIT.txt</p>")
A("</div></body></html>")

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\n".join(H))
print("written:", OUT, os.path.getsize(OUT))
