# -*- coding: utf-8 -*-
"""Consolidated evidence ledger for paper A, after L1-L3.

One page that a referee could read to decide whether the claims are supported.
Every number here is read from a file produced by the analysis scripts, or is
a verbatim constant from those reports -- nothing is retyped by hand where it
can be avoided.
"""

import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = r"D:\2026codes\outputs"
OUT = os.path.join(OUT_DIR, "论文A_最终证据链.html")

CSS = """
body{font-family:-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;
background:#f7f8fa;color:#1b1f24;margin:0;padding:32px;line-height:1.65}
.wrap{max-width:1080px;margin:0 auto}
h1{font-size:26px;margin:0 0 6px}
h2{font-size:19px;margin:34px 0 12px;padding-bottom:6px;
border-bottom:2px solid #e3e6ea}
h3{font-size:15px;margin:20px 0 8px;color:#37414d}
p{margin:8px 0}
.sub{color:#6a7480;font-size:13px}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13.5px;
background:#fff}
th,td{border:1px solid #dfe3e8;padding:7px 9px;text-align:left}
th{background:#eef1f5;font-weight:600}
td.num{text-align:right;font-variant-numeric:tabular-nums}
.card{background:#fff;border:1px solid #dfe3e8;border-radius:8px;
padding:16px 18px;margin:14px 0}
.kpis{display:flex;gap:12px;flex-wrap:wrap;margin:14px 0}
.kpi{flex:1;min-width:180px;background:#fff;border:1px solid #dfe3e8;
border-radius:8px;padding:14px}
.kpi .v{font-size:24px;font-weight:700;color:#1b4d8f}
.kpi .l{font-size:12px;color:#6a7480;margin-top:4px}
.ok{color:#1a7f37;font-weight:600}
.bad{color:#b42318;font-weight:600}
.warn{color:#b54708;font-weight:600}
.tag{display:inline-block;padding:1px 7px;border-radius:4px;font-size:12px;
font-weight:600}
.t-ok{background:#e6f4ea;color:#1a7f37}
.t-bad{background:#fdecea;color:#b42318}
.t-warn{background:#fef3e2;color:#b54708}
.quote{border-left:3px solid #1b4d8f;background:#f2f6fc;padding:10px 14px;
margin:12px 0;font-size:14px}
ul{margin:8px 0;padding-left:22px}
li{margin:4px 0}
code{background:#eef1f5;padding:1px 5px;border-radius:3px;font-size:12.5px}
"""


def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def main():
    H = []
    A = H.append
    A("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
    A("<title>论文 A 最终证据链</title>")
    A("<style>{}</style></head><body><div class='wrap'>".format(CSS))

    A("<h1>论文 A · 最终证据链</h1>")
    A("<p class='sub'>离散 Ollivier–Ricci 流的曲率均质化：从「观察到均质化」"
      "到「残余由什么决定」。L1–L3 全部执行完毕后的合并审计。"
      "</p>")

    # ---- KPI
    A("<div class='kpis'>")
    for v, l in (("31×", "固定度序列下，传递性 0.12→0.46 带来的饱和平台上升"
                         "（3/3 基图，中位）"),
                 ("74%", "该效应解释的合成图–真实网络差距（log₁₀ 口径）"),
                 ("22/22", "sham 对照中方差不降反升的配对数（p=2.4e-7）"),
                 ("6.4×", "加 surgery 后残余中位下降倍数")):
        A("<div class='kpi'><div class='v'>{}</div><div class='l'>{}</div>"
          "</div>".format(v, esc(l)))
    A("</div>")

    # ---- 1. 主张
    A("<h2>1. 现在可以写成什么主张</h2>")
    A("<div class='quote'>拓扑保持的离散 Ollivier–Ricci 流只能<b>部分</b>均质化"
      "曲率：流收敛到一个由网络结构决定的平台。在度序列固定的前提下，"
      "该平台主要由<b>局部聚类（传递性）</b>设定——把传递性从 0.12 提到 0.46"
      "（真实网络中位数）使饱和平台上升约 31 倍，解释了随机图与实证网络之间"
      "残余差距的约 74%。均质化需要真实的曲率—边映射（sham 对照成立），"
      "而残余的大部分可通过 surgery 移除；因此该平台度量的是"
      "「在不改变边集的前提下离常数曲率还有多远」。</div>")
    A("<p class='sub'>与 Ma–Yang 2025（带 surgery 的流收敛到常数曲率）"
      "互补而非重复：他们删边，我们量化不删边要付多少代价。</p>")

    # ---- 2. 证据表
    A("<h2>2. 证据账本</h2>")
    A("<table><tr><th>主张</th><th>判据</th><th>结果</th><th>状态</th></tr>")
    rows = [
        ("均质化不是更新规则的机械副产品",
         "sham 流：每步随机置换曲率—边的对应，其余一字不改，两臂步数配对",
         "22/22 可辩护配对中 sham 更高（符号检验 p=2.4e-7）；"
         "家族层面 8/8（p=0.0039）；real 中位 0.103 vs sham 1.159",
         "ok"),
        ("均质化是部分的，不是「还没跑够」",
         "长程跑到饱和",
         "原 10 步的 0.781 是瞬态；饱和平台中位 0.204，中位 130 步",
         "ok"),
        ("残余主要是「没删边」造成的",
         "surgery 配对（tol=0.10/0.30）",
         "删边 ≥6 且跑完的 3/3 全部下降：6.1× / 131.8× / 6.6×；"
         "符号检验 p=0.062（n=4，靠效应量不靠 p）",
         "warn"),
        ("同配性因果地抬高残余",
         "实验 B：固定度序列，双边交换推向目标 r；基图固定效应",
         "Spearman +0.796 (p=3.0e-14)；基图固定效应 β=+0.341 "
         "(p=2.7e-10)，传递性 p=0.765 消失",
         "ok"),
        ("社区结构抬高残余（幅度很小）",
         "实验 C：固定度序列，推向目标混合参数 μ；4 基图 × 3 水平",
         "4/4 基图内严格单调，ρ=+0.902 (p=6e-5)；"
         "<b>但跑到饱和后仅 1.35×</b>——匹配时程比较高估了效应",
         "bad"),
        ("<b>传递性是主导因子</b>",
         "实验 D：固定度序列，直接操控传递性 0.12→0.46",
         "3/3 基图上升，25.4× / 30.9× / 37.0×（中位 30.9×）；"
         "400 步检验在 176 步饱和于 0.0525，与 170 步的 0.0530 一致"
         "——<b>不是收敛速度伪迹</b>",
         "ok"),
        ("该效应不是同配性漂移",
         "用实验 B 标定的 β=+0.341 归因",
         "同配性漂移 Δr≈+0.12，只能预测 1.10×；"
         "占观测效应的 2.5%–3.4%",
         "ok"),
        ("桥边构成解析层面的障碍",
         "桥边上界引理（对一切正边权成立）+ 数值枚举",
         "引理成立（9 条桥中 5 条上界为负）；"
         "但桥密度不预测平台（ρ=−0.063）",
         "warn"),
        ("边瓶颈异质性假说",
         "合并真实+合成 97 图，边介数 CV",
         "合并 ρ=+0.672，真实/合成相差 3.3×；"
         "<b>但真实网络内部 ρ=+0.115 (ns)</b> → 假说被否",
         "bad"),
        ("跨家族可预测",
         "留一家族交叉验证（25 网络 / 7 家族）",
         "基线 oos_R²=−0.169；传递性 −0.330、同配性 +0.039、"
         "四变量 −0.106。<b>无一跑赢基线</b>",
         "bad"),
        ("家族内部可预测",
         "ENZYMES 内部（n=11，唯一 n≥5 的家族）",
         "传递性 ρ=+0.782 (p=0.0045，Bonferroni 0.031)，"
         "jackknife +0.709..+0.879；<b>但控制规模后偏相关 +0.281 "
         "(p=0.40)</b>，与规模（ρ=−0.879）不可分离",
         "warn"),
    ]
    for claim, test, result, state in rows:
        tag = {"ok": ("t-ok", "成立"), "bad": ("t-bad", "否定"),
               "warn": ("t-warn", "受限")}[state]
        A("<tr><td>{}</td><td class='sub'>{}</td><td>{}</td>"
          "<td><span class='tag {}'>{}</span></td></tr>".format(
              claim, esc(test), result, tag[0], tag[1]))
    A("</table>")

    # ---- 3. 必须写进论文的负面结果
    A("<h2>3. 必须主动写进论文的负面结果</h2>")
    A("<div class='card'><ul>")
    for item in (
        "<b>社区结构不是答案。</b>实验 C 在匹配时程下给出 4.6× 的漂亮效应，"
        "但跑到饱和只剩 1.35×。若不披露，这是审稿人最容易抓的把柄。",
        "<b>同配性不能当预测因子。</b>家族层面 ρ=+0.857 (p=0.014) 是"
        "生态学相关；留一家族 oos_R²=+0.039 仅勉强超过基线，"
        "且 ENZYMES 内部方向相反（ρ=−0.345）。",
        "<b>跨家族外推不存在。</b>所有模型的 oos_R² 都不高于"
        "「用留出家族自己的均值」。结构性→残余的映射是家族特有的。",
        "<b>家族内传递性效应与规模共线。</b>ENZYMES 内控制规模后偏相关"
        "塌到 +0.281 (p=0.40)，规模本身 ρ=−0.879。"
        "只能说「方向一致、相互印证」，不能说「确认」。",
        "<b>仍有 3.4× 差距未解释。</b>把传递性匹配到真实网络中位数后，"
        "合成图残余 0.0587 vs 真实 0.1967。",
    ):
        A("<li>{}</li>".format(item))
    A("</ul></div>")

    # ---- 4. 剩余缺口
    A("<h2>4. 挡在更高层级前的缺口</h2>")
    A("<table><tr><th>缺口</th><th>为什么挡路</th><th>难度</th></tr>")
    for gap, why, diff in (
        ("残余方差的解析下界",
         "目前只有桥边的上界引理（单条边）。要变成方差下界，"
         "需要对非桥边曲率给出匹配的下界。这是从 PRE 档迈向 "
         "PRL/NatComm 档唯一真正的杠杆。",
         "高（纯数学）"),
        ("surgery 配对样本量",
         "目前只有 4 对可辩护配对，符号检验 p=0.062。补到 8–10 对"
         "才能让「残余主要是没删边造成的」这句话站得住。",
         "中（纯算力）"),
        ("实验 D 的独立基图数",
         "3 个完整配对基图。方向无异议，但补到 6 个能显著增强。",
         "中（纯算力）"),
        ("区分传递性与规模",
         "真实网络里两者共线。需要固定规模、只改聚类的真实数据子集，"
         "或在合成图上做规模梯度。",
         "中"),
    ):
        A("<tr><td>{}</td><td>{}</td><td>{}</td></tr>".format(
            gap, esc(why), esc(diff)))
    A("</table>")

    # ---- 5. 期刊档位
    A("<h2>5. 期刊档位判断</h2>")
    A("<div class='card'>")
    A("<p><span class='ok'>现实目标</span>："
      "<b>PRE / Communications Physics / Network Science / EPJ Data Science</b>"
      "——已具备「对照 + 因果操控 + 定量幅度」三件套。"
      "实验 D 给出的是一条定量规律（平台随传递性上升约 31 倍，"
      "解释 74% 的差距），而不是一条相关。</p>")
    A("<p><span class='warn'>不建议</span>：NatComm / PNAS / PRL。"
      "这类期刊要的是解析结果（下界或收敛率），目前只有一条单边的桥边引理，"
      "没有下界，也没有应用侧验证。</p>")
    A("</div>")

    # ---- 6. 文件索引
    A("<h2>6. 结果文件索引</h2>")
    A("<table><tr><th>内容</th><th>路径</th></tr>")
    files = [
        ("实验 D（传递性操控）报告",
         r"D:\2026codes\2026821100data\synthetic_results\CLUSTERING_REPORT.txt"),
        ("实验 D 合并数据",
         r"D:\2026codes\2026821100data\synthetic_results\clustering_merged.csv"),
        ("实验 C（社区）报告",
         r"D:\2026codes\2026821100data\synthetic_results\COMMUNITY_REPORT.txt"),
        ("预测性审计 2（传递性）",
         r"D:\2026codes\2026821100data\ricci_flow_longrun_results\analysis\PREDICTIVE_AUDIT2.txt"),
        ("家族内稳健性检验",
         r"D:\2026codes\2026821100data\ricci_flow_longrun_results\analysis\WITHIN_FAMILY_ROBUSTNESS.txt"),
        ("桥边引理数值验证",
         r"D:\2026codes\2026821100data\ricci_flow_longrun_results\analysis\BRIDGE_REPORT.txt"),
        ("瓶颈异质性检验",
         r"D:\2026codes\2026821100data\ricci_flow_longrun_results\analysis\BOTTLENECK_REPORT.txt"),
        ("surgery 配对报告",
         r"D:\2026codes\2026821100data\surgery_results\SURGERY_REPORT.txt"),
        ("sham 对照报告",
         r"D:\2026codes\2026821100data\sham_results\analysis\ROBUST_REPORT.txt"),
    ]
    for label, path in files:
        A("<tr><td>{}</td><td><code>{}</code></td></tr>".format(
            esc(label), esc(path)))
    A("</table>")

    A("</div></body></html>")

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(H))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
