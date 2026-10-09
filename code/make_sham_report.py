"""Render the sham-flow control result as a standalone HTML report.

Reads sham_results/analysis/{03_robust_table.csv,04_family_table.csv,
ROBUST_REPORT.txt} and writes outputs/sham_control_report.html.
"""

import os

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BASE)
OUT_DIR = os.path.join(ROOT, "outputs")
ANALYSIS = os.path.join(BASE, "sham_results", "analysis")

CSS = """
:root{--bg:#f7f8fa;--card:#ffffff;--ink:#1a1d21;--muted:#5c6470;
--line:#e2e6ec;--accent:#2f6fed;--good:#1a7f4b;--bad:#b4342b;--warn:#9a6b00;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
font-family:-apple-system,"Segoe UI","Microsoft YaHei",Helvetica,Arial,sans-serif;
line-height:1.65;font-size:15px}
.wrap{max-width:1040px;margin:0 auto;padding:36px 24px 72px}
h1{font-size:27px;margin:0 0 6px;letter-spacing:-.2px}
h2{font-size:19px;margin:38px 0 12px;padding-bottom:7px;border-bottom:1px solid var(--line)}
h3{font-size:15.5px;margin:22px 0 8px;color:var(--muted)}
p{margin:10px 0}
.sub{color:var(--muted);font-size:13.5px;margin-bottom:20px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;
padding:18px 20px;margin:14px 0}
.kpis{display:flex;gap:14px;flex-wrap:wrap;margin:18px 0}
.kpi{flex:1 1 190px;background:var(--card);border:1px solid var(--line);
border-radius:10px;padding:14px 16px}
.kpi .v{font-size:25px;font-weight:650;letter-spacing:-.5px}
.kpi .l{font-size:12.5px;color:var(--muted);margin-top:3px}
.kpi.ok .v{color:var(--good)} .kpi.bad .v{color:var(--bad)}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:10px 0}
th,td{border-bottom:1px solid var(--line);padding:6px 9px;text-align:right;
font-variant-numeric:tabular-nums}
th:first-child,td:first-child{text-align:left}
th{background:#eef1f5;font-weight:600;color:var(--muted);font-size:12.5px}
tr:hover td{background:#f4f6f9}
.tag{display:inline-block;padding:1px 8px;border-radius:999px;font-size:12px;
font-weight:600}
.tag.ok{background:#e6f4ec;color:var(--good)}
.tag.bad{background:#fbeae8;color:var(--bad)}
.tag.warn{background:#fdf3e0;color:var(--warn)}
.note{background:#f2f5fa;border-left:3px solid var(--accent);padding:12px 15px;
margin:14px 0;border-radius:0 8px 8px 0;font-size:14px}
.warnbox{background:#fdf6e8;border-left:3px solid var(--warn);padding:12px 15px;
margin:14px 0;border-radius:0 8px 8px 0;font-size:14px}
pre{background:#f4f6f9;border:1px solid var(--line);border-radius:8px;
padding:14px;overflow-x:auto;font-size:12.5px;line-height:1.5}
code{font-family:"SFMono-Regular",Consolas,monospace}
svg{max-width:100%}
"""


def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;"))


def slope_chart(pairs, width=980, height=300):
    """Two-column slope chart: real plateau vs sham plateau, one line each."""
    pad_l, pad_r, pad_t, pad_b = 62, 62, 22, 34
    usable = height - pad_t - pad_b
    lo, hi = 0.0, 1.6

    def y(v):
        v = min(max(float(v), lo), hi)
        return pad_t + (hi - v) / (hi - lo) * usable

    x1, x2 = pad_l + 130, width - pad_r - 130
    parts = ['<svg viewBox="0 0 {} {}" xmlns="http://www.w3.org/2000/svg">'.format(
        width, height)]
    parts.append('<rect width="{}" height="{}" fill="#ffffff"/>'.format(
        width, height))
    # reference line at 1.0 (no change)
    parts.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="#c8ced8" '
                 'stroke-dasharray="4 4"/>'.format(pad_l, y(1.0),
                                                   width - pad_r, y(1.0)))
    parts.append('<text x="{}" y="{}" font-size="11" fill="#8a929e">'
                 'ratio 1.0 = 无变化</text>'.format(pad_l + 2, y(1.0) - 5))
    for gv in [0.0, 0.4, 0.8, 1.2, 1.6]:
        parts.append('<text x="{}" y="{}" font-size="11" fill="#8a929e" '
                     'text-anchor="end">{:.1f}</text>'.format(
                         pad_l - 8, y(gv) + 4, gv))
    parts.append('<text x="{}" y="{}" font-size="13" font-weight="600" '
                 'fill="#1a1d21" text-anchor="middle">真实 Ricci 流</text>'.format(
                     x1, pad_t - 6))
    parts.append('<text x="{}" y="{}" font-size="13" font-weight="600" '
                 'fill="#1a1d21" text-anchor="middle">Sham 流（每步随机置换）</text>'.format(
                     x2, pad_t - 6))

    for real, sham in pairs:
        colour = "#1a7f4b" if sham > real else "#b4342b"
        parts.append('<line x1="{}" y1="{}" x2="{}" y2="{}" stroke="{}" '
                     'stroke-width="1.6" opacity="0.55"/>'.format(
                         x1, y(real), x2, y(sham), colour))
        parts.append('<circle cx="{}" cy="{}" r="2.6" fill="{}"/>'.format(
            x1, y(real), colour))
        parts.append('<circle cx="{}" cy="{}" r="2.6" fill="{}"/>'.format(
            x2, y(sham), colour))
    parts.append("</svg>")
    return "\n".join(parts)


def main():
    robust = pd.read_csv(os.path.join(ANALYSIS, "03_robust_table.csv"))
    family = pd.read_csv(os.path.join(ANALYSIS, "04_family_table.csv"))
    with open(os.path.join(ANALYSIS, "ROBUST_REPORT.txt"),
              encoding="utf-8") as fh:
        report = fh.read()

    keep = robust[~robust["sham_diverged"] & ~robust["real_inert"]].copy()
    pairs = list(zip(keep["ratio_real"], keep["ratio_sham"]))

    html = []
    html.append("<!DOCTYPE html><html lang='zh-CN'><head><meta charset='utf-8'>")
    html.append("<meta name='viewport' content='width=device-width,initial-scale=1'>")
    html.append("<title>Sham 流对照实验报告</title>")
    html.append("<style>{}</style></head><body><div class='wrap'>".format(CSS))

    html.append("<h1>Sham 流对照实验：曲率均质化是结构性的，还是更新规则的机械后果？</h1>")
    html.append("<div class='sub'>离散 Ollivier–Ricci 流 · 长程跑到饱和 · "
                "每步随机置换曲率—边对应关系 · 步数严格配对</div>")

    # ---- verdict banner
    html.append("<div class='card'>")
    html.append("<span class='tag ok'>结论：H1 成立</span> "
                "在 <b>22/22</b> 个可辩护配对上，破坏曲率—边映射后方差几乎不再下降，"
                "而真实流降到初始值的约 10%。")
    html.append("</div>")

    med_real = float(np.median(keep["ratio_real"]))
    med_sham = float(np.median(keep["ratio_sham"]))
    html.append("<div class='kpis'>")
    html.append("<div class='kpi ok'><div class='v'>{:.3f}</div>"
                "<div class='l'>真实流饱和平台（方差比中位数）</div></div>".format(med_real))
    html.append("<div class='kpi bad'><div class='v'>{:.3f}</div>"
                "<div class='l'>Sham 流平台（同一步数）</div></div>".format(med_sham))
    html.append("<div class='kpi'><div class='v'>22 / 22</div>"
                "<div class='l'>sham 高于 real 的数据集</div></div>")
    html.append("<div class='kpi'><div class='v'>2.4e-07</div>"
                "<div class='l'>符号检验 / Wilcoxon p 值</div></div>")
    html.append("<div class='kpi'><div class='v'>8 / 8</div>"
                "<div class='l'>家族层面一致（p = 0.0039）</div></div>")
    html.append("</div>")

    # ---- why this experiment
    html.append("<h2>1. 这个实验为什么非做不可</h2>")
    html.append("<p>长程实验已经证明：把流从原来的 10 步（t=0.5）一直跑到饱和，"
                "曲率方差比从 0.78 继续降到约 0.20，并且<b>停在一个平台而非降到 0</b>。"
                "对此有两种互斥解释：</p>")
    html.append("<div class='card'>")
    html.append("<p><b>H1（结构性）</b>：残余方差由网络几何决定。曲率—边的对应关系"
                "是真实几何信息，流只是在读取它。</p>")
    html.append("<p><b>H2（机械性）</b>：残余只是<b>没跑够</b>或<b>更新规则本身</b>造成的。"
                "权重更新 <code>w ← w·exp(−dt·κ)</code> 加上每步归一化，"
                "无论哪条边拿到哪个 κ，都可能机械地压缩权重离散度，"
                "那么“残余是结构性的”这个说法就必须撤回。</p>")
    html.append("</div>")
    html.append("<p><b>判别设计</b>：Sham 流每一步都把计算出的曲率值在边之间"
                "随机置换，然后再做完全相同的权重更新。它保留了 κ 的取值集合、"
                "步长、归一化与停止准则，<b>只破坏 κ 与边的对应关系</b>。"
                "两臂步数严格相同（25/26 完全一致）。</p>")

    # ---- slope chart
    html.append("<h2>2. 主结果</h2>")
    html.append("<p>每条线是一个网络：左端是真实流平台，右端是 sham 流平台。"
                "真实流几乎全部落到 0.3 以下，sham 流几乎全部停在 1.0 附近或更高"
                "（即方差不但没降，还略有上升）。</p>")
    html.append("<div class='card'>{}</div>".format(slope_chart(pairs)))

    html.append("<div class='note'><b>读图要点</b>：sham 平台中位数 {:.3f} &gt; 1，"
                "意味着随机化后的流不仅没有均质化，反而让方差略微发散。"
                "这说明方差下降<b>不是</b>更新规则的机械副产品。</div>".format(med_sham))

    # ---- tiered results
    html.append("<h3>三重稳健性：从原始到最保守</h3>")
    html.append("<table><thead><tr><th>分析层级</th><th>n</th><th>real 中位</th>"
                "<th>sham 中位</th><th>sham 更高</th><th>符号检验 p</th>"
                "<th>Wilcoxon p</th></tr></thead><tbody>")
    html.append("<tr><td>全部配对（原始）</td><td>26</td><td>0.150</td><td>1.138</td>"
                "<td>25</td><td>4.0e-07</td><td>4.5e-08</td></tr>")
    html.append("<tr><td>剔除发散 + 惰性臂</td><td>22</td><td>0.103</td><td>1.159</td>"
                "<td>22</td><td>2.4e-07</td><td>2.4e-07</td></tr>")
    html.append("<tr><td>严格匹配步数的 sham 端点</td><td>22</td><td>0.103</td>"
                "<td>1.159</td><td>22</td><td>2.4e-07</td><td>2.4e-07</td></tr>")
    html.append("<tr><td><b>家族层面（独立单元）</b></td><td>8</td>"
                "<td>0.070</td><td>1.160</td><td>8</td><td>0.0039</td>"
                "<td>0.0039</td></tr>")
    html.append("</tbody></table>")
    html.append("<div class='warnbox'><b>为什么必须做家族层面</b>：ENZYMES_g* 是同一"
                "生成器的 12 次抽样，scc_rt_* 是同一语料的 4 个子网络。"
                "把它们当作 12 个独立样本是伪重复。按家族聚合后 n 只剩 8，"
                "但方向仍然 8/8 一致，p = 0.0039 —— 这是最保守也最可信的一行。</div>")

    # ---- family table
    html.append("<h2>3. 家族层面明细</h2>")
    html.append("<table><thead><tr><th>家族</th><th>网络数</th><th>real 中位</th>"
                "<th>sham 中位</th><th>差值</th></tr></thead><tbody>")
    for _, row in family.iterrows():
        html.append("<tr><td>{}</td><td>{:.0f}</td><td>{:.4f}</td>"
                    "<td>{:.4f}</td><td>{:+.4f}</td></tr>".format(
                        esc(row["family"]), row["n"], row["real"],
                        row["sham"], row["diff"]))
    html.append("</tbody></table>")

    # ---- exclusions
    html.append("<h2>4. 排除项与被排除的理由</h2>")
    html.append("<table><thead><tr><th>数据集</th><th>排除原因</th>"
                "<th>real</th><th>sham</th></tr></thead><tbody>")
    for _, row in robust[robust["sham_diverged"] | robust["real_inert"]].iterrows():
        reason = ("sham 臂发散（权重坍缩 → κ 爆炸）" if row["sham_diverged"]
                  else "真实流本身也未均质化（平台 ≈ 1）")
        html.append("<tr><td>{}</td><td>{}</td><td>{:.4f}</td>"
                    "<td>{:.4f}</td></tr>".format(
                        esc(row["dataset"]), reason, row["ratio_real"],
                        row["ratio_sham"]))
    html.append("</tbody></table>")
    html.append("<p class='sub'>eco-mangwet 的 sham 臂在 7 步后发散"
                "（exponent ≈ 1092），剔除后结论不变；"
                "scc_rt_ksa、scc_rt_dash、bn-macaque-rhesus_interareal 三个网络"
                "真实流也几乎没动，属于无信息配对，同样剔除。</p>")

    # ---- what this unlocks
    html.append("<h2>5. 这个结果解锁了什么</h2>")
    html.append("<div class='card'>")
    html.append("<p><b>1. 撤掉“没跑够”这个替代解释。</b>两臂步数完全相同，"
                "sham 多跑也回不来。残余方差不是时间不足。</p>")
    html.append("<p><b>2. 论文 A 的核心主张可以成立。</b>“部分均质化 + 由结构决定的残余”"
                "不再是一个观察，而是一个有对照支撑的结论。"
                "主主张应表述为：<i>均质化需要真实的曲率—边映射；"
                "残余水平随网络几何（同配性）变化</i>。</p>")
    html.append("<p><b>3. 但距离“机制已验证”还差一步。</b>真实网络里同配性与"
                "聚类、密度、模块度纠缠在一起，相关性无法分离因果。"
                "正在进行的合成图实验（固定度序列、只改连线方式）"
                "是补齐这一步的下一项工作。</p>")
    html.append("</div>")

    html.append("<h2>6. 原始统计输出</h2>")
    html.append("<pre>{}</pre>".format(esc(report)))

    html.append("</div></body></html>")

    if not os.path.isdir(OUT_DIR):
        os.makedirs(OUT_DIR)
    path = os.path.join(OUT_DIR, "sham_control_report.html")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(html))
    print("written:", path)


if __name__ == "__main__":
    main()
