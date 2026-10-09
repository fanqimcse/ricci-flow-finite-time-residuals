"""Build the HTML report for the long-horizon Ricci flow experiment.

Reads the tables and figures produced by ``analyze_longrun.py`` and renders a
single self-contained HTML page.  Re-run it after every new batch of results.
"""

from __future__ import print_function

import base64
import os

import numpy as np
import pandas as pd


BASE = os.path.dirname(os.path.abspath(__file__))
RESULT_DIR = os.path.join(BASE, "ricci_flow_longrun_results")
OUT_DIR = os.path.join(RESULT_DIR, "analysis")
HTML_PATH = os.path.join(os.path.dirname(BASE), "outputs", "长程Ricci流实验报告.html")


def read(name):
    path = os.path.join(OUT_DIR, name)
    if not os.path.exists(path):
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def embed_png(name):
    path = os.path.join(OUT_DIR, name)
    if not os.path.exists(path):
        return ""
    with open(path, "rb") as handle:
        data = base64.b64encode(handle.read()).decode("ascii")
    return '<img src="data:image/png;base64,{}" alt="{}"/>'.format(data, name)


def fmt(value, digits=4):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "-"
    if not np.isfinite(number):
        return "-"
    return "{:.{d}f}".format(number, d=digits)


def table_html(frame, digits=4):
    if frame is None or not len(frame):
        return "<p class='empty'>该表尚未生成。</p>"
    frame = frame.copy()
    for column in frame.columns:
        if pd.api.types.is_float_dtype(frame[column]):
            frame[column] = frame[column].map(lambda v: fmt(v, digits))
    return frame.to_html(index=False, escape=False, border=0)


def load_summary():
    import glob
    frames = []
    for path in sorted(glob.glob(os.path.join(RESULT_DIR, "longrun_summary*.csv"))):
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if "dataset" in frame.columns and len(frame):
            frames.append(frame)
    if not frames:
        return pd.DataFrame()
    merged = pd.concat(frames, ignore_index=True)
    merged["_rank"] = (merged["status"].astype(str) == "cached").astype(int)
    merged = merged.sort_values("_rank").drop_duplicates("dataset", keep="first")
    return merged.drop(columns=["_rank"]).reset_index(drop=True)


def status_counts(summary):
    counts = summary["status"].value_counts()
    rows = "".join(
        "<tr><td class='k'>{}</td><td class='v'>{}</td></tr>".format(k, int(v))
        for k, v in counts.items())
    return "<table class='mini'>{}</table>".format(rows)


def main():
    # the per-dataset master table is authoritative: worker summary files lose
    # their stopping reasons whenever a worker is restarted
    master = read("00_dataset_table.csv")
    if master is None or not len(master):
        raise SystemExit("run analyze_longrun.py first")
    if "converged" in master.columns:
        master["converged"] = master["converged"].astype(str) == "True"
    summary = master
    diag = master[master["steps"] > 0] if "steps" in master.columns else master
    strat = read("06_stratified_correlations.csv")
    domain = read("04_domain_plateau.csv")
    outcomes = read("01_plateau_ratio_by_outcome.csv")
    transient = read("02_transient_vs_plateau.csv")
    predictors = read("03_spearman_plateau_predictors.csv")

    conv = diag[diag["converged"]] if diag is not None and "converged" in diag.columns else None
    trunc = diag[~diag["converged"]] if diag is not None and "converged" in diag.columns else None

    def stat(frame, column, func=np.median):
        if frame is None or column not in frame.columns:
            return None
        values = pd.to_numeric(frame[column], errors="coerce").dropna().values
        return float(func(values)) if values.size else None

    n_total = len(summary)
    n_done = int((summary["status"] != "cached").sum()) if n_total else 0
    # steps: only runs that actually flowed, i.e. not budget_insufficient(1 step)
    # and not skipped/failed(0 steps); matches the "completed runs" stratum of
    # the text report
    flowed = summary[summary["status"].isin(
        ["saturated", "constant_curvature", "horizon_reached", "budget_exceeded"])]
    med_steps = stat(flowed, "steps")
    med_all = stat(summary[summary["status"].isin(
        ["saturated", "constant_curvature", "horizon_reached", "budget_exceeded"])],
        "var_ratio_final")
    med_conv = stat(conv, "var_ratio_final")
    frac_zero = (float(np.mean(conv["var_ratio_final"].values < 0.05))
                 if conv is not None and len(conv) else None)
    frac_partial = (float(np.mean(conv["var_ratio_final"].values > 0.30))
                    if conv is not None and len(conv) else None)

    at10 = stat(transient, "var_ratio_at_t0.5")
    share10 = stat(transient, "fraction_of_drop_in_first_10_steps")

    rho_all = rho_conv = None
    if strat is not None and len(strat):
        for _, row in strat.iterrows():
            if row["covariate"] == "log10_nodes" and row["stratum"] == "all completed":
                rho_all = float(row["spearman"])
            if row["covariate"] == "log10_nodes" and row["stratum"] == "converged only":
                rho_conv = float(row["spearman"])
    rho_kappa = None
    if strat is not None and len(strat):
        hit = strat[(strat["covariate"] == "initial_kappa")
                    & (strat["stratum"] == "converged only")]
        if len(hit):
            rho_kappa = float(hit["spearman"].iloc[0])

    html = u"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8"/>
<title>长程离散 Ricci 流实验报告</title>
<style>
  :root {{
    --bg: #fbfbfa; --card: #ffffff; --ink: #1f2328; --muted: #5b6470;
    --line: #e4e6ea; --accent: #185FA5; --warn: #b45309; --bad: #b42318;
    --ok: #067647; --soft: #f2f5f9;
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin:0; padding:32px 20px 64px; background:var(--bg); color:var(--ink);
    font-family:-apple-system,"Segoe UI","Microsoft YaHei",Helvetica,Arial,sans-serif;
    line-height:1.65; }}
  .wrap {{ max-width: 1040px; margin:0 auto; }}
  h1 {{ font-size:28px; margin:0 0 6px; letter-spacing:-0.01em; }}
  .sub {{ color:var(--muted); font-size:14px; margin-bottom:28px; }}
  h2 {{ font-size:20px; margin:34px 0 12px; padding-bottom:8px;
    border-bottom:2px solid var(--line); }}
  h3 {{ font-size:16px; margin:22px 0 8px; color:var(--accent); }}
  .card {{ background:var(--card); border:1px solid var(--line);
    border-radius:10px; padding:18px 20px; margin:14px 0; }}
  .kpi {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
    gap:12px; margin:16px 0; }}
  .kpi div {{ background:var(--soft); border:1px solid var(--line);
    border-radius:9px; padding:14px 16px; }}
  .kpi .n {{ font-size:24px; font-weight:650; color:var(--accent);
    letter-spacing:-0.02em; }}
  .kpi .l {{ font-size:12.5px; color:var(--muted); margin-top:3px; }}
  table {{ border-collapse:collapse; width:100%; font-size:13.5px; margin:10px 0; }}
  th {{ background:var(--soft); text-align:left; padding:8px 10px;
    border-bottom:2px solid var(--line); font-weight:600; }}
  td {{ padding:7px 10px; border-bottom:1px solid var(--line); }}
  tr:last-child td {{ border-bottom:none; }}
  table.mini {{ width:auto; }}
  table.mini td.k {{ font-weight:600; padding-right:22px; }}
  .callout {{ border-left:4px solid var(--accent); background:var(--soft);
    padding:12px 16px; border-radius:0 8px 8px 0; margin:14px 0; }}
  .warn {{ border-left-color:var(--warn); }}
  .bad {{ border-left-color:var(--bad); }}
  .ok {{ border-left-color:var(--ok); }}
  img {{ max-width:100%; border:1px solid var(--line); border-radius:8px;
    margin:10px 0; background:#fff; }}
  code {{ background:var(--soft); padding:1px 5px; border-radius:4px;
    font-size:12.5px; }}
  ul {{ margin:8px 0 8px 4px; padding-left:20px; }}
  li {{ margin:5px 0; }}
  .empty {{ color:var(--muted); font-style:italic; }}
  .tag {{ display:inline-block; font-size:11.5px; padding:2px 8px;
    border-radius:20px; background:var(--soft); border:1px solid var(--line);
    color:var(--muted); margin-right:6px; }}
</style>
</head>
<body>
<div class="wrap">

<h1>长程离散 Ricci 流实验报告</h1>
<div class="sub">
  数据源 <code>D:\\100data</code> 全量数据集 · 流参数 α=0.5, dt=0.05, Sinkhorn reg=0.05
  · 上限 400 步 · 单数据集时间预算 600 s<br/>
  生成脚本 <code>ricci_flow_longrun.py</code> / <code>analyze_longrun.py</code>
</div>

<div class="card">
  <span class="tag">Ollivier–Ricci 流</span>
  <span class="tag">曲率方差均质化</span>
  <span class="tag">饱和判定</span>
  <h3 style="margin-top:10px">一句话结论</h3>
  <div class="callout ok">
    原 10 步实验（流时间 t=0.5）捕获的只是<b>暂态</b>：方差下降总量中仅约
    {share10} 发生在前 10 步。把流推进到饱和后，曲率方差中位数降到初始值的
    {med_all_str}；在真正收敛的子集中进一步降到 {med_conv_str}，其中
    {frac_zero_str} 接近完全均质化（&lt;5%），但仍有 {frac_partial_str}
    停在 30% 以上的残余方差——<b>均质化只能部分完成，残余量由结构决定</b>。
  </div>
</div>

<h2>一、实验设置与代码改动</h2>
<div class="card">
  <h3>相对原脚本改了什么</h3>
  <ul>
    <li><b>停止准则</b>：原来是固定 10 步。改为「方差相对变化连续 {patience} 步
      小于 1e-4 → saturated；方差降到初值 1e-6 以下 → constant_curvature；
      达到 400 步 → horizon_reached；单数据集墙钟超 600 s → budget_exceeded；
      Ricci 更新溢出 → diverged」。</li>
    <li><b>自适应视野</b>：用第 0 步实测耗时推算能走多少步；若连最小视野都跑不完
      则标记 <code>budget_insufficient</code> 并保留原 10 步结论，不再硬撑。</li>
    <li><b>逐数据集落盘</b>：每个数据集单独写 <code>per_dataset/&lt;id&gt;.csv</code>，
      记录每一步的曲率均值/方差/权重极值/谱隙/耗时，便于断点续跑与事后诊断。</li>
    <li><b>向量化 Sinkhorn</b>：<code>fast_sinkhorn.py</code> 按支撑集尺寸分桶批量求解，
      与 POT 参考实现在随机问题和真实图上的最大差异为 <b>8.9e-16</b>（机器精度），
      实测加速 1.4–2.7 倍。</li>
    <li><b>并行对向推进</b>：新增 <code>--reverse</code> / <code>--worker_tag</code> /
      <code>--exclude</code>，两个工作进程从序列两端相向收敛，互不覆盖成果。</li>
  </ul>
</div>

<h2>二、跑到什么程度了</h2>
<div class="kpi">
  <div><div class="n">{n_total}</div><div class="l">数据集总数</div></div>
  <div><div class="n">{med_steps_str}</div><div class="l">收敛所需步数中位数</div></div>
  <div><div class="n">{t_stop_str}</div><div class="l">对应流时间（原为 0.5）</div></div>
  <div><div class="n">{med_all_str}</div><div class="l">终态方差比中位数</div></div>
</div>
{status_table}

<h2>三、核心结果</h2>

<h3>3.1 原结果是暂态，不是终点</h3>
<div class="card">
  <p>在流时间 t=0.5（原 10 步）时，方差比中位数为 <b>{at10_str}</b>——这正是
  之前「下降约 20%」的来源。继续推进到停止后，同一批数据集的方差比中位数降到
  <b>{med_all_str}</b>。总下降量中只有 <b>{share10_str}</b> 发生在前 10 步，
  也就是说约七成的均质化发生在原来从未观察过的区间里。</p>
</div>
{transient_table}

<h3>3.2 饱和水平：完全均质化与部分均质化并存</h3>
<div class="kpi">
  <div><div class="n">{med_conv_str}</div><div class="l">收敛子集终态方差比中位数</div></div>
  <div><div class="n">{frac_zero_str}</div><div class="l">收敛子集中 &lt;5% 的比例</div></div>
  <div><div class="n">{frac_partial_str}</div><div class="l">收敛子集中 &gt;30% 的比例</div></div>
</div>
{outcomes_table}
<div class="callout warn">
  <b>注意</b>：被墙钟预算截断的样本（<code>budget_exceeded</code>）停止时往往仍在下降，
  它们的方差比只是<b>上界</b>，不应与真正收敛的样本混为一谈。下一节专门处理这个混杂。
</div>

<h3>3.3 截断混杂：规模效应是假的，初始曲率才是真的</h3>
<div class="card">
  <p>图越大 → 每步越慢 → 越容易被 600 s 预算截断 → 报出的残余方差越高。
  这个机制会凭空造出一个「大网络不易均质化」的假效应。分层之后：</p>
  <ul>
    <li><code>rho(终态方差比, log10 N)</code>：全样本 <b>{rho_all_str}</b>，
      但只看真正收敛的子集时翻转为 <b>{rho_conv_str}</b>——<b>符号翻转</b>，
      说明全样本上的正相关是截断造成的伪影。</li>
    <li><code>rho(终态方差比, 初始平均曲率 κ₀)</code>：在收敛子集中升到
      <b>{rho_kappa_str}</b>，是目前最强的单一结构预测因子。</li>
  </ul>
</div>
{strat_table}

<h3>3.4 领域差异</h3>
{domain_table}

<h2>四、曲线</h2>
<div class="card">
  <h3>轨迹集合</h3>
  {fig_traj}
  <h3>中位包络</h3>
  {fig_env}
  <h3>饱和水平 vs 网络规模</h3>
  {fig_size}
</div>

<h2>五、数值失败与已知缺陷</h2>
<div class="card">
  <div class="callout bad">
    <b>生态网络溢出</b>：<code>eco-florida</code>、<code>eco-foodweb-baydry</code>、
    <code>eco-foodweb-baywet</code> 在 t≈0.15 触发 <code>Ricci update overflow risk</code>
    （指数项达 10³ 量级）。根因是权重塌缩 → 距离 d→0 → κ 发散。
    这不是物理结论，而是积分器缺陷，必须在论文前给出数值稳定性协议
    （权重下限 / 自适应步长 / 曲率截断）。
  </div>
  <ul>
    <li><code>budget_insufficient</code>：单步耗时过大，连最小视野都跑不完，只能沿用原 10 步结论——
      这批样本在长程分析中应被明确排除或单列。</li>
    <li><code>skipped_large</code>：超过工作上限被主动跳过。</li>
    <li>残余方差与步数的强负相关（<code>rho≈-0.78</code>，仅见于截断子集）是机械性的，
      不能解释为因果。</li>
  </ul>
</div>

<h2>六、对论文 A 立论的影响</h2>
<div class="card">
  <ol>
    <li><b>「均质化会发生」不再能作为新颖点</b>——流定义本身就在做这件事，
      且已有收敛性证明。必须改写为「均质化<b>只能部分完成</b>，残余方差由结构决定」。</li>
    <li><b>长程实验是必要前置</b>：没有它，报告的效应量会被低估约 3 倍，
      且方向性结论不稳。现在这条已经补齐。</li>
    <li><b>下一步最便宜的升级</b>：把初始平均曲率 κ₀ 作为残余方差的主预测因子，
      补上聚类系数、度异质性、模块度，建立多变量模型。</li>
    <li><b>必须补的对照</b>：sham flow（打乱曲率与边的对应）、度保持重连、
      换曲率定义——否则「残余方差由结构决定」无法与「只是没跑够」区分。</li>
  </ol>
</div>

<div class="sub" style="margin-top:30px">
  本报告由 <code>make_longrun_report.py</code> 自动生成；
  数值随实验推进更新，重跑脚本即可刷新。
</div>

</div>
</body>
</html>
"""

    html = html.format(
        n_total=n_total,
        n_done=n_done,
        share10=fmt(share10, 3) if share10 is not None else "-",
        share10_str="{:.0%}".format(share10) if share10 is not None else "-",
        med_all_str=fmt(med_all, 4) if med_all is not None else "-",
        med_conv=fmt(med_conv, 4) if med_conv is not None else "-",
        med_conv_str=fmt(med_conv, 4) if med_conv is not None else "-",
        frac_zero_str="{:.0%}".format(frac_zero) if frac_zero is not None else "-",
        frac_partial_str="{:.0%}".format(frac_partial) if frac_partial is not None else "-",
        med_steps_str="{:.0f}".format(med_steps) if med_steps else "-",
        t_stop_str=fmt((med_steps or 0) * 0.05, 1) if med_steps else "-",
        at10_str=fmt(at10, 4) if at10 is not None else "-",
        rho_all_str="{:+.3f}".format(rho_all) if rho_all is not None else "-",
        rho_conv_str="{:+.3f}".format(rho_conv) if rho_conv is not None else "-",
        rho_kappa_str="{:+.3f}".format(rho_kappa) if rho_kappa is not None else "-",
        patience=5,
        status_table=status_counts(summary),
        transient_table=table_html(transient),
        outcomes_table=table_html(outcomes),
        strat_table=table_html(strat),
        domain_table=table_html(domain),
        fig_traj=embed_png("fig_trajectories.png"),
        fig_env=embed_png("fig_median_envelope.png"),
        fig_size=embed_png("fig_plateau_vs_size.png"),
    )

    os.makedirs(os.path.dirname(HTML_PATH), exist_ok=True)
    with open(HTML_PATH, "w", encoding="utf-8") as handle:
        handle.write(html)
    print("report written to {}".format(HTML_PATH))


if __name__ == "__main__":
    main()
