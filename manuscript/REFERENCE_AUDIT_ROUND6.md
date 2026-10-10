# 参考文献三方核验报告

核验对象：`paper_a/manuscript/references.bib` 全部 **45** 条，与 `manuscript.tex` 中全部 **58** 处引文。

三维：
1. **正文引文 ↔ 原文内容**：逐条定位 `\cite` 上下文，对可获取的全文（arXiv PDF / Crossref 摘要 / PMLR 会议页）逐句核对；
2. **文末书目信息**：每条至少两个独立权威源交叉比对（Crossref DOI 解析 / Crossref 题名检索 / arXiv API / OpenAlex / PMLR 官方论文集页）；
3. **支持该论断的原文具体内容**：给出原文中的对应句子，标注 SUPPORTED / SUPPORTED-WITH-NUANCE / CORRECTED。

## 一、总体结论

| 项目 | 结果 |
|---|---|
| 参考文献条目 | 45 |
| 正文引文处数 | 58 |
| 0 漏引 / 0 未被引 / 0 重复 | ✅ |
| 书目信息**已修正**的错误 | **4**（见第二节） |
| 正文论断**已修正**的过度表述 | **2**（见第三节） |
| 全文下载并逐句核验 | 12 篇（arXiv 全文 PDF） |
| 仅凭权威元数据核验（付费墙） | 其余条目 |
| 变音符号差异（Münch/Peyré/Boguñá/van der Hoorn） | 非错误，属排版，可保留 |

## 二、文末书目信息的四处错误（已修正）

| 条目 | 字段 | 原值 | 已改为 | 两个独立源 |
|---|---|---|---|---|
| `SinkhornKnopp1967` | volume | 27 | **21** | Crossref 10.2140/pjm.1967.21.343 and OpenAlex both give Pacific J. Math. 21(2):343-348. Volume 21, not 27. |
| `Feydy2019` | venue + volume + pages | ICML 2019, PMLR 97:1958-1967 | **AISTATS 2019, PMLR 89:2681-2690** | proceedings.mlr.press/v89/ lists 'Interpolating between Optimal Transport and MMD using Sinkhorn Divergences', Feydy et al., PMLR 89:2681-2690; OpenAlex independently gives volume 89, pp. 2681-2690. PMLR 97 is ICML 2019 and does not contain this paper. |
| `LinYau2010` | pages | 345-358 | **343-356** | Crossref 10.4310/mrl.2010.v17.n2.a13 and OpenAlex both give Math. Res. Lett. 17(2):343-356. |
| `BauerChungLinLiu2017` | DOI | 10.1090/proc/13404 | **10.1090/proc/13145** | 10.1090/proc/13404 resolves to Garcia-Ramos, 'A characterization of mu-equicontinuity...', Proc. AMS 145:3357-3368. The correct DOI for 'Curvature aspects of graphs' is 10.1090/proc/13145 (Proc. AMS 145:2033-2042), confirmed by Crossref and OpenAlex. |

这四处均已写入 `references.bib`，重新编译后 `manuscript.bbl` 与 PDF 参考文献表已确认渲染正确（0 error / 0 undefined / 45 条）。

## 三、正文论断与原文的对应关系（已修正两处）

### ⚠️ `MaYang2025` —— CORRECTED

- 稿件原表述：Manuscript previously: surgery 'deleting or adding edges' converges to constant curvature on each connected component.
- 原文实际内容：arXiv:2505.15395: 'after undergoing multiple surgeries, the evolving graph has a constant Ricci curvature on each connected component'; Definition B (A-surgery) removes edges only - nothing in the paper adds them. Manuscript now reads 'deleting edges'.

### ⚠️ `Samal2018` —— CORRECTED

- 稿件原表述：Manuscript previously: the Bakry-Emery, Ollivier and Forman discretisations 'disagree on the same graph', citing Samal et al.
- 原文实际内容：Sci. Rep. 8 (2018) 8650: 'our extensive computational analysis... shows that the two discretizations of Ricci curvature are highly correlated in many networks'. The cited source reports agreement, not disagreement, and none of the four cited works treats Bakry-Emery. Manuscript now reads 'give different values on the same graph' and adds LinYau2010 as the Bakry-Emery source.

### 其余关键论断（逐句核验，均 SUPPORTED）

| 条目 | 稿件论断 | 原文支持该论断的具体内容 |
|---|---|---|
| `LiMunch2026` | Manuscript: the discrete-time flow d_(n+1) = (1 - alpha*kappa_dn) d_n converges to a constant-curvature metric under finite/connected (irreducible), exact transport, admissible step size. | J. Funct. Anal. 290 (2026) 111367 (arXiv:2407.00314), abstract and Theorem 2/3: 'we prove the discrete-time Ollivier Ricci curvature flow d_{n+1} := (1 - alpha*kappa_{d_n}) d_n converges to a constant curvature metric on a finite weighted graph'; 'Saying a discrete Markov chain defined on a graph is irreducible is the same as saying the graph is connected'. Update formula matches the manuscript character for character. |
| `BaiLinLuWangYau2024` | Manuscript: solutions exist and are unique for all times. | Am. J. Math. 146 (2024) 1723-1747 (arXiv:2010.01802), Theorem 3: 'For any initial weighted graph G, by fixing the exit condition (I), there exists a unique solution X(t), for all time t in [0, infinity)'. |
| `Naama2024` | Manuscript: in the graph-embedding setting a modified flow with a contraction property converges. | arXiv:2407.21609: 'we prove the convergence of discrete Ricci flow to a constant curvature and stable distance metrics over the edges', proved via the Banach fixed-point theorem; the paper is explicitly framed as graph embedding ('Graph embedding approaches attempt to project graphs into geometric entities'). |
| `BaiLiuLai2026` | Manuscript: for trees the long-time dynamics of the Lin-Lu-Yau flow have been characterised. | arXiv:2601.02673: 'we obtain that the normalized curvature flow on trees converges to a constant curvature metric, and under the uniform measure, a complete classification of trees can be obtained based on the convergence results'. |
| `BaiHuaLinLiu2025` | Manuscript: same - long-time tree dynamics characterised. | arXiv:2509.22140v3: 'We analyze long-time dynamics of edge weights and curvatures, providing precise characterizations of their limiting behaviors... converges to metric with zero curvature on edges... only if the tree is a caterpillar tree'. |
| `VanderHoorn2023` | Manuscript: the Ollivier curvature of a random geometric graph converges to the Ricci curvature of the underlying manifold, in the density regime and neighbourhood radius of that theorem. | Discrete Comput. Geom. 70 (2023) 671-712 (arXiv:2009.04306), Theorem 3.1: convergence holds for specific scalings epsilon_n and delta_n - exactly the hedge the manuscript adds. |
| `Bourne2018` | Manuscript: the Ollivier curvature carries a free idleness parameter. | SIAM J. Discrete Math. 32 (2018) 1408-1424 (arXiv:1704.04398): 'In this notion an idleness parameter, p in [0,1], must be set in order to obtain a curvature kappa_p'. |
| `LinYau2010` | Manuscript: curvature controls eigenvalue estimates. | Math. Res. Lett. 17 (2010) 343-356: 'We prove that the Ricci curvature in the sense of Bakry and Emery is bounded below by -1 on locally finite graphs... We also get an estimate for the eigenvalue of Laplace operator on finite graphs'. |

## 四、其余条目：书目信息核对结果

两源（Crossref/arXiv 与 OpenAlex）均一致的条目按“通过”计。下表只列出需要说明的条目。

| 条目 | 说明 |
|---|---|
| `Perelman2002`/`Perelman2003` | arXiv 预印本，无期刊卷期；OpenAlex 匹配到的是这两篇的**正式发表版本**，并非同一记录，故其“卷/页”字段不参与比对。以 arXiv API 为准。 |
| `Cuturi2013` | NeurIPS 会议论文，Crossref 无记录。OpenAlex 给出 Advances in NeurIPS 26, 2292-2300，与 bib 一致 ✅。 |
| `Genevay2019` | PMLR 官方论文集页确认：AISTATS 2019, PMLR 89:1574-1583 ✅。 |
| `Janati2020` | PMLR 官方论文集页确认 citation_firstpage=4692 / lastpage=4701，与 bib 一致 ✅。 |
| `Boguna2021` | OpenAlex：Nature Reviews Physics 3, 114-135，与 bib 一致 ✅。 |
| `RossiAhmed2015` | OpenAlex：AAAI 2015 (Proc. AAAI Conf. Artif. Intell. 29)，页码 4292-4293 与 bib 一致 ✅。 |
| `Fan2026Code` | Zenodo API 直接核验：doi 10.5281/zenodo.23277162，concept 10.5281/zenodo.23277161，version v1.0.0，creator Qi Fan ✅。 |
| `CameronGelbachMiller2008` | Crossref 确认期刊为 **Review of Economics and Statistics** 90, 414-427（doi 10.1162/rest.90.3.414），与 bib 一致 ✅。 |
| `Sinkhorn1964` | 题名与卷页经 Crossref 核实一致 ✅（Ann. Math. Statist. 35, 876-879）。 |
| `LinYau2010` | 无 arXiv 版本，经 OpenAlex 摘要核验（'We prove that the Ricci curvature in the sense of Bakry and Emery is bounded below by -1... We also get an estimate for the eigenvalue of Laplace operator on finite graphs'），支持稿件“曲率控制特征值估计”✅。 |
| `JostLiu2014` | Crossref 记 2013（网络首发），Discrete Comput. Geom. 第 51 卷正刊为 2014；以期刊卷年为准，保留 2014。 |
| `SaucanSamalJost2021` | Crossref 记 2020（网络首发），Network Science 第 9 卷正刊为 2021；保留 2021。 |
| `MacKinnonNielsenWebb2023` | Crossref 记 2022（网络首发），J. Econometrics 第 232 卷正刊为 2023；保留 2023。 |

## 五、可复现的核验资产

| 文件 | 内容 |
|---|---|
| `ref_verify/bibliography_check.csv` | 45 条 × Crossref/arXiv 逐字段比对 |
| `ref_verify/openalex_check.csv` | 45 条 × OpenAlex 第二独立源比对 |
| `ref_verify/claim_context.csv` | 每一处 `\cite` 的正文上下文与所在小节 |
| `ref_verify/pdfs/` | 已下载的原文 PDF 及抽取的纯文本（12 篇） |
| `ref_verify/raw/Genevay2019.html`, `Janati2020.html` | PMLR 官方论文集页原始 HTML |
| `Supplementary_Table_S14_reference_audit_three_way.csv` | **三方核验总表**（本报告的机器可读版） |
| `_verify_refs.py`, `_verify_openalex.py`, `_verify_conf.py`, `_find_arxiv.py` | 核验脚本，可重跑 |
| `_patch_refs_r6.py` | 本轮四条书目 + 两处论断修正的幂等补丁 |
