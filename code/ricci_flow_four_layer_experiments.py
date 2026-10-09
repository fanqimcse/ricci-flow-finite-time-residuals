#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ricci-flow validation suite:
  Layer 1: real networks vs topology-preserving null audit
  Layer 2: synthetic controlled networks
  Layer 3: trajectory relaxation / Lyapunov analysis
  Layer 4: leave-one-domain-out prediction
  Experiment 5.1: curvature-energy evolution

No TensorFlow. No PyTorch. No graph-learning framework.

IMPORTANT:
raw_results.csv contains aggregate curvature statistics only. Therefore
Experiment 5.1 can test L_t = Var(kappa_t) and Delta L_t <= 0, but cannot
test the stronger edge-level dissipation inequality without saving the
edge-level curvature vector at every flow step.

Default inputs are the paths supplied in the conversation.
"""

from pathlib import Path
import argparse
import random
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

SEED = 20260830
BOOT = 5000
NULL = 5000
RAW_DEFAULT = Path(r"D:\2026codes\2026821100data\results_correct_fast\raw_results.csv")
SUMMARY_DEFAULT = Path(r"D:\2026codes\2026821100data\results_correct_fast\dataset_summary.csv")
OUT_DEFAULT = Path("ricci_flow_four_layer_results")


def finite(x):
    a = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(float)
    return a[np.isfinite(a)]


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    m = len(order)
    adj = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    q[order] = np.clip(adj, 0, 1)
    return q


def bootstrap(x, fn, seed):
    x = finite(x)
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    point = float(fn(x))
    if len(x) == 1:
        return point, point, point
    rng = np.random.default_rng(seed)
    z = x[rng.integers(0, len(x), size=(BOOT, len(x)))]
    v = fn(z, axis=1)
    return point, float(np.quantile(v, .025)), float(np.quantile(v, .975))


def wilcoxon(x):
    x = finite(x)
    x = x[np.abs(x) > 0]
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    try:
        r = stats.wilcoxon(x, zero_method="wilcox")
        ranks = stats.rankdata(np.abs(x))
        rb = (ranks[x > 0].sum() - ranks[x < 0].sum()) / ranks.sum()
        return float(r.statistic), float(r.pvalue), float(rb)
    except Exception:
        return np.nan, np.nan, np.nan


def load(raw_path, summary_path):
    raw = pd.read_csv(raw_path)
    req = {"dataset","step","mean_curvature","curvature_variance","lambda1","nodes","edges"}
    miss = sorted(req - set(raw.columns))
    if miss:
        raise ValueError("Missing raw columns: " + ", ".join(miss))
    raw["dataset"] = raw["dataset"].astype(str)
    raw["step"] = pd.to_numeric(raw["step"], errors="coerce")
    for c in req - {"dataset","step"}:
        raw[c] = pd.to_numeric(raw[c], errors="coerce")
    raw = raw.dropna(subset=["dataset","step"]).copy()
    raw["step"] = raw["step"].astype(int)
    raw = raw.sort_values(["dataset","step"]).reset_index(drop=True)

    if summary_path.exists():
        sm = pd.read_csv(summary_path)
        sm["dataset"] = sm["dataset"].astype(str)
        if "domain" not in sm:
            sm["domain"] = "Unknown"
        sm = sm[["dataset","domain"]].drop_duplicates("dataset")
    else:
        sm = pd.DataFrame({"dataset": raw.dataset.unique(), "domain":"Unknown"})
    raw = raw.merge(sm, on="dataset", how="left")
    raw["domain"] = raw["domain"].fillna("Unknown")
    return raw, sm


def dataset_energy(raw):
    rows = []
    for ds, g in raw.groupby("dataset", sort=True):
        g = g.sort_values("step")
        if len(g) < 2:
            continue
        L = g.curvature_variance.to_numpy(float)
        dL = np.diff(L)
        rows.append({
            "dataset": ds, "domain": g.domain.iloc[0],
            "nodes": g.nodes.iloc[0], "edges": g.edges.iloc[0],
            "n_steps": len(g),
            "initial_energy_L0": g.curvature_variance.iloc[0],
            "final_energy_LT": g.curvature_variance.iloc[-1],
            "energy_change_delta_L": g.curvature_variance.iloc[-1]-g.curvature_variance.iloc[0],
            "energy_reduction": g.curvature_variance.iloc[0]-g.curvature_variance.iloc[-1],
            "relative_energy_reduction": (
                (g.curvature_variance.iloc[0]-g.curvature_variance.iloc[-1]) /
                abs(g.curvature_variance.iloc[0])
                if abs(g.curvature_variance.iloc[0]) > 1e-15 else np.nan),
            "monotonic_L": bool(np.all(dL <= 1e-12)),
            "monotonic_fraction": float(np.mean(dL <= 1e-12)),
            "n_increasing_steps": int(np.sum(dL > 1e-12)),
            "max_step_increase": float(np.max(dL)),
            "initial_kappa": g.mean_curvature.iloc[0],
            "final_kappa": g.mean_curvature.iloc[-1],
            "lambda_initial": g.lambda1.iloc[0],
            "lambda_final": g.lambda1.iloc[-1]
        })
    return pd.DataFrame(rows)


def step_energy(raw):
    rows=[]
    for ds,g in raw.groupby("dataset", sort=True):
        g=g.sort_values("step")
        prev=None
        for _,r in g.iterrows():
            L=float(r.curvature_variance)
            d=np.nan if prev is None else L-prev
            rows.append({
                "dataset":ds,"domain":r.domain,"step":int(r.step),
                "L_t":L,"delta_L":d,
                "relative_delta_L": d/abs(prev) if prev is not None and abs(prev)>1e-15 else np.nan,
                "dissipative_observation": bool(d<=1e-12) if prev is not None else np.nan
            })
            prev=L
    return pd.DataFrame(rows)


def energy_stats(de):
    specs=[
        ("energy_change_delta_L", de.energy_change_delta_L),
        ("energy_reduction", de.energy_reduction),
        ("relative_energy_reduction", de.relative_energy_reduction),
        ("monotonic_fraction", de.monotonic_fraction)]
    rows=[]
    for name,x in specs:
        x=finite(x)
        if len(x)==0: continue
        mean=float(np.mean(x)); med=float(np.median(x))
        ci=stats.t.interval(.95,len(x)-1,loc=mean,scale=stats.sem(x)) if len(x)>1 else (mean,mean)
        bm,blo,bhi=bootstrap(x,np.mean,SEED+11)
        bd,blo2,bhi2=bootstrap(x,np.median,SEED+12)
        w,p,rb=wilcoxon(x)
        rows.append({
            "endpoint":name,"n":len(x),"mean":mean,
            "sd":float(np.std(x,ddof=1)) if len(x)>1 else np.nan,
            "median":med,"mean_ci95_low_t":ci[0],"mean_ci95_high_t":ci[1],
            "bootstrap_mean":bm,"bootstrap_mean_ci95_low":blo,"bootstrap_mean_ci95_high":bhi,
            "bootstrap_median":bd,"bootstrap_median_ci95_low":blo2,"bootstrap_median_ci95_high":bhi2,
            "wilcoxon_statistic":w,"wilcoxon_p":p,"rank_biserial":rb})
    out=pd.DataFrame(rows)
    if len(out): out["wilcoxon_q_BH"]=bh(out.wilcoxon_p)
    return out


def signflip(de):
    rows=[]
    rng=np.random.default_rng(SEED+100)
    for j,col in enumerate(["energy_change_delta_L","energy_reduction","relative_energy_reduction"]):
        x=finite(de[col])
        if len(x)==0: continue
        s=rng.choice([-1,1],size=(NULL,len(x)))
        means=(s*x).mean(1); meds=np.median(s*x,axis=1)
        om=float(x.mean()); od=float(np.median(x))
        rows.append({
            "endpoint":col,"n":len(x),"observed_median":od,
            "median_sign_flip_p":(np.sum(np.abs(meds)>=abs(od))+1)/(NULL+1),
            "observed_mean":om,
            "mean_sign_flip_p":(np.sum(np.abs(means)>=abs(om))+1)/(NULL+1)})
    out=pd.DataFrame(rows)
    if len(out):
        out["median_q_BH"]=bh(out.median_sign_flip_p)
        out["mean_q_BH"]=bh(out.mean_sign_flip_p)
    return out


# ---------- Layer 1: topology-preserving null audit ----------

def read_edges(path):
    es=[]
    with open(path,"r",encoding="utf-8",errors="ignore") as f:
        for line in f:
            p=line.replace(","," ").split()
            if len(p)<2 or p[0].startswith("#"): continue
            try: u=int(float(p[0])); v=int(float(p[1]))
            except: continue
            if u!=v: es.append(tuple(sorted((u,v))))
    return np.unique(np.asarray(es,dtype=int),axis=0) if es else None


def degree_seq(n,e):
    d=np.zeros(n,dtype=int)
    for u,v in e: d[u]+=1; d[v]+=1
    return d


def components(n,e):
    p=np.arange(n)
    def f(a):
        while p[a]!=a:
            p[a]=p[p[a]]; a=p[a]
        return a
    for a,b in e:
        ra,rb=f(a),f(b)
        if ra!=rb:p[rb]=ra
    return len({f(i) for i in range(n)})


def swap(e,seed,attempts_factor=30):
    rng=np.random.default_rng(seed)
    e=[tuple(x) for x in e]; S=set(e)
    target=max(10,5*len(e)); acc=att=0
    while acc<target and att<target*attempts_factor:
        att+=1
        i,j=rng.choice(len(e),2,replace=False)
        a,b=e[i]; c,d=e[j]
        if len({a,b,c,d})<4: continue
        if rng.random()<.5: x,y=(min(a,d),max(a,d)),(min(c,b),max(c,b))
        else: x,y=(min(a,c),max(a,c)),(min(b,d),max(b,d))
        if x[0]==x[1] or y[0]==y[1] or x==y or x in S or y in S: continue
        S.remove(e[i]); S.remove(e[j]); S.add(x); S.add(y)
        e[i],e[j]=x,y; acc+=1
    return np.asarray(e,dtype=int)


def layer1(edge_dir,de,nrep):
    rows=[]
    for _,r in de.iterrows():
        ds=str(r.dataset)
        p=None if edge_dir is None else next(
            (q for q in [edge_dir/f"{ds}.edges",edge_dir/f"{ds}.edge",
                         edge_dir/f"{ds}.txt",edge_dir/ds] if q.exists()),None)
        if p is None:
            rows.append({"dataset":ds,"status":"edge_file_not_found_or_not_supplied",
                         "real_nodes":r.nodes,"real_edges":r.edges})
            continue
        e=read_edges(p)
        if e is None:
            rows.append({"dataset":ds,"status":"unreadable","real_nodes":r.nodes,"real_edges":r.edges})
            continue
        nodes=np.unique(e); mp={int(v):i for i,v in enumerate(nodes)}
        e=np.asarray([[mp[int(a)],mp[int(b)]] for a,b in e])
        n=len(nodes); rd=np.sort(degree_seq(n,e)); rc=components(n,e)
        dok=[]; cok=[]
        for k in range(nrep):
            z=swap(e,SEED+k+abs(hash(ds))%10000)
            dok.append(np.array_equal(rd,np.sort(degree_seq(n,z))))
            cok.append(rc==components(n,z))
        rows.append({"dataset":ds,"status":"ok","edge_file":str(p),
                     "real_nodes":n,"real_edges":len(e),"null_replicates":nrep,
                     "degree_sequence_preserved_rate":np.mean(dok),
                     "component_count_preserved_rate":np.mean(cok)})
    return pd.DataFrame(rows)


# ---------- Layer 2: controlled synthetic topology ----------

def er(n,p,rng):
    ij=np.triu_indices(n,1); keep=rng.random(len(ij[0]))<p
    return np.column_stack((ij[0][keep],ij[1][keep]))


def ws(n,k,beta,rng):
    k=max(2,min(n-1,k-(k%2))); E=set()
    for i in range(n):
        for d in range(1,k//2+1):
            j=(i+d)%n; E.add(tuple(sorted((i,j))))
    for i in range(n):
        for d in range(1,k//2+1):
            j=(i+d)%n; old=tuple(sorted((i,j)))
            if rng.random()<beta:
                E.discard(old)
                cand=[x for x in range(n) if x!=i and tuple(sorted((i,x))) not in E]
                if cand:
                    x=int(rng.choice(cand)); E.add(tuple(sorted((i,x))))
                else:E.add(old)
    return np.asarray(sorted(E),int)


def ba(n,m,rng):
    m=max(1,min(m,n-1)); E=set()
    for i in range(m+1):
        for j in range(i+1,m+1):E.add((i,j))
    deg=np.zeros(n)
    for a,b in E:deg[a]+=1;deg[b]+=1
    for new in range(m+1,n):
        pr=deg[:new]; pr=pr/pr.sum()
        targets=rng.choice(new,m,replace=False,p=pr)
        for t in targets:
            E.add((min(new,int(t)),max(new,int(t))))
            deg[new]+=1;deg[int(t)]+=1
    return np.asarray(sorted(E),int)


def topo(e):
    n=int(e.max()+1) if len(e) else 0; m=len(e)
    d=degree_seq(n,e) if n else np.array([])
    du=d[e[:,0]]; dv=d[e[:,1]] if len(e) else np.array([])
    return {"nodes":n,"edges":m,"mean_degree":d.mean() if n else np.nan,
            "degree_sd":d.std(ddof=1) if n>1 else 0,
            "density":2*m/(n*(n-1)) if n>1 else np.nan,
            "components":components(n,e) if n else 0,
            "degree_assortativity":np.corrcoef(du,dv)[0,1]
            if len(e)>1 and du.std()>0 and dv.std()>0 else np.nan}


def layer2():
    rng=np.random.default_rng(SEED); rows=[]
    settings=[("ER",100,.03,0),("ER",100,.08,0),("ER",250,.02,0),
              ("WS",100,6,.05),("WS",100,6,.30),("WS",250,10,.10),
              ("BA",100,3,0),("BA",250,4,0)]
    for model,n,a,b in settings:
        for rep in range(20):
            e=er(n,a,rng) if model=="ER" else ws(n,int(a),b,rng) if model=="WS" else ba(n,int(a),rng)
            z=topo(e); z.update(model=model,parameter_1=a,parameter_2=b,replicate=rep)
            rows.append(z)
    return pd.DataFrame(rows)


# ---------- Layer 4: leave-one-domain-out prediction ----------

def ridge_predict(Xtr,ytr,Xte,alpha=1.0):
    mu=np.nanmean(Xtr,0); sd=np.nanstd(Xtr,0); sd[sd<1e-12]=1
    a=(Xtr-mu)/sd; b=(Xte-mu)/sd
    a=np.c_[np.ones(len(a)),a]; b=np.c_[np.ones(len(b)),b]
    P=np.eye(a.shape[1])*alpha; P[0,0]=0
    beta=np.linalg.solve(a.T@a+P,a.T@ytr)
    return b@beta


def layer4(de):
    feats=["nodes","edges","initial_kappa","initial_energy_L0"]
    x=de[feats+["energy_reduction","domain"]].replace([np.inf,-np.inf],np.nan).dropna()
    pred=[]; summ=[]
    for dom in sorted(x.domain.unique()):
        tr=x[x.domain!=dom]; te=x[x.domain==dom]
        if len(tr)<10 or len(te)<1:continue
        p=ridge_predict(tr[feats].to_numpy(float),tr.energy_reduction.to_numpy(float),
                        te[feats].to_numpy(float))
        y=te.energy_reduction.to_numpy(float)
        for a,b in zip(y,p):pred.append({"held_out_domain":dom,"actual":a,"predicted":b,"residual":a-b})
        rho,pp=(stats.spearmanr(y,p) if len(y)>1 and np.std(y)>0 and np.std(p)>0 else (np.nan,np.nan))
        summ.append({"held_out_domain":dom,"n_test":len(te),"n_train":len(tr),
                     "spearman_rho":rho,"spearman_p":pp,
                     "rmse":np.sqrt(np.mean((y-p)**2)),"mae":np.mean(np.abs(y-p))})
    return pd.DataFrame(pred),pd.DataFrame(summ)


# ---------- figures ----------

def figures(raw,de,ood,out):
    fig,ax=plt.subplots(figsize=(10,7))
    A=[]
    for _,g in raw.groupby("dataset"):
        g=g.sort_values("step"); ax.plot(g.step,g.curvature_variance,alpha=.10,lw=.8)
        A.append(g.curvature_variance.to_numpy(float))
    if A:
        A=np.asarray(A); k=A.shape[1]; s=np.arange(k)
        ax.plot(s,np.nanmedian(A,0),lw=2.5,label="Dataset median")
        ax.fill_between(s,np.nanquantile(A,.025,0),np.nanquantile(A,.975,0),alpha=.2,label="95% envelope")
    ax.set(xlabel="Ricci-flow step",ylabel=r"Curvature energy $L_t=Var(\kappa_t)$",
           title="Experiment 5.1 — Curvature Energy Evolution"); ax.legend()
    fig.tight_layout(); fig.savefig(out/"11_energy_trajectories.png",dpi=300); plt.close(fig)

    fig,ax=plt.subplots(figsize=(9,6)); x=finite(de.relative_energy_reduction)
    ax.hist(x,bins=25,alpha=.75); ax.axvline(np.median(x),lw=2,label="Median"); ax.axvline(0,ls="--",lw=1.5,label="No reduction")
    ax.set(xlabel="Relative curvature-energy reduction",ylabel="Datasets",title="Dataset-level energy relaxation"); ax.legend()
    fig.tight_layout(); fig.savefig(out/"12_energy_reduction_distribution.png",dpi=300); plt.close(fig)

    d=de.groupby("domain").relative_energy_reduction.median().sort_values()
    fig,ax=plt.subplots(figsize=(10,7)); ax.barh(np.arange(len(d)),d.values); ax.set_yticks(np.arange(len(d))); ax.set_yticklabels(d.index); ax.axvline(0,ls="--")
    ax.set(xlabel="Median relative energy reduction",title="Cross-domain curvature-energy relaxation")
    fig.tight_layout(); fig.savefig(out/"13_domain_energy.png",dpi=300); plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,7))
    if len(ood):
        ax.scatter(ood.actual,ood.predicted,alpha=.65)
        lo=min(ood.actual.min(),ood.predicted.min()); hi=max(ood.actual.max(),ood.predicted.max())
        ax.plot([lo,hi],[lo,hi],ls="--")
    ax.set(xlabel="Observed energy reduction",ylabel="OOD predicted energy reduction",title="Layer 4 — Leave-one-domain-out prediction")
    fig.tight_layout(); fig.savefig(out/"14_ood_prediction.png",dpi=300); plt.close(fig)


def report(out,raw,sm,de,es,mono,nul,l1,syn,ood_s):
    lines=["RICCI FLOW FOUR-LAYER + EXPERIMENT 5.1 REPORT","="*78,"",
           f"Raw records: {len(raw)}",f"Raw datasets: {raw.dataset.nunique()}",
           f"Summary datasets: {sm.dataset.nunique()}",f"Domains: {sm.domain.nunique()}","",
           "EXPERIMENT 5.1 — ENERGY EVOLUTION","-"*78,
           f"Dataset-level n = {len(de)}",
           f"Monotone trajectories: {int(de.monotonic_L.sum())}/{len(de)} ({de.monotonic_L.mean():.4f})",
           f"Median relative energy reduction: {de.relative_energy_reduction.median():.8g}","",
           es.to_string(index=False),"","SIGN-FLIP NULL","-"*78,nul.to_string(index=False),"",
           "LAYER 1","-"*78,l1.status.value_counts(dropna=False).to_string(),
           "Topology rewiring is an audit here; null curvature trajectories require rerunning the identical flow engine on rewired graphs.","",
           "LAYER 2","-"*78,f"Synthetic networks: {len(syn)}",syn.groupby("model").size().to_string(),
           "Synthetic topology is generated here; Ricci-flow endpoints must be computed by the same curvature engine before making a flow-effect claim.","",
           "LAYER 3","-"*78,
           "Aggregate test: Delta L_t <= 0. The stronger edge-level dissipation inequality cannot be evaluated from raw_results.csv.","",
           "LAYER 4","-"*78,ood_s.to_string(index=False) if len(ood_s) else "No valid OOD folds.","",
           "SCIENTIFIC BOUNDARY","-"*78,
           "Dataset is the independent experimental unit; ten steps are repeated trajectory observations, not ten independent replicates.",
           "No TensorFlow is used."]
    (out/"analysis_report.txt").write_text("\n".join(lines),encoding="utf-8")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--raw",type=Path,default=RAW_DEFAULT)
    ap.add_argument("--summary",type=Path,default=SUMMARY_DEFAULT)
    ap.add_argument("--out",type=Path,default=OUT_DEFAULT)
    ap.add_argument("--edge-dir",type=Path,default=None)
    ap.add_argument("--null-replicates",type=int,default=10)
    a=ap.parse_args()
    random.seed(SEED); np.random.seed(SEED); a.out.mkdir(parents=True,exist_ok=True)

    print("Loading:",a.raw)
    raw,sm=load(a.raw,a.summary)
    de=dataset_energy(raw); se=step_energy(raw); es=energy_stats(de); nul=signflip(de)
    mono=pd.DataFrame([{
        "datasets":len(de),"monotone_n":int(de.monotonic_L.sum()),
        "monotone_rate":de.monotonic_L.mean(),
        "mean_stepwise_dissipation_fraction":de.monotonic_fraction.mean(),
        "median_stepwise_dissipation_fraction":de.monotonic_fraction.median(),
        "total_transitions":int((de.n_steps-1).sum()),
        "increasing_transitions":int(de.n_increasing_steps.sum())}])
    dom=de.groupby("domain").agg(n_datasets=("dataset","size"),median_nodes=("nodes","median"),
        median_edges=("edges","median"),energy_reduction_mean=("energy_reduction","mean"),
        energy_reduction_median=("energy_reduction","median"),
        relative_reduction_mean=("relative_energy_reduction","mean"),
        relative_reduction_median=("relative_energy_reduction","median"),
        monotone_rate=("monotonic_L","mean")).reset_index()
    l1=layer1(a.edge_dir,de,a.null_replicates)
    syn=layer2(); ood,oods=layer4(de)

    de.to_csv(a.out/"01_energy_dataset_level.csv",index=False)
    se.to_csv(a.out/"02_energy_stepwise.csv",index=False)
    es.to_csv(a.out/"03_energy_statistics.csv",index=False)
    mono.to_csv(a.out/"04_energy_monotonicity.csv",index=False)
    nul.to_csv(a.out/"05_sign_flip_null.csv",index=False)
    dom.to_csv(a.out/"06_domain_energy_summary.csv",index=False)
    l1.to_csv(a.out/"07_layer1_topology_null_audit.csv",index=False)
    syn.to_csv(a.out/"08_layer2_synthetic_networks.csv",index=False)
    ood.to_csv(a.out/"09_layer4_ood_predictions.csv",index=False)
    oods.to_csv(a.out/"10_layer4_ood_summary.csv",index=False)
    figures(raw,de,ood,a.out)
    report(a.out,raw,sm,de,es,mono,nul,l1,syn,oods)

    print("Finished:",a.out.resolve())
    print("TensorFlow: NOT USED")
    print("Experiment 5.1 + Layers 1–4 outputs generated.")


if __name__=="__main__":
    main()
