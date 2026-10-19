#!/usr/bin/env python3
"""The control claims, as tables rather than console output.

    .venv/bin/python analysis/paper_statistics/control_tables.py

Four claims that the revision rests on and that had no deposited table, so a
reader could not check them without re-running a figure script and reading its
stdout:

  5  WHAT THE INJECTED NOISE BUYS. Task input x injected noise, on the corrected
     rule. The covariance the method uses is noise-induced: removing the noise
     costs a great deal, removing the task input costs nothing.

  6  WHAT A GLOBAL FACTOR BUYS, BY MASK. The same two derived constants applied
     to the noise score, to magnitude and to the random mask. A global factor is
     not uniformly useless -- it lifts the unbiased mask and destroys both
     score-selected ones -- so this is the table behind the narrowed claim.

  7  WHY DM1 AND ANTI ARE EXCLUDED. Not that pruning leaves them intact, but
     that the methods stop separating on them: the spread between best and worst
     of five methods is about twenty times smaller than on the retained tasks.

  8  SAMPLING AGAINST A GLOBAL FACTOR, WITHIN A SCORE. The contrast that was
     previously confounded, because it compared noise-prune WITH sampling
     against MAGNITUDE with a global factor.

Every effect is the mean of the 8 paired task-level differences, every interval
is a percentile of 10,000 two-stage hierarchical resamples (tasks, then seeds
within task), and every p is the exact sign-flip test over all 256 assignments.
A p at 0.0078 is the floor, marked *, and means the data are as extreme as this
test can register -- read the effect size and sign count alongside it.
"""
from __future__ import annotations

import glob
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_stats import (Cube, GRID, N_BOOT, auc, contrast,  # noqa: E402
                         make_surrogates)

REPO = str(pathlib.Path(__file__).resolve().parents[2])
R = f"{REPO}/results"
B = f"{R}/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_"
NEED = {"strategy", "amount", "post_acc_sequence", "pre_acc_sequence", "task",
        "source_network_seed"}
EXTRA = ("prune_gain_restore_mode", "prune_sigma_factor", "prune_zero_task_input",
         "sim_np_sigma_factor", "sim_np_zero_task_input", "pruning_seed")
_POOL: dict = {}


def pool(pat: str) -> pd.DataFrame:
    if pat in _POOL:
        return _POOL[pat]
    fr = []
    for f in sorted(glob.glob(pat)):
        cols = set(pd.read_csv(f, nrows=0).columns)
        if not NEED <= cols:
            continue
        fr.append(pd.read_csv(f, usecols=list(NEED | {e for e in EXTRA if e in cols}),
                              low_memory=False))
    d = pd.concat(fr, ignore_index=True)
    d = d[d.strategy != "none"].copy()
    d["t"] = d.task.str.replace("modcog:", "", regex=False)
    d["pct"] = (d.amount * 100).round().astype(int)
    d["ret"] = d.post_acc_sequence / d.pre_acc_sequence
    _POOL[pat] = d
    return d


def make_cube(arms: dict, qs) -> Cube:
    qs = np.array(sorted(qs))
    tasks = sorted(set.intersection(*[set(df.t) for df in arms.values()]))
    seeds = sorted(set.intersection(*[set(df.source_network_seed) for df in arms.values()]))
    names = list(arms)
    Rm = np.full((len(tasks), len(seeds), len(names), len(qs)), np.nan)
    D = np.zeros_like(Rm)
    for mi, n in enumerate(names):
        g = arms[n].groupby(["t", "source_network_seed", "pct"]).ret.agg(["mean", "size"])
        for (t, s, q), row in g.iterrows():
            if t in tasks and s in seeds and q in qs:
                Rm[tasks.index(t), seeds.index(s), mi, list(qs).index(q)] = row["mean"]
                D[tasks.index(t), seeds.index(s), mi, list(qs).index(q)] = row["size"]
    if np.isnan(Rm).any():
        raise SystemExit(f"{int(np.isnan(Rm).sum())} empty cells in {names}")
    return Cube(R=Rm, tasks=tasks, seeds=seeds, methods=names, qs=qs, draws=D)


def rows_for(cube, sur, a, b, label, qs):
    out = []
    tl = cube.task_level()
    ia, ib = cube.m(a), cube.m(b)
    for q in list(qs) + [None]:
        r = contrast(cube, sur, a, b, q=q)
        if q is None:
            ma, mb = (auc(tl, cube.qs)[:, i].mean() for i in (ia, ib))
        else:
            k = list(cube.qs).index(q)
            ma, mb = tl[:, ia, k].mean(), tl[:, ib, k].mean()
        out.append(dict(contrast=label,
                        sparsity=("AUC 10-80%" if q is None else f"{q}%"),
                        mean_a=round(float(ma), 4), mean_b=round(float(mb), 4),
                        effect=round(r["effect"], 4),
                        CI95=f"[{r['lo95']:+.3f}, {r['hi95']:+.3f}]",
                        tasks=f"{r['n_positive']}/{r['n_tasks']}",
                        p=f"{r['p_two']:.4f}" + ("*" if np.isclose(r["p_two"], r["p_floor"]) else "")))
    return out


def main() -> None:
    out = {}

    # ---- 5: what the injected noise buys -------------------------------
    d = pool(f"{B}cov2x2_p10_80/*.csv")
    sf = d["prune_sigma_factor"] if "prune_sigma_factor" in d else d["sim_np_sigma_factor"]
    zi = d["prune_zero_task_input"] if "prune_zero_task_input" in d else d["sim_np_zero_task_input"]
    d = d.assign(noise_on=sf.astype(float) > 0.5,
                 task_on=~zi.astype(str).str.lower().isin(["true", "1"]))
    REF, TASKONLY, NOISEONLY = "task+noise", "task only (no noise)", "noise only (no task input)"
    arms = {REF: d[d.task_on & d.noise_on],
            TASKONLY: d[d.task_on & ~d.noise_on],
            NOISEONLY: d[~d.task_on & d.noise_on]}
    cube = make_cube(arms, GRID)
    sur = make_surrogates(cube, n_boot=N_BOOT)
    T5 = pd.DataFrame(rows_for(cube, sur, REF, TASKONLY, "cost of removing the injected noise", GRID)
                      + rows_for(cube, sur, REF, NOISEONLY, "cost of removing the task input", GRID))
    print("=" * 140)
    print("TABLE 5  What the injected noise buys (sigma -> 0 is approached at 1e-3 x sigma_nat; exactly 0 is undefined)")
    print("  the covariance the method uses is noise-induced: removing the noise costs a great deal, removing the task input costs nothing")
    print("=" * 140)
    print(T5.to_string(index=False))
    out["control_table5_noise_vs_task_covariance"] = T5

    # ---- 6 + 8: what a global factor buys, by mask ---------------------
    sm = pool(f"{B}score_matched_gain_p10_80/*.csv")
    syn = pool(f"{B}synaptic_p50_80/*.csv")
    rg = pool(f"{B}random_gain_p10_80/*.csv")
    sc = pool(f"{B}score_control_full_p50_80/*.csv")
    base = pool(f"{B}renorm_p10_80/*.csv")
    mag = pool(f"{B}fig3_lowsparsity_p10_40/*.csv")

    def gain(frames, strategy, mode):
        x = pd.concat([f[(f.strategy == strategy) & (f.get("prune_gain_restore_mode") == mode)]
                       for f in frames], ignore_index=True)
        return x.drop_duplicates(subset=["strategy", "amount", "task", "source_network_seed",
                                         "pruning_seed", "post_acc_sequence"])

    SETS = [
        ("noise score", base[base.strategy == "simulation_noise_prune_rescale"],
         gain([sm, syn], "simulation_noise_prune_gain", "match_l1"),
         gain([sm, syn], "simulation_noise_prune_gain", "inv_density")),
        # the magnitude gain arms span 10-80% across TWO suites: 10-40% in the
        # low-sparsity suite and 50-80% in the score-control suite. Reading only
        # the first silently truncated this row to 10-40% and made its AUC
        # incomparable with the noise score's.
        ("magnitude score", base[base.strategy == "simulation_noise_prune_magnitude_rescale"],
         gain([mag, sc], "l1_unstructured_gain", "match_l1"),
         gain([mag, sc], "l1_unstructured_gain", "inv_density")),
        ("random mask", None,
         gain([rg], "random_unstructured_gain", "match_l1"),
         gain([rg, sc], "random_unstructured_gain", "inv_density")),
    ]
    rows = []
    for name, sampled, l1, dens in SETS:
        arms = {"|W|-matched": l1, "1/density": dens}
        if sampled is not None:
            arms["sample-and-rescale"] = sampled
        qs = sorted(set.intersection(*[set(a.pct) for a in arms.values()]))
        cu = make_cube(arms, qs)
        su = make_surrogates(cu, n_boot=N_BOOT)
        tl = cu.task_level()
        for lab in ("|W|-matched", "1/density"):
            k = cu.m(lab)
            r = {"mask": name, "rescale": lab}
            for qi, q in enumerate(cu.qs):
                r[f"{q}%"] = round(float(tl[:, k, qi].mean()), 4)
            if sampled is not None:
                c = contrast(cu, su, "sample-and-rescale", lab,
                             lo=min(qs), hi=max(qs))
                r["sampling minus this (AUC)"] = round(c["effect"], 4)
                r["CI95"] = f"[{c['lo95']:+.3f}, {c['hi95']:+.3f}]"
                r["tasks"] = f"{c['n_positive']}/8"
                r["p"] = f"{c['p_two']:.4f}" + ("*" if np.isclose(c["p_two"], c["p_floor"]) else "")
            else:
                r.update({"sampling minus this (AUC)": "-- (no sampled arm)",
                          "CI95": "--", "tasks": "--", "p": "--"})
            rows.append(r)
    T6 = pd.DataFrame(rows)
    print("\n" + "=" * 150)
    print("TABLE 6  A global factor against per-connection sample-and-rescale, WITHIN each mask")
    print("  the random mask has no sampled counterpart: sampling IS the informative mask, and a random mask sampled at")
    print("  uniform probability is the uniform-p control, which appears in Table 8 of the main statistics instead")
    print("  NOTE the sign: a global factor LIFTS the unbiased random mask and destroys both score-selected masks")
    print("=" * 150)
    print(T6.to_string(index=False))
    out["control_table6_global_factor_by_mask"] = T6

    # ---- 7: why dm1 and anti are excluded ------------------------------
    ex = pool(f"{B.replace('revised8_12k', 'dm1anti_12k')}pruning_p10_80/*.csv") \
        if glob.glob(f"{B.replace('revised8_12k', 'dm1anti_12k')}pruning_p10_80/*.csv") \
        else pool(f"{R}/task_preservation_tanh_h512_modcog_dm1anti_12k_seqbest_pruning_p10_80/*.csv")
    ARMS5 = {"noise-prune + capped": "simulation_noise_prune_capped_rescale",
             "noise-prune + rescale": "simulation_noise_prune_rescale",
             "OBS": "obs_compensated", "magnitude": "l1_unstructured",
             "random": "random_unstructured"}
    rows = []
    for task in sorted(ex.t.unique()):
        sub = ex[ex.t == task]
        cur = {lab: sub[sub.strategy == st].groupby(["pct", "source_network_seed"]).ret.mean()
               .reset_index().groupby("pct").ret.mean() for lab, st in ARMS5.items()}
        tbl = pd.DataFrame(cur)
        for pct in tbl.index:
            rows.append(dict(task=task, sparsity=f"{pct}%",
                             **{k: round(float(v), 4) for k, v in tbl.loc[pct].items()},
                             spread=round(float(tbl.loc[pct].max() - tbl.loc[pct].min()), 4)))
    # The retained-task reference comes from paper_stats.build_cube, which is
    # the only place the eight arms are routed to their correct pools -- an
    # earlier version read them out of whichever suite happened to be loaded
    # here and silently produced NaN for the three deterministic arms.
    from paper_stats import ARMS as PAPER_ARMS, build_cube as paper_cube
    pc = paper_cube()
    ptl = pc.task_level()
    NAME = {"noise-prune + capped": "noise-prune + capped rescale",
            "noise-prune + rescale": "noise-prune + rescale",
            "OBS": "OBS", "magnitude": "magnitude", "random": "random"}
    rt = pd.DataFrame({lab: {int(q): float(ptl[:, pc.m(full), qi].mean())
                             for qi, q in enumerate(pc.qs)}
                       for lab, full in NAME.items()})
    for pct in rt.index:
        rows.append(dict(task="8 retained tasks (pooled)", sparsity=f"{pct}%",
                         **{k: round(float(v), 4) for k, v in rt.loc[pct].items()},
                         spread=round(float(rt.loc[pct].max() - rt.loc[pct].min()), 4)))
    T7 = pd.DataFrame(rows)
    print("\n" + "=" * 140)
    print("TABLE 7  Why dm1 and anti are excluded: the methods stop separating on them")
    print("  'spread' is best minus worst of the five methods -- the benchmark's resolving power on that task")
    print("=" * 140)
    print(T7.to_string(index=False))
    out["control_table7_excluded_tasks"] = T7
    print("\n  max spread:  " + ",  ".join(
        f"{t} {T7[T7.task == t].spread.max():.3f}" for t in T7.task.unique()))

    for name, df in out.items():
        df.to_csv(HERE / f"{name}.csv", index=False)
    print(f"\nwritten to {HERE}/control_table*.csv")


if __name__ == "__main__":
    main()
