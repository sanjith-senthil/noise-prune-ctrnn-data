#!/usr/bin/env python3
"""Each score at its OWN constant cap -- the matched comparison.

    .venv/bin/python analysis/paper_statistics/matched_cap_comparison.py

Figure 5b caps both scores at 3.695. That number is the mean over 50-80% of each
network's median candidate-edge 1/p, measured on the NOISE score. The identical
construction applied to the MAGNITUDE score gives 3.584. Holding the *number*
fixed asks whether the noise score wins at one particular cap; holding the
*rule* fixed and letting each score set its own constant asks whether the noise
score wins when each is given its own best-calibrated version. The second is the
matched question, and it is the one a referee will ask.

Because the cap-dependent arm set runs both arms, the 3.584 pass also produced
noise-prune at 3.584. With the existing 3.695 rows that completes a 2x2 over
{score} x {cap}, so the score effect and the cap effect can be read apart
instead of being confounded.

Everything uses the paper's procedure: pruning draws averaged within a network,
then the 3 seeds within a task, then 8 tasks; intervals from 10,000 hierarchical
resamples shared across every contrast; p from the exact sign-flip over all 256
assignments, floor 2/256 = 0.0078 (marked *).
"""
from __future__ import annotations

import glob
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
pd.set_option("display.width", 260)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_stats import (CAP_CONSTANT, Cube, GRID, N_BOOT, auc, contrast,  # noqa: E402
                         make_surrogates)

REPO = str(pathlib.Path(__file__).resolve().parents[2])
R = f"{REPO}/results"
MAG_CAP = 3.584
MAG_OPT = 4.549      # magnitude's own sweep optimum (pooled q60)

SUITE = {
    "recap":  f"{R}/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_renorm_recap_p10_80/*.csv",
    "magcap": f"{R}/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_renorm_magcap_p10_80/*.csv",
    "magopt": f"{R}/task_preservation_tanh_h512_modcog_revised8_12k_seqbest_renorm_magopt_p10_80/*.csv",
}
NEED = ["strategy", "amount", "post_acc_sequence", "pre_acc_sequence", "task",
        "source_network_seed", "prune_rescale_cap_value"]
NOISE = "simulation_noise_prune_capped_rescale"
MAGN = "simulation_noise_prune_capped_magnitude_rescale"


def pool(name: str) -> pd.DataFrame:
    fr = []
    for f in sorted(glob.glob(SUITE[name])):
        cols = set(pd.read_csv(f, nrows=0).columns)
        fr.append(pd.read_csv(f, usecols=[c for c in NEED if c in cols], low_memory=False))
    d = pd.concat(fr, ignore_index=True)
    d = d[d.strategy != "none"].copy()
    d["t"] = d.task.str.replace("modcog:", "", regex=False)
    d["pct"] = (d.amount * 100).round().astype(int)
    d["ret"] = d.post_acc_sequence / d.pre_acc_sequence
    return d


def pick(name, strategy, cap):
    d = pool(name)
    s = d[(d.strategy == strategy) & np.isclose(d.prune_rescale_cap_value.astype(float), cap)]
    if s.empty:
        raise SystemExit(f"no rows: {name}/{strategy}/cap={cap}")
    return s


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
    miss = np.argwhere(np.isnan(Rm))
    if len(miss):
        raise SystemExit(f"{len(miss)} empty cells, e.g. "
                         f"{[f'{tasks[a]}/{names[c]}/{qs[e]}' for a, b, c, e in miss[:4]]}")
    return Cube(R=Rm, tasks=tasks, seeds=seeds, methods=names, qs=qs, draws=D)


NP695, NP584 = f"noise @ {CAP_CONSTANT}", f"noise @ {MAG_CAP}"
MG695, MG584 = f"magnitude @ {CAP_CONSTANT}", f"magnitude @ {MAG_CAP}"
NPOPT, MGOPT = f"noise @ {MAG_OPT}", f"magnitude @ {MAG_OPT}"

arms = {
    NP695: pick("recap", NOISE, CAP_CONSTANT),
    NP584: pick("magcap", NOISE, MAG_CAP),
    NPOPT: pick("magopt", NOISE, MAG_OPT),
    MG695: pick("recap", MAGN, CAP_CONSTANT),
    MG584: pick("magcap", MAGN, MAG_CAP),
    MGOPT: pick("magopt", MAGN, MAG_OPT),
}
cube = make_cube(arms, GRID)
sur = make_surrogates(cube, n_boot=N_BOOT)


def row(a, b, q=None):
    r = contrast(cube, sur, a, b, q=q)
    return dict(comparison=f"{a}  -  {b}",
                at=("AUC 10-80%" if q is None else f"{q}%"),
                effect=round(r["effect"], 4),
                CI95=f"[{r['lo95']:+.3f}, {r['hi95']:+.3f}]",
                tasks=f"{r['n_positive']}/{r['n_tasks']}",
                p=f"{r['p_two']:.4f}" + ("*" if np.isclose(r["p_two"], r["p_floor"]) else ""))


obs = cube.task_level()
A = auc(obs, cube.qs)
print("\n" + "=" * 96)
print("ARM MEANS, normalised AUC 10-80%")
print("=" * 96)
for i, n in enumerate(cube.methods):
    print(f"  {n:<26} {A[:, i].mean():.4f}")

print("\n" + "=" * 96)
print("EACH SCORE AT ITS OWN SWEEP OPTIMUM  (noise q50 -> 3.695;  magnitude q60 -> 4.549)")
print("=" * 96)
print(pd.DataFrame([row(NP695, MGOPT)] + [row(NP695, MGOPT, q) for q in GRID]).to_string(index=False))

print("\n" + "=" * 96)
print("AGAINST THE CURRENT FIGURE -- both at 3.695")
print("=" * 96)
print(pd.DataFrame([row(NP695, MG695)] + [row(NP695, MG695, q) for q in GRID]).to_string(index=False))

print("\n" + "=" * 96)
print("WITHIN-SCORE: what each cap does to each score")
print("=" * 96)
print(pd.DataFrame([row(MGOPT, MG695), row(MGOPT, MG584),
                    row(NP695, NPOPT), row(NP695, NP584)]).to_string(index=False))

print("\nper-task AUC at each score's own optimum (noise @ 3.695 - magnitude @ 4.549)")
ia, ib = cube.m(NP695), cube.m(MGOPT)
per = pd.DataFrame({"task": cube.tasks,
                    "noise": A[:, ia].round(4),
                    "magnitude": A[:, ib].round(4),
                    "difference": (A[:, ia] - A[:, ib]).round(4)})
print(per.to_string(index=False))
print(f"{'MEAN':>16} {A[:, ia].mean():.4f} {A[:, ib].mean():.4f} {(A[:, ia] - A[:, ib]).mean():+.4f}")
