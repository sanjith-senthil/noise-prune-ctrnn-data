#!/usr/bin/env python3
"""The same contrast table at H = 1024, which reviewer 1's comment 4 asks for.

    .venv/bin/python analysis/paper_statistics/h1024_contrasts.py

Every other statistics file in this directory is H = 512. The size result
existed only as a figure, so there were no effect sizes, intervals or p-values
for it anywhere -- a reader could see that the curves separate but not whether
any separation survives inference. This produces the identical table at the
larger size: 8 paired task-level differences, 10,000 two-stage hierarchical
resamples (tasks, then seeds within task), exact sign-flip over all 256
assignments, normalised AUC over 10-80%.

TWO THINGS DIFFER FROM H = 512, AND BOTH ARE STATED IN THE OUTPUT.

The cap is 3.695 for BOTH scores here. That constant was derived from the noise
score's sweep at H = 512, and no sweep was run at H = 1024, so neither arm is at
a cap derived at this size. The H = 512 tables give each score its own swept
constant; these do not, and the capped contrasts should be read as "both scores
at the H = 512 noise-score constant" rather than as a matched comparison.

The sampling arms come from the renormalisation re-run and the deterministic
arms from the original size suite. That suite also holds PRE-correction copies
of the sampling arms, so each arm is pinned to the pool it is allowed to come
from -- reading both and de-duplicating on value silently averages corrected and
uncorrected rows, which moved noise-prune's 40% retention from 0.934 to 0.843
when this was first attempted in the size figure.
"""
from __future__ import annotations

import glob
import itertools
import pathlib
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from paper_stats import (CAP_CONSTANT, Cube, GRID, N_BOOT, auc, contrast,  # noqa: E402
                         make_surrogates)

REPO = pathlib.Path(__file__).resolve().parents[2]
NEED = {"strategy", "amount", "post_acc_sequence", "pre_acc_sequence", "task",
        "source_network_seed"}

SAMPLING = {"simulation_noise_prune_rescale",
            "simulation_noise_prune_capped_rescale",
            "simulation_noise_prune_magnitude_rescale",
            "simulation_noise_prune_capped_magnitude_rescale"}
RENORM = ["results/task_preservation_tanh_h1024_modcog_revised8_12k_seqbest_renorm_p10_80/*.csv",
          "results/task_preservation_tanh_h1024_modcog_revised8_12k_seqbest_renorm_rest_p10_80/*.csv"]
LEGACY = ["results/task_preservation_tanh_h1024_modcog_revised8_12k_seqbest_revised_task_only_p10_40/*.csv",
          "results/task_preservation_tanh_h1024_modcog_revised8_12k_seqbest_revised_task_only_p50_80/*.csv"]

ARMS = {
    "noise-prune + capped rescale": "simulation_noise_prune_capped_rescale",
    "noise-prune + rescale": "simulation_noise_prune_rescale",
    "noise-prune score only": "simulation_noise_prune_mask_only",
    "magnitude + capped rescale": "simulation_noise_prune_capped_magnitude_rescale",
    "magnitude + rescale": "simulation_noise_prune_magnitude_rescale",
    "magnitude": "l1_unstructured",
    "OBS": "obs_compensated",
    "random": "random_unstructured",
}


def _pool(pats):
    fr = []
    for pat in pats:
        for f in sorted(glob.glob(str(REPO / pat))):
            cols = set(pd.read_csv(f, nrows=0).columns)
            if not NEED <= cols:
                continue
            extra = {c for c in ("prune_rescale_cap_mode", "prune_rescale_cap_value",
                                 "pruning_seed") if c in cols}
            fr.append(pd.read_csv(f, usecols=list(NEED | extra), low_memory=False))
    d = pd.concat(fr, ignore_index=True)
    d = d[d.strategy != "none"].copy()
    d["t"] = d.task.str.replace("modcog:", "", regex=False)
    d["pct"] = (d.amount * 100).round().astype(int)
    d["ret"] = d.post_acc_sequence / d.pre_acc_sequence
    return d


def build() -> Cube:
    fresh = _pool(RENORM)
    fresh = fresh[fresh.strategy.isin(SAMPLING)]
    legacy = _pool(LEGACY)
    legacy = legacy[~legacy.strategy.isin(SAMPLING)]
    d = pd.concat([fresh, legacy], ignore_index=True)
    d = d.drop_duplicates(subset=["strategy", "amount", "task", "source_network_seed",
                                  "post_acc_sequence"])
    qs = np.array(sorted(GRID))
    frames = {}
    for name, strat in ARMS.items():
        s = d[d.strategy == strat]
        if strat in {"simulation_noise_prune_capped_rescale",
                     "simulation_noise_prune_capped_magnitude_rescale"}:
            s = s[np.isclose(s.prune_rescale_cap_value.astype(float), CAP_CONSTANT)]
        if s.empty:
            raise SystemExit(f"no rows for {name!r} ({strat})")
        frames[name] = s
    tasks = sorted(set.intersection(*[set(f.t) for f in frames.values()]))
    seeds = sorted(set.intersection(*[set(f.source_network_seed) for f in frames.values()]))
    names = list(frames)
    R = np.full((len(tasks), len(seeds), len(names), len(qs)), np.nan)
    D = np.zeros_like(R)
    for mi, n in enumerate(names):
        g = frames[n].groupby(["t", "source_network_seed", "pct"]).ret.agg(["mean", "size"])
        for (t, s, q), row in g.iterrows():
            if t in tasks and s in seeds and q in qs:
                R[tasks.index(t), seeds.index(s), mi, list(qs).index(q)] = row["mean"]
                D[tasks.index(t), seeds.index(s), mi, list(qs).index(q)] = row["size"]
    miss = np.argwhere(np.isnan(R))
    if len(miss):
        raise SystemExit(f"{len(miss)} empty cells, e.g. "
                         f"{[f'{tasks[a]}/{names[c]}/{qs[e]}' for a, b, c, e in miss[:4]]}")
    return Cube(R=R, tasks=tasks, seeds=seeds, methods=names, qs=qs, draws=D)


def main() -> None:
    cube = build()
    sur = make_surrogates(cube, n_boot=N_BOOT)
    tl = cube.task_level()
    A = auc(tl, cube.qs)

    per_method = pd.DataFrame([
        dict(method=n, AUC=round(float(A[:, cube.m(n)].mean()), 6),
             sd_across_tasks=round(float(A[:, cube.m(n)].std(ddof=1)), 6),
             n_tasks=len(cube.tasks))
        for n in cube.methods])

    rows = []
    for a, b in itertools.combinations(cube.methods, 2):
        for q in list(cube.qs) + [None]:
            r = contrast(cube, sur, a, b, q=q)
            rows.append(dict(
                arm_a=a, arm_b=b,
                measure=("AUC 10-80%" if q is None else f"{q}% sparsity"),
                sparsity=("AUC" if q is None else int(q)),
                effect=round(r["effect"], 6),
                lo95=round(r["lo95"], 6), hi95=round(r["hi95"], 6),
                lo68=round(r["lo68"], 6), hi68=round(r["hi68"], 6),
                boot_sd=round(r["boot_sd"], 6),
                n_tasks=r["n_tasks"], n_positive=r["n_positive"],
                signs=f"{r['n_positive']}/{r['n_tasks']}",
                p_two=round(r["p_two"], 6), p_one=round(r["p_one"], 6),
                p_floor=round(r["p_floor"], 6),
                at_floor=bool(np.isclose(r["p_two"], r["p_floor"])),
                CI95=f"[{r['lo95']:+.3f}, {r['hi95']:+.3f}]",
                sig=bool(r["lo95"] > 0 or r["hi95"] < 0)))
    con = pd.DataFrame(rows)

    print("=" * 130)
    print("H = 1024 -- AUC by method (24 networks: 8 tasks x 3 seeds)")
    print("  BOTH capped arms are at 3.695, the constant derived from the NOISE score's sweep at H = 512.")
    print("  No cap sweep was run at this size, so neither score is at a cap derived here.")
    print("=" * 130)
    print(per_method.to_string(index=False))

    head = con[(con.sparsity == "AUC") & (con.arm_a == "noise-prune + capped rescale")]
    print("\n" + "=" * 130)
    print("H = 1024 -- the capped rule against every other arm, normalised AUC 10-80%")
    print("=" * 130)
    print(head[["arm_b", "effect", "CI95", "signs", "p_two", "sig"]].to_string(index=False))

    key = [("noise-prune + capped rescale", "noise-prune + rescale", "what capping buys"),
           ("noise-prune + capped rescale", "magnitude + capped rescale", "the covariance term, both capped"),
           ("noise-prune + rescale", "magnitude + rescale", "the covariance term, neither capped")]
    print("\n" + "=" * 130)
    print("H = 1024 -- the three contrasts the size claim rests on, by sparsity")
    print("=" * 130)
    for a, b, lab in key:
        sub = con[(con.arm_a == a) & (con.arm_b == b)]
        print(f"\n{lab}:  {a}  minus  {b}")
        print(sub[["measure", "effect", "CI95", "signs", "p_two"]].to_string(index=False))

    per_method.to_csv(HERE / "h1024_auc_by_method.csv", index=False)
    con.to_csv(HERE / "h1024_contrasts.csv", index=False)
    print(f"\nwritten to {HERE}/h1024_auc_by_method.csv and h1024_contrasts.csv")


if __name__ == "__main__":
    main()
