#!/usr/bin/env python3
"""Two-stage random-effects meta-analysis for the within-host contrast (P1).

Stage 1 (per study): paired effect d_i = mean over animals of [ALR_hindgut - ALR_foregut],
where ALR = log((acsB + c) / (mcrA + c)), with variance s_i^2 / n_i.
Stage 2: random-effects pooling (DerSimonian-Laird and REML tau^2), Hartung-Knapp-Sidik-Jonkman
CI, heterogeneity (Q, I^2, tau^2), 95% prediction interval, and leave-one-out.

This mirrors metafor::rma(method="REML", test="knha") so results can be cross-checked in R.
Run `python random_effects.py --selftest` to exercise the code on SIMULATED data (no biology).
"""
import argparse, json, sys
import numpy as np
import pandas as pd
from scipy import stats, optimize


def alr(num, den, c_num, c_den):
    return np.log((np.asarray(num) + c_num) / (np.asarray(den) + c_den))


def half_min_nonzero(x):
    x = np.asarray(x, float); nz = x[x > 0]
    return 0.5 * nz.min() if nz.size else 1e-6


def study_effects(df):
    """df columns: study, animal, site in {foregut, hindgut}, acsB, mcrA (copies per GE)."""
    c_a, c_m = half_min_nonzero(df["acsB"]), half_min_nonzero(df["mcrA"])
    df = df.assign(alr=alr(df["acsB"], df["mcrA"], c_a, c_m))
    w = df.pivot_table(index=["study", "animal"], columns="site", values="alr").dropna()
    w["diff"] = w["hindgut"] - w["foregut"]
    g = w.groupby(level="study")["diff"]
    out = pd.DataFrame({"yi": g.mean(), "n": g.size(), "sd": g.std(ddof=1)})
    out["vi"] = out["sd"] ** 2 / out["n"]
    return out.reset_index()


def tau2_dl(y, v):
    w = 1 / v; ybar = np.sum(w * y) / np.sum(w)
    Q = np.sum(w * (y - ybar) ** 2); k = len(y)
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, (Q - (k - 1)) / C), Q


def tau2_reml(y, v):
    def nll(t2):
        w = 1 / (v + t2); mu = np.sum(w * y) / np.sum(w)
        return 0.5 * (np.sum(np.log(v + t2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))
    upper = max(10 * np.var(y, ddof=1), 1e-6)
    r = optimize.minimize_scalar(nll, bounds=(0, upper), method="bounded")
    return float(r.x)


def pool(y, v, method="REML", knha=True):
    y, v = np.asarray(y, float), np.asarray(v, float); k = len(y)
    t2_dl, Q = tau2_dl(y, v)
    t2 = tau2_reml(y, v) if method == "REML" else t2_dl
    w = 1 / (v + t2); mu = np.sum(w * y) / np.sum(w)
    if knha and k > 1:
        q = np.sum(w * (y - mu) ** 2) / (k - 1)
        se = np.sqrt(max(q, 1e-12) / np.sum(w)); crit = stats.t.ppf(0.975, k - 1)
        p = 2 * stats.t.sf(abs(mu / se), k - 1)
    else:
        se = np.sqrt(1 / np.sum(w)); crit = stats.norm.ppf(0.975); p = 2 * stats.norm.sf(abs(mu / se))
    vbar = (k - 1) / (np.sum(1 / v) - np.sum(1 / v ** 2) / np.sum(1 / v)) if k > 1 else np.nan
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    pi_crit = stats.t.ppf(0.975, max(k - 2, 1))
    pi = (mu - pi_crit * np.sqrt(t2 + se ** 2), mu + pi_crit * np.sqrt(t2 + se ** 2))
    return {"k": k, "mu": mu, "se": se, "ci_lo": mu - crit * se, "ci_hi": mu + crit * se, "p": p,
            "tau2": t2, "Q": Q, "Q_p": stats.chi2.sf(Q, k - 1) if k > 1 else np.nan, "I2": I2,
            "pi_lo": pi[0], "pi_hi": pi[1], "typical_v": vbar}


def leave_one_out(es):
    rows = []
    for s in es["study"]:
        sub = es[es["study"] != s]
        r = pool(sub["yi"], sub["vi"]); r["omitted"] = s; rows.append(r)
    return pd.DataFrame(rows)


def simulate(k=8, n_range=(6, 20), mu=1.0, tau=0.4, seed=1):
    """SIMULATED paired data with a known true effect mu (log-ratio units). For code testing only."""
    rng = np.random.default_rng(seed); rows = []
    for i in range(k):
        theta = rng.normal(mu, tau); n = rng.integers(*n_range)
        for a in range(n):
            base_m = rng.lognormal(0, 1); base_a = rng.lognormal(-1, 1)
            eff = np.exp(rng.normal(theta, 0.8))
            rows += [{"study": f"S{i}", "animal": a, "site": "foregut", "acsB": base_a, "mcrA": base_m},
                     {"study": f"S{i}", "animal": a, "site": "hindgut", "acsB": base_a * eff, "mcrA": base_m}]
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", help="TSV with study, animal, site, acsB, mcrA")
    ap.add_argument("--out", default="p1_meta.json")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        recov = []
        for seed in range(200):
            es = study_effects(simulate(seed=seed)); r = pool(es["yi"], es["vi"])
            recov.append((r["mu"], r["ci_lo"] <= 1.0 <= r["ci_hi"]))
        mus = np.array([m for m, _ in recov]); cover = np.mean([c for _, c in recov])
        es = study_effects(simulate(seed=0)); r = pool(es["yi"], es["vi"])
        print("SELFTEST on SIMULATED data (true mu = 1.0, tau = 0.4, k = 8) - not a biological result")
        print(f"mean estimated mu over 200 sims = {mus.mean():.3f}; HKSJ 95% CI coverage = {cover:.3f}")
        print("example pooled fit:", json.dumps({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}))
        print("leave-one-out mu range:", tuple(round(x, 3) for x in leave_one_out(es)["mu"].agg(["min", "max"])))
        ok = abs(mus.mean() - 1.0) < 0.1 and 0.90 <= cover <= 0.99
        print("SELFTEST", "PASSED" if ok else "FAILED"); sys.exit(0 if ok else 1)
    df = pd.read_csv(a.input, sep="\t"); es = study_effects(df)
    res = {"pooled": pool(es["yi"], es["vi"]), "pooled_DL": pool(es["yi"], es["vi"], method="DL"),
           "study_effects": es.to_dict("records"), "leave_one_out": leave_one_out(es).to_dict("records")}
    json.dump(res, open(a.out, "w"), indent=2, default=float)


if __name__ == "__main__":
    main()
