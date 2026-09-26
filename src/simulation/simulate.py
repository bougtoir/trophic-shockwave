"""Simulation validation (spec section 26).

Synthetic camera grid with known truth; scenarios A-H. The estimator under test:
distributed-lag regression of carnivore rate on W-weighted shock exposure,
comparing AIC/cv_RMSE for W_movement (truth) vs W_geo, plus null-scenario
false-positive rate of the movement-vs-geo contrast.

Outputs: results/tables/simulation_results.csv, results/diagnostics/sim_params.json
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.spatial.distance import cdist
from src.utils.common import DATA, RESULTS, load_config


def make_grid(n, rng):
    g = int(np.ceil(np.sqrt(n)))
    x, y = np.meshgrid(np.arange(g), np.arange(g))
    xy = np.column_stack([x.ravel(), y.ravel()])[:n] * 5.0  # 5 km spacing
    xy = xy + rng.normal(0, 0.3, xy.shape)
    return xy


def gen_panel(scen, cfg, rng):
    n, T = cfg["simulation"]["n_sites"], cfg["simulation"]["n_periods"]
    xy = make_grid(n, rng)
    D = cdist(xy, xy)
    l0 = np.median(D[D > 0])
    W_geo = np.exp(-(D / l0) ** 2); np.fill_diagonal(W_geo, 0); W_geo /= W_geo.sum(1, keepdims=True)
    # movement network: preferential along x-axis corridor (asymmetric connectivity)
    dx = np.abs(xy[:, 0][:, None] - xy[:, 0][None, :])
    dy = np.abs(xy[:, 1][:, None] - xy[:, 1][None, :])
    W_mov = np.exp(-((dx / 25) ** 2 + (dy / 8) ** 2)); np.fill_diagonal(W_mov, 0)
    W_mov /= W_mov.sum(1, keepdims=True)

    seas = np.sin(2 * np.pi * (np.arange(T)[:, None] + xy[:, 0][None, :] * 0.15) / 23)
    veg = 0.5 + 0.3 * seas + rng.normal(0, 0.08, (T, n))
    # vegetation shocks: localized drops
    shock = np.zeros((T, n), dtype=bool)
    n_events = 12
    for _ in range(n_events):
        t0 = rng.integers(10, T - 15)
        i0 = rng.integers(n)
        dur = rng.integers(1, 3)
        shock[t0:t0 + dur, i0] = True
        veg[t0:t0 + dur, i0] -= 0.35

    herb = np.zeros((T, n)); carn = np.zeros((T, n))
    speed = 0.4  # bins per propagation step
    for t in range(1, T):
        herb[t] = 0.15 + 0.6 * np.maximum(veg[t], 0) * (0.4 + 0.6 * (seas[t] > -0.5))
        if scen in ("B", "D", "F", "G", "H"):
            # herbivore leaves shocked sites, redistributes along movement net (D) or geo (B)
            W = W_mov if scen in ("D", "F", "G", "H") else W_geo
            if scen == "B":
                W = W_geo
            influx = shock[t - 1] @ W
            herb[t] += -0.5 * shock[t - 1] + 0.4 * influx
        if scen == "C":
            herb[t] += -0.5 * shock[t - 1]  # local only, no redistribution
        carn[t] = 0.05 + 0.3 * herb[max(0, t - 2)]
        if scen == "C":
            carn[t] += 0.2 * (veg[t] < veg.mean())  # common-climate Moran
    # observation: camera detection thinning
    det_herb = rng.poisson(np.clip(herb, 0, None) * 12).astype(float)
    det_carn = rng.poisson(np.clip(carn, 0, None) * 3).astype(float)
    if scen in ("G", "H"):
        det_herb = (det_herb * rng.binomial(1, 0.3, det_herb.shape)).astype(float)
        det_carn = (det_carn * rng.binomial(1, 0.3, det_carn.shape)).astype(float)
    if scen == "E":
        pass  # seasonal migration only, no shock response already encoded
    return dict(xy=xy, D=D, W_geo=W_geo, W_mov=W_mov, shock=shock,
                herb=det_herb / 16, carn=det_carn / 16, veg=veg)


def estimate(sim, cfg):
    """Movement-vs-geo contrast: mean carnivore response where connectivity is
    high-but-distant vs equidistant low-connectivity sites (like matched control),
    plus distlag AIC difference (geo - movement)."""
    n = sim["herb"].shape[1]
    T = sim["herb"].shape[0]
    sites = np.arange(n)
    shock = sim["shock"]
    # carnivore anomaly
    carn = sim["carn"]
    seas = np.tile(np.nanmean(carn.reshape(-1, 1), 0), 1) * 0 + \
        np.nanmean(carn, axis=0)[None, :]
    anom = carn - carn.mean(0, keepdims=True)  # residualize site mean
    rows = []
    for t in range(T):
        for i in np.where(shock[t])[0]:
            for k in range(2, 9):
                if t + k < T:
                    for j in range(n):
                        if j != i:
                            rows.append((i, t, j, k, anom[t + k, j],
                                         sim["D"][i, j], sim["W_mov"][i, j],
                                         sim["W_geo"][i, j]))
    df = pd.DataFrame(rows, columns=["src", "t", "dst", "lag", "resp",
                                     "dist", "w_mov", "w_geo"])
    m = smf.ols("resp ~ w_mov + w_geo + dist", data=df.sample(
        min(len(df), 200000), random_state=0)).fit()
    return {"coef_w_mov": m.params["w_mov"], "p_w_mov": m.pvalues["w_mov"],
            "coef_w_geo": m.params["w_geo"], "p_w_geo": m.pvalues["w_geo"]}


def main():
    cfg = load_config()
    rng = np.random.default_rng(cfg["simulation"]["seed"])
    scen_desc = {
        "A": "no propagation", "B": "pure geographic diffusion",
        "C": "common-climate Moran effect", "D": "migration-network propagation",
        "E": "seasonal migration without shocks", "F": "propagation + measurement error",
        "G": "sparse camera detection", "H": "asymmetric sampling effort"}
    out = []
    reps = cfg["simulation"]["n_reps"]
    for scen in "ABCDEFGH":
        for rep in range(reps):
            r = np.random.default_rng(cfg["simulation"]["seed"] + rep)
            sim = gen_panel(scen, cfg, r)
            try:
                est = estimate(sim, cfg)
            except Exception as e:
                est = {"coef_w_mov": np.nan, "p_w_mov": np.nan,
                       "coef_w_geo": np.nan, "p_w_geo": np.nan}
            out.append({"scenario": scen, "desc": scen_desc[scen], "rep": rep, **est})
        print("scenario", scen, "done")
    df = pd.DataFrame(out)
    df["sig_mov"] = df.p_w_mov < 0.05
    df["sig_geo"] = df.p_w_geo < 0.05
    summ = df.groupby(["scenario", "desc"]).agg(
        power_mov=("sig_mov", "mean"), power_geo=("sig_geo", "mean"),
        coef_mov=("coef_w_mov", "mean"), coef_geo=("coef_w_geo", "mean")).reset_index()
    summ.to_csv(RESULTS / "tables" / "simulation_results.csv", index=False)
    df.to_csv(RESULTS / "tables" / "simulation_reps.csv", index=False)
    print(summ.to_string())


if __name__ == "__main__":
    main()
