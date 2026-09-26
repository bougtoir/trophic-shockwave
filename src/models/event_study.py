"""Event-study models (spec sections 13, 20).

Test A: vegetation shock -> herbivore event-time response at shocked site.
Test C: remote carnivore response at connected recipient sites vs controls.
Y_it = sum_k beta_k * I(event time = k) + site FE + bin FE + eps, cluster by site.

Outputs results/tables/event_study_{species}.csv with coef/CI per lag,
plus results/models/testB_redistribution.csv.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from pathlib import Path
from src.utils.common import DATA, RESULTS, load_config

OUT_T = RESULTS / "tables"
OUT_M = RESULTS / "models"


def local_events():
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    return ev[(ev["index"] == "evi") & (~ev.regionwide)]


def treated_bins(ev):
    t = {}
    for _, r in ev.iterrows():
        t.setdefault(r.site, set()).add(r.bin_start)
    return t


def event_study(panel, species, ev, lag_min, lag_max, ref=-1):
    tb = treated_bins(ev)
    p = panel[(panel.species == species) & (panel.cam_days > 0)].copy()
    p = p.dropna(subset=["encounter_rate"])
    evmap = {}
    for _, r in ev.iterrows():
        evmap.setdefault(r.site, r.bin_start)
    p["ev_start"] = p.site.map(evmap)
    tre = p[p.ev_start.notna()].copy()
    tre["etime"] = tre["bin"] - tre.ev_start
    for k in range(lag_min, lag_max + 1):
        if k == ref:
            continue
        tre[f"L{k}"] = ((tre.etime == k)).astype(float)
    dummies = [f"L{k}" for k in range(lag_min, lag_max + 1) if k != ref]
    rhs = " + ".join(dummies) + " + C(site) + C(bin)"
    m = smf.ols(f"encounter_rate ~ {rhs}", data=tre).fit(
        cov_type="cluster", cov_kwds={"groups": tre.site})
    rows = []
    for k in range(lag_min, lag_max + 1):
        if k == ref:
            rows.append({"lag": k, "coef": 0.0, "lo": 0.0, "hi": 0.0, "p": np.nan})
            continue
        c = m.params.get(f"L{k}", np.nan)
        ci = m.conf_int().loc[f"L{k}"] if f"L{k}" in m.params.index else [np.nan, np.nan]
        rows.append({"lag": k, "coef": c, "lo": ci[0], "hi": ci[1],
                     "p": m.pvalues.get(f"L{k}", np.nan)})
    out = pd.DataFrame(rows)
    out["species"] = species
    out["n_site"] = tre.site.nunique()
    out["n_obs"] = len(tre)
    out["n_events"] = len(ev)
    return out


def test_B(panel, ev, cfg):
    """Herbivore redistribution: recipient response vs connectivity & distance."""
    sites = np.load(DATA / "processed" / "sites_order.npy")
    sidx = {s: i for i, s in enumerate(sites)}
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    Wmov = np.load(DATA / "processed" / "W_movement.npy")
    species = cfg["species"]["herbivores_primary"]
    p = panel[panel.species.isin(species)]
    rate = p.groupby(["site", "bin"]).encounter_rate.mean().unstack("site")
    seas = rate.groupby(rate.index % 23).transform("mean")  # crude seasonal baseline
    anom = rate - seas
    rows = []
    for _, r in ev.iterrows():
        i = r.site
        if i not in sidx:
            continue
        for j in sites:
            if j == i:
                continue
            post = anom.loc[r.bin_start + 2: r.bin_start + 8, j].mean()
            rows.append({"event_site": i, "recipient": j, "bin_start": r.bin_start,
                         "dH": post, "dist_km": D.loc[i, j],
                         "w_movement": Wmov[sidx[i], sidx[j]]})
    df = pd.DataFrame(rows)
    OUT_T.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_T / "testB_redistribution.csv", index=False)
    m = smf.ols("dH ~ dist_km + w_movement", data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df.event_site})
    m.save(OUT_M / "testB.pkl")
    return df, m


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    ev = local_events()
    OUT_T.mkdir(parents=True, exist_ok=True)
    OUT_M.mkdir(parents=True, exist_ok=True)
    res = []
    for sp in (cfg["species"]["herbivores_primary"] + cfg["species"]["herbivores_secondary"]
               + cfg["species"]["carnivores_primary"] + cfg["species"]["carnivores_secondary"]
               + cfg["species"]["negative_control_species"]):
        try:
            res.append(event_study(panel, sp, ev, cfg["temporal"]["lag_min"],
                                   cfg["temporal"]["lag_max"],
                                   cfg["models"]["event_study_reference_lag"]))
            print("event study:", sp)
        except Exception as e:
            print("event study FAILED:", sp, repr(e)[:150])
    if res:
        pd.concat(res).to_csv(OUT_T / "event_study_coefficients.csv", index=False)
    try:
        test_B(panel, ev, cfg)
        print("test B done")
    except Exception as e:
        print("test B FAILED:", repr(e)[:150])


if __name__ == "__main__":
    main()
