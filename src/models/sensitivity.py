"""Sensitivity + negative controls (revised spec).

- leave-one-year-out (2010/2011/2012) and leave-one-shock-bin-out for the
  key contrasts (wildebeest event-study post-mean lags 2-8; eco-model eco_sum).
- Permuted connectivity: permute site labels of W_movement (preserves degree
  sequence), refit carnivore eco model -> null distribution of eco_sum.
- Redistribution (Test B): recipient herbivore anomaly on geo exposure,
  cross-fitted movement exposure, and matched-distance negative control.
- Temporal placebo: event-study with fake shock dates shifted earlier.
Output: results/tables/{sensitivity_loyo_loso.csv, permuted_network.csv,
redistribution.csv, placebo.csv}
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from src.utils.common import DATA, RESULTS, load_config
from src.models.estimation import (load_all, stack_events, fit_es, carnivore_model,
                                   expo_matrix, K_EXPO, LAG_MIN, LAG_MAX, TAB)

POST = (2, 8)


def post_mean(panel, evs, species):
    d = stack_events(panel, evs, species)
    c, m = fit_es(d)
    return c[c.lag.between(*POST)].coef.mean()


def main():
    cfg = load_config()
    panel, ev, veg, clim, shock = load_all()
    ev = ev[~ev.regionwide].reset_index(drop=True)
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    Wg = np.load(DATA / "processed" / "W_geo.npy")
    # cross-fit connectivity: split time at median shock bin
    from src.preprocess.connectivity import movement_W
    split_bin = ev.bin_start.median()
    WA = movement_W(panel, sites, mask_bins=range(int(split_bin), 76))
    WB = movement_W(panel, sites, mask_bins=range(0, int(split_bin)))
    rng = np.random.default_rng(cfg["seed"])

    def set_shock(evs):
        sub = shock.copy()
        ok = set()
        for r in evs.itertuples():
            for b in range(r.bin_start, r.bin_end + 1):
                ok.add((r.site, b))
        sub["shock"] = sub.shock & sub.apply(lambda rr: (rr.site, rr.bin) in ok, axis=1)
        return sub
    TAB.mkdir(parents=True, exist_ok=True)
    sens = []

    # ---------- leave-one-year-out ----------
    for yr in sorted(ev.year.unique()):
        sub = ev[ev.year != yr]
        for sp in cfg["species"]["herbivores_primary"]:
            try:
                v = post_mean(panel, sub, sp)
                sens.append({"analysis": f"LOYO_{yr}", "species": sp,
                             "stat": "post_mean_2_8", "value": v})
            except Exception as e:
                sens.append({"analysis": f"LOYO_{yr}", "species": sp,
                             "stat": "error", "value": np.nan})
        # eco model (cross-fitted)
        try:
            res, m, d = carnivore_model(panel, set_shock(sub), veg, clim, None, sites,
                                        cfg["species"]["carnivores_primary"], "eco",
                                        ev=sub, WA=WA, WB=WB, split_bin=split_bin)
            sens.append({"analysis": f"LOYO_{yr}", "species": "carnivores",
                         "stat": "eco_sum", "value": res["eco_sum"]})
        except Exception as e:
            print("loyo eco fail", yr, repr(e)[:120])

    # ---------- leave-one-shock-bin-out ----------
    for g in sorted(ev.bin_start.unique()):
        sub = ev[ev.bin_start != g]
        try:
            v = post_mean(panel, sub, "wildebeest")
            sens.append({"analysis": f"LOSBO_bin{g}", "species": "wildebeest",
                         "stat": "post_mean_2_8", "value": v})
        except Exception:
            pass
    pd.DataFrame(sens).to_csv(TAB / "sensitivity_loyo_loso.csv", index=False)
    print("sensitivity done", len(sens))

    # ---------- permuted network null ----------
    out = []
    sub_shock = set_shock(ev)
    We_full = movement_W(panel, sites)
    try:
        res_obs, _, _ = carnivore_model(panel, sub_shock, veg, clim, We_full, sites,
                                        cfg["species"]["carnivores_primary"], "eco")
        obs = res_obs["eco_sum"]
    except Exception as e:
        obs = np.nan; print("obs eco fail", repr(e)[:120])
    perm_sums = []
    for b in range(cfg["models"]["permutations_network"]):
        P = rng.permutation(len(sites))
        Wp = We_full[P][:, P]  # relabel sites: preserves degree distribution
        try:
            rp, _, _ = carnivore_model(panel, sub_shock, veg, clim, Wp, sites,
                                       cfg["species"]["carnivores_primary"], "perm")
            perm_sums.append(rp["eco_sum"])
        except Exception:
            pass
    perm_sums = np.array(perm_sums)
    out.append({"stat": "observed_eco_sum", "value": obs})
    out.append({"stat": "perm_mean", "value": np.nanmean(perm_sums)})
    out.append({"stat": "perm_sd", "value": np.nanstd(perm_sums)})
    out.append({"stat": "perm_p_gt_obs", "value": float(np.nanmean(perm_sums >= obs))})
    pd.DataFrame(out).to_csv(TAB / "permuted_network.csv", index=False)
    print("permutation done")

    # ---------- redistribution (Test B) ----------
    foc = cfg["species"]["herbivores_primary"]
    rate = panel[panel.species.isin(foc)].groupby(["site", "bin"]).encounter_rate.mean().unstack("site")
    seas = rate.groupby(rate.index % 23).transform("mean")
    anom = rate - seas
    sidx = {s: i for i, s in enumerate(sites)}
    rows = []
    for r in ev.itertuples():
        i = r.site
        if i not in sidx:
            continue
        for j in sites:
            if j == i:
                continue
            dH = anom.loc[r.bin_start + POST[0]: r.bin_start + POST[1], j].mean()
            rows.append({"src": i, "dst": j, "bin_start": r.bin_start, "grp": r.bin_start,
                         "year": r.year, "dH": dH,
                         "dist_km": D.loc[i, j],
                         "w_geo": Wg[sidx[i], sidx[j]],
                         "w_eco": We_full[sidx[i], sidx[j]],
                         "isolated": r.isolated})
    df = pd.DataFrame(rows)
    # matched-distance negative control: for each (src,dst) pair, the expected eco
    # weight among equidistant sites -> residual eco effect
    df["w_eco"] = [We_full[sidx[r.src], sidx[r.dst]] for r in df.itertuples()]
    df["dist_bin"] = pd.cut(df.dist_km, np.quantile(df.dist_km, np.linspace(0, 1, 21)))
    df["w_eco_dm"] = df.groupby("dist_bin", observed=True).w_eco.transform(
        lambda x: x - x.mean())
    df = df.dropna(subset=["dH", "dist_km", "w_eco", "w_eco_dm"])
    m = smf.ols("dH ~ dist_km + w_eco + w_eco_dm", data=df).fit(
        cov_type="cluster", cov_kwds={"groups": df.grp})
    df.to_csv(TAB / "redistribution.csv", index=False)
    pd.DataFrame([{"term": t, "coef": m.params[t], "se": m.bse[t], "p": m.pvalues[t]}
                  for t in m.params.index]).to_csv(TAB / "redistribution_fit.csv", index=False)
    print("redistribution done", m.params.to_dict())

    # ---------- temporal placebo ----------
    prow = []
    for shift in cfg["models"]["placebo_shifts"]:
        evp = ev.copy(); evp["bin_start"] = evp.bin_start + shift
        evp = evp[evp.bin_start.between(-evp.bin_start.min() + 8, 60)]
        try:
            v = post_mean(panel, evp, "wildebeest")
            prow.append({"shift_bins": shift, "species": "wildebeest",
                         "post_mean_2_8": v})
        except Exception as e:
            print("placebo fail", shift, repr(e)[:100])
    pd.DataFrame(prow).to_csv(TAB / "placebo.csv", index=False)


if __name__ == "__main__":
    main()
