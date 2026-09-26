"""Distributed-lag spatial model (spec section 14).

For each connectivity matrix W in {geo, neighbor, movement, habitat}:
  exposure_jt = sum_i W_ij * shock_i,t-k  for lags k in [0,12]
  Y_jt (carnivore encounter_rate) ~ sum_k b_k exposure_k + ndvi_jt + precip + t2m
                                  + C(site) + C(bin), cluster by site.
Model comparison via AIC and 5-fold blocked (by site) out-of-sample RMSE.
Output: results/tables/connectivity_model_comparison.csv, distlag_{W}.csv
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from src.utils.common import DATA, RESULTS, load_config

K_MAX = 12


def build_exposure(shock, W, sites):
    sidx = {s: i for i, s in enumerate(sites)}
    S = shock.set_index(["site", "bin"]).shock.unstack("site").reindex(columns=sites).fillna(0)
    bins = S.index.values
    Sarr = S.values  # bins x sites
    exp = {}
    for k in range(1, K_MAX + 1):
        Ek = Sarr @ W.T if False else (Sarr * 1.0)  # placeholder
        # exposure_j,t = sum_i W_ij S_i,t-k -> (S_{t-k} dot W) row? W_ij: i->j.
        Ek = np.zeros_like(Sarr)
        Ek[k:, :] = Sarr[:-k, :] @ W  # bins x sites: sum_i S_i,t-k * W_i,j
        exp[k] = Ek
    return bins, exp


def fit_W(name, W, panel, shock, veg, clim, species, sites):
    bins, exp = build_exposure(shock, W, sites)
    expo = pd.DataFrame({"bin": np.repeat(bins, len(sites)),
                         "site": list(sites) * len(bins)})
    for k in range(1, K_MAX + 1):
        expo[f"E{k}"] = exp[k].ravel()
    p = panel[panel.species.isin(species)].groupby(["site", "bin"]).encounter_rate.mean().reset_index()
    d = p.merge(expo, on=["site", "bin"]).merge(
        veg[["site", "bin", "evi_z"]], on=["site", "bin"], how="left").merge(
        clim, on=["site", "bin"], how="left")
    d = d.dropna(subset=["encounter_rate"])
    rhs = " + ".join([f"E{k}" for k in range(1, K_MAX + 1)]) + \
          " + evi_z + precip_sum + t2m_mean + C(site) + C(bin)"
    m = smf.ols(f"encounter_rate ~ {rhs}", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d.site})
    # blocked CV by site
    rng = np.random.default_rng(0)
    site_ids = np.array(sorted(d.site.unique()))
    rng.shuffle(site_ids)
    folds = np.array_split(site_ids, 5)
    rmses = []
    for f in folds:
        tr, te = d[~d.site.isin(f)], d[d.site.isin(f)]
        try:
            mm = smf.ols(f"encounter_rate ~ {rhs}", data=tr).fit()
            pred = mm.predict(te)
            rmses.append(np.sqrt(np.mean((te.encounter_rate - pred) ** 2)))
        except Exception:
            rmses.append(np.nan)
    rows = [{"lag": k, "coef": m.params.get(f"E{k}", np.nan),
             "se": m.bse.get(f"E{k}", np.nan)} for k in range(1, K_MAX + 1)]
    coef = pd.DataFrame(rows); coef["W"] = name
    return {"W": name, "aic": m.aic, "bic": m.bic,
            "cv_rmse": np.nanmean(rmses), "sum_coef": np.nansum(coef.coef),
            "n": len(d)}, coef


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    shock = pd.read_parquet(DATA / "processed" / "site_shock.parquet")
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    # climate 16-day aggregates
    clim = pd.read_parquet(DATA / "interim" / "site_climate.parquet")
    clim["bin"] = ((clim.time - pd.Timestamp("2010-01-01")).dt.days //
                   cfg["temporal"]["primary_bin_days"])
    clim = clim.groupby(["site", "bin"]).agg(precip_sum=("precip", "sum"),
                                           t2m_mean=("t2m", "mean")).reset_index()
    species = cfg["species"]["carnivores_primary"]
    results, coefs = [], []
    for name in ["geo", "neighbor", "movement", "habitat"]:
        W = np.load(DATA / "processed" / f"W_{name}.npy")
        try:
            r, c = fit_W(name, W, panel, shock, veg, clim, species, sites)
            results.append(r); coefs.append(c)
            print(name, "aic=%.1f cv_rmse=%.5f" % (r["aic"], r["cv_rmse"]))
        except Exception as e:
            print("fit", name, "FAILED:", repr(e)[:150])
    if results:
        pd.DataFrame(results).to_csv(RESULTS / "tables" / "connectivity_model_comparison.csv",
                                     index=False)
        pd.concat(coefs).to_csv(RESULTS / "tables" / "distlag_coefficients.csv", index=False)


if __name__ == "__main__":
    main()
