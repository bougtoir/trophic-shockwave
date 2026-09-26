"""Residual spatial dependence diagnostics (spec section 25).

Moran's I on residuals of the primary carnivore event-study model,
plus a semivariogram summary (binned). Output results/diagnostics/.
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from src.utils.common import DATA, RESULTS, load_config


def morans_I(resid_by_site, D, band_km=10):
    sites = resid_by_site.index.values
    W = ((D.loc[sites, sites] > 0) & (D.loc[sites, sites] <= band_km)).astype(float).values
    x = resid_by_site.values - resid_by_site.mean()
    num = (W * np.outer(x, x)).sum()
    den = (x ** 2).sum()
    S0 = W.sum()
    n = len(x)
    return (n / S0) * num / den if den > 0 and S0 > 0 else np.nan


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    out = []
    for sp in cfg["species"]["carnivores_primary"]:
        p = panel[panel.species == sp].dropna(subset=["encounter_rate"])
        if len(p) < 500:
            continue
        m = smf.ols("encounter_rate ~ C(site) + C(bin)", data=p).fit()
        p["resid"] = m.resid
        rsite = p.groupby("site").resid.mean()
        for band in [6, 10, 15]:
            out.append({"species": sp, "band_km": band,
                        "morans_I": morans_I(rsite.reindex(D.index).dropna(), D, band)})
        # binned semivariogram
        r = rsite.reindex(D.index).dropna()
        dd = D.loc[r.index, r.index].values
        iu = np.triu_indices(len(r), 1)
        dist = dd[iu]
        gam = 0.5 * (r.values[iu[0]] - r.values[iu[1]]) ** 2
        bins = np.linspace(0, dist.max(), 15)
        vg = pd.DataFrame({"dist": dist, "gamma": gam})
        vg["band"] = pd.cut(vg.dist, bins)
        vg.groupby("band", observed=True).agg(d=("dist", "mean"),
                                              gamma=("gamma", "mean")).to_csv(
            RESULTS / "diagnostics" / f"variogram_{sp}.csv")
    pd.DataFrame(out).to_csv(RESULTS / "diagnostics" / "morans_I.csv", index=False)
    print(out)


if __name__ == "__main__":
    RESULTS_DIAG = RESULTS / "diagnostics"
    RESULTS_DIAG.mkdir(parents=True, exist_ok=True)
    main()
