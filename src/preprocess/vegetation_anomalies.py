"""Vegetation anomalies per spec section 11.

Seasonal expectation = site-specific mean of the index for that day-of-year band
(bin mapped to day-of-year). Z = (V - seasonal mean) / seasonal SD per site.
Output: data/processed/vegetation.parquet  (site x bin, ndvi/evi + z scores)
"""
import numpy as np
import pandas as pd
from src.utils.common import DATA, load_config, bin_index, doy_of_bin, bin_start_date


def main():
    cfg = load_config()
    veg = pd.read_parquet(DATA / "interim" / "site_vegetation.parquet")
    bin_days = cfg["temporal"]["primary_bin_days"]
    prim_buf = cfg["spatial"]["primary_buffer_m"]
    veg = veg[veg.buffer_m == prim_buf].copy()
    veg["bin"] = bin_index(veg["date"], bin_days, ORIGIN := "2010-01-01")
    # combine terra/aqua: mean per site-bin (same 16-day grid offset by ~8d)
    g = veg.groupby(["site", "bin"]).agg(
        ndvi=("ndvi", "mean"), evi=("evi", "mean"), n_pix=("n_pix", "sum")).reset_index()
    bins_per_year = int(np.ceil(366 / bin_days))  # 23
    g["biny"] = g["bin"] % bins_per_year
    # Seasonal climatology: per-site harmonic fit V = a + b*cos + c*sin + d*cos2 + e*sin2
    # (raw site x bin-of-year means have only ~3 obs -> z bounded by sqrt(n-1); harmonic
    # smoothing gives a usable baseline, see DECISION_LOG)
    def harmonic_fit(df, idx):
        x = 2 * np.pi * df.biny / bins_per_year
        X = np.column_stack([np.ones(len(df)), np.cos(x), np.sin(x),
                             np.cos(2 * x), np.sin(2 * x)])
        y = df[idx].values
        ok = ~np.isnan(y)
        if ok.sum() < 8:
            return pd.Series(np.nan, index=df.index)
        beta, *_ = np.linalg.lstsq(X[ok], y[ok], rcond=None)
        return pd.Series(X @ beta, index=df.index)

    out_rows = []
    for site, sv in g.groupby("site"):
        sv = sv.copy()
        for idx in ["ndvi", "evi"]:
            sv[f"{idx}_seas_mean"] = harmonic_fit(sv, idx)
            resid = sv[idx] - sv[f"{idx}_seas_mean"]
            sv[f"{idx}_seas_sd"] = resid.std()
            sv[f"{idx}_anom"] = resid
            sv[f"{idx}_z"] = resid / resid.std() if resid.std() > 0 else np.nan
        out_rows.append(sv)
    g = pd.concat(out_rows)
    out = DATA / "processed" / "vegetation.parquet"
    g.to_parquet(out, index=False)
    print("vegetation:", g.shape, "sites", g.site.nunique())


if __name__ == "__main__":
    main()
