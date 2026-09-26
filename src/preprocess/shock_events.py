"""Shock-event construction (spec section 12).

Event = site x run of consecutive bins with evi_z <= threshold (primary index EVI;
NDVI for robustness). Region-wide events (>= regionwide_frac of sites shocked in
the same bin) are flagged separately and excluded from the principal analysis.
Output: data/processed/shock_events.parquet + site_shock indicator panel.
"""
import numpy as np
import pandas as pd
from src.utils.common import DATA, load_config


def build_events(veg, index, thresh, cfg):
    v = veg[["site", "bin", f"{index}_z"]].dropna().copy()
    v["shocked"] = v[f"{index}_z"] <= thresh
    n_sites = v.site.nunique()
    # per site runs
    evs = []
    for site, sv in v.sort_values("bin").groupby("site"):
        run_start = None
        prev = None
        for _, r in sv.iterrows():
            if r.shocked and run_start is None:
                run_start = r.bin
            if not r.shocked and run_start is not None:
                evs.append((site, run_start, prev))
                run_start = None
            prev = r.bin
        if run_start is not None:
            evs.append((site, run_start, prev))
    rows = []
    for site, b0, b1 in evs:
        seg = v[(v.site == site) & (v["bin"] >= b0) & (v["bin"] <= b1)]
        rows.append({"site": site, "bin_start": b0, "bin_end": b1,
                     "duration_bins": b1 - b0 + 1,
                     "min_z": seg[f"{index}_z"].min(),
                     "frac_sites_bin": np.nan})
    ev = pd.DataFrame(rows)
    if len(ev):
        frac = (v[v.shocked].groupby("bin").site.nunique() / n_sites)
        ev["frac_sites_bin"] = ev.bin_start.map(frac).fillna(0)
        ev["regionwide"] = ev.frac_sites_bin >= cfg["vegetation"]["regionwide_frac"]
    return ev


def main():
    cfg = load_config()
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    thr = cfg["vegetation"]["shock_z_threshold"]
    evs = []
    for index in ["evi", "ndvi"]:
        e = build_events(veg, index, thr, cfg)
        e["index"] = index
        evs.append(e)
    ev = pd.concat(evs, ignore_index=True)
    out = DATA / "processed" / "shock_events.parquet"
    ev.to_parquet(out, index=False)
    # site x bin shock indicator (primary index = evi, local events only)
    veg["shock"] = veg.evi_z <= thr
    veg[["site", "bin", "shock"]].to_parquet(
        DATA / "processed" / "site_shock.parquet", index=False)
    print("events:", len(ev), "| regionwide:", int(ev.regionwide.sum()) if len(ev) else 0,
          "| local:", int((~ev.regionwide).sum()) if len(ev) else 0)


if __name__ == "__main__":
    main()
