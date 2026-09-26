"""Spatiotemporal impulse-response surface (spec section 15).

For each local shock event, collect recipient-site response anomalies (herbivore
and carnivore guilds) by distance band x lag. Bootstrap CIs over events.
Output: results/models/irf_{herb,carn}.parquet (dist_band, lag, mean, lo, hi)
"""
import numpy as np
import pandas as pd
from src.utils.common import DATA, RESULTS, load_config

DIST_EDGES = [0, 3, 6, 10, 15, 22, 35, 60]


def guild_anom(panel, species):
    p = panel[panel.species.isin(species)].groupby(["site", "bin"]).encounter_rate.mean()
    p = p.unstack("site")
    seas = p.groupby(p.index % 23).transform("mean")
    return p - seas


def collect(ev, anom, D, sites, lag_min, lag_max):
    recs = []
    for eid, r in enumerate(ev.itertuples()):
        i = r.site
        if i not in D.index:
            continue
        for j in sites:
            if j == i or j not in anom.columns or j not in D.columns:
                continue
            d = D.loc[i, j]
            for k in range(lag_min, lag_max + 1):
                b = r.bin_start + k
                if b in anom.index and not np.isnan(anom.loc[b, j]):
                    recs.append({"event": eid, "dist": d, "lag": k,
                                 "resp": anom.loc[b, j]})
    df = pd.DataFrame(recs)
    if len(df) == 0:
        return df
    df["dist_band"] = pd.cut(df.dist, DIST_EDGES)
    return df


def summarize(df, ev, seed=0):
    grid = df.groupby(["dist_band", "lag"], observed=True).resp.mean().rename("mean")
    rng = np.random.default_rng(seed)
    # block bootstrap over shock bins (events sharing a bin start resampled together)
    evgrp = df[["event"]].merge(ev[["bin_start"]].reset_index().rename(
        columns={"index": "event"}), on="event")
    grps = evgrp.bin_start.unique()
    boots = []
    for _ in range(200):
        pick = set(rng.choice(grps, size=len(grps), replace=True))
        sel = evgrp[evgrp.bin_start.isin(pick)].event.unique()
        sub = df[df.event.isin(sel)]
        boots.append(sub.groupby(["dist_band", "lag"], observed=True).resp.mean())
    B = pd.concat(boots, axis=1)
    out = grid.to_frame()
    out["lo"] = B.quantile(0.025, axis=1)
    out["hi"] = B.quantile(0.975, axis=1)
    out["n"] = df.groupby(["dist_band", "lag"], observed=True).size()
    out["sig"] = (out.lo > 0) | (out.hi < 0)
    return out.reset_index()


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    ev = ev[(ev["index"] == "evi") & (~ev.regionwide)].reset_index(drop=True)
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    lag_min, lag_max = -6, 12  # spec window
    for name, spp in [("herb", cfg["species"]["herbivores_primary"]),
                      ("carn", cfg["species"]["carnivores_primary"])]:
        anom = guild_anom(panel, spp)
        df = collect(ev, anom, D, sites, lag_min, lag_max)
        if len(df) == 0:
            print("irf", name, "empty")
            continue
        out = summarize(df, ev, seed=cfg["seed"])
        out["dist_band"] = out.dist_band.astype(str)
        out.to_parquet(RESULTS / "models" / f"irf_{name}.parquet", index=False)
        out.to_csv(RESULTS / "models" / f"irf_{name}.csv", index=False)
        print("irf", name, out.shape)


if __name__ == "__main__":
    main()
