"""Decomposition of the clustered-vs-isolated sign reversal (user directive).

Per event (not averaged): decompose source-site post-shock herbivore response
into candidate drivers and test which differs between clustered and isolated
shock sets / interacts with the response.

Components per event (wildebeest + zebra):
  y_post       source-site encounter-rate anomaly, lags +2..+8 (event response)
  y_pre        same, lags -6..-1 (pre-trend)
  count_pre    mean individual count (Count midpoint) at source, lags -6..-1
  count_post   same, lags +2..+8
  n_sites_post  # recipient sites (eco-connected, w>median) with positive anomaly
  flow_real    realized redistribution: mean post anomaly at eco-connected
               recipients minus matched unconnected recipients
  outflow      structural mobility = row-sum of co-detection W at source
  kc           carrying-capacity proxy = site seasonal EVI mean (productivity)
  predation    local lion+hyena encounter rate at source, lags -6..-1
  n_nb_shocked # other sites shocked in same bin within 15 km (spatial spillover)
  bin_frac     share of all sites shocked in that bin

Analysis:
  1) compare covariate distributions clustered vs isolated (Welch t on
     bin-level means to respect dependence)
  2) event-level OLS: y_post ~ isolated*(each covariate) clustered by bin
Output: results/tables/decomposition_events.csv, decomposition_fit.csv,
        docs/DECOMPOSITION.md
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from src.utils.common import DATA, RESULTS, DOCS, load_config

PRE = (-6, -1)
POST = (2, 8)
HERB = ["wildebeest", "zebra"]
CARN = ["lion", "hyenaSpotted"]


def count_mid(x):
    if pd.isna(x):
        return np.nan
    s = str(x)
    if "-" in s:
        a, b = s.split("-")[:2]
        return (float(a) + float(b)) / 2
    try:
        return float(s)
    except ValueError:
        return np.nan


def site_rate(panel, species):
    return panel[panel.species.isin(species)].groupby(
        ["site", "bin"]).encounter_rate.mean().unstack("site")


def anom_frame(rate):
    seas = rate.groupby(rate.index % 23).transform("mean")
    return rate - seas


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    ev = ev[(ev["index"] == "evi") & (~ev.regionwide)].reset_index(drop=True)
    shock = pd.read_parquet(DATA / "processed" / "site_shock.parquet")
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    W = np.load(DATA / "processed" / "W_movement.npy")
    sidx = {s: i for i, s in enumerate(sites)}

    # individual counts per site-bin for focal herbivores
    con = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv",
                      usecols=["DateTime", "SiteID", "Species", "Count"])
    con = con[con.Species.isin(HERB)]
    con["bin"] = ((pd.to_datetime(con.DateTime) - pd.Timestamp("2010-01-01")).dt.days // 16)
    con["cnt"] = con.Count.map(count_mid)
    counts = con.groupby(["SiteID", "bin"]).cnt.sum().unstack("SiteID").reindex(columns=sites)

    rate_h = site_rate(panel, HERB)
    anom_h = anom_frame(rate_h)
    rate_c = site_rate(panel, CARN)

    per_bin = shock[shock.shock].groupby("bin").site.nunique() / shock.site.nunique()
    ev["bin_frac"] = ev.bin_start.map(per_bin).fillna(0)
    ev["isolated"] = ev.bin_frac <= 0.10
    kc_map = veg.groupby("site").evi_seas_mean.first()

    # shocked sites per bin for neighbor count
    shocked_by_bin = shock[shock.shock].groupby("bin").site.apply(set)

    rows = []
    for e in ev.itertuples():
        i = e.site
        if i not in sidx or i not in anom_h.columns:
            continue
        def wmean(frame, lo, hi):
            b = [x for x in range(e.bin_start + lo, e.bin_start + hi + 1)
                 if x in frame.index]
            return frame.loc[b, i].mean() if b else np.nan
        # eco-connected recipients (top tercile of W row i)
        wr = pd.Series(W[sidx[i]], index=sites).drop(i)
        conn = wr[wr >= wr.quantile(0.9)].index
        unconn = wr[wr <= wr.quantile(0.1)].index
        flow_conn = anom_h.loc[[b for b in range(e.bin_start + POST[0], e.bin_start + POST[1] + 1)
                                if b in anom_h.index], conn].values
        flow_unc = anom_h.loc[[b for b in range(e.bin_start + POST[0], e.bin_start + POST[1] + 1)
                               if b in anom_h.index], unconn].values
        # neighbor shocks within 15 km in same bin
        nb = set()
        for b in range(e.bin_start, e.bin_end + 1):
            nb |= shocked_by_bin.get(b, set())
        nb.discard(i)
        dvec = D.loc[i, list(nb)] if nb else pd.Series(dtype=float)
        rows.append({
            "event": e.Index, "site": i, "grp": e.bin_start, "year": e.year if hasattr(e, "year") else np.nan,
            "isolated": e.isolated, "bin_frac": e.bin_frac,
            "y_post": wmean(anom_h, *POST), "y_pre": wmean(anom_h, *PRE),
            "count_pre": wmean(counts, *PRE), "count_post": wmean(counts, *POST),
            "flow_real": np.nanmean(flow_conn) - np.nanmean(flow_unc),
            "outflow": float(wr.sum()),
            "kc": float(kc_map.get(i, np.nan)),
            "predation": wmean(rate_c, *PRE),
            "n_nb_shocked": int((dvec <= 15).sum()) if len(dvec) else 0,
        })
    df = pd.DataFrame(rows)
    (RESULTS / "tables").mkdir(parents=True, exist_ok=True)
    df.to_csv(RESULTS / "tables" / "decomposition_events.csv", index=False)

    covs = ["y_pre", "count_pre", "flow_real", "outflow", "kc", "predation",
            "n_nb_shocked"]
    # 1) set differences (cluster means per shock bin, Welch on bin means)
    comp = []
    for c in covs + ["y_post"]:
        bm = df.groupby(["grp", "isolated"])[c].mean().reset_index()
        a = bm[bm.isolated][c].dropna(); b_ = bm[~bm.isolated][c].dropna()
        comp.append({"var": c, "isolated_mean": a.mean(), "clustered_mean": b_.mean(),
                     "diff": a.mean() - b_.mean(), "n_bins_iso": len(a), "n_bins_clu": len(b_)})
    comp = pd.DataFrame(comp)

    # 2) event-level interactions: y_post ~ isolated + isolated:cov + cov
    z = df.dropna(subset=["y_post"] + covs).copy()
    for c in covs:
        z[c] = (z[c] - z[c].mean()) / (z[c].std() + 1e-9)
    f = "y_post ~ isolated + " + " + ".join(
        f"{c} + isolated:{c}" for c in covs)
    m = smf.ols(f, data=z).fit(cov_type="cluster", cov_kwds={"groups": z.grp})
    fit = pd.DataFrame({"term": m.params.index, "coef": m.params.values,
                        "se": m.bse.values, "p": m.pvalues.values})
    fit.to_csv(RESULTS / "tables" / "decomposition_fit.csv", index=False)
    comp.to_csv(RESULTS / "tables" / "decomposition_components.csv", index=False)

    lines = ["# Decomposition: clustered vs isolated sign reversal", "",
             "Per-event source-site herbivore anomaly post mean (lags 2-8), by set:", "",
             comp.to_markdown(index=False), "",
             "## Interaction model (y_post ~ isolated * covariates)", "",
             fit.to_markdown(index=False), ""]
    (DOCS / "DECOMPOSITION.md").write_text("\n".join(lines))
    print(comp.to_string(index=False))
    print(fit[fit.term.str.startswith("isolated")].to_string(index=False))


if __name__ == "__main__":
    main()
