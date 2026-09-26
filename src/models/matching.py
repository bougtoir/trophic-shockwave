"""Exchangeability tests + matched re-estimation (user directive).

First hypothesis: clustered and isolated shocks are NOT exchangeable —
they may draw from different event populations (e.g. regional dry pulses vs
ongoing local declines). This module:

  1) builds per-event covariates: season (bin-of-year), magnitude (min_z),
     duration, baseline pre-trend (y_pre), carrying capacity (kc);
  2) propensity-style overlap check: logistic P(clustered) on covariates,
     and whether clustered events lie inside the isolated covariate support
     (convex hull / nearest-neighbor distance);
  3) matches each clustered event to its nearest isolated events on the
     standardized covariate vector (season, magnitude, duration, y_pre, kc)
     with a caliper; reports how many clustered events are matchable at all;
  4) re-estimates post-mean herbivore response on:
       a. all events minus 2010 events (2010 dry pulse excluded)
       b. isolated events with flat pre-trend (|y_pre| < 1 MAD)
       c. matched isolated vs matched clustered subsets;
  5) leave-one-shock-bin-out post-means (uses estimation.post_mean).

Output: results/tables/exchangeability.csv, matched_event_study.csv,
        docs/MATCHING.md
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.spatial.distance import cdist
from src.utils.common import DATA, RESULTS, DOCS, load_config
from src.models.estimation import load_all, stack_events, fit_es, LAG_MIN, LAG_MAX, TAB

POST = (2, 8)


def post_mean(panel, evs, species="wildebeest"):
    if len(evs) < 4:
        return np.nan
    d = stack_events(panel, evs, species)
    c, _ = fit_es(d)
    return float(c[c.lag.between(*POST)].coef.mean())


def main():
    cfg = load_config()
    panel, ev, veg, clim, shock = load_all()
    ev = ev[~ev.regionwide].reset_index(drop=True)
    de = pd.read_csv(RESULTS / "tables" / "decomposition_events.csv")
    df = ev.merge(de[["event", "y_post", "y_pre", "count_pre", "kc"]],
                  left_index=True, right_on="event", how="left")
    df["season"] = df.bin_start % 23
    df["magnitude"] = -df.min_z  # larger = stronger shock
    covs = ["season", "magnitude", "duration_bins", "y_pre", "kc"]
    df = df.dropna(subset=covs + ["y_post"]).reset_index(drop=True)
    TAB.mkdir(parents=True, exist_ok=True)

    z = df.copy()
    z["isolated"] = z.isolated.astype(int)
    for c in covs:
        z[c] = (z[c] - z[c].mean()) / (z[c].std() + 1e-9)

    # --- propensity overlap ---
    pm_ = smf.logit("isolated ~ " + " + ".join(covs), data=z).fit(disp=0)
    z["ps"] = pm_.predict(z)
    ps_iso = z[z.isolated==1].ps; ps_clu = z[z.isolated==0].ps
    overlap_lo = max(ps_iso.min(), ps_clu.min())
    overlap_hi = min(ps_iso.max(), ps_clu.max())
    n_clu_in_support = int(((ps_clu >= overlap_lo) & (ps_clu <= overlap_hi)).sum())

    X = z[covs].values
    clu = z[z.isolated==0]; iso = z[z.isolated==1]
    dmat = cdist(clu[covs].values, iso[covs].values)
    nn = dmat.min(1)
    # caliper: median iso-iso NN distance * 1.5
    dii = cdist(iso[covs].values, iso[covs].values)
    np.fill_diagonal(dii, np.inf)
    caliper = np.median(dii.min(1)) * 1.5
    matchable = int((nn <= caliper).sum())

    exch = {"n_events": len(df), "n_isolated": int(df.isolated.sum()),
            "n_clustered": int((~df.isolated).sum()),
            "ps_overlap_range": [float(overlap_lo), float(overlap_hi)],
            "n_clustered_in_ps_support": n_clu_in_support,
            "nn_caliper": float(caliper),
            "n_clustered_matchable": matchable,
            "frac_clustered_matchable": matchable / max(1, len(clu))}
    pd.DataFrame([exch]).to_csv(TAB / "exchangeability.csv", index=False)

    # --- re-estimations ---
    rows = []
    ev_all = ev.merge(df[["event", "y_pre"]], left_index=True, right_on="event")
    rows.append({"analysis": "all_pooled", "n": len(ev_all),
                 "post_mean": post_mean(panel, ev_all)})
    rows.append({"analysis": "exclude_2010", "n": int((ev_all.year != 2010).sum()),
                 "post_mean": post_mean(panel, ev_all[ev_all.year != 2010])})
    iso_ev = ev_all[ev_all.index.isin(df[df.isolated].event)]
    # pre-trend filter: |y_pre| within 1.5 MAD of isolated events' median
    yp = df[df.isolated].y_pre
    med, mad = yp.median(), (yp - yp.median()).abs().median() + 1e-9
    clean_iso = df[df.isolated & ((df.y_pre - med).abs() <= 1.5 * mad)]
    rows.append({"analysis": "isolated_no_pretrend", "n": len(clean_iso),
                 "post_mean": post_mean(panel, ev.loc[ev.index.isin(
                     clean_iso.event)]) if len(clean_iso) >= 4 else np.nan})
    # matched sets
    match_idx = np.where(nn <= caliper)[0]
    mrows = []
    if len(match_idx) >= 3:
        picked_iso = set()
        for ci in match_idx:
            for j in np.argsort(dmat[ci])[:2]:
                picked_iso.add(j)
        mclu = ev.loc[ev.index.isin(clu.iloc[match_idx].event)]
        miso = ev.loc[ev.index.isin(iso.iloc[list(picked_iso)].event)]
        rows.append({"analysis": "clustered_matched_only", "n": len(mclu),
                     "post_mean": post_mean(panel, mclu)})
        rows.append({"analysis": "isolated_matched", "n": len(miso),
                     "post_mean": post_mean(panel, miso)})
    else:
        rows.append({"analysis": "matched_set", "n": 0,
                     "post_mean": np.nan})  # match failure IS the finding

    # LOSBO
    for g in sorted(ev.bin_start.unique()):
        sub = ev[ev.bin_start != g]
        rows.append({"analysis": f"LOSBO_drop_bin{g}", "n": len(sub),
                     "post_mean": post_mean(panel, sub)})
    pd.DataFrame(rows).to_csv(TAB / "matched_event_study.csv", index=False)

    lines = ["# Exchangeability: clustered vs isolated", "",
             "```", pd.DataFrame([exch]).T.to_string(), "```", "",
             "## Propensity model (isolated ~ covariates)", "",
             pm_.summary2().tables[1].to_string(), "",
             "## Re-estimation results", "",
             pd.DataFrame(rows).to_markdown(index=False), ""]
    (DOCS / "MATCHING.md").write_text("\n".join(lines))
    print(exch)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
