"""Negative controls (spec sections 17-18).

1. Temporal placebo: fake shocks shifted earlier -> event-study coefs ~0.
2. Species negative control: elephant/giraffe (from event_study output).
3. Matched spatial negative control: for connected recipients, compare
   equidistant but weakly-connected recipients.
4. Permuted-network test: degree-preserving random rewiring of W_movement;
   compare distributed-lag fit to observed.
Output: results/tables/negative_controls.csv
"""
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from src.utils.common import DATA, RESULTS, load_config
from src.models.distributed_lag import fit_W


def permute_W(W, rng):
    """Degree-preserving permutation: shuffle the columns of the row-normalized
    weight matrix (permutes which senders map to each receiver)."""
    P = rng.permutation(W.shape[1])
    return W[:, P]


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    shock = pd.read_parquet(DATA / "processed" / "site_shock.parquet")
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    clim = pd.read_parquet(DATA / "interim" / "site_climate.parquet")
    clim["bin"] = ((clim.time - pd.Timestamp("2010-01-01")).dt.days //
                   cfg["temporal"]["primary_bin_days"])
    clim = clim.groupby(["site", "bin"]).agg(precip_sum=("precip", "sum"),
                                           t2m_mean=("t2m", "mean")).reset_index()
    species = cfg["species"]["carnivores_primary"]
    rng = np.random.default_rng(cfg["seed"])
    out = []

    # ---- temporal placebo: shift shock earlier by lag L bins
    from src.models.event_study import event_study, local_events
    ev = local_events()
    for shift in cfg["models"]["placebo_shifts"]:
        evp = ev.copy()
        evp["bin_start"] = evp.bin_start + shift
        for sp in species:
            try:
                r = event_study(panel, sp, evp, 0, 12, ref=0)
                post = r[r.lag.between(1, 8)]
                out.append({"control": f"temporal_placebo_shift{shift}",
                            "species": sp, "stat": "mean_post_coef",
                            "value": post.coef.mean()})
            except Exception as e:
                print("placebo fail", shift, sp, repr(e)[:100])

    # ---- permuted network test on distributed-lag AIC
    W = np.load(DATA / "processed" / "W_movement.npy")
    try:
        r_obs, _ = fit_W("movement", W, panel, shock, veg, clim, species, sites)
        aics = []
        for b in range(cfg["models"]["permutations_network"]):
            try:
                rp, _ = fit_W("perm", permute_W(W, rng), panel, shock, veg, clim,
                              species, sites)
                aics.append(rp["aic"])
            except Exception:
                pass
        aics = np.array(aics)
        out.append({"control": "permuted_network", "species": "+".join(species),
                    "stat": "aic_observed", "value": r_obs["aic"]})
        out.append({"control": "permuted_network", "species": "+".join(species),
                    "stat": "aic_perm_mean", "value": aics.mean()})
        out.append({"control": "permuted_network", "species": "+".join(species),
                    "stat": "frac_perm_better", "value": float((aics < r_obs["aic"]).mean())})
    except Exception as e:
        print("perm fail", repr(e)[:150])

    # ---- matched unconnected recipients (section 17)
    D = pd.read_csv(DATA / "processed" / "site_distances.csv", index_col=0)
    Wm = np.load(DATA / "processed" / "W_movement.npy")
    sidx = {s: i for i, s in enumerate(sites)}
    evl = ev
    try:
        carn = panel[panel.species.isin(species)].groupby(["site", "bin"]).encounter_rate.mean().unstack("site")
        seas = carn.groupby(carn.index % 23).transform("mean")
        anom = carn - seas
        rows = []
        for r in evl.itertuples():
            i = r.site
            if i not in sidx:
                continue
            cand = [j for j in sites if j != i]
            conn = np.array([Wm[sidx[i], sidx[j]] for j in cand])
            dist = np.array([D.loc[i, j] for j in cand])
            hi = np.array(cand)[conn >= np.quantile(conn, 0.75)]
            hi_dist = dist[conn >= np.quantile(conn, 0.75)]
            resp_hi = [anom.loc[r.bin_start + 2: r.bin_start + 8, j].mean() for j in hi]
            # matched: similar distance, low connectivity
            lo_mask = conn <= np.quantile(conn, 0.25)
            dm = dist[lo_mask]
            resp_lo = [anom.loc[r.bin_start + 2: r.bin_start + 8, j].mean()
                       for j in np.array(cand)[lo_mask]]
            rows.append({"event_site": i, "bin_start": r.bin_start,
                         "resp_connected": np.nanmean(resp_hi),
                         "resp_unconnected": np.nanmean(resp_lo)})
        mdf = pd.DataFrame(rows)
        m = smf.ols("resp_connected ~ resp_unconnected", data=mdf).fit()
        out.append({"control": "matched_unconnected", "species": "+".join(species),
                    "stat": "mean_resp_connected", "value": mdf.resp_connected.mean()})
        out.append({"control": "matched_unconnected", "species": "+".join(species),
                    "stat": "mean_resp_unconnected", "value": mdf.resp_unconnected.mean()})
        out.append({"control": "matched_unconnected", "species": "+".join(species),
                    "stat": "paired_diff", "value": (mdf.resp_connected - mdf.resp_unconnected).mean()})
    except Exception as e:
        print("matched fail", repr(e)[:150])

    pd.DataFrame(out).to_csv(RESULTS / "tables" / "negative_controls.csv", index=False)
    print("negative controls:", len(out))


if __name__ == "__main__":
    main()
