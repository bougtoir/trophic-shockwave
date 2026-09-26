"""Core estimation engine (per revised spec).

Design
------
Shock sets:
  - "all": all 149 site-level EVI shock events (z<=-2, non-regionwide)
  - "clustered": events in bins where >10% of sites shocked (temporal clusters)
  - "isolated": events in bins where <=10% of sites shocked

Stacked event-study: for each event (site i, shock bin b0), include all sites'
rows for bins in [b0+lag_min, b0+lag_max]; treated site gets event-time dummies,
all other sites contribute as controls via C(bin) FE (plus event-block FE to
absorb calendar-time differences across events). SEs clustered by shock bin;
block bootstrap over shock bins is the mandatory uncertainty validation.

Carnivore remote response: per site-bin, exposures
  geo_jt  = sum_i Wgeo_ij * shock_i,t-?  (distributed lags 1..K)
  eco_jt  = sum_i Weco_ij * shock_i,t-?
fit geo-only / eco-only / combined with local EVI + precip + temp + siteFE + binFE.

Effective replication diagnostics: ICC of post-shock response within shock bin,
ESS approximation, bootstrap uncertainty width.
"""
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from pathlib import Path
from src.utils.common import DATA, RESULTS, load_config

LAG_MIN, LAG_MAX = -6, 12
K_EXPO = 8  # exposure lags for carnivore remote model
TAB = RESULTS / "tables"
MOD = RESULTS / "models"


# ---------------- data builders ----------------

def load_all():
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    ev = ev[ev["index"] == "evi"].reset_index(drop=True)
    # shock-bin share of sites shocked at the event's onset bin
    shock = pd.read_parquet(DATA / "processed" / "site_shock.parquet")
    per_bin = shock[shock.shock].groupby("bin").site.nunique() / shock.site.nunique()
    ev["bin_frac"] = ev.bin_start.map(per_bin).fillna(0)
    ev["isolated"] = ev.bin_frac <= 0.10
    ev["shock_grp"] = ev.bin_start.astype(int)  # temporal cluster = shock bin
    ev["year"] = (pd.Timestamp("2010-01-01") + pd.to_timedelta(ev.bin_start * 16, "D")).dt.year
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    clim = pd.read_parquet(DATA / "interim" / "site_climate.parquet")
    clim["bin"] = ((clim.time - pd.Timestamp("2010-01-01")).dt.days // 16)
    clim = clim.groupby(["site", "bin"]).agg(precip_sum=("precip", "sum"),
                                           t2m_mean=("t2m", "mean")).reset_index()
    return panel, ev, veg, clim, shock


def stack_events(panel, ev_sub, species, lag_min=LAG_MIN, lag_max=LAG_MAX):
    """Stacked event-study frame. Returns DataFrame with y, etime, grp, site, bin."""
    p = panel[panel.species == species][["site", "bin", "encounter_rate", "cam_days"]]
    p = p[p.cam_days > 0]
    rate = p.set_index(["site", "bin"]).encounter_rate
    frames = []
    for eid, e in enumerate(ev_sub.itertuples()):
        b0 = e.bin_start
        bmin, bmax = b0 + lag_min, b0 + lag_max
        sub = p[(p.bin >= bmin) & (p.bin <= bmax)].copy()
        sub["etime"] = np.where(sub.site == e.site, sub.bin - b0, np.nan)
        sub["event_id"] = eid
        sub["grp"] = e.bin_start
        sub["year"] = e.year
        frames.append(sub)
    d = pd.concat(frames, ignore_index=True)
    d["block"] = d.event_id.astype(str) + "_" + d.site.astype(str)
    return d


def fit_es(d, lag_min=LAG_MIN, lag_max=LAG_MAX, ref=-1):
    """Fit event-study on stacked frame; returns per-lag stats."""
    dd = d.copy()
    for k in range(lag_min, lag_max + 1):
        if k == ref:
            continue
        dd[f"Xm{-k}" if k<0 else f"X{k}"] = (dd.etime == k).astype(float)
    dums = [(f"Xm{-k}" if k<0 else f"X{k}") for k in range(lag_min, lag_max + 1) if k != ref]
    f = "encounter_rate ~ " + " + ".join(dums) + " + C(site) + C(bin)"
    m = smf.ols(f, data=dd).fit(cov_type="cluster",
                               cov_kwds={"groups": dd.grp})
    # per-lag contributing shocks/sites
    tr = dd[dd.etime.notna()]
    rows = []
    for k in range(lag_min, lag_max + 1):
        if k == ref:
            rows.append(dict(lag=k, coef=0.0, lo=np.nan, hi=np.nan,
                             boot_lo=np.nan, boot_hi=np.nan,
                             n_shocks=np.nan, n_sites=np.nan))
            continue
        contrib = tr[tr.etime == k]
        c = m.params.get(f"Xm{-k}" if k<0 else f"X{k}", np.nan)
        nm = f"Xm{-k}" if k<0 else f"X{k}"; ci = m.conf_int().loc[nm] if nm in m.params.index else (np.nan, np.nan)
        rows.append(dict(lag=k, coef=c, lo=ci[0], hi=ci[1],
                         n_shocks=contrib.event_id.nunique(),
                         n_sites=contrib.site.nunique()))
    return pd.DataFrame(rows), m


def fast_es_coefs(d, lag_min=LAG_MIN, lag_max=LAG_MAX, ref=-1, iters=8):
    """Low-memory event-study coefs: absorb C(site)+C(bin) by alternating
    demeaning, OLS via lstsq. Used inside bootstrap loops."""
    y = d.encounter_rate.values.astype(np.float64)
    klist = [k for k in range(lag_min, lag_max + 1) if k != ref]
    X = np.column_stack([(d.etime.values == k).astype(np.float64) for k in klist])
    site = pd.factorize(d.site)[0]
    binc = pd.factorize(d.bin)[0]
    def absorb(v):
        r = v.copy()
        ns, nb = np.bincount(site), np.bincount(binc)
        for _ in range(iters):
            r -= (np.bincount(site, r) / ns)[site]
            r -= (np.bincount(binc, r) / nb)[binc]
        return r
    ya = absorb(y)
    Xa = np.column_stack([absorb(X[:, j]) for j in range(X.shape[1])])
    coef, *_ = np.linalg.lstsq(Xa, ya, rcond=None)
    return dict(zip(klist, coef))


def block_boot(panel, ev_sub, species, coefs_of_interest=None, nboot=300, seed=0):
    """Resample shock bins (grp) with replacement; refit; return per-lag boot CI.
    Stack is built once; each replicate subsets/duplicates rows of the sampled
    shock bins (standard block bootstrap over clusters)."""
    rng = np.random.default_rng(seed)
    d0 = stack_events(panel, ev_sub, species)
    by_grp = {g: d0[d0.grp == g] for g in d0.grp.unique()}
    grps = np.array(list(by_grp.keys()))
    boots = []
    lags = [k for k in range(LAG_MIN, LAG_MAX + 1) if k != -1]
    for b in range(nboot):
        pick = rng.choice(grps, size=len(grps), replace=True)
        d = pd.concat([by_grp[g] for g in pick], ignore_index=True)
        try:
            cc = fast_es_coefs(d)
            boots.append([cc[k] for k in lags])
        except Exception:
            boots.append([np.nan] * len(lags))
    B = np.array(boots)
    return pd.DataFrame({"lag": lags,
                         "boot_lo": np.nanpercentile(B, 2.5, axis=0),
                         "boot_hi": np.nanpercentile(B, 97.5, axis=0)})


# ---------------- carnivore remote-response ----------------

def expo_matrix(shock, W, sites, k):
    """exposure_j,t = sum_i W_ij shock_i,t-k  -> DataFrame site,bin,E"""
    sidx = {s: i for i, s in enumerate(sites)}
    S = shock.pivot_table(index="bin", columns="site", values="shock",
                          aggfunc="max").reindex(columns=sites).fillna(0).values
    bins = np.sort(shock.bin.unique())
    E = np.zeros(S.shape, dtype=float)
    E[k:, :] = S[:-k, :] @ W  # W_ij maps source i -> recipient j
    return pd.DataFrame({"bin": np.repeat(bins, len(sites)),
                         "site": list(sites) * len(bins),
                         "E": E.ravel()})


def expo_matrix_xf(shock, ev, sites, k, WA, WB, split_bin):
    """Cross-fitted exposure: shock events starting before split_bin are
    propagated through WB (connectivity estimated on the later period) and
    vice versa, so W can never encode the outcomes it is tested on."""
    em = ev[["site", "bin_start", "bin_end"]].drop_duplicates()
    j = shock.merge(em, on="site")
    j = j[(j.bin >= j.bin_start) & (j.bin <= j.bin_end)]
    inA = set(zip(j.site[j.bin_start < split_bin], j.bin[j.bin_start < split_bin]))
    sA = shock.copy(); sA["shock"] = sA.shock & sA.apply(lambda r: (r.site, r.bin) in inA, axis=1)
    sB = shock.copy(); sB["shock"] = sB.shock & ~sB.apply(lambda r: (r.site, r.bin) in inA, axis=1)
    eA = expo_matrix(sA, WB, sites, k)
    eB = expo_matrix(sB, WA, sites, k)
    out = eA.merge(eB, on=["site", "bin"], suffixes=("_a", "_b"))
    out["E"] = out.E_a + out.E_b
    return out[["site", "bin", "E"]]


def carnivore_model(panel, shock, veg, clim, W, sites, species, label,
                    combined_with=None, ev=None, WA=None, WB=None, split_bin=None):
    d = panel[panel.species.isin(species)].groupby(
        ["site", "bin"]).encounter_rate.mean().reset_index()
    expo = (lambda k: expo_matrix_xf(shock, ev, sites, k, WA, WB, split_bin)
           ) if WA is not None else (lambda k: expo_matrix(shock, W, sites, k))
    for k in range(1, K_EXPO + 1):
        e = expo(k).rename(columns={"E": f"E{k}"})
        d = d.merge(e, on=["site", "bin"], how="left")
    if combined_with is not None:
        for k in range(1, K_EXPO + 1):
            e = expo_matrix(shock, combined_with, sites, k).rename(columns={"E": f"G{k}"})
            d = d.merge(e, on=["site", "bin"], how="left")
    d = d.merge(veg[["site", "bin", "evi_z", "ndvi_z"]], on=["site", "bin"], how="left") \
         .merge(clim, on=["site", "bin"], how="left")
    d = d.dropna(subset=["encounter_rate"])
    eco = " + ".join(f"E{k}" for k in range(1, K_EXPO + 1))
    rhs = eco + " + evi_z + precip_sum + t2m_mean + C(site) + C(bin)"
    if combined_with is not None:
        rhs = " + ".join(f"G{k}" for k in range(1, K_EXPO + 1)) + " + " + rhs
    m = smf.ols(f"encounter_rate ~ {rhs}", data=d).fit(
        cov_type="cluster", cov_kwds={"groups": d.bin})
    eco_coefs = [m.params.get(f"E{k}", np.nan) for k in range(1, K_EXPO + 1)]
    geo_coefs = [m.params.get(f"G{k}", np.nan) for k in range(1, K_EXPO + 1)] \
        if combined_with is not None else [np.nan] * K_EXPO
    return {"model": label, "aic": m.aic, "bic": m.bic, "n": len(d),
            # 'expo_sum' = sum of distributed-lag coefs on the model's primary
            # exposure matrix (E cols); 'expo2_sum' = the secondary (G cols)
            "expo_sum": np.nansum(eco_coefs), "expo2_sum": np.nansum(geo_coefs),
            "eco_sum": np.nansum(eco_coefs), "geo_sum": np.nansum(geo_coefs),
            "coefs": {"expo": eco_coefs, "expo2": geo_coefs}}, m, d


def cv_rmse(d, nfold=5, seed=0):
    """Site-blocked CV RMSE on FE-absorbed design (site + bin demeaned).
    numpy implementation: predictors E1..E8 + evi_z + precip_sum + t2m_mean."""
    cols = [f"E{k}" for k in range(1, K_EXPO + 1)] + ["evi_z", "precip_sum", "t2m_mean"]
    rng = np.random.default_rng(seed)
    dd = d.dropna(subset=["encounter_rate"] + cols)
    ids = np.array(sorted(dd.site.unique()))
    rng.shuffle(ids)
    rmses = []
    site_c = pd.factorize(dd.site)[0]
    bin_c = pd.factorize(dd.bin)[0]
    def demean(v, mask_tr):
        r = v.copy()
        for _ in range(6):
            r[mask_tr] -= np.bincount(site_c[mask_tr], r[mask_tr],
                                      minlength=site_c.max() + 1)[site_c[mask_tr]] / \
                          np.maximum(1, np.bincount(site_c[mask_tr],
                                      minlength=site_c.max() + 1))[site_c[mask_tr]]
            r[mask_tr] -= np.bincount(bin_c[mask_tr], r[mask_tr],
                                      minlength=bin_c.max() + 1)[bin_c[mask_tr]] / \
                          np.maximum(1, np.bincount(bin_c[mask_tr],
                                      minlength=bin_c.max() + 1))[bin_c[mask_tr]]
        return r
    X = dd[cols].values.astype(float)
    y = dd.encounter_rate.values.astype(float)
    for f in np.array_split(ids, nfold):
        tr = ~dd.site.isin(f).values
        rX = X.copy(); ry = y.copy()
        for j in range(rX.shape[1]):
            rX[tr, j] = demean(rX[:, j], tr)[tr]
        ry[tr] = demean(y, tr)[tr]
        try:
            coef, *_ = np.linalg.lstsq(rX[tr], ry[tr], rcond=None)
            pred = X[~tr] @ coef
            rmses.append(np.sqrt(np.mean((y[~tr] - pred) ** 2)))
        except Exception:
            pass
    return float(np.nanmean(rmses)) if rmses else np.nan


# ---------------- main ----------------

def main():
    cfg = load_config()
    panel, ev, veg, clim, shock = load_all()
    TAB.mkdir(parents=True, exist_ok=True)
    MOD.mkdir(parents=True, exist_ok=True)
    ev = ev[~ev.regionwide]
    sets = {"all": ev, "clustered": ev[~ev.isolated], "isolated": ev[ev.isolated]}
    summary_rows = []

    # ---- effective replication table
    per_bin_sites = ev.groupby("bin_start").site.nunique()
    er = {"total_site_shocks": len(ev), "unique_shock_bins": ev.bin_start.nunique(),
          "isolated_shock_events": int(ev.isolated.sum()),
          "clustered_shock_events": int((~ev.isolated).sum()),
          "sites_per_bin_mean": float(per_bin_sites.mean()),
          "sites_per_bin_median": float(per_bin_sites.median()),
          "events_by_year": ev.groupby("year").size().to_dict()}
    # ICC of post-shock herbivore response within shock bin
    try:
        d0 = stack_events(panel, ev, "wildebeest")
        tr = d0[d0.etime.between(2, 8)]
        grp_mean = tr.groupby("grp").encounter_rate.mean()
        resid = tr.encounter_rate - tr.etime.map(tr.groupby("etime").encounter_rate.mean())
        gm = tr.assign(r=resid).groupby("grp").r.mean()
        # ICC via one-way ANOVA
        dfa = tr.assign(r=resid).groupby("grp").r
        grand = resid.mean()
        ssb = ((dfa.mean() - grand) ** 2 * dfa.size()).sum()
        ssw = ((resid - tr.grp.map(gm)) ** 2).sum()
        k_, n_ = dfa.size().shape[0], len(resid)
        m0 = dfa.size().mean()
        msw = ssw / (n_ - k_); msb = ssb / (k_ - 1)
        icc = max(0.0, (msb - msw) / (msb + (m0 - 1) * msw))
        er["icc_shock_bin_wildebeest_post"] = icc
        er["ess_approx"] = len(ev) / (1 + (ev.groupby("bin_start").size().mean() - 1) * icc)
    except Exception as e:
        er["icc_error"] = repr(e)
    pd.DataFrame([er]).to_csv(TAB / "effective_replication.csv", index=False)
    print("effective replication:", er)

    # ---- herbivore event studies
    herb = cfg["species"]["herbivores_primary"]
    boot_species = set(herb + cfg["species"]["carnivores_primary"])
    for sname, evs in sets.items():
        for sp in herb + cfg["species"]["carnivores_primary"] + cfg["species"]["negative_control_species"]:
            if len(evs) < 8:
                continue
            try:
                d = stack_events(panel, evs, sp)
                coefs, m = fit_es(d)
                # mandatory block bootstrap on primary species x all/isolated
                if sp in boot_species and sname in ("all", "isolated"):
                    bb = block_boot(panel, evs, sp, nboot=200,
                                    seed=cfg["seed"])
                    coefs = coefs.merge(bb, on="lag", how="left")
                coefs["species"] = sp
                coefs["shock_set"] = sname
                coefs["n_events"] = len(evs)
                coefs.to_csv(TAB / f"es_{sname}_{sp}.csv", index=False)
                summary_rows.append({"set": sname, "species": sp,
                                     "post_mean_2_8": coefs[coefs.lag.between(2, 8)].coef.mean(),
                                     "n_events": len(evs)})
                print("es", sname, sp, "done")
            except Exception as e:
                print("es fail", sname, sp, repr(e)[:150])

    # ---- carnivore remote model: geo vs eco vs combined
    sites = np.load(DATA / "processed" / "sites_order.npy", allow_pickle=True)
    Wg = np.load(DATA / "processed" / "W_geo.npy")
    # cross-fit connectivity: split time at median shock bin; W for evaluating
    # events in one half is estimated on co-detections in the other half only
    from src.preprocess.connectivity import movement_W
    split_bin = ev.bin_start.median()
    WA = movement_W(panel, sites, mask_bins=range(int(split_bin), 76))
    WB = movement_W(panel, sites, mask_bins=range(0, int(split_bin)))
    comp_rows = []
    for sname, evs in sets.items():
        sub_shock = shock.copy()
        ok = set()
        for r in evs.itertuples():
            for b in range(r.bin_start, r.bin_end + 1):
                ok.add((r.site, b))
        sub_shock["shock"] = sub_shock.shock & sub_shock.apply(
            lambda r: (r.site, r.bin) in ok, axis=1)
        evs_x = evs.assign(bin_end=evs.bin_end)
        for (lbl, W1, W2) in [("geo", Wg, None), ("eco", None, None),
                              ("combined", None, Wg)]:
            try:
                res, m, d = carnivore_model(panel, sub_shock, veg, clim, W1, sites,
                                            cfg["species"]["carnivores_primary"], lbl,
                                            combined_with=W2, ev=evs,
                                            WA=(WA if W1 is None else None),
                                            WB=WB, split_bin=split_bin)
                res["shock_set"] = sname
                res["connectivity"] = "xfit"
                # blocked CV
                eco = " + ".join(f"E{k}" for k in range(1, K_EXPO + 1))
                rhs = f"encounter_rate ~ {eco} + evi_z + precip_sum + t2m_mean + C(site) + C(bin)"
                res["cv_rmse"] = cv_rmse(d, rhs, seed=cfg["seed"])
                comp_rows.append(res)
                with open(MOD / f"carn_{sname}_{lbl}.json", "w") as f:
                    json.dump({k: v for k, v in res.items() if k != "coefs"} |
                              res["coefs"], f, indent=2)
                print("carn", sname, lbl, "aic=%.0f" % res["aic"])
            except Exception as e:
                print("carn fail", sname, lbl, repr(e)[:200])
    pd.DataFrame(comp_rows).to_csv(TAB / "connectivity_model_comparison.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(TAB / "es_summary.csv", index=False)


if __name__ == "__main__":
    main()
