"""Figures A-H and tables 1-7 (spec sections 29-30). All scripted; no manual edits.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from src.utils.common import DATA, RESULTS, DOCS, load_config

FIG = RESULTS / "figures"
TAB = RESULTS / "tables"


def fig_A():
    sites = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv",
                        usecols=["SiteID", "LocationX", "LocationY"]).drop_duplicates("SiteID")
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    snap = veg.groupby("site").evi_z.min()
    d = sites.set_index("SiteID").join(snap.rename("min_z"))
    fig, ax = plt.subplots(figsize=(7, 6))
    sc = ax.scatter(d.LocationX / 1000, d.LocationY / 1000, c=d.min_z,
                    cmap="RdYlGn", vmin=-3, vmax=0, s=25)
    ax.set_xlabel("UTM36S x (km)"); ax.set_ylabel("UTM36S y (km)")
    ax.set_title("Figure A. Camera sites, min EVI z-score (2010-2013)")
    fig.colorbar(sc, label="min EVI z"); fig.tight_layout()
    fig.savefig(FIG / "figA_map.png", dpi=200); plt.close(fig)


def fig_BD():
    es = pd.read_csv(TAB / "event_study_coefficients.csv")
    for name, spp, ttl in [("figB_herb_eventstudy.png",
                            ["wildebeest", "zebra"], "Figure B. Herbivore event-study"),
                           ("figD_carn_eventstudy.png",
                            ["lion", "hyenaSpotted"], "Figure D. Carnivore event-study")]:
        fig, ax = plt.subplots(figsize=(8, 5))
        for sp in spp:
            d = es[es.species == sp].sort_values("lag")
            ax.errorbar(d.lag, d.coef, yerr=[d.coef - d.lo, d.hi - d.coef],
                        label=sp, marker="o", ms=3, capsize=2)
        ax.axhline(0, color="k", lw=0.7); ax.axvline(0, color="r", ls="--", lw=0.8)
        ax.set_xlabel("event time (16-day bins)"); ax.set_ylabel("encounter-rate effect")
        ax.legend(); ax.set_title(ttl); fig.tight_layout()
        fig.savefig(FIG / name, dpi=200); plt.close(fig)


def fig_C():
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    ev = ev[(ev["index"] == "evi") & (~ev.regionwide)].head(6)
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    sites = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv",
                        usecols=["SiteID", "LocationX", "LocationY"]).drop_duplicates("SiteID").set_index("SiteID")
    herb = panel[panel.species == "wildebeest"].groupby(["site", "bin"]).encounter_rate.mean().unstack()
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharex=True, sharey=True)
    for ax, (_, r) in zip(axes.ravel(), ev.iterrows()):
        for k, alpha in [(0, 0.4), (4, 0.7), (8, 1.0)]:
            b = r.bin_start + k
            if b in herb.index:
                ax.scatter(sites.LocationX / 1000, sites.LocationY / 1000,
                           c=herb.loc[b].reindex(sites.index), cmap="Blues", s=8,
                           alpha=alpha)
        i = r.site
        ax.scatter(sites.loc[i].LocationX / 1000, sites.loc[i].LocationY / 1000,
                   marker="*", s=120, c="red")
        ax.set_title(f"event {i} bin {r.bin_start}")
    fig.suptitle("Figure C. Wildebeest encounter rate maps around shocks (* = source)")
    fig.tight_layout(); fig.savefig(FIG / "figC_redistribution.png", dpi=200); plt.close(fig)


def fig_E():
    for name, ttl in [("herb", "Figure E1. Herbivore IRF"), ("carn", "Figure E2. Carnivore IRF")]:
        p = RESULTS / "models" / f"irf_{name}.parquet"
        if not p.exists():
            continue
        d = pd.read_parquet(p)
        piv = d.pivot(index="dist_band", columns="lag", values="mean")
        fig, ax = plt.subplots(figsize=(9, 5))
        im = ax.imshow(piv.values, aspect="auto", cmap="RdBu_r",
                       vmin=-np.nanmax(np.abs(piv.values)),
                       vmax=np.nanmax(np.abs(piv.values)),
                       extent=[piv.columns.min(), piv.columns.max(), len(piv) - .5, -.5])
        ax.set_yticks(range(len(piv))); ax.set_yticklabels([str(x) for x in piv.index])
        ax.set_xlabel("lag (16-day bins)"); ax.set_ylabel("distance band (km)")
        ax.axvline(0, color="k", ls="--")
        ax.set_title(ttl); fig.colorbar(im, label="mean anomaly")
        fig.tight_layout(); fig.savefig(FIG / f"figE_irf_{name}.png", dpi=200); plt.close(fig)


def fig_FG():
    p = TAB / "connectivity_model_comparison.csv"
    if p.exists():
        d = pd.read_csv(p).sort_values("aic")
        fig, ax = plt.subplots(figsize=(7, 4))
        x = np.arange(len(d))
        ax.bar(x - 0.2, d.aic - d.aic.min(), 0.35, label="AIC (rel)")
        ax.bar(x + 0.2, (d.cv_rmse - d.cv_rmse.min()) / (d.cv_rmse.max() - d.cv_rmse.min() + 1e-12),
               0.35, label="CV RMSE (norm)")
        ax.set_xticks(x); ax.set_xticklabels(d.W)
        ax.set_title("Figure F. Connectivity model comparison"); ax.legend()
        fig.tight_layout(); fig.savefig(FIG / "figF_connectivity.png", dpi=200); plt.close(fig)
    nc = TAB / "negative_controls.csv"
    if nc.exists():
        d = pd.read_csv(nc)
        mc = d[d.control == "matched_unconnected"]
        if len(mc):
            fig, ax = plt.subplots(figsize=(6, 4))
            v = mc.set_index("stat")["value"]
            ax.bar(["connected", "unconnected"], [v.get("mean_resp_connected", 0),
                   v.get("mean_resp_unconnected", 0)])
            ax.set_title("Figure G. Connected vs matched-unconnected recipients")
            fig.tight_layout(); fig.savefig(FIG / "figG_matched.png", dpi=200); plt.close(fig)


def fig_H():
    p = TAB / "simulation_results.csv"
    if p.exists():
        d = pd.read_csv(p)
        fig, ax = plt.subplots(figsize=(8, 4))
        x = np.arange(len(d))
        ax.bar(x - 0.2, d.power_mov, 0.35, label="movement coef p<.05")
        ax.bar(x + 0.2, d.power_geo, 0.35, label="geo coef p<.05")
        ax.set_xticks(x); ax.set_xticklabels(d.scenario)
        ax.set_ylabel("fraction of reps"); ax.axhline(0.05, ls="--", c="gray")
        ax.set_title("Figure H. Simulation validation")
        ax.legend(); fig.tight_layout()
        fig.savefig(FIG / "figH_simulation.png", dpi=200); plt.close(fig)


def tables():
    cfg = load_config()
    cons = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv")
    t2 = cons.groupby("Species").CaptureEventID.nunique().sort_values(ascending=False)
    eff = pd.read_parquet(DATA / "processed" / "panel.parquet")
    eff_sp = eff.groupby("species").agg(cam_days=("cam_days", "sum"),
                                        events=("encounters", "sum")).reset_index()
    t2.to_csv(TAB / "table2_species_counts.csv")
    eff_sp.to_csv(TAB / "table2_effort.csv", index=False)
    ev = pd.read_parquet(DATA / "processed" / "shock_events.parquet")
    ev[ev["index"] == "evi"].to_csv(TAB / "table3_shock_events.csv", index=False)


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    for f in [fig_A, fig_C, fig_E, fig_FG, fig_H, tables]:
        try:
            f(); print("ok", f.__name__)
        except Exception as e:
            print("FAIL", f.__name__, repr(e)[:150])
    p = TAB / "event_study_coefficients.csv"
    if p.exists():
        fig_BD(); print("ok fig_BD")


if __name__ == "__main__":
    main()
