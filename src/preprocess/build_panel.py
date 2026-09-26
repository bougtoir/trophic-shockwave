"""Build master panel: site x 16-day bin.

Responses (spec section 10):
- detection: any consensus event of species in bin (0/1)
- encounters: number of capture events of species in bin
- encounter_rate: encounters / active camera-days in bin  (PRIMARY)
Camera effort = days within search_effort deployment windows intersecting bin.
"""
import pandas as pd
import numpy as np
from src.utils.common import DATA, RESULTS, load_config, bin_index

ORIGIN = "2010-01-01"


def load_consensus():
    df = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv")
    df["DateTime"] = pd.to_datetime(df["DateTime"])
    return df


def effort_windows():
    e = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "search_effort.csv")
    e.columns = ["site", "start", "end"]
    e["start"], e["end"] = pd.to_datetime(e["start"]), pd.to_datetime(e["end"])
    return e


def camera_days(eff, site, t0, t1):
    w = eff[(eff.site == site) & (eff.end >= t0) & (eff.start <= t1)]
    tot = 0.0
    for _, r in w.iterrows():
        lo, hi = max(r.start, t0), min(r.end, t1)
        tot += max(0.0, (hi - lo).days + 1)
    return tot


def build(bin_days=16):
    cfg = load_config()
    cons = load_consensus()
    eff = effort_windows()
    cons["bin"] = bin_index(cons["DateTime"], bin_days, ORIGIN)
    bins = np.arange(cons["bin"].min(), cons["bin"].max() + 1)
    sites = sorted(cons.SiteID.unique())
    species_groups = (cfg["species"]["herbivores_primary"] + cfg["species"]["herbivores_secondary"]
                      + cfg["species"]["carnivores_primary"] + cfg["species"]["carnivores_secondary"]
                      + cfg["species"]["negative_control_species"])
    # lionFemale+lionMale -> lion
    cons["Species2"] = cons.Species.replace({"lionFemale": "lion", "lionMale": "lion"})

    # effort per site-bin
    eff_rows = []
    for s in sites:
        for b in bins:
            t0 = pd.Timestamp(ORIGIN) + pd.Timedelta(days=int(b) * bin_days)
            t1 = t0 + pd.Timedelta(days=bin_days - 1)
            eff_rows.append({"site": s, "bin": b, "cam_days": camera_days(eff, s, t0, t1)})
    effort = pd.DataFrame(eff_rows)

    # observations per species
    obs = cons[cons.Species2.isin(species_groups)].groupby(
        ["SiteID", "bin", "Species2"]).agg(events=("CaptureEventID", "nunique")).reset_index()
    obs.columns = ["site", "bin", "species", "events"]

    panel = effort.merge(obs, on=["site", "bin"], how="left")
    panel["species"] = panel.species.fillna("")
    recs = []
    for sp in species_groups:
        p = effort.copy()
        sub = obs[obs.species == sp][["site", "bin", "events"]]
        p = p.merge(sub, on=["site", "bin"], how="left").fillna({"events": 0})
        p["species"] = sp
        p["detection"] = (p.events > 0).astype(int)
        p["encounters"] = p.events
        p["encounter_rate"] = np.where(p.cam_days > 0, p.events / p.cam_days, np.nan)
        recs.append(p[["site", "bin", "species", "cam_days", "detection",
                       "encounters", "encounter_rate"]])
    panel = pd.concat(recs, ignore_index=True)
    out = DATA / "processed" / "panel.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    panel.to_parquet(out, index=False)
    print("panel:", len(panel), "rows", panel.site.nunique(), "sites",
          panel["bin"].nunique(), "bins")
    return panel


def main():
    build()


if __name__ == "__main__":
    main()
