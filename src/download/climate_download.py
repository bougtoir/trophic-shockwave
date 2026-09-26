"""Download ERA5-Land hourly temperature + precipitation per camera site via
Open-Meteo archive API (no key; CC-BY). Writes data/interim/site_climate.parquet
and data/interim/climate_manifest.json.
"""
import json, time
from datetime import datetime, timezone

import pandas as pd
import requests

from src.utils.common import DATA, load_config
from src.download.modis_extract import site_table

API = "https://archive-api.open-meteo.com/v1/archive"


def fetch_site(lat, lon, start="2010-07-01", end="2013-05-31", retries=5):
    for a in range(retries):
        r = requests.get(API, params={
        "latitude": lat, "longitude": lon,
        "start_date": start, "end_date": end,
        "hourly": "temperature_2m,precipitation",
        # era5_land lacks precipitation on Open-Meteo; era5_seamless blends ERA5-Land/ERA5
        "models": "era5_seamless"}, timeout=120)
        if r.status_code == 429:
            time.sleep(20 * (a + 1))
            continue
        r.raise_for_status()
        break
    else:
        r.raise_for_status()
    h = r.json()["hourly"]
    return pd.DataFrame({"time": pd.to_datetime(h["time"]),
                         "t2m": h["temperature_2m"], "precip": h["precipitation"]})


def main():
    cfg = load_config()
    sites = site_table(cfg)
    existing = DATA / "interim" / "site_climate.parquet"
    frames = []
    if existing.exists():
        done = pd.read_parquet(existing)
        frames.append(done)
        sites = sites[~sites.SiteID.isin(done.site.unique())]
        print("resuming:", len(sites), "sites remaining")
    for i, s in sites.iterrows():
        try:
            df = fetch_site(s.lat, s.lon)
            df["site"] = s.SiteID
            frames.append(df)
        except Exception as e:
            print("site", s.SiteID, "failed:", repr(e)[:120])
        if i % 25 == 0:
            print("site", i, "of", len(sites))
        time.sleep(0.2)
    out = pd.concat(frames, ignore_index=True)
    p = DATA / "interim" / "site_climate.parquet"
    p.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(p, index=False)
    (DATA / "interim" / "climate_manifest.json").write_text(json.dumps({
        "api": API, "model": "era5_seamless", "variables": ["temperature_2m", "precipitation"],
        "downloaded_utc": datetime.now(timezone.utc).isoformat()}, indent=2))
    print("wrote", p, len(out))


if __name__ == "__main__":
    main()
