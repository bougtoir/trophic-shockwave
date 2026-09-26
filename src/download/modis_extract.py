"""Extract MOD13Q1/MYD13Q1 v6.1 NDVI+EVI at camera sites via Planetary Computer STAC.

Reads only the window covering the camera-grid bbox from each scene COG
(vsicurl range requests), computes buffer means per site, and writes
data/interim/site_vegetation.parquet with one row per site x composite.
Quality: pixels masked when VI_Quality indicates non-useful data
(MODIS_Quality 0b11 = pixel not produced/other, or 0b10 = snow/ice)
or pixel_reliability > 2.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.windows import from_bounds

from src.utils.common import DATA, load_config

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def site_table(cfg):
    cons = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv",
                       usecols=["SiteID", "LocationX", "LocationY"])
    sites = cons.drop_duplicates("SiteID").set_index("SiteID")[["LocationX", "LocationY"]]
    # Arc1960 UTM36S -> WGS84
    tr = Transformer.from_crs(f"EPSG:{cfg['study_system']['site_crs_epsg']}", "EPSG:4326",
                              always_xy=True)
    lon, lat = tr.transform(sites.LocationX.values, sites.LocationY.values)
    sites["lon"], sites["lat"] = lon, lat
    return sites.reset_index()


def extract(cfg):
    sites = site_table(cfg)
    bbox = cfg["study_system"]["bbox"]
    catalog = pystac_client.Client.open(STAC)
    items = list(catalog.search(
        collections=["modis-13Q1-061"],  # MODIS Collection 6.1 only
        bbox=bbox,
        datetime="2010-07-01/2013-05-31",
    ).items())
    print(f"{len(items)} MODIS scenes found")
    rows = []
    tr_ll2sin = Transformer.from_crs("EPSG:4326",
        "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs",
        always_xy=True)
    for it in sorted(items, key=lambda i: i.id):
        signed = pc.sign(it)
        platform = "terra" if it.id.startswith("MOD") else "aqua"
        # doy encoded in id: MOD13Q1.A2011033...
        year, doy = int(it.id.split(".")[1][1:5]), int(it.id.split(".")[1][5:8])
        date = pd.Timestamp(year=year, month=1, day=1) + pd.Timedelta(days=doy - 1)
        try:
            arrays = {}
            for asset_key, band_name in [("250m_16_days_NDVI", "NDVI"),
                                         ("250m_16_days_EVI", "EVI"),
                                         ("250m_16_days_pixel_reliability", "rel")]:
                href = signed.assets[asset_key].href
                with rasterio.open(href) as ds:
                    xs, ys = tr_ll2sin.transform(
                        [bbox[0], bbox[0], bbox[2], bbox[2]],
                        [bbox[1], bbox[3], bbox[1], bbox[3]])
                    sb = (min(xs), min(ys), max(xs), max(ys))
                    w0 = from_bounds(*sb, transform=ds.transform)
                    win = rasterio.windows.Window(
                        int(np.floor(w0.col_off)), int(np.floor(w0.row_off)),
                        int(np.ceil(w0.width)), int(np.ceil(w0.height)))
                    if win.width <= 0 or win.height <= 0:
                        continue
                    arr = ds.read(1, window=win)
                    wtr = ds.window_transform(win)
                    arrays[band_name] = (arr, wtr)
            if not arrays:
                continue
            arr_ndvi, wtr = arrays["NDVI"]
            arr_evi = arrays["EVI"][0]
            arr_rel = arrays["rel"][0]
            good = (arr_rel <= 2) & (arr_ndvi > -3000)
            for _, s in sites.iterrows():
                sx, sy = tr_ll2sin.transform(s.lon, s.lat)
                for rad in cfg["spatial"]["sensitivity_buffers_m"] + [cfg["spatial"]["primary_buffer_m"]]:
                    col, row = ~wtr * (sx, sy)
                    px = rad / abs(wtr.a)
                    c0, c1 = int(col - px), int(col + px) + 1
                    r0, r1 = int(row - px), int(row + px) + 1
                    c0, r0 = max(c0, 0), max(r0, 0)
                    c1, r1 = min(c1, arr_ndvi.shape[1]), min(r1, arr_ndvi.shape[0])
                    yy, xx = np.mgrid[r0:r1, c0:c1]
                    # circle mask
                    cx = (xx * wtr.a + wtr.c + wtr.a / 2 - sx) ** 2 + \
                         (yy * wtr.e + wtr.f + wtr.e / 2 - sy) ** 2 <= rad ** 2
                    mask = good[r0:r1, c0:c1] & cx
                    rec = {"site": s.SiteID, "date": date, "platform": platform,
                           "buffer_m": rad, "n_pix": int(mask.sum())}
                    rec["ndvi"] = float(np.nanmean(arr_ndvi[r0:r1, c0:c1][mask])) * 1e-4 if mask.any() else np.nan
                    rec["evi"] = float(np.nanmean(arr_evi[r0:r1, c0:c1][mask])) * 1e-4 if mask.any() else np.nan
                    rows.append(rec)
        except Exception as e:
            import traceback
            print("skip", it.id, repr(e)[:120])
            if os.environ.get("DEBUG"):
                traceback.print_exc()
    df = pd.DataFrame(rows)
    out = DATA / "interim" / "site_vegetation.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    (DATA / "interim" / "modis_manifest.json").write_text(json.dumps({
        "collection": "modis-13Q1-061", "stac": STAC,
        "n_scenes": len(items), "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "item_ids": [i.id for i in items]}, indent=2))
    print("wrote", out, len(df))


def main():
    extract(load_config())


if __name__ == "__main__":
    main()
