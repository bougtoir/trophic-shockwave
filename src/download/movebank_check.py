"""Movement data acquisition check (spec section 3D, 23).

- Movebank live API (direct-read) requires a Movebank user account -> NOT used;
  documented as access limitation.
- Movebank Data Repository is public DSpace: we fetch the wildebeest-Mara
  reference dataset (Stabach et al. 2020, doi:10.5441/001/1.h0t27719/3) which
  includes Mara-population tracks, used for general corridor context only.
Writes data/external/movebank/ + ACCESS_NOTES.json.
"""
import json
from datetime import datetime, timezone
from pathlib import Path
import requests

from src.utils.common import DATA

DSPACE = "https://datarepository.movebank.org/server/api"
ITEM = "5b6706c8-e7e5-46e4-82ba-da5a82324298"  # study item: wildebeest Kenya (3 pops incl. Mara)
KEEP = {"White-bearded wildebeest in Kenya-gps.csv", "README.txt",
        "White-bearded wildebeest in Kenya-reference-data.csv"}


def main():
    outdir = DATA / "external" / "movebank"
    outdir.mkdir(parents=True, exist_ok=True)
    notes = {
        "movebank_direct_read": "requires authenticated Movebank account; not available in this environment",
        "movebank_data_repository": "public; no auth needed for CC0/CC-BY study downloads",
        "serengeti_overlap": "No contemporaneous Serengeti-NP GPS tracking dataset is public on the Movebank Data Repository as of access date. The Mara wildebeest dataset (Stabach et al. 2020) partially overlaps the wider Serengeti-Mara ecosystem (Kenya side) and post-dates the camera study; per spec section 23 it can inform general corridor structure only, never contemporaneous trajectories.",
        "accessed_utc": datetime.now(timezone.utc).isoformat(),
        "downloaded": [],
    }
    try:
        bundles = requests.get(f"{DSPACE}/core/items/{ITEM}/bundles", timeout=60).json()
        for b in bundles["_embedded"]["bundles"]:
            bss = requests.get(b["_links"]["bitstreams"]["href"], timeout=60).json()
            for bs in bss["_embedded"]["bitstreams"]:
                name = bs["name"]
                if name in KEEP:
                    r = requests.get(bs["_links"]["content"]["href"], timeout=300)
                    p = outdir / name
                    p.write_bytes(r.content)
                    notes["downloaded"].append({"file": name, "bytes": len(r.content)})
                    print("downloaded", name, len(r.content))
    except Exception as e:
        notes["error"] = repr(e)
        print("movebank fetch failed:", repr(e)[:200])
    (outdir / "ACCESS_NOTES.json").write_text(json.dumps(notes, indent=2))


if __name__ == "__main__":
    main()
