"""Download Snapshot Serengeti files from Dryad (doi:10.5061/dryad.5pt92).

Dryad currently gates file downloads behind an Anubis proof-of-work challenge
(sha256(randomData + nonce) with `difficulty` leading zero nibbles). This module
solves the challenge once per session and then downloads the listed files.
Records sha256 + timestamp into data/raw/snapshot_serengeti/MANIFEST.json.
"""
import hashlib, json, re, sys, time
from datetime import datetime, timezone
from pathlib import Path
import requests

from src.utils.common import DATA, load_sources

UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")
BASE = "https://datadryad.org"


def _solve_pow(random_data: str, difficulty: int):
    nonce = 0
    c, odd = difficulty // 2, difficulty % 2
    while True:
        h = hashlib.sha256(f"{random_data}{nonce}".encode()).digest()
        if all(b == 0 for b in h[:c]) and (not odd or h[c] >> 4 == 0):
            return h.hex(), nonce
        nonce += 1


def _pass_anubis(s: requests.Session, url: str):
    r = s.get(url, headers={"User-Agent": UA})
    if not r.content.startswith(b"<!"):
        return r
    m = re.search(r'<script id="anubis_challenge" type="application/json">(.*?)</script>',
                  r.text, re.S)
    if not m:
        return r
    ch = json.loads(m.group(1))
    t0 = time.time()
    hash_hex, nonce = _solve_pow(ch["challenge"]["randomData"], ch["rules"]["difficulty"])
    elapsed = int((time.time() - t0) * 1000)
    s.get(f"{BASE}/.within.website/x/cmd/anubis/api/pass-challenge",
          params={"id": ch["challenge"]["id"], "response": hash_hex,
                  "nonce": nonce, "redir": url, "elapsedTime": elapsed},
          headers={"User-Agent": UA, "Referer": url})
    return s.get(url, headers={"User-Agent": UA, "Referer": BASE})


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    src = load_sources()["data_sources"]["snapshot_serengeti"]
    outdir = DATA / "raw" / "snapshot_serengeti"
    outdir.mkdir(parents=True, exist_ok=True)
    s = requests.Session()
    manifest = {"source": src["doi"], "license": src["license"],
                "downloaded_utc": datetime.now(timezone.utc).isoformat(), "files": {}}
    for key, spec in src["files"].items():
        out = outdir / spec["name"]
        if out.exists() and out.stat().st_size > 1000:
            manifest["files"][spec["name"]] = {"file_id": spec["file_id"],
                                               "sha256": sha256_file(out), "cached": True}
            print(spec["name"], "cached")
            continue
        url = f"{BASE}/downloads/file_stream/{spec['file_id']}"
        r = _pass_anubis(s, url)
        out.write_bytes(r.content)
        manifest["files"][spec["name"]] = {"file_id": spec["file_id"],
                                           "sha256": sha256_file(out), "cached": False}
        print(spec["name"], len(r.content))
    (outdir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
