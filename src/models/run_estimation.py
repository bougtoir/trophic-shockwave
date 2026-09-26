"""Driver: runs each (shock_set, species) event study in a fresh subprocess to
bound memory, then carnivore model comparison. Skips combos whose output CSV
already exists (resumable)."""
import subprocess, sys, os

COMBOS = []
SETS = ["isolated", "clustered", "all"]
SPECIES = ["wildebeest", "zebra", "lion", "hyenaSpotted",
           "cheetah", "leopard", "elephant", "giraffe"]

SNIPPET = """
import sys
sys.path.insert(0, '.')
from src.models.estimation import *
cfg = load_config()
panel, ev, veg, clim, shock = load_all()
ev = ev[~ev.regionwide]
sets = {"all": ev, "clustered": ev[~ev.isolated], "isolated": ev[ev.isolated]}
evs = sets[sys.argv[1]]
sp = sys.argv[2]
if len(evs) >= 8:
    d = stack_events(panel, evs, sp)
    coefs, m = fit_es(d)
    if sp in ('wildebeest','zebra','lion','hyenaSpotted') and sys.argv[1] in ('all','isolated'):
        bb = block_boot(panel, evs, sp, nboot=150, seed=cfg['seed'])
        coefs = coefs.merge(bb, on='lag', how='left')
    coefs['species'] = sp; coefs['shock_set'] = sys.argv[1]; coefs['n_events'] = len(evs)
    TAB.mkdir(parents=True, exist_ok=True)
    coefs.to_csv(TAB / f'es_{sys.argv[1]}_{sp}.csv', index=False)
    print('post_mean_2_8', coefs[coefs.lag.between(2,8)].coef.mean())
"""

for s in SETS:
    for sp in SPECIES:
        out = f"results/tables/es_{s}_{sp}.csv"
        if os.path.exists(out):
            print("skip", out)
            continue
        print("RUN", s, sp, flush=True)
        r = subprocess.run([sys.executable, "-c", SNIPPET, s, sp],
                           capture_output=True, text=True)
        print(r.stdout[-300:], r.stderr[-300:], flush=True)

# carnivore model comparison + effective replication + summary in-process
print("RUN carnivore comparison", flush=True)
r = subprocess.run([sys.executable, "-c", """
import sys; sys.path.insert(0,'.')
import json, numpy as np, pandas as pd
from src.models.estimation import *
from src.preprocess.connectivity import movement_W
cfg = load_config()
panel, ev, veg, clim, shock = load_all()
ev = ev[~ev.regionwide]
sets = {"all": ev, "clustered": ev[~ev.isolated], "isolated": ev[ev.isolated]}
sites = np.load(DATA/"processed"/"sites_order.npy", allow_pickle=True)
Wg = np.load(DATA/"processed"/"W_geo.npy")
split_bin = ev.bin_start.median()
WA = movement_W(panel, sites, mask_bins=range(int(split_bin), 76))
WB = movement_W(panel, sites, mask_bins=range(0, int(split_bin)))
comp_rows = []
MOD.mkdir(parents=True, exist_ok=True)
for sname, evs in sets.items():
    sub_shock = shock.copy()
    ok = set()
    for r_ in evs.itertuples():
        for b in range(r_.bin_start, r_.bin_end + 1):
            ok.add((r_.site, b))
    sub_shock["shock"] = sub_shock.shock & sub_shock.apply(lambda rr: (rr.site, rr.bin) in ok, axis=1)
    for (lbl, W1, W2) in [("geo", Wg, None), ("eco", None, None), ("combined", None, Wg)]:
        try:
            res, m, d = carnivore_model(panel, sub_shock, veg, clim, W1, sites,
                                        cfg["species"]["carnivores_primary"], lbl,
                                        combined_with=W2, ev=evs,
                                        WA=(WA if W1 is None else None),
                                        WB=WB, split_bin=split_bin)
            res["shock_set"] = sname; res["connectivity"] = "xfit"
            res["cv_rmse"] = cv_rmse(d, seed=cfg["seed"])
            comp_rows.append(res)
            with open(MOD / f"carn_{sname}_{lbl}.json", "w") as f:
                json.dump({k: v for k, v in res.items() if k != "coefs"} | res["coefs"], f, indent=2)
            print("carn", sname, lbl, "aic=%.0f" % res["aic"], flush=True)
        except Exception as e:
            print("carn fail", sname, lbl, repr(e)[:200], flush=True)
pd.DataFrame(comp_rows).to_csv(TAB/"connectivity_model_comparison.csv", index=False)

# es_summary from written files
rows = []
import glob
for f in glob.glob(str(TAB/"es_*.csv")):
    c = pd.read_csv(f)
    rows.append({"set": c.shock_set.iloc[0], "species": c.species.iloc[0],
                 "post_mean_2_8": c[c.lag.between(2,8)].coef.mean(),
                 "n_events": c.n_events.iloc[0]})
pd.DataFrame(rows).to_csv(TAB/"es_summary.csv", index=False)
print("done")
"""], capture_output=True, text=True)
print(r.stdout[-2000:], r.stderr[-500:], flush=True)
