# Trophic Shockwave — Phase 1 Decision Memo (Serengeti)

All numbers below are produced by code under `src/` into `results/tables/`
and `results/models/`; regenerate with `make analysis report`.

## The eight questions

**1. Temporal ordering (veg shock → herbivore response)?**
No clean ordering. Stacked event studies (lags -6..+12, site+bin FE,
clustered by shock bin): wildebeest post-mean(2–8) = +0.03 (all),
+0.27 (isolated); zebra +0.60 / +0.99 — positive, not declines.
Pre-shock coefficients are non-zero in several sets (zebra all: -0.36),
and 0–1 of 7 post lags clear the shock-bin block-bootstrap CI.
→ Criterion A: NOT supported as specified (no local decline; direction
opposite or null under clustered uncertainty).

**2. Spatial redistribution (herbivore flow to connected sites)?**
Weak/none. Source→recipient anomaly regression
(`redistribution_fit.csv`): distance decay ≈ 0 (-0.004/km, ns);
eco-connectivity coefficient not distinguishable once distance-matched
(`w_eco_dm` ns). Moran's I of carnivore residuals ≈ 0 at 6–15 km.
→ B: not demonstrated.

**3. Remote carnivore response with delay?**
Equivocal at best. Distributed-lag exposure models (K=8):
geo exposure sum positive (all: +2.49; clustered: +3.94); eco exposure
sum positive but weaker (all: +0.90 xf). Sign is positive (carnivore
activity *rises* at exposed sites), magnitudes unstable across LOYO
(eco_sum 0.10–1.39). AIC: geo beats eco in all/isolated; combined wins
only for clustered (the 2-bin set). Site-blocked CV RMSE ≈0.30 both.
→ C: weak, sign unstable, geo-favoured.

**4. Ecological connectivity > geography?**
No. Cross-fitted eco matrix adds little; combined-model eco coefs small;
permutation null (degree-preserving site relabel, 200 perms): observed
eco_sum=1.00 vs perm mean -0.16 (sd 1.42), p(perm>=obs)=0.21 —
within the null. → D: rejected.

**5. Survives clustered uncertainty?**
No. Shock-bin block bootstrap CIs are wide; effective temporal
replication = 14 shock bins (ESS ≈130 site-events but only 14 time
slots). → E: fails.

**6. Isolated subset consistent?**
No — and per user's hypothesis, clustered vs isolated are NOT
exchangeable: clustered events live in ~2 shock bins (n_nb_shocked
≈20 vs 2.6), isolated carry a strong pre-existing decline
(y_pre -0.84 vs -0.14). Propensity-matched subsets still diverge
(clustered-matched -0.24 vs isolated-matched +0.52), but only ~70%
of clustered events are matchable and the matchable subset is drawn
from the same 2 bins — residual confounding by bin dominates.
→ G: fails; interpreted as different event populations, not a
spatial-clustering mechanism.

**7. Dependent on 2010/2011?**
Yes for the negative herbivore signal (clustered decline vanishes
outside those bins); LOSBO post-means mostly +0.02–0.08 with two bins
(31, 32) driving deviations. LOYO: unstable sign/magnitude.
→ F: fails.

**8. Strongest falsification result?**
Temporal placebos: shocks shifted earlier by 4–8 bins produce
post-means +0.27 to +0.45 — LARGER than the observed all-set
wildebeest estimate (+0.03). Negative-control species
(elephant +0.006, giraffe +0.023) show the same order of "response"
as focal carnivores (lion +0.007, hyena +0.024). The pipeline finds
the same signature everywhere it should NOT exist.
→ H: fails.

## Verdict

**NO GO** (on the "spatiotemporal propagation signature" hypothesis as
specified — this is the falsification outcome the design allows).

- Temporal ordering: absent/unstable
- Spatial redistribution: not detected
- Remote carnivore delay: geo-exposure positive but eco adds nothing
  beyond permutation null; sign positive (aggregation, not flight)
- Eco > geo: rejected (p≈0.21 under degree-preserving permutation)
- Clustered uncertainty: fails
- Isolated subset: opposite direction; clustered/isolated shown to be
  non-exchangeable event populations (2-bin concentration + pre-trend)
- Dependence on 2010/2011: yes
- Strongest falsification: placebo shifts reproduce/exceed the effect

## Caveats for the record

- `encounter_rate` at the shocked site rising post-shock is consistent
  with herbivore concentration on remaining green patches OR detection
  artefacts during sparse effort; the data cannot separate these.
- Movebank contemporaneous GPS was unavailable; W_movement is
  camera co-detection only (external GPS used as structural info only).
- Camera grid covers a small window of wildebeest seasonal range —
  redistributive exits may occur outside the grid (unobservable here).

## If pursued further (not done — per stop condition)

Richer systems with GPS + cameras (Serengeti full extent via Snapshot
Safari, or fenced reserves), individual-level data, or GRAD-based
indices rather than EVI anomalies.
