# ANALYSIS_PLAN (prespecified)

Primary unit: site x 16-day bin (MOD13Q1 grid). Primary buffer: 1 km.
Primary shock: EVI z <= -2 vs site x day-of-year climatology.
Primary response: encounter_rate = capture events / camera-days (effort from
search_effort.csv deployment windows).
Primary species: herbivores = wildebeest, zebra; carnivores = lion (F+M pooled),
spotted hyena. Negative-control species: elephant, giraffe.

Tests:
- A: event study, site+bin FE, lags -12..+24, ref lag -1, site-clustered SEs.
- B: recipient herbivore anomaly (lags 2-8) on distance + movement connectivity,
  clustered by event.
- C: remote carnivore event-study + distributed-lag model comparing
  W_geo / W_neighbor / W_movement / W_habitat by AIC and blocked-by-site CV RMSE.
- Negative controls: temporal placebos (shock shifted -8/-6/-4 bins),
  species controls, matched equidistant-unconnected recipients,
  column-permuted movement network (200 reps).
- Diagnostics: Moran's I and binned variograms of residuals.
- Simulation: scenarios A-H as coded in src/simulation/simulate.py.

Decisions taken before seeing inferential outcomes:
- EVI chosen as primary index alongside NDVI; both reported (no cherry-pick).
- Region-wide shocks (>=80% of sites in a bin) excluded from the principal
  causal analysis; reported separately.
