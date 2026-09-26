# Trophic Shockwave

**Do local vegetation shocks propagate through mobile terrestrial food webs with a delay?**

Phase-1 falsification-oriented GO/NO-GO analysis on the Serengeti ecosystem
(Snapshot Serengeti camera traps × MODIS MOD13Q1 v6.1 NDVI/EVI × ERA5-Land).

## Reproduce

```bash
conda env create -f environment.yml && conda activate trophic-shockwave
make data        # download + preprocess
make analysis    # event studies, distributed-lag W-comparison, negative controls
make simulate    # synthetic-data validation of the estimator
make figures     # figures A-H and tables 1-7
make report      # DATA_AUDIT / GO_NO_GO rendering
```

## Layout

- `config/` — prespecified analysis parameters and data sources
- `src/download/` — data acquisition (Dryad, Planetary Computer MODIS, Open-Meteo ERA5-Land, Movebank)
- `src/preprocess/` — camera panel, vegetation anomalies, shock events, connectivity matrices
- `src/models/` — event studies, distributed-lag spatial models, IRF, negative controls, autocorrelation
- `src/simulation/` — synthetic validation scenarios A–H
- `src/figures/` — figure/table generation
- `docs/` — DATA_AUDIT, ANALYSIS_PLAN, ASSUMPTIONS, DECISION_LOG, GO_NO_GO, LITERATURE_MATRIX
- `results/` — generated only; never hand-edited
- `manuscript/outline.md` — manuscript skeleton only

## Data

See `config/data_sources.yaml` and `docs/DATA_AUDIT.md`. Primary: Snapshot
Serengeti (Dryad doi:10.5061/dryad.5pt92, CC0), MODIS MOD13Q1/MYD13Q1 Collection
6.1 (Planetary Computer, NASA public domain), ERA5-Land via Open-Meteo (CC-BY),
Movebank Data Repository wildebeest-Mara (Stabach et al. 2020, CC0).
