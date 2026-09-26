# DECISION_LOG

| date | old spec | new spec | reason | outcomes known? |
|---|---|---|---|---|
| 2026-09-20 | CDS ERA5-Land | Open-Meteo archive API (ERA5-Land) | no CDS key in environment; identical product | no |
| 2026-09-20 | Movebank live API | Movebank Data Repository + documentation of restriction | live API requires Movebank account | no |
| 2026-09-20 | elephant/giraffe as species negative controls | kept | not preyed on by focal carnivores; weak expected response | no |
| 2026-09-20 | site x day-of-year vegetation climatology | per-site harmonic (2 harmonics over 16-day bin-of-year) | single-obs-per-cell climatology gave degenerate z-scores bounded ~+/-1.15, zero shocks | no |
| 2026-09-20 | treat 149 site shocks as independent | shock-bin clustered SEs + mandatory shock-bin block bootstrap + LOYO/LOSBO | shocks concentrated in 14/64 bins; user directive | no |
| 2026-09-20 | full-sample co-detection connectivity | split-time cross-fitted W (W for events in half A estimated on half B) | prevents connectivity matrix encoding the tested post-shock signal | no |
| 2026-09-20 | event bootstrap for IRF | shock-bin block bootstrap | same dependence structure | no |
