# GO_NO_GO — Phase 1 checkpoint (pre-inferential)

Written after data audit and shock-event construction, before inferential results.

## Independent local vegetation shocks (primary definition: EVI z <= -2)

- Total shock events (EVI): **149**
- Local (spatially heterogeneous, usable for causal test): **149**
- Region-wide (excluded/separately classified): **0**
- Sites with >=1 shock: **115** of 225
- Median event duration: 1 x 16-day bins
- Events by year: {2010: 84, 2011: 61, 2012: 4}

### Rule (prespecified)
- Strong GO requires >1 shock/event/species; >= ~20 independent local events gives reasonable event-study power.
- If independent local shocks < ~10, expand ecosystems before heavy modeling.

---

# Post-estimation evaluation (generated from results/tables)

## Effective replication

```
 total_site_shocks  unique_shock_bins  isolated_shock_events  clustered_shock_events  sites_per_bin_mean  sites_per_bin_median                events_by_year  icc_shock_bin_wildebeest_post  ess_approx
               149                 14                     51                      98           10.642857                   2.5 {2010: 84, 2011: 61, 2012: 4}                       0.015187  129.967329
```

## Criteria

- A (veg shock precedes herbivore response): post-mean(2-8) wildebeest=0.033, zebra=0.598 -> weak/no
- B (structured redistribution): see redistribution_fit.csv
- C (remote carnivore delay): eco exposure sum=0.905
- D (eco > geo): eco_sum=0.905 vs geo_sum=0.000; permutation null in permuted_network.csv
- E (survives shock-bin block bootstrap): see es_<set>_<species>.csv boot_lo/boot_hi columns
- F (not one year/bin): see sensitivity_loyo_loso.csv
- G (isolated subset consistent): isolated post-mean=0.274
- H (negative controls): placebo.csv, permuted_network.csv, elephant/giraffe es_* files