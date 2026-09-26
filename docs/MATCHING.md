# Exchangeability: clustered vs isolated

```
                                                                   0
n_events                                                          89
n_isolated                                                        26
n_clustered                                                       63
ps_overlap_range           [0.08141808516632819, 0.6733822751042782]
n_clustered_in_ps_support                                         57
nn_caliper                                                  1.240695
n_clustered_matchable                                             44
frac_clustered_matchable                                    0.698413
```

## Propensity model (isolated ~ covariates)

                  Coef.  Std.Err.         z     P>|z|    [0.025    0.975]
Intercept     -1.040686  0.274977 -3.784634  0.000154 -1.579630 -0.501741
season        -0.823338  0.293168 -2.808416  0.004979 -1.397938 -0.248739
magnitude     -0.471668  0.306150 -1.540644  0.123403 -1.071712  0.128375
duration_bins -0.044197  0.323490 -0.136625  0.891327 -0.678225  0.589832
y_pre         -0.610257  0.359956 -1.695362  0.090007 -1.315758  0.095245
kc            -0.337680  0.277832 -1.215410  0.224210 -0.882221  0.206861

## Re-estimation results

| analysis               |   n |   post_mean |
|:-----------------------|----:|------------:|
| all_pooled             |  89 | -0.00994834 |
| exclude_2010           |  41 | -0.0612573  |
| isolated_no_pretrend   |  14 |  0.691199   |
| clustered_matched_only |  44 | -0.236485   |
| isolated_matched       |  19 |  0.522955   |
| LOSBO_drop_bin10       | 147 |  0.0214199  |
| LOSBO_drop_bin20       |  79 |  0.0801136  |
| LOSBO_drop_bin21       | 138 |  0.00990836 |
| LOSBO_drop_bin22       | 148 |  0.032908   |
| LOSBO_drop_bin25       | 146 |  0.0524722  |
| LOSBO_drop_bin26       | 139 | -0.0237551  |
| LOSBO_drop_bin29       | 148 |  0.0329324  |
| LOSBO_drop_bin30       | 144 |  0.0651701  |
| LOSBO_drop_bin31       | 121 |  0.270084   |
| LOSBO_drop_bin32       | 137 | -0.0735979  |
| LOSBO_drop_bin33       | 147 |  0.0194471  |
| LOSBO_drop_bin48       | 148 |  0.0466717  |
| LOSBO_drop_bin50       | 147 |  0.0225318  |
| LOSBO_drop_bin66       | 148 |  0.0164466  |
