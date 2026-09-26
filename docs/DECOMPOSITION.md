# Decomposition: clustered vs isolated sign reversal

Per-event source-site herbivore anomaly post mean (lags 2-8), by set:

| var          |   isolated_mean |   clustered_mean |         diff |   n_bins_iso |   n_bins_clu |
|:-------------|----------------:|-----------------:|-------------:|-------------:|-------------:|
| y_pre        |      -0.840202  |       -0.139661  |  -0.700541   |            9 |            2 |
| count_pre    |     173.107     |      126.896     |  46.2116     |            8 |            2 |
| flow_real    |      -0.100975  |       -0.0584422 |  -0.0425328  |           12 |            2 |
| outflow      |       0.992424  |        0.992857  |  -0.0004329  |           12 |            2 |
| kc           |       0.251211  |        0.255313  |  -0.00410175 |           12 |            2 |
| predation    |       0.0687884 |        0.0651916 |   0.00359682 |            9 |            2 |
| n_nb_shocked |       2.60076   |       20.2143    | -17.6135     |           12 |            2 |
| y_post       |      -0.0704922 |       -0.358797  |   0.288305   |           10 |            2 |

## Interaction model (y_post ~ isolated * covariates)

| term                          |              coef |          se |           p |
|:------------------------------|------------------:|------------:|------------:|
| Intercept                     |      -0.299472    | 0.022099    | 7.77674e-42 |
| isolated[T.True]              |       0.310898    | 0.19715     | 0.114806    |
| y_pre                         |      -0.187562    | 0.259213    | 0.469322    |
| isolated[T.True]:y_pre        |       0.350674    | 0.278514    | 0.207999    |
| count_pre                     |      -0.0232867   | 0.0150063   | 0.120711    |
| isolated[T.True]:count_pre    |       0.158475    | 0.236178    | 0.502223    |
| flow_real                     |       0.463737    | 0.205023    | 0.023705    |
| isolated[T.True]:flow_real    |      -0.297268    | 0.325211    | 0.360676    |
| outflow                       |       2.00378e+06 | 1.68299e+06 | 0.233809    |
| isolated[T.True]:outflow      | -522169           | 3.74913e+06 | 0.889231    |
| kc                            |       0.0325944   | 0.0452306   | 0.47114     |
| isolated[T.True]:kc           |       0.230099    | 0.197562    | 0.244143    |
| predation                     |      -0.0226044   | 0.047564    | 0.634614    |
| isolated[T.True]:predation    |      -0.0700397   | 0.071734    | 0.328876    |
| n_nb_shocked                  |       0.000394913 | 0.311284    | 0.998988    |
| isolated[T.True]:n_nb_shocked |      -0.0958332   | 0.406413    | 0.813586    |
