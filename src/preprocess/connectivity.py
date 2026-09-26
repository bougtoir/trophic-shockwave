"""Connectivity matrices (spec section 14/23).

Outputs data/processed/W_{geo,neighbor,movement,habitat}.npy (row-normalized)
plus data/processed/site_distances.csv.

- W_geo: Gaussian kernel exp(-(d/l0)^2), l0 = median pairwise distance.
- W_neighbor: rook adjacency on the 5km grid (sites ~5km apart; threshold = 1.5x min dist).
- W_movement: camera-derived directional co-occurrence of focal herbivores:
  fraction of bins where herbivore presence at i in bin t predicts presence at j in t+1
  (symmetrized), beyond site x season baseline. Fallback = distance decay.
- W_habitat: correlation of site NDVI seasonal profiles (shared habitat suitability).
"""
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from src.utils.common import DATA, load_config

FOCAL = ["wildebeest", "zebra"]


def movement_W(panel, sites, mask_bins=None):
    """Camera-derived connectivity from time-lagged co-detection of focal herbivores.
    mask_bins: iterable of bins to EXCLUDE (cross-fitting: drop shock bins or a
    time subset so the matrix can't encode the post-shock response signal)."""
    foc = panel[panel.species.isin(FOCAL)].groupby(["site", "bin"])["detection"].max().unstack("site")
    foc = foc.reindex(columns=sites).fillna(0)
    if mask_bins is not None:
        foc = foc.loc[~foc.index.isin(list(mask_bins))]
    foc = foc.values
    M = np.zeros((len(sites), len(sites)))
    a = foc[:-1, :]; b = foc[1:, :]
    num = a.T @ b / a.shape[0]
    pij = num / (a.mean(0)[:, None] + 1e-9)
    M = np.clip(pij - b.mean(0)[None, :], 0, None)
    M = (M + M.T) / 2
    return row_norm(M)


def row_norm(W):
    W = W.copy().astype(float)
    np.fill_diagonal(W, 0.0)
    rs = W.sum(1, keepdims=True)
    rs[rs == 0] = 1.0
    return W / rs


def main():
    cfg = load_config()
    panel = pd.read_parquet(DATA / "processed" / "panel.parquet")
    veg = pd.read_parquet(DATA / "processed" / "vegetation.parquet")
    sites = pd.read_parquet(DATA / "processed" / "panel.parquet").site.unique()
    coords = pd.read_csv(DATA / "raw" / "snapshot_serengeti" / "consensus_data.csv",
                         usecols=["SiteID", "LocationX", "LocationY"]).drop_duplicates("SiteID")
    coords = coords.set_index("SiteID").loc[sites]
    D = cdist(coords[["LocationX", "LocationY"]], coords[["LocationX", "LocationY"]]) / 1000.0
    pd.DataFrame(D, index=sites, columns=sites).to_csv(DATA / "processed" / "site_distances.csv")

    np.save(DATA / "processed" / "sites_order.npy", sites)
    l0 = np.median(D[D > 0])
    W_geo = row_norm(np.exp(-(D / l0) ** 2))
    dmin = D[D > 0].min()
    W_neighbor = row_norm((D <= dmin * 1.5).astype(float))

    # W_movement: lagged co-detection asymmetry for focal herbivores
    foc = panel[panel.species.isin(FOCAL)].groupby(["site", "bin"])["detection"].max().unstack("site")
    foc = foc.reindex(columns=sites).fillna(0).values  # bins x sites
    M = np.zeros((len(sites), len(sites)))
    for lag in (1,):
        a = foc[:-lag, :]  # bins x sites presence at t
        b = foc[lag:, :]   # bins x sites presence at t+lag
        num = a.T @ b / a.shape[0]                    # E[S_i,t * S_j,t+lag] -> (i,j)
        pij = num / (a.mean(0)[:, None] + 1e-9)       # P(j active | i active)
        M += np.clip(pij - b.mean(0)[None, :], 0, None)
    M = (M + M.T) / 2
    W_movement = row_norm(M)

    # W_habitat: correlation of seasonal ndvi profile
    seas = veg.groupby(["site", "biny"]).ndvi.mean().unstack("biny").reindex(sites)
    seas = seas.T.fillna(seas.mean(axis=1)).T
    C = np.corrcoef(seas.values)
    W_habitat = row_norm(np.clip(C, 0, None))

    np.save(DATA / "processed" / "W_geo.npy", W_geo)
    np.save(DATA / "processed" / "W_neighbor.npy", W_neighbor)
    np.save(DATA / "processed" / "W_movement.npy", W_movement)
    np.save(DATA / "processed" / "W_habitat.npy", W_habitat)
    print("connectivity done: l0=%.1f km, min dist=%.2f km" % (l0, dmin))


if __name__ == "__main__":
    main()
