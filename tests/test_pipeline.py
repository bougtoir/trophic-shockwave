import pandas as pd
import numpy as np
from pathlib import Path
from src.utils.common import DATA, bin_index


def test_bin_index_monotonic():
    d = pd.Series(pd.date_range("2010-01-01", "2010-03-01"))
    b = bin_index(d, 16)
    assert b.is_monotonic_increasing
    assert b.iloc[-1] - b.iloc[0] == 3


def test_panel_structure():
    p = DATA / "processed" / "panel.parquet"
    if not p.exists():
        return
    df = pd.read_parquet(p)
    assert {"site", "bin", "species", "cam_days", "detection",
            "encounters", "encounter_rate"} <= set(df.columns)
    assert (df.encounter_rate[df.cam_days > 0] >= 0).all()


def test_shock_indicator_binary():
    p = DATA / "processed" / "site_shock.parquet"
    if not p.exists():
        return
    s = pd.read_parquet(p)
    assert set(s.shock.unique()) <= {True, False}


def test_connectivity_normalized():
    for name in ["geo", "neighbor", "movement", "habitat"]:
        f = DATA / "processed" / f"W_{name}.npy"
        if not f.exists():
            continue
        W = np.load(f)
        np.testing.assert_allclose(W.sum(1), 1.0, atol=1e-6)
        assert (np.diag(W) == 0).all()
