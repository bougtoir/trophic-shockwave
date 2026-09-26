"""Shared helpers: config loading, paths, binning."""
from pathlib import Path
import yaml
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"


def load_config():
    return yaml.safe_load(open(ROOT / "config" / "config.yaml"))


def load_sources():
    return yaml.safe_load(open(ROOT / "config" / "data_sources.yaml"))


def bin_index(dates: pd.Series, bin_days: int, origin: str = "2010-01-01") -> pd.Series:
    """16-day (or bin_days) bin index anchored at MODIS-compatible origin."""
    o = pd.Timestamp(origin)
    return ((pd.to_datetime(dates) - o).dt.days // bin_days).astype(int)


def bin_start_date(idx: pd.Series, bin_days: int, origin: str = "2010-01-01") -> pd.Series:
    o = pd.Timestamp(origin)
    return o + pd.to_timedelta(idx * bin_days, unit="D")


def doy_of_bin(idx: pd.Series, bin_days: int, origin: str = "2010-01-01") -> pd.Series:
    return bin_start_date(idx, bin_days, origin).dt.dayofyear
