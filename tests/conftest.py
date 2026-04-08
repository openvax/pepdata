"""Pytest fixtures: patch IEDB local_path to use bundled test data instead of downloading."""

import os
import pytest

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(autouse=True)
def patch_iedb_paths(monkeypatch):
    """Redirect all IEDB modules to use small bundled CSVs."""
    monkeypatch.setattr(
        "pepdata.iedb.tcell.local_path",
        lambda auto_download=True: os.path.join(FIXTURES_DIR, "tcell_test_fixture.csv"),
    )
    monkeypatch.setattr(
        "pepdata.iedb.mhc.local_path",
        lambda auto_download=True: os.path.join(FIXTURES_DIR, "mhc_test_fixture.csv"),
    )
    # Clear memoize caches so each test gets fresh data from fixtures
    from pepdata.iedb.tcell import load_dataframe as tcell_load
    from pepdata.iedb.mhc import load_dataframe as mhc_load
    if hasattr(tcell_load, '__wrapped__'):
        pass  # functools.wraps doesn't expose the cache
    # The memoize decorator stores in a closure dict — simplest to just
    # not worry about it since tests use nrows= which varies the cache key
