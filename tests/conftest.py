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
