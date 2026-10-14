"""Use synthetic IEDB-format fixtures instead of live downloads in unit tests."""

import os
import pytest

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(autouse=True)
def patch_iedb_paths(monkeypatch):
    """Redirect assay tables and the allele catalogue to bundled fixtures."""
    monkeypatch.setattr(
        "pepdata.iedb.alleles.local_path",
        lambda force_download=False: os.path.join(FIXTURES_DIR, "alleles_test_fixture.xml"),
    )
    monkeypatch.setattr(
        "pepdata.iedb.tcell.local_path",
        lambda auto_download=True: os.path.join(FIXTURES_DIR, "tcell_test_fixture.csv"),
    )
    monkeypatch.setattr(
        "pepdata.iedb.mhc.local_path",
        lambda auto_download=True: os.path.join(FIXTURES_DIR, "mhc_test_fixture.csv"),
    )
