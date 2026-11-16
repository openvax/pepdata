"""The two assay loaders agree about file freshness and ownership."""

from pathlib import Path
import shutil

import pytest

from pepdata.iedb import alleles, mhc, tcell


@pytest.mark.parametrize("module", [mhc, tcell, alleles])
def test_delete_never_downloads_and_is_idempotent(module, monkeypatch, tmp_path):
    path = tmp_path / "assays.csv"
    monkeypatch.setattr(module.cache, "local_path", lambda **kwargs: str(path))

    def fail_download(**kwargs):
        pytest.fail("deleting a local export must never download it")

    monkeypatch.setattr(module.cache, "fetch", fail_download)
    module.delete()
    path.write_text("local export")
    module.delete()
    assert not path.exists()
    module.delete()


@pytest.mark.parametrize("module", [mhc, tcell])
def test_load_returns_independent_frames_and_observes_replaced_file(module, monkeypatch, tmp_path):
    source = tmp_path / "assays.csv"
    fixture = Path(module.local_path())
    shutil.copyfile(fixture, source)
    monkeypatch.setattr(module, "local_path", lambda auto_download=True: str(source))

    first = module.load_dataframe()
    count = len(first)
    assert count > 0
    first.drop(first.index, inplace=True)
    assert len(module.load_dataframe()) == count

    lines = source.read_text().splitlines(keepends=True)
    source.write_text("".join(lines + lines[2:]))
    assert len(module.load_dataframe()) == count * 2


def test_allele_loaders_return_independent_values():
    first = alleles.load_alleles()
    assert first
    first[0].synonyms.add("caller-added-synonym")
    first.clear()
    second = alleles.load_alleles()
    assert second
    assert "caller-added-synonym" not in second[0].synonyms

    first_dict = alleles.load_alleles_dict()
    assert first_dict
    first_dict.clear()
    assert alleles.load_alleles_dict()
