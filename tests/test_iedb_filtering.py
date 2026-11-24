# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Offline parity tests for unified IEDB MHC and T-cell table filtering."""

from pathlib import Path

import pandas as pd
import pytest

from pepdata import iedb
from pepdata.iedb import alleles, columns, common, mhc, tcell


# Bound at import, before conftest's autouse fixture replaces it with a fixture path.
original_tcell_local_path = tcell.local_path

PARITY_CSV = """\
Epitope,MHC Restriction,MHC Restriction,Host,Assay
Name,Name,Class,Name,Method
acdefghik,HLA-A*02:01,I,Homo sapiens,ELISPOT
ACDEFGHIL,HLA-A2,I,Homo sapiens (human),purified MHC/competitive
ACDEFGHIM,HLA-A*02:01,I,Mus musculus,ELISPOT
ACDEFGHIN,H2-Kb,I,Homo sapiens,ELISA
ACDEFGHIPQ,HLA-DRB1*01:01,II,Homo sapiens,
AcDeFgHiKlMnOpQ,HLA-DRB1*04:01,II,Homo sapiens,ELISA
ACDEFGHV,H2-IAb,II,Mus musculus,ELISPOT
ACDEFXGHK,HLA-A*02:01,I,Homo sapiens,ELISPOT
acdefughk,HLA-A*02:01,I,Homo sapiens,ELISPOT
ACDEFBGHK,HLA-A*02:01,I,Homo sapiens,ELISPOT
ACDEFJGHK,HLA-A*02:01,I,Homo sapiens,ELISPOT
ACDEFZGHK,HLA-A*02:01,I,Homo sapiens,ELISPOT
,HLA-A*02:01,I,Homo sapiens,ELISPOT
SIINFEKL,,I,Homo sapiens,ELISPOT
"""

LEGACY_HEADERS_CSV = """\
Epitope,MHC,MHC,Host,Assay
Description,Allele Name,MHC allele class,Name,Method/Technique
siinfekl,HLA-A*02:01,I,Homo sapiens,ELISPOT
ACDEFGHIKLMNPQR,HLA-DRB1*01:01,II,Mus musculus,ELISA
"""


@pytest.fixture
def parity_csv_path(monkeypatch, tmp_path) -> Path:
    csv_path = tmp_path / "iedb_parity.csv"
    csv_path.write_text(PARITY_CSV)
    monkeypatch.setattr(mhc, "local_path", lambda auto_download=True: str(csv_path))
    monkeypatch.setattr(tcell, "local_path", lambda auto_download=True: str(csv_path))
    return csv_path


def assert_loaders_match(parity_csv_path, **kwargs) -> tuple[pd.DataFrame, pd.DataFrame]:
    df_mhc = mhc.load_dataframe(**kwargs)
    df_tcell = tcell.load_dataframe(**kwargs)
    pd.testing.assert_frame_equal(df_mhc, df_tcell)

    raw_df = common.read_assay_csv(parity_csv_path)
    df_direct = common.filter_assay_table(raw_df, **kwargs)
    pd.testing.assert_frame_equal(df_mhc, df_direct)
    return df_mhc, df_tcell


def test_epitopes_are_uppercased_and_nulls_dropped_consistently(parity_csv_path):
    df_mhc, _ = assert_loaders_match(parity_csv_path)
    epitopes = list(columns.get_epitope_name(df_mhc))
    assert epitopes == [
        "ACDEFGHIK",
        "ACDEFGHIL",
        "ACDEFGHIM",
        "ACDEFGHIN",
        "ACDEFGHIPQ",
        "ACDEFGHIKLMNOPQ",
        "ACDEFGHV",
        "SIINFEKL",
    ]
    assert all(seq == seq.upper() for seq in epitopes)


def test_only_standard_amino_acids_parity(parity_csv_path):
    df_strict, _ = assert_loaders_match(parity_csv_path, only_standard_amino_acids=True)
    df_all, _ = assert_loaders_match(parity_csv_path, only_standard_amino_acids=False)

    strict_epitopes = set(columns.get_epitope_name(df_strict))
    all_epitopes = set(columns.get_epitope_name(df_all))

    # Ambiguous/non-standard residues (including lowercase 'u') are dropped when strict
    # and retained (upper-cased) when only_standard_amino_acids=False.
    nonstandard = {"ACDEFXGHK", "ACDEFUGHK", "ACDEFBGHK", "ACDEFJGHK", "ACDEFZGHK"}
    assert strict_epitopes.isdisjoint(nonstandard)
    assert nonstandard.issubset(all_epitopes)
    # Null epitope rows are always dropped.
    assert columns.get_epitope_name(df_all).isna().sum() == 0


def test_human_only_filters_on_host_organism_not_allele_prefix(parity_csv_path):
    df_human, _ = assert_loaders_match(parity_csv_path, human_only=True)
    epitopes = set(columns.get_epitope_name(df_human))
    # Includes Homo sapiens + Homo sapiens (human), even with non-HLA allele H2-Kb
    assert "ACDEFGHIK" in epitopes
    assert "ACDEFGHIL" in epitopes
    assert "ACDEFGHIN" in epitopes
    # Excludes HLA-transgenic Mus musculus row (HLA-A*02:01 in mouse host)
    assert "ACDEFGHIM" not in epitopes
    assert "ACDEFGHV" not in epitopes


@pytest.mark.parametrize(
    ("class_arg", "expected_epitopes"),
    [
        (1, {"ACDEFGHIK", "ACDEFGHIL", "ACDEFGHIM", "ACDEFGHIN", "SIINFEKL"}),
        ("I", {"ACDEFGHIK", "ACDEFGHIL", "ACDEFGHIM", "ACDEFGHIN", "SIINFEKL"}),
        (2, {"ACDEFGHIPQ", "ACDEFGHIKLMNOPQ", "ACDEFGHV"}),
        ("II", {"ACDEFGHIPQ", "ACDEFGHIKLMNOPQ", "ACDEFGHV"}),
    ],
)
def test_mhc_class_parity_and_no_allele_xml_access(
    parity_csv_path, monkeypatch, class_arg, expected_epitopes
):
    def fail_alleles(*args, **kwargs):
        pytest.fail("mhc_class filtering must not load or download allele XML")

    monkeypatch.setattr(alleles, "local_path", fail_alleles)
    monkeypatch.setattr(alleles, "load_alleles", fail_alleles)
    monkeypatch.setattr(alleles, "load_alleles_dict", fail_alleles)

    df_mhc, _ = assert_loaders_match(parity_csv_path, mhc_class=class_arg)
    assert set(columns.get_epitope_name(df_mhc)) == expected_epitopes


@pytest.mark.parametrize("bad_class", [0, 3, "1", "2", "i", "ii", "III", "", True, False])
@pytest.mark.parametrize("loader", [mhc.load_dataframe, tcell.load_dataframe])
def test_invalid_mhc_class_raises_value_error(parity_csv_path, loader, bad_class):
    with pytest.raises(ValueError, match="Invalid MHC class"):
        loader(mhc_class=bad_class)


def test_assay_method_handles_na_and_matches_regex_consistently(parity_csv_path):
    df_elispot, _ = assert_loaders_match(parity_csv_path, assay_method="ELISPOT")
    assert set(columns.get_epitope_name(df_elispot)) == {
        "ACDEFGHIK",
        "ACDEFGHIM",
        "ACDEFGHV",
        "SIINFEKL",
    }

    df_elisa_or_comp, _ = assert_loaders_match(
        parity_csv_path, assay_method="ELISA|competitive"
    )
    assert set(columns.get_epitope_name(df_elisa_or_comp)) == {
        "ACDEFGHIL",
        "ACDEFGHIN",
        "ACDEFGHIKLMNOPQ",
    }
    # Row with missing assay method (ACDEFGHIPQ) never matches and does not raise.
    assert "ACDEFGHIPQ" not in set(columns.get_epitope_name(df_elisa_or_comp))


def test_hla_and_exclude_hla_parity(parity_csv_path):
    df_hla, _ = assert_loaders_match(parity_csv_path, hla=r"HLA-A\*02:01")
    assert set(columns.get_epitope_name(df_hla)) == {"ACDEFGHIK", "ACDEFGHIM"}

    df_exclude, _ = assert_loaders_match(parity_csv_path, exclude_hla=r"HLA-A\*02:01")
    excluded_epitopes = set(columns.get_epitope_name(df_exclude))
    assert "ACDEFGHIK" not in excluded_epitopes
    assert "ACDEFGHIM" not in excluded_epitopes
    # Row with missing allele name is kept by exclude_hla (na=False)
    assert "SIINFEKL" in excluded_epitopes

    df_both, _ = assert_loaders_match(
        parity_csv_path, hla="HLA-DRB1", exclude_hla=r"04:01"
    )
    assert set(columns.get_epitope_name(df_both)) == {"ACDEFGHIPQ"}


def test_peptide_length_parity_and_validation(parity_csv_path):
    df_8, _ = assert_loaders_match(parity_csv_path, peptide_length=8)
    assert set(columns.get_epitope_name(df_8)) == {"ACDEFGHV", "SIINFEKL"}

    df_15, _ = assert_loaders_match(parity_csv_path, peptide_length=15)
    assert set(columns.get_epitope_name(df_15)) == {"ACDEFGHIKLMNOPQ"}

    for bad_length in [0, -1, True]:
        with pytest.raises(ValueError, match="peptide_length"):
            mhc.load_dataframe(peptide_length=bad_length)
        with pytest.raises(ValueError, match="peptide_length"):
            tcell.load_dataframe(peptide_length=bad_length)


def test_legacy_column_headers_supported_by_both_loaders(monkeypatch, tmp_path):
    csv_path = tmp_path / "legacy_iedb.csv"
    csv_path.write_text(LEGACY_HEADERS_CSV)
    monkeypatch.setattr(mhc, "local_path", lambda auto_download=True: str(csv_path))
    monkeypatch.setattr(tcell, "local_path", lambda auto_download=True: str(csv_path))

    df_mhc = mhc.load_dataframe(mhc_class=1, human_only=True, assay_method="ELISPOT")
    df_tcell = tcell.load_dataframe(mhc_class="I", human_only=True, assay_method="ELISPOT")
    pd.testing.assert_frame_equal(df_mhc, df_tcell)
    assert list(columns.get_epitope_name(df_mhc)) == ["SIINFEKL"]


def test_filter_assay_table_does_not_mutate_input_dataframe(parity_csv_path):
    raw_df = common.read_assay_csv(parity_csv_path)
    original_first_epitope = columns.get_epitope_name(raw_df).iloc[0]
    assert original_first_epitope == "acdefghik"

    filtered = common.filter_assay_table(raw_df, mhc_class=1)
    assert columns.get_epitope_name(filtered).iloc[0] == "ACDEFGHIK"
    assert columns.get_epitope_name(raw_df).iloc[0] == "acdefghik"


def test_trailing_comma_and_all_nan_assay_method_column(tmp_path):
    csv_path = tmp_path / "trailing_comma.csv"
    csv_path.write_text(
        "Epitope,MHC Restriction,MHC Restriction,Host,Assay,\n"
        "Name,Name,Class,Name,Method,\n"
        "ACDEFGHIK,HLA-A*02:01,I,Homo sapiens,,\n"
    )
    df = common.read_assay_csv(csv_path)
    assert ("Assay", "Method") in df.columns
    assert all("Unnamed:" not in str(part) for col in df.columns for part in col)

    filtered = common.filter_assay_table(df, assay_method="ELISA")
    assert len(filtered) == 0


def test_missing_required_filter_columns_raise_value_error():
    minimal_df = pd.DataFrame({("Epitope", "Name"): ["ACDEFGHIK"]})
    # Works when optional filters are not requested
    assert len(common.filter_assay_table(minimal_df)) == 1

    with pytest.raises(ValueError, match="host organism"):
        common.filter_assay_table(minimal_df, human_only=True)
    with pytest.raises(ValueError, match="MHC class"):
        common.filter_assay_table(minimal_df, mhc_class=1)
    with pytest.raises(ValueError, match="MHC allele"):
        common.filter_assay_table(minimal_df, hla="HLA-A2")
    with pytest.raises(ValueError, match="assay method"):
        common.filter_assay_table(minimal_df, assay_method="ELISPOT")

    no_epitope_df = pd.DataFrame({("Host", "Name"): ["Homo sapiens"]})
    with pytest.raises(ValueError, match="epitope name"):
        common.filter_assay_table(no_epitope_df)


def test_unused_get_mhc_assay_removed():
    assert not hasattr(columns, "get_mhc_assay")
    assert hasattr(columns, "get_assay_method")


def test_tcell_missing_file_error_names_the_download_function(monkeypatch, tmp_path):
    monkeypatch.setattr(
        iedb.tcell.cache, "local_path", lambda **kwargs: str(tmp_path / "missing.csv")
    )
    with pytest.raises(ValueError, match=r"pepdata\.iedb\.tcell\.download\(\)"):
        original_tcell_local_path(auto_download=False)
