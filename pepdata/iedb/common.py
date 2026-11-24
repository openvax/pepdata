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

from __future__ import annotations

import logging
import os
from pathlib import Path

import datacache
import pandas as pd

from .columns import (
    get_assay_method,
    get_epitope_name,
    get_host_name,
    get_mhc_allele,
    get_mhc_class,
)

cache = datacache.Cache("pepdata")

bad_amino_acids = "U|X|J|B|Z"


def download_export(
        filename: str,
        url: str,
        decompress: bool = True,
        force: bool = False,
        data_cache: datacache.Cache | None = None) -> str:
    """Download an IEDB export into the local pepdata cache and return its path."""
    active_cache = cache if data_cache is None else data_cache
    return active_cache.fetch(
        filename=filename,
        url=url,
        decompress=decompress,
        force=force,
    )


def local_export_path(
        filename: str,
        url: str,
        decompress: bool = True,
        auto_download: bool = True,
        download_instruction: str | None = None,
        data_cache: datacache.Cache | None = None) -> str:
    """Return the local path to a cached IEDB export, optionally downloading it."""
    active_cache = cache if data_cache is None else data_cache
    path = active_cache.local_path(
        filename=filename,
        url=url,
        decompress=decompress,
    )
    if not os.path.exists(path):
        if auto_download:
            return download_export(
                filename=filename,
                url=url,
                decompress=decompress,
                data_cache=active_cache,
            )
        hint = download_instruction or "pepdata.iedb.common.download_export()"
        raise ValueError(
            f"IEDB data file {path} does not exist locally, "
            f"call {hint} to get a copy from IEDB"
        )
    return path


def delete_export(
        filename: str,
        url: str,
        decompress: bool = True,
        data_cache: datacache.Cache | None = None) -> None:
    """Remove a local cached IEDB export without downloading; missing files are ignored."""
    active_cache = cache if data_cache is None else data_cache
    path = active_cache.local_path(
        filename=filename,
        url=url,
        decompress=decompress,
    )
    Path(path).unlink(missing_ok=True)


def normalize_mhc_class(mhc_class: int | str) -> str:
    """Normalize an MHC class selector (1, 2, 'I', or 'II') to 'I' or 'II'."""
    if isinstance(mhc_class, bool):
        raise ValueError(
            f"Invalid MHC class: {mhc_class!r}. Expected 1, 2, 'I', or 'II'."
        )
    if mhc_class == 1 or mhc_class == "I":
        return "I"
    if mhc_class == 2 or mhc_class == "II":
        return "II"
    raise ValueError(
        f"Invalid MHC class: {mhc_class!r}. Expected 1, 2, 'I', or 'II'."
    )


def read_assay_csv(
        path: str | os.PathLike[str],
        nrows: int | None = None,
        warn_bad_lines: bool = True) -> pd.DataFrame:
    """Read an IEDB two-header-row CSV export and drop unnamed all-NaN columns."""
    df = pd.read_csv(
        path,
        header=[0, 1],
        skipinitialspace=True,
        nrows=nrows,
        low_memory=False,
        on_bad_lines="warn" if warn_bad_lines else "skip",
        encoding="latin-1",
    )
    # Sometimes IEDB exports include a trailing comma in the header line,
    # creating an unnamed column of NaNs. Drop unnamed all-NaN columns without
    # dropping real named columns that happen to be empty in a slice.
    unnamed_all_na_columns = [
        col for col in df.columns
        if isinstance(col, tuple)
        and any(str(part).startswith("Unnamed:") for part in col)
        and df[col].isna().all()
    ]
    if unnamed_all_na_columns:
        df = df.drop(columns=unnamed_all_na_columns)
    return df


def filter_assay_table(
        df: pd.DataFrame,
        mhc_class: int | str | None = None,
        hla: str | None = None,
        exclude_hla: str | None = None,
        human_only: bool = False,
        peptide_length: int | None = None,
        assay_method: str | None = None,
        only_standard_amino_acids: bool = True) -> pd.DataFrame:
    """
    Filter an IEDB assay DataFrame (MHC-binding or T-cell) using consistent semantics.

    Columns are resolved via :mod:`pepdata.iedb.columns` so both legacy and
    current IEDB two-level column headers are supported.

    Parameters
    ----------
    df : pandas.DataFrame
        Raw IEDB assay table with two-level column headers.

    mhc_class : {1, 2, "I", "II", None}, optional
        Restrict to MHC Class I or Class II using the table's own MHC class
        column (without downloading the IEDB allele XML catalogue). Raises
        ``ValueError`` for any other non-``None`` value.

    hla : str, optional
        Regex pattern restricting results to matching MHC allele names.

    exclude_hla : str, optional
        Regex pattern excluding matching MHC allele names.

    human_only : bool, default False
        Restrict to rows whose host organism starts with ``"Homo sapiens"``.

    peptide_length : int, optional
        Restrict epitopes to amino acid strings of this exact positive length.

    assay_method : str, optional
        Regex pattern restricting results to matching assay methods. Rows with
        missing assay methods do not match.

    only_standard_amino_acids : bool, default True
        Drop sequences containing ambiguous or non-standard residues
        (``U``, ``X``, ``J``, ``B``, ``Z``).

    Returns
    -------
    pandas.DataFrame
        Filtered copy of ``df`` with upper-cased epitope sequences.
    """
    normalized_class = (
        normalize_mhc_class(mhc_class) if mhc_class is not None else None
    )

    if peptide_length is not None:
        if (
            isinstance(peptide_length, bool)
            or not isinstance(peptide_length, int)
            or peptide_length <= 0
        ):
            raise ValueError(
                f"peptide_length must be a positive integer, got {peptide_length!r}"
            )

    epitopes = get_epitope_name(df)
    if epitopes is None:
        raise ValueError(
            "Could not find epitope name column in IEDB assay data. "
            f"Available columns: {list(df.columns)}"
        )

    n = len(df)
    normalized_epitopes = epitopes.astype("string").str.upper()

    null_epitope_seq = epitopes.isnull()
    n_null = int(null_epitope_seq.sum())
    if n_null > 0:
        logging.info("Dropping %d null sequences", n_null)

    mask = (~null_epitope_seq).astype(bool)

    if only_standard_amino_acids:
        bad_epitope_seq = (
            normalized_epitopes.str.contains(bad_amino_acids, na=False).astype(bool)
        )
        n_bad = int(bad_epitope_seq.sum())
        if n_bad > 0:
            logging.info("Dropping %d bad sequences", n_bad)
        mask &= ~bad_epitope_seq

    if human_only:
        organism = get_host_name(df)
        if organism is None:
            raise ValueError(
                "Could not find host organism column in IEDB assay data. "
                f"Available columns: {list(df.columns)}"
            )
        mask &= (
            organism.astype("string")
            .str.startswith("Homo sapiens", na=False)
            .astype(bool)
        )

    if normalized_class is not None:
        mhc_class_series = get_mhc_class(df)
        if mhc_class_series is None:
            raise ValueError(
                "Could not find MHC class column in IEDB assay data. "
                f"Available columns: {list(df.columns)}"
            )
        mask &= (
            (mhc_class_series.astype("string").str.strip() == normalized_class)
            .fillna(False)
            .astype(bool)
        )

    if hla or exclude_hla:
        mhc_alleles = get_mhc_allele(df)
        if mhc_alleles is None:
            raise ValueError(
                "Could not find MHC allele column in IEDB assay data. "
                f"Available columns: {list(df.columns)}"
            )
        mhc_allele_strings = mhc_alleles.astype("string")
        if hla:
            mask &= mhc_allele_strings.str.contains(hla, na=False).astype(bool)
        if exclude_hla:
            mask &= ~(
                mhc_allele_strings.str.contains(exclude_hla, na=False).astype(bool)
            )

    if assay_method:
        assay_method_series = get_assay_method(df)
        if assay_method_series is None:
            raise ValueError(
                "Could not find assay method column in IEDB assay data. "
                f"Available columns: {list(df.columns)}"
            )
        mask &= (
            assay_method_series.astype("string")
            .str.contains(assay_method, na=False)
            .astype(bool)
        )

    if peptide_length is not None:
        mask &= (
            (normalized_epitopes.str.len() == peptide_length)
            .fillna(False)
            .astype(bool)
        )

    result = df.loc[mask].copy()
    if len(result) > 0:
        result[epitopes.name] = result[epitopes.name].str.upper()

    logging.info("Returning %d / %d entries after filtering", len(result), n)
    return result
