from __future__ import annotations

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

import pandas as pd

from .common import (
    cache,
    delete_export,
    download_export,
    filter_assay_table,
    local_export_path,
    read_assay_csv,
)


MHC_URL = "https://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full_single_file.zip"
MHC_LOCAL_FILENAME = "mhc_ligand_full.csv"
MHC_DECOMPRESS = True


def download(force: bool = False) -> str:
    return download_export(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS,
        force=force,
        data_cache=cache,
    )


def local_path(auto_download: bool = True) -> str:
    return local_export_path(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS,
        auto_download=auto_download,
        download_instruction="pepdata.iedb.mhc.download()",
        data_cache=cache,
    )


def delete() -> None:
    """Remove the local assay export without downloading; missing files are ignored."""
    delete_export(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS,
        data_cache=cache,
    )


def load_dataframe(
        mhc_class: int | str | None = None,
        hla: str | None = None,
        exclude_hla: str | None = None,
        human_only: bool = False,
        peptide_length: int | None = None,
        assay_method: str | None = None,
        only_standard_amino_acids: bool = True,
        warn_bad_lines: bool = True,
        nrows: int | None = None) -> pd.DataFrame:
    """
    Load IEDB MHC data without aggregating multiple entries for the same epitope.

    Parameters
    ----------
    mhc_class : {1, 2, "I", "II", None}, optional
        Restrict to MHC Class I or Class II using the table's MHC class column.

    hla : str, optional
        Restrict results to specific MHC alleles matching this regex pattern.

    exclude_hla : str, optional
        Regex pattern to exclude matching MHC alleles.

    human_only : bool, default False
        Restrict to rows whose host organism starts with ``"Homo sapiens"``.

    peptide_length : int, optional
        Restrict epitopes to amino acid strings of given positive length.

    assay_method : str, optional
        Regex pattern limiting results to matching assay methods (null values
        do not match).

    only_standard_amino_acids : bool, default True
        Drop sequences which use non-standard amino acids, anything outside
        the core 20, such as X or U.

    warn_bad_lines : bool, default True
        The full MHC ligand dataset can contain lines with too many fields.
        When True, emit warnings from pandas; when False, skip them silently.

    nrows : int, optional
        Don't load the full IEDB dataset but instead read only the first nrows.
    """
    df = read_assay_csv(
        local_path(),
        nrows=nrows,
        warn_bad_lines=warn_bad_lines,
    )
    return filter_assay_table(
        df,
        mhc_class=mhc_class,
        hla=hla,
        exclude_hla=exclude_hla,
        human_only=human_only,
        peptide_length=peptide_length,
        assay_method=assay_method,
        only_standard_amino_acids=only_standard_amino_acids,
    )
