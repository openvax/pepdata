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

from __future__ import print_function, division, absolute_import
import logging
import os

import pandas as pd

from .memoize import  memoize
from .common import bad_amino_acids, cache


MHC_URL = "http://www.iedb.org/downloader.php?file_name=doc/mhc_ligand_full.zip"
MHC_LOCAL_FILENAME = "mhc_ligand_full.csv"
MHC_DECOMPRESS = True

def download(force=False):
    return cache.fetch(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS,
        force=force)

def local_path(auto_download=True):
    path = cache.local_path(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS)
    if not os.path.exists(path):
        if auto_download:
            return download()
        raise ValueError(
            ("MHC data file %s does not exist locally,"
             " call pepdata.mhc.download() to get a copy from IEDB") % path)
    return path

def delete():
    os.remove(local_path())

@memoize
def load_dataframe(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        human_only=False,
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        only_standard_amino_acids=True,
        reduced_alphabet=None,  # 20 letter AA strings -> simpler alphabet
        warn_bad_lines=True,
        nrows=None):
    """
    Load IEDB MHC data without aggregating multiple entries for the same epitope

    Parameters
    ----------
    mhc_class : {None, 1, 2}
        Restrict to MHC Class I or Class II (or None for neither)

    hla : regex pattern, optional
        Restrict results to specific HLA type used in assay

    exclude_hla : regex pattern, optional
        Exclude certain HLA types

    human_only : bool
        Restrict to human samples (default False)

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    assay_method : string, optional
        Limit to assay methods which contain the given string

    assay_group : string, optional
        Limit to assay groups which contain the given string

    only_standard_amino_acids : bool, optional
        Drop sequences which use non-standard amino acids, anything outside
        the core 20, such as X or U (default = True)

    reduced_alphabet : dictionary, optional
        Remap amino acid letters to some other alphabet

    warn_bad_lines : bool, optional
        The full MHC ligand dataset seems to contain several dozen lines with
        too many fields. This currently results in a lot of warning messages
        from Pandas, which you can turn off with this option (default = True)

    nrows : int, optional
        Don't load the full IEDB dataset but instead read only the first nrows
    """
    df = pd.read_csv(
            local_path(),
            header=[0, 1],
            skipinitialspace=True,
            nrows=nrows,
            low_memory=False,
            error_bad_lines=False,
            encoding="latin-1",
            warn_bad_lines=warn_bad_lines)

    # Sometimes the IEDB seems to put in an extra comma in the
    # header line, which creates an unnamed column of NaNs.
    # To deal with this, drop any columns which are all NaN
    df = df.dropna(axis=1, how="all")

    n = len(df)

    epitope_column_key = ("Epitope", "Description")
    mhc_allele_column_key = ("MHC", "Allele Name")

    epitopes = df[epitope_column_key] = df[epitope_column_key].str.upper()

    null_epitope_seq = epitopes.isnull()
    n_null = null_epitope_seq.sum()
    if n_null > 0:
        logging.info("Dropping %d null sequences", n_null)

    mask = ~null_epitope_seq

    if only_standard_amino_acids:
        # if have rare or unknown amino acids, drop the sequence
        bad_epitope_seq = \
            epitopes.str.contains(bad_amino_acids, na=False).astype("bool")
        n_bad = bad_epitope_seq.sum()
        if n_bad > 0:
            logging.info("Dropping %d bad sequences", n_bad)

        mask &= ~bad_epitope_seq

    if human_only:
        mask &= df[mhc_allele_column_key].str.startswith("HLA").astype("bool")

    if mhc_class == 1:
        mask &= df["MHC"]["MHC allele class"] == "I"
    elif mhc_class == 2:
        mask &= df["MHC"]["MHC allele class"] == "II"

    if hla:
        mask &= df[mhc_allele_column_key].str.contains(hla, na=False)

    if exclude_hla:
        mask &= ~(df[mhc_allele_column_key].str.contains(exclude_hla, na=False))

    if assay_group:
        mask &= df["Assay"]["Assay Group"].str.contains(assay_group)

    if assay_method:
        mask &= df["Assay"]["Method/Technique"].str.contains(assay_method)

    if peptide_length:
        assert peptide_length > 0
        mask &= df[epitope_column_key].str.len() == peptide_length

    df = df[mask].copy()

    logging.info("Returning %d / %d entries after filtering", len(df), n)

    return df
