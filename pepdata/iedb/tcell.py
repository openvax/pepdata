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

import numpy as np
import pandas as pd


from .alleles import load_alleles_dict
from .memoize import  memoize
from .common import  bad_amino_acids, cache
from .columns import (
    get_assay_method,
    get_assay_num_tested,
    get_assay_response_measured,
    get_assay_units,
    get_host_name,
    get_mhc_allele,
    get_mhc_assay,
    get_mhc_class,
    get_epitope_source_organism,
    get_epitope_type,
    get_epitope_name,
    
)

TCELL_COMPACT_FILENAME = "tcell_full.csv"
TCELL_COMPACT_URL = "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
TCELL_COMPACT_DECOMPRESS = True

def download(force=False):
    return cache.fetch(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS,
        force=force)

def local_path(auto_download=True):
    path = cache.local_path(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS)
    if not os.path.exists(path):
        if auto_download:
            return download()
        raise ValueError(
            ("Local file %s does not exist, call"
            " pepdata.iedb.tcell.download()") % path)
    return path

def delete():
    os.remove(local_path())

@memoize
def load_dataframe(
        mhc_class : str | None = None,  # 1, 2, or None for neither
        mhc_pattern : str | None  = None,
        exclude_mhc : str | None  = None,
        human_only : bool =False,
        peptide_length : int | None = None,
        assay_method : str | None = None,
        only_standard_amino_acids : bool = True,
        reduced_alphabet : dict | None = None,  # 20 letter AA strings -> simpler alphabet
        nrows : int | None = None):
    """
    Load IEDB T-cell data without aggregating multiple entries for same epitope

    Parameters
    ----------
    mhc_class: {None, 1, 2}
        Restrict to MHC Class I or Class II (or None for neither)

    mhc_pattern: regex pattern, optional
        Restrict results to specific MHC used in assay

    exclude_mhc: regex pattern, optional
        Exclude certain MHC allele patterns 

    human_only: bool
        Restrict to human samples (default False)

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    assay_method string, optional
        Only collect results with assay methods containing the given string

    only_standard_amino_acids : bool, optional
        Drop sequences which use non-standard amino acids, anything outside
        the core 20, such as X or U (default = True)

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet

    nrows: int, optional
        Don't load the full IEDB dataset but instead read only the first nrows
    """
    path = local_path()
    df = pd.read_csv(
            path,
            header=[0, 1],
            skipinitialspace=True,
            nrows=nrows,
            low_memory=False,
            on_bad_lines='warn',
            encoding="latin-1")
    
    mhc = get_mhc_allele(df)
    mhc_class = get_mhc_class(df)
    epitopes = get_epitope_name(df)
    organism = get_host_name(df)
    assay_method = get_assay_method(df)


    # Sometimes the IEDB seems to put in an extra comma in the
    # header line, which creates an unnamed column of NaNs.
    # To deal with this, drop any columns which are all NaN
    df = df.dropna(axis=1, how="all")

    n = len(df)
    
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
        mask &= organism.str.startswith('Homo sapiens', na=False).astype('bool')


    if mhc_class is not None:
        # since MHC classes can be specified as either strings ("I") or integers
        # standard them to be strings
        if mhc_class == 1:
            mhc_class = "I"
        elif mhc_class == 2:
            mhc_class = "II"
        if mhc_class not in {"I", "II"}:
            raise ValueError("Invalid MHC class: %s" % mhc_class)
        allele_dict = load_alleles_dict()
        mhc_class_mask = [False] * len(df)
        for i, allele_name in enumerate(mhc):
            allele_object = allele_dict.get(allele_name)
            if allele_object and allele_object.mhc_class == mhc_class:
                mhc_class_mask[i] = True
        mask &= np.array(mhc_class_mask)

    # Match known alleles such as "HLA-A*02:01",
    # broader groupings such as "HLA-A2"
    # and unknown alleles of the MHC-1 listed either as
    #  "HLA-Class I,allele undetermined"
    #  or
    #  "Class I,allele undetermined"
]

    if hla:
        mask &= df[mhc_allele_column_key].str.contains(hla, na=False)

    if exclude_hla:
        mask &= ~(df[mhc_allele_column_key].str.contains(exclude_hla, na=False))

    if assay_group:
        mask &= df[assay_group_column_key].str.contains(assay_group)

    if assay_method:
        mask &= df[assay_method_column_key].str.contains(assay_method)

    if peptide_length:
        assert peptide_length > 0
        mask &= df[epitope_column_key].str.len() == peptide_length

    df = df[mask]

    logging.info("Returning %d / %d entries after filtering", len(df), n)
    return df
