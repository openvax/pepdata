# Copyright (c) 2014. Mount Sinai School of Medicine
#
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

import os

import pandas as pd

from ..common import split_classes, bad_amino_acids, cache
from ..reduced_alphabet import make_alphabet_transformer
from .common import group_peptides

TCELL_COMPACT_FILENAME = "tcell_compact.csv"
TCELL_COMPACT_URL = "http://www.iedb.org/doc/tcell_compact.zip"
TCELL_COMPACT_DECOMPRESS = True

def download_compact(force=False):
    return cache.fetch(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS,
        force=force)

def compact_local_path():
    path = cache.local_path(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS)
    if not os.exists(path):
        raise ValueError(
            ("Local file %s does not exist, call"
            " pepdata.iedb.tcell.download_compact()") % path)
    return path

def load_compact(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        human=True,
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        reduced_alphabet=None,  # 20 letter AA strings -> simpler alphabet
        nrows=None,
        verbose=False):
    """
    Load IEDB T-cell data without aggregating multiple entries for same epitope

    Parameters
    ----------
    mhc_class: {None, 1, 2}
        Restrict to MHC Class I or Class II (or None for neither)

    hla: regex pattern, optional
        Restrict results to specific HLA type used in assay

    exclude_hla: regex pattern, optional
        Exclude certain HLA types

    human: bool
        Restrict to human samples (default True)

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    assay_method string, optional
        Only collect results with assay methods containing the given string

    assay_group: string, optional
        Only collect results with assay groups containing the given string

    nrows: int, optional
        Don't load the full IEDB dataset but instead read only the first nrows

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet

    verbose: bool
        Print debug output
    """
    df = pd.read_csv(
            compact_local_path(),
            header=[0, 1],
            skipinitialspace=True,
            nrows=nrows,
            low_memory=False,
            error_bad_lines=False)

    # Sometimes the IEDB seems to put in an extra comma in the
    # header line, which creates an unnamed column of NaNs.
    # To deal with this, drop any columns which are all NaN
    df = df.dropna(axis=1, how='all')

    n = len(df)

    epitopes = df['Epitope']['Description'].str.upper()

    null_epitope_seq = epitopes.isnull()

    # if have rare or unknown amino acids, drop the sequence
    bad_epitope_seq = \
        epitopes.str.contains(bad_amino_acids, na=False).astype('bool')

    if verbose:
        print "Dropping %d null sequences" % null_epitope_seq.sum()
        print "Dropping %d bad sequences" % bad_epitope_seq.sum()

    mask = ~(bad_epitope_seq | null_epitope_seq)

    if human:
        mask &= df['MHC']['Allele Name'].str.startswith('HLA').astype('bool')

    # Match known alleles such as 'HLA-A*02:01',
    # broader groupings such as 'HLA-A2'
    # and unknown alleles of the MHC-1 listed either as
    #  'HLA-Class I,allele undetermined'
    #  or
    #  'Class I,allele undetermined'

    if mhc_class == 1:
        mask &= df['MHC']['MHC allele class'] == 'I'
    elif mhc_class == 2:
        mask &= df['MHC']['MHC allele class'] == 'II'

    if hla:
        mask &= df['MHC']['MHC Allele'].str.contains(hla, na=False)

    if exclude_hla:
        mask &= ~(df['MHC']['MHC Allele'].str.contains(exclude_hla, na=False))

    if assay_group:
        mask &= df['Assay']['Assay Group'].str.contains(assay_group)

    if assay_method:
        mask &= df['Assay']['Method/Technique'].str.contains(assay_method)

    if peptide_length:
        assert peptide_length > 0
        mask &= df['Epitope Linear Sequence'].str.len() == peptide_length

    df = df[mask]

    if verbose:
        print "Returning %d / %d entries after filtering" % (len(df), n)

    if reduced_alphabet:
        epitopes = df['Epitope']['Description']
        df['Epitope']['ReducedAlphabet'] = \
            epitopes.map(make_alphabet_transformer(reduced_alphabet))
    return df