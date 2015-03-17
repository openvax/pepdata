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

import logging
import os

import pandas as pd

from ..common import split_classes, bad_amino_acids, cache, memoize
from ..features import make_ngram_dataset_from_args
from ..reduced_alphabet import make_alphabet_transformer
from .common import group_peptides

TCELL_COMPACT_FILENAME = "tcell_compact.csv"
TCELL_COMPACT_URL = "http://www.iedb.org/doc/tcell_compact.zip"
TCELL_COMPACT_DECOMPRESS = True

def download(force=False):
    return cache.fetch(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS,
        force=force)

def local_path():
    path = cache.local_path(
        filename=TCELL_COMPACT_FILENAME,
        url=TCELL_COMPACT_URL,
        decompress=TCELL_COMPACT_DECOMPRESS)
    if not os.path.exists(path):
        raise ValueError(
            ("Local file %s does not exist, call"
            " pepdata.iedb.tcell.download()") % path)
    return path

def delete():
    os.remove(local_path())

@memoize
def load_dataframe(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        human=True,
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        reduced_alphabet=None,  # 20 letter AA strings -> simpler alphabet
        nrows=None):
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
        Don"t load the full IEDB dataset but instead read only the first nrows

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet
    """
    df = pd.read_csv(
            local_path(),
            skipinitialspace=True,
            nrows=nrows,
            low_memory=False,
            error_bad_lines=True)

    # Sometimes the IEDB seems to put in an extra comma in the
    # header line, which creates an unnamed column of NaNs.
    # To deal with this, drop any columns which are all NaN
    df = df.dropna(axis=1, how="all")

    n = len(df)

    epitopes = df["Epitope Linear Sequence"].str.upper()

    null_epitope_seq = epitopes.isnull()
    n_null = null_epitope_seq.sum()

    if n_null > 0:
        logging.info("Dropping %d null sequences", n_null)

    # if have rare or unknown amino acids, drop the sequence
    bad_epitope_seq = \
        epitopes.str.contains(bad_amino_acids, na=False).astype("bool")
    n_bad = bad_epitope_seq.sum()
    if n_bad > 0:
        logging.info("Dropping %d bad sequences", n_bad)

    mask = ~(bad_epitope_seq | null_epitope_seq)

    if human:
        organism = df['Host Organism Name']
        mask &= organism.str.startswith('Homo sapiens', na=False).astype('bool')

    # Match known alleles such as "HLA-A*02:01",
    # broader groupings such as "HLA-A2"
    # and unknown alleles of the MHC-1 listed either as
    #  "HLA-Class I,allele undetermined"
    #  or
    #  "Class I,allele undetermined"
    mhc = df['MHC Allele Name']

    if mhc_class == 1:
        mhc1_pattern = 'Class I,|HLA-[A-C](?:[0-9]|\*)'
        mask &= mhc.str.contains(mhc1_pattern, na=False).astype('bool')
    elif mhc_class == 2:
        mhc2_pattern = "Class II,|HLA-D(?:P|M|O|Q|R)"
        mask &= mhc.str.contains(mhc2_pattern, na=False).astype('bool')

    if mhc_class == 1:
        mask &= df["MHC"]["MHC allele class"] == "I"
    elif mhc_class == 2:
        mask &= df["MHC"]["MHC allele class"] == "II"

    if hla:
        mask &= df["MHC"]["MHC Allele"].str.contains(hla, na=False)

    if exclude_hla:
        mask &= ~(df["MHC"]["MHC Allele"].str.contains(exclude_hla, na=False))

    if assay_group:
        mask &= df["Assay"]["Assay Group"].str.contains(assay_group)

    if assay_method:
        mask &= df["Assay"]["Method/Technique"].str.contains(assay_method)

    if peptide_length:
        assert peptide_length > 0
        mask &= df["Epitope Linear Sequence"].str.len() == peptide_length

    df = df[mask]

    logging.info("Returning %d / %d entries after filtering", len(df), n)

    if reduced_alphabet:
        epitopes = df["Epitope"]["Description"]
        df["Epitope"]["ReducedAlphabet"] = \
            epitopes.map(make_alphabet_transformer(reduced_alphabet))
    return df

@memoize
def load_groups(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        human=True,
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        reduced_alphabet=None,  # 20 letter AA strings -> simpler alphabet
        nrows=None,
        min_count=0,
        group_by_allele=False,
        verbose=False):
    """
    Load the T-cell response data from IEDB, collect into a dataframe mapping
    epitopes to percentage positive results.

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

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet

    nrows: int, optional
        Don"t load the full IEDB dataset but instead read only the first nrows

    group_by_allele:
        Don"t combine epitopes across multiple HLA types

    min_count: int, optional
        Exclude epitopes which appear fewer times than min_count

    verbose: bool
        Print debug output
    """

    df = load_dataframe(
        mhc_class=mhc_class,
        hla=hla,
        human=human,
        exclude_hla=exclude_hla,
        peptide_length=peptide_length,
        assay_method=assay_method,
        assay_group=assay_group,
        reduced_alphabet=reduced_alphabet,
        nrows=nrows,
        verbose=verbose)

    return group_peptides(
            df,
            group_by_allele=group_by_allele,
            min_count=min_count,
            verbose=verbose)

@memoize
def load_classes(*args, **kwargs):
    """
    Split the T-cell assay results into positive and negative sets.

    Parameters
    ----------
    noisy_labels : {"majority", "negative", "positive"}
        Which class do we assign an epitope with contradictory labels?

    *args, **kwargs : same as "load_tcell"
    """
    noisy_labels = kwargs.pop("noisy_labels", "majority")
    verbose = kwargs.get("verbose")
    tcell_values = load_groups(*args, **kwargs)
    return split_classes(
        tcell_values.value,
        noisy_labels=noisy_labels,
        verbose=verbose)

@memoize
def load_ngrams(*args, **kwargs):
    """
    Construct n-gram input features X and output labels Y for T-cell responses

    Parameters:
    ----------
    max_ngram : int
        Order of n-grams to consider when constructing X.
        For example, when max_ngram = 1, the vector space is the individual
        frequencies of letters in the amino acid strings.

    normalize_row : bool, optional
        If True (default), then return frequencies, else raw counts.

    subsample_bigger_class: bool, optional
        When the number of samples in both classes are unbalanced,
        randomly drop some samples from the larger class (default = False).

    return_transformer: bool
        Return `X, Y, f` instead of just `X, Y`,
        where f can be used to transform new amino acid strings into
        the same space as the training data.


    *args, **kwargs : same as `load_tcell_classes`
    """
    kwargs["training_already_reduced"] = True
    return make_ngram_dataset_from_args(load_classes, *args, **kwargs)
