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

from ..reduced_alphabet import make_alphabet_transformer
from ..features import make_ngram_dataset_from_args
from ..common import split_classes, bad_amino_acids, cache
from .common import group_peptides


MHC_URL = "http://www.iedb.org/doc/mhc_ligand_full.zip"
MHC_LOCAL_FILENAME = "mhc_ligand_full.csv"
MHC_DECOMPRESS = True

def download(force=False):
    return cache.fetch(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS,
        force=force)

def local_path():
    path = cache.local_path(
        filename=MHC_LOCAL_FILENAME,
        url=MHC_URL,
        decompress=MHC_DECOMPRESS)
    if not os.exists(path):
        raise ValueError(
            ("MHC data file %s does not exist locally,"
             " call pepdata.mhc.download() to get a copy from IEDB") % path)
    return path

def delete():
    os.remove(local_path())

def _load_dataframe(
        filename,
        human=True,
        mhc_class=None,  # 1, 2, or None for both
        hla=None,  # regex pattern i.e. '(HLA-A2)|(HLA-A\*02)'
        exclude_hla=None,  # regex pattern i.e. '(HLA-A2)|(HLA-A\*02)'
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        nrows=None,
        reduced_alphabet=None,
        verbose=False):
    """
    Load a CSV of IEDB's full MHC data into a Pandas DataFrame and filter
    using the criteria given as function arguments
    """
    df = pd.read_csv(
            filename,
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


def load_full(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        human=True,
        peptide_length=None,
        assay_method=None,
        assay_group=None,
        reduced_alphabet=None,  # 20 letter AA strings -> simpler alphabet
        nrows=None,
        verbose=False,
        cache_download=True):
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

    human : bool
        Restrict to human samples (default True)

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    assay_method : string, optional
        Limit to assay methods which contain the given string

    assay_group : string, optional
        Limit to assay groups which contain the given string

    reduced_alphabet : dictionary, optional
        Remap amino acid letters to some other alphabet

    nrows : int, optional
        Don't load the full IEDB dataset but instead read only the first nrows

    verbose : bool
        Print debug output
    """
    return _load_dataframe(
                local_path(),
                mhc_class=mhc_class,
                hla=hla,
                exclude_hla=exclude_hla,
                human=human,
                peptide_length=peptide_length,
                assay_method=assay_method,
                assay_group=assay_group,
                reduced_alphabet=reduced_alphabet,
                nrows=nrows,
                verbose=verbose)

def _group_mhc_peptides(
        df,
        unique_sequences=True,
        min_count=0,
        group_by_allele=False,
        verbose=False):
    """
    Given a dataframe of epitopes and qualitative measures,
    group the epitope strings (optionally also grouping by allele),
    and associate each group with its percentage of Positive
    Qualitative Measure results.
    """
    epitopes = df['Epitope']['Description']
    measure = df['Assay']['Qualitative Measure']
    pos_mask = measure.str.startswith('Positive').astype('bool')
    mhc = df['MHC']['Allele Name']
    return group_peptides(
        epitopes,
        mhc,
        pos_mask,
        group_by_allele=group_by_allele,
        min_count=min_count)

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
        group_by_allele=False,
        min_count=0,
        verbose=False):
    """
    Load the MHC binding results from IEDB, collect into a dataframe mapping
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
        Don't load the full IEDB dataset but instead read only the first nrows

    group_by_allele:
        Don't combine epitopes across multiple HLA types

    min_count: int, optional
        Exclude epitopes which appear fewer times than min_count

    verbose: bool
        Print debug output
    """
    df = load_full(
        mhc_class=mhc_class,
        hla=hla,
        exclude_hla=exclude_hla,
        human=human,
        peptide_length=peptide_length,
        assay_method=assay_method,
        assay_group=assay_group,
        reduced_alphabet=reduced_alphabet,
        nrows=nrows,
        verbose=verbose)

    return _group_mhc_peptides(
            df,
            group_by_allele=group_by_allele,
            min_count=min_count,
            verbose=verbose)
