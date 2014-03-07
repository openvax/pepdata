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


from os.path import join

import numpy as np
import pandas as pd

from base import DATA_DIR
from features import make_ngram_dataset

TCELL_CSV = join(DATA_DIR, 'tcell_compact.csv')
MHC_CSV = join(DATA_DIR, 'elution_compact.csv')

def load_dataframe(
        filename = TCELL_CSV,
        human = True,
        mhc_class = None, # 1, 2, or None for both
        hla_type = None, # regex pattern i.e. '(HLA-A2)|(HLA-A\*02)'
        exclude_hla_type = None, # regex pattern i.e. '(HLA-A2)|(HLA-A\*02)'
        peptide_length = None,
        assay_group=None,
        nrows = None,
        verbose= True):
    """
    Load an IEDB csv into a pandas dataframe and filter using the
    criteria given as function arguments
    """
    df = pd.read_csv(filename, skipinitialspace=True, nrows = nrows)
    mhc = df['MHC Allele Name']

    #
    # Match known alleles such as 'HLA-A*02:01',
    # broader groupings such as 'HLA-A2'
    # and unknown alleles of the MHC-1 listed either as
    #  'HLA-Class I,allele undetermined'
    #  or
    #  'Class I,allele undetermined'
    mhc1_pattern = 'Class I,|HLA-[A-C]([0-9]|\*)'
    mhc1_mask = mhc.str.contains(mhc1_pattern, na=False).astype('bool')

    mhc2_pattern = "Class II,|HLA-D(P|M|O|Q|R)"
    mhc2_mask = mhc.str.contains(mhc2_pattern, na=False).astype('bool')

    # just in case any of the results were from mice or other species,
    # restrict to humans
    organism = df['Host Organism Name']
    human_organism_mask = \
        organism.str.startswith('Homo sapiens', na=False).astype('bool')
    human_hla_mask = mhc.str.startswith("HLA", na=False).astype('bool')
    human_mask = human_organism_mask | human_hla_mask

    if verbose:
        print "Human entries", human_mask.sum()
        print "Class I MHC Entries", mhc1_mask.sum()
        print "Class II MHC Entries", mhc2_mask.sum()
        print "Human Class I MHCs", (human_mask & mhc1_mask).sum()
        print "Human Class II MHCs", (human_mask & mhc2_mask).sum()

    epitopes = df['Epitope Linear Sequence'].str.upper()
    df['Epitope Linear Sequence'] = epitopes
    null_epitope_seq = epitopes.isnull()

    bad_amino_acids = 'U|X|J|B|Z'
    # if have rare or unknown amino acids, drop the sequence
    bad_epitope_seq = \
        epitopes.str.contains(bad_amino_acids, na=False).astype('bool')

    if verbose:
        print "Dropping %d null sequences" % null_epitope_seq.sum()
        print "Dropping %d bad sequences" % bad_epitope_seq.sum()

    has_epitope_seq = ~(bad_epitope_seq | null_epitope_seq)

    mask = has_epitope_seq

    if human:
        mask &= human_mask
    if mhc_class == 1:
        mask &= mhc1_mask
    if mhc_class == 2:
        mask &= mhc2_mask

    if assay_group:
        mask &= df['Assay Group'] == assay_group

    if hla_type:
        mask &= mhc.str.contains(hla_type, na=False)

    if exclude_hla_type:
        mask &= ~mhc.str.contains(exclude_hla_type, na=False)

    if peptide_length:
        assert peptide_length > 0
        mask &=  epitopes.str.len() == peptide_length

    if verbose:
        print "Filtered sequences epitope sequences", mask.sum()

    return df[mask]

def group_epitopes(
        df,
        unique_sequences = True,
        min_count = 0,
        group_by_allele = False,
        reduced_alphabet = None,
        verbose = True):
    """
    Given a dataframe of epitopes and qualitative measures,
    group the epitope strings (optionally also grouping by allele),
    and associate each group with its percentage of Positive
    Qualitative Measure results.
    """
    epitopes = df['Epitope Linear Sequence'].str.upper()
    if reduced_alphabet:
        def transform(s):
            return ''.join([chr(48 + reduced_alphabet[char]) for char in s])
        epitopes = epitopes.map(transform)
        df['Epitope Linear Sequence'] = epitopes


    measure = df['Qualitative Measure']
    mhc = df['MHC Allele Name']
    pos_mask = measure.str.startswith('Positive').astype('bool')

    if group_by_allele:
        groups = pos_mask.groupby([epitopes, mhc])
    else:
        groups = pos_mask.groupby(epitopes)

    values = groups.mean()

    if min_count:
        counts = groups.count()
        values = values[counts >= min_count]

    return values


def split_classes(
        values,
        noisy_labels = 'majority',
        unique_sequences = True,
        verbose = True):
    """
    - majority = epitope is Positive if majority of its results are Positive
    - negative = epitope is Negative if any result is Negative
    - positive = epitope is Positive if any result is Positive
    - drop = remove any epitopes with contradictory results
    - keep = leave contradictory results in both positive and negative sets
    """
    if noisy_labels == 'majority':
        pos_mask = values >= 0.5
    elif noisy_labels == 'positive':
        pos_mask = values > 0
    else:
        pos_mask = values == 1.0

    neg_mask = ~pos_mask
    pos = pos_mask.index[pos_mask]
    neg = neg_mask.index[neg_mask]

    pos_set = set(pos)
    neg_set = set(neg)

    if verbose:
        print "# positive sequences", len(pos)
        print "# negative sequences", len(neg)

    noisy_set = pos_set.intersection(neg_set)

    if verbose:
        print "# unique positive", len(pos_set)
        print "# unique negative", len(neg_set)
        print "# overlap %d (%0.4f)" % (len(noisy_set), \
          float(len(noisy_set)) / len(pos_set))

    if noisy_labels != 'majority':
        if (noisy_labels == 'drop') or (noisy_labels == 'negative'):
            pos_set = pos_set.difference(noisy_set)
        if (noisy_labels == 'drop') or (noisy_labels == 'positive'):
            neg_set = neg_set.difference(noisy_set)
    if unique_sequences:
        return pos_set, neg_set
    else:
        pos = [epitope for epitope in pos if epitope not in pos_set]
        neg = [epitope for epitope in neg if epitope not in neg_set]
        return pos, neg


def load_tcell(
        mhc_class = None, # 1, 2, or None for neither
        hla_type = None,
        exclude_hla_type = None,
        peptide_length = None,
        min_count = 0,
        assay_group=None,
        nrows = None,
        group_by_allele = False,
        reduced_alphabet = None, # 20 letter AA strings -> simpler alphabet
        verbose= True):
    """
    Load the T-cell response data from IEDB, collect into a dataframe mapping
    epitopes to percentage positive results.

    Parameters
    ----------
    mhc_class: {None, 1, 2}
        Restrict to MHC Class I or Class II (or None for neither)

    hla_type: regex pattern, optional
        Restrict results to specific HLA type used in assay

    exclude_hla_type: regex pattern, optional
        Exclude certain HLA types

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    min_count: int, optional
        Exclude epitopes which appear fewer times than min_count

    assay_group: string, optional
        Only collect results from assays of the given type

    nrows: int, optional
        Don't load the full IEDB dataset but instead read only the first nrows

    group_by_allele:
        Don't combine epitopes across multiple HLA types

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet

    verbose: bool
        Print debug output
    """


    df = load_dataframe(
            TCELL_CSV,
            assay_group = assay_group,
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            nrows = nrows,
            verbose = verbose)

    return group_epitopes(
            df,
            group_by_allele = group_by_allele,
            min_count = min_count,
            reduced_alphabet = reduced_alphabet,
            verbose = verbose)

def load_tcell_classes(*args, **kwargs):
    """
    Split the T-cell assay results into positive and negative sets.

    Parameters
    ----------
    noisy_labels : {'majority', 'negative', 'positive'}
        Which class do we assign an epitope with contradictory labels?

    *args, **kwargs : same as 'load_tcell'
    """
    noisy_labels = kwargs.pop('noisy_labels', None)
    verbose = kwargs.get('verbose')
    tcell_values = load_tcell(*args, **kwargs)
    return split_classes(
        tcell_values,
        noisy_labels = noisy_labels,
        verbose = verbose)

def load_tcell_ngrams(*args, **kwargs):
    """
    Construct n-gram input features X and output labels Y for T-cell responses

    Parameters:
    ----------
    ngram : int
        Order of n-grams to consider when constructing X.
        For example, when ngram = 1, the vector space is the individual
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

    ngram = kwargs.pop('ngram', 1)
    normalize_row = kwargs.pop('normalize_row', True)
    subsample_bigger_class = kwargs.pop('subsample_bigger_class', False)
    return_transformer = kwargs.pop('return_transformer', False)
    verbose = kwargs.get('verbose')

    pos, neg = load_tcell_classes(*args, **kwargs)
    return make_ngram_dataset(
        pos,
        neg,
        max_ngram = ngram,
        normalize_row = normalize_row,
        subsample_bigger_class = subsample_bigger_class,
        return_transformer = return_transformer)

def load_mhc(
        mhc_class = None, # 1, 2, or None for neither
        hla_type = None,
        exclude_hla_type = None,
        peptide_length = None,
        min_count = 0,
        assay_group=None,
        nrows = None,
        group_by_allele = False,
        reduced_alphabet = None, # 20 letter AA strings -> simpler alphabet
        verbose= True):
    """
    Load the MHC binding results from IEDB, collect into a dataframe mapping
    epitopes to percentage positive results.

    Parameters
    ----------
    mhc_class: {None, 1, 2}
        Restrict to MHC Class I or Class II (or None for neither)

    hla_type: regex pattern, optional
        Restrict results to specific HLA type used in assay

    exclude_hla_type: regex pattern, optional
        Exclude certain HLA types

    peptide_length: int, optional
        Restrict epitopes to amino acid strings of given length

    min_count: int, optional
        Exclude epitopes which appear fewer times than min_count

    assay_group: string, optional
        Only collect results from assays of the given type

    nrows: int, optional
        Don't load the full IEDB dataset but instead read only the first nrows

    group_by_allele:
        Don't combine epitopes across multiple HLA types

    reduced_alphabet: dictionary, optional
        Remap amino acid letters to some other alphabet

    verbose: bool
        Print debug output
    """
    df = load_dataframe(
            MHC_CSV,
            assay_group = assay_group,
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            nrows = nrows,
            verbose = verbose)

    return group_epitopes(
            df,
            group_by_allele = group_by_allele,
            min_count = min_count,
            reduced_alphabet = reduced_alphabet,
            verbose = verbose)

def load_mhc_classes(*args, **kwargs):
    """
    Split the MHC binding assay results into positive and negative sets.

    Parameters
    ----------
    noisy_labels : 'majority' | 'negative' | 'positive'
        Which class do we assign an epitope with contradictory labels?

    *args, **kwargs : same as 'load_tcell'
    """
    noisy_labels = kwargs.pop('noisy_labels', None)
    verbose = kwargs.get('verbose')
    mhc_values = load_mhc(*args, **kwargs)
    return split_classes(
        mhc_values,
        noisy_labels = noisy_labels,
        verbose = verbose)

def load_mhc_ngrams(*args, **kwargs):
    """
    Construct n-gram input features X and output labels Y for MHC binding

    Parameters:
    ----------
    ngram : int
        Order of n-grams to consider when constructing X.
        For example, when ngram = 1, the vector space is the individual
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

    ngram = kwargs.pop('ngram', 1)
    normalize_row = kwargs.pop('normalize_row', True)
    subsample_bigger_class = kwargs.pop('subsample_bigger_class', False)
    return_transformer = kwargs.pop('return_transformer', False)
    verbose = kwargs.get('verbose')

    pos, neg = load_mhc_classes(*args, **kwargs)
    return make_ngram_dataset(
        pos,
        neg,
        max_ngram = ngram,
        normalize_row = normalize_row,
        subsample_bigger_class = subsample_bigger_class,
        return_transformer = return_transformer)


def load_tcell_vs_mhc(
        mhc_class = None, # 1, 2, or None for neither
        hla_type = None,
        exclude_hla_type = None,
        peptide_length = None,
        min_count = 0,
        assay_group=None,
        nrows = None,
        group_by_allele = False,
        verbose= True):
    """
    Percentage positive results for both T-cell response assays
    and MHC binding assays (keyed by epitopes for which we have data
    for both)
    """
    mhc = load_mhc(
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            assay_group=assay_group,
            nrows = nrows,
            min_count = min_count,
            group_by_allele = group_by_allele,
            verbose = verbose)
    tcell = load_tcell(
                mhc_class = mhc_class,
                hla_type = hla_type,
                exclude_hla_type = exclude_hla_type,
                assay_group=assay_group,
                peptide_length = peptide_length,
                nrows = nrows,
                min_count = min_count,
                group_by_allele = group_by_allele,
                verbose = verbose)
    df_combined = pd.DataFrame({'mhc':mhc, 'tcell':tcell})
    both = ~(df_combined.mhc.isnull() | df_combined.tcell.isnull())
    return df_combined[both]

