import numpy as np
import pandas as pd

def load_dataframe(
        filename = 'data/tcell_compact.csv',
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

    df = df[mask]
    return df

def group_epitopes(
        df,
        unique_sequences = True,
        noisy_labels = 'majority', # 'majority', 'percent', 'negative', 'positive'
        min_count = 0,
        group_by_allele = False,
        verbose = True):
    """
    Given a dataframe of epitopes and qualitative measures,
    group the epitope strings (optionally also grouping by allele),
    and combine the response using the mechanism specified by the argument 'noisy_labels':
        - majority = Split into two sets, epitope is Positive if majority of its results are Positive
        - negative = epitope is Negative if any result is Negative
        - positive = epitope is Positive if any result is Positive
        - percent = instead of returning two sets, return a value in [0,1] indicating the percentage of positive results
    """
    measure = df['Qualitative Measure']
    pos_mask = measure.str.startswith('Positive').astype('bool')

    mhc = df['MHC Allele Name']
    epitopes = df['Epitope Linear Sequence'].str.upper()

    if group_by_allele:
        groups = pos_mask.groupby([epitopes, mhc])
    else:
        groups = pos_mask.groupby(epitopes)

    values = groups.mean()

    if min_count:
        counts = groups.count()
        values = values[counts > min_count]

    if noisy_labels == 'percent':
        return values
    elif noisy_labels == 'majority':
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


def load_tcell_prct(
        mhc_class = None, # 1, 2, or None for neither
        hla_type = None,
        exclude_hla_type = None,
        peptide_length = None,
        min_count = 0,
        assay_group=None,
        nrows = None,
        group_by_allele = False,
        verbose= True):

    df = load_dataframe(
            'data/tcell_compact.csv',
            assay_group = assay_group,
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            nrows = nrows,
            verbose = verbose)
    return group_epitopes(
            df,
            noisy_labels = 'percent',
            group_by_allele = group_by_allele,
            min_count = min_count,
            verbose = verbose)

def load_mhc_prct(
        mhc_class = None, # 1, 2, or None for neither
        hla_type = None,
        exclude_hla_type = None,
        peptide_length = None,
        min_count = 0,
        assay_group=None,
        nrows = None,
        group_by_allele = False,
        verbose= True):
    df = load_dataframe(
            'data/elution_compact.csv',
            assay_group = assay_group,
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            nrows = nrows,
            verbose = verbose)

    return group_epitopes(
            df,
            noisy_labels = 'percent',
            group_by_allele = group_by_allele,
            min_count = min_count,
            verbose = verbose)

def load_tcell_vs_mhc_prct(
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
    mhc = load_mhc_prct(
            mhc_class = mhc_class,
            hla_type = hla_type,
            exclude_hla_type = exclude_hla_type,
            peptide_length = peptide_length,
            assay_group=assay_group,
            nrows = nrows,
            min_count = min_count,
            group_by_allele = group_by_allele,
            verbose = verbose)
    tcell = load_tcell_prct(
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


"""
import numpy as np
import amino_acid
from amino_acid import letter_to_index

fns = [amino_acid.hydropathy,
       amino_acid.volume,
       amino_acid.pK_side_chain,
       amino_acid.polarity,
       amino_acid.prct_exposed_residues,
       amino_acid.hydrophilicity,
       amino_acid.accessible_surface_area,
       amino_acid.local_flexibility,
       amino_acid.accessible_surface_area_folded,
       amino_acid.refractivity
       ]


import data
def load_dataset(filename = 'tcell_compact.csv',
                 assay_group=None,
                 unique_sequences = True,
                 # 'drop' | 'keep' | 'positive' | 'negative'
                 noisy_labels = 'majority',
                 human = True,
                 hla_type1 = True,
                 exclude_hla_a2 = False,
                 only_hla_a2 = False,
                 max_ngram = 1,
                         normalize_row = True,
                         reduced_alphabet = None,
                         rebalance = False,
                         return_transformer = False):
    assert noisy_labels in ('majority', 'drop', 'keep', 'positive', 'negative'), \
      "Invalid option: %s" % noisy_labels


    imm, non = load_csv(filename,
       assay_group,
       unique_sequences,
       noisy_labels,
       human,
       hla_type1,
       exclude_hla_a2,
       only_hla_a2)
    return data.make_ngram_dataset(imm,non, max_ngram, normalize_row, reduced_alphabet, rebalance, return_transformer)
"""