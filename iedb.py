import numpy as np
import pandas as pd

def load_csv(filename = 'data/tcell_compact.csv',
             assay_group=None,
             unique_sequences = True,
             noisy_labels = 'majority', # 'majority', 'percent', 'negative', 'positive'
             human = True,
             hla_type = None, # 1, 2, or None for neither
             hla_type1 = False,  # deprecated, kept for old scripts
             exclude_hla_a2 = False,
             only_hla_a2 = False,
             peptide_length = None,
             nrows = None,
             min_count = 0,
             key_by_allele = False,
             return_dataframe = False,
             verbose = True):

    df = pd.read_csv(filename, skipinitialspace=True, nrows = nrows)
    mhc = df['MHC Allele Name']

    #
    # Match known alleles such as 'HLA-A*02:01',
    # broader groupings such as 'HLA-A2'
    # and unknown alleles of the MHC-1 listed either as
    #  'HLA-Class I,allele undetermined'
    #  or
    #  'Class I,allele undetermined'
    class_1_mhc_mask = mhc.str.contains('Class I,|HLA-[A-C]([0-9]|\*)', na=False).astype('bool')
    class_2_mhc_mask = mhc.str.contains("Class II,|HLA-D(P|M|O|Q|R)", na=False).astype('bool')

    if verbose:
        print "Class I MHC Entries", class_1_mhc_mask.sum()
        print "Class II MHC Entries", class_2_mhc_mask.sum()

    # just in case any of the results were from mice or other species,
    # restrict to humans
    human_mask = df['Host Organism Name'].str.startswith('Homo sapiens', na=False).astype('bool') | df["MHC Allele Name"].str.startswith("HLA", na=False).astype('bool')

    if verbose:
        print "Human entries", human_mask.sum()
        print "Human Class I MHCs", (human_mask & class_1_mhc_mask).sum()
        print "Human Class II MHCs", (human_mask & class_2_mhc_mask).sum()

    null_epitope_seq = df['Epitope Linear Sequence'].isnull()
    if verbose:
        print "Dropping %d null sequences" % null_epitope_seq.sum()

    # if have rare or unknown amino acids, drop the sequence
    bad_epitope_seq = df['Epitope Linear Sequence'].str.contains('u|x|j|b|z|U|X|J|B|Z', na=False).astype('bool')
    if verbose:
        print "Dropping %d bad sequences" % bad_epitope_seq.sum()
    has_epitope_seq = ~(bad_epitope_seq | null_epitope_seq)

    mask = has_epitope_seq
    if human:
        mask &= human_mask
    if hla_type1 or hla_type == 1:
        mask &= class_1_mhc_mask
    if hla_type == 2:
        mask &= class_2_mhc_mask

    if assay_group:
        mask &= df['Assay Group'] == assay_group

    hla_a2_mask = (mhc == 'HLA-A2') | mhc.str.startswith('HLA-A\*02', na=False)
    if exclude_hla_a2:
        mask &= ~hla_a2_mask

    if only_hla_a2:
        mask &= hla_a2_mask

    epitopes = df['Epitope Linear Sequence'].str.upper()

    if peptide_length:
        assert peptide_length > 0
        mask &=  epitopes.str.len() == peptide_length

    if verbose:
        print "Filtered sequences epitope sequences", mask.sum()

    df = df[mask]

    if return_dataframe:
        return df

    pos_mask = df['Qualitative Measure'].str.startswith('Positive').astype('bool')

    if key_by_allele:
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

def load_tcell(assay_group=None,
               hla_type = None, # 1, 2, or None for neither
               peptide_length = None,
               nrows = None,
               min_count = 0,
               key_by_allele = False,
               return_dataframe = False,
               verbose = True):
    return load_csv('data/tcell_compact.csv',
                    assay_group = assay_group,
                    noisy_labels = 'percent',
                    hla_type = hla_type,
                    peptide_length = peptide_length,
                    nrows = nrows,
                    min_count = min_count,
                    key_by_allele = key_by_allele,
                    return_dataframe = return_dataframe,
                    verbose = verbose)

def load_mhc(assay_group=None,
             hla_type = None, # 1, 2, or None for neither
             peptide_length = None,
             nrows = None,
             min_count = 0,
             key_by_allele = False,
             return_dataframe = False,
             verbose = True):
    return load_csv('data/elution_compact.csv',
                    assay_group = assay_group,
                    noisy_labels = 'percent',
                    hla_type = hla_type,
                    peptide_length = peptide_length,
                    nrows = nrows,
                    min_count = min_count,
                    key_by_allele = key_by_allele,
                    return_dataframe = return_dataframe,
                    verbose = verbose)

def load_tcell_vs_mhc(assay_group=None,
                      hla_type = None, # 1, 2, or None for neither
                      peptide_length = None,
                      nrows = None,
                      min_count = 0,
                      key_by_allele = False,
                      verbose= True):
    mhc = load_mhc(assay_group=assay_group,
                       hla_type = hla_type,
                       peptide_length = peptide_length,
                       nrows = nrows,
                       min_count = min_count,
                       key_by_allele = key_by_allele,
                       verbose = verbose)
    tcell = load_tcell(assay_group=assay_group,
                       hla_type = hla_type,
                       peptide_length = peptide_length,
                       nrows = nrows,
                       min_count = min_count,
                       key_by_allele = key_by_allele)

    df_combined = pd.DataFrame({'mhc':mhc, 'tcell':tcell})
    both = ~(df_combined.mhc.isnull() | df_combined.tcell.isnull())
    return df_combined[both]

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
