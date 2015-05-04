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

"""
TANTIGEN: Tumor T-cell Antigen Database from Dana Farber CVC
http://cvc.dfci.harvard.edu/tadb/index.html
"""


from __future__ import print_function, division, absolute_import
from os.path import join

import pandas as pd

from .common import bad_amino_acids
from .features import make_unlabeled_ngram_dataset_from_args
from .reduced_alphabet import make_alphabet_transformer
from .static_data import DATA_DIR

def _load_dataframe(
        path,
        epitope_column_name,
        mhc_class=None,
        hla_type=None,
        exclude_hla_type=None,  # regex pattern i.e. '(HLA-A2)|(HLA-A\*02)'
        peptide_length=None,
        reduced_alphabet=None,
        nrows=None):
    df = pd.read_csv(path, skipinitialspace=True, nrows=nrows)
    epitopes = df[epitope_column_name]
    hla = df['HLA allele']
    mask = ~(epitopes.str.contains(bad_amino_acids, na=False).astype('bool'))
    if mhc_class == 1:
        a = hla.str.startswith('A')
        b = hla.str.startswith('B')
        c = hla.str.startswith('C')
        mask &= (a | b | c)
    elif mhc_class == 2:
        mask &= hla.str.startswith('D')
    if hla_type:
        mask &= hla.str.contains(hla_type, na=False).astype('bool')
    if exclude_hla_type:
        mask &= ~(hla.str.contains(exclude_hla_type, na=True).astype('bool'))
    if peptide_length:
        mask &= epitopes.str.len() == peptide_length
    df = df[mask]
    if reduced_alphabet:
        epitopes = df[epitope_column_name]
        df[epitope_column_name] = \
            epitopes.map(make_alphabet_transformer(reduced_alphabet))
    return df

def load_tcell(*args, **kwargs):
    """
    Return a dataframe with accession IDs, peptide sequence, and MHC allele
    for T-cell responsive tumor antigens
    """
    tcell_path = join(DATA_DIR, 'tantigen_tcell.csv')
    return _load_dataframe(tcell_path, 'Epitope sequence', *args, **kwargs)

def load_tcell_set(*args, **kwargs):
    df = load_tcell(*args, **kwargs)
    return set(df['Epitope sequence'])

def load_tcell_ngrams(*args, **kwargs):
    return make_unlabeled_ngram_dataset_from_args(
        load_tcell_set, *args, **kwargs)

def load_mhc(*args, **kwargs):
    """
    Return a dataframe with accession IDs, peptide sequence, and MHC allele
    for MHC-binding tumor antigens
    """
    mhc_path = join(DATA_DIR, 'tantigen_mhc.csv')
    return _load_dataframe(mhc_path, 'Ligand sequence', *args, **kwargs)

def load_mhc_set(*args, **kwargs):
    df = load_mhc(*args, **kwargs)
    return set(df['Ligand sequence'])

def load_mhc_ngrams(*args, **kwargs):
    return make_unlabeled_ngram_dataset_from_args(load_mhc_set, *args, **kwargs)