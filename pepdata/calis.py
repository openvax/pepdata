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
Datasets from Calis, et al. "Properties of MHC Class I Presented Peptides That Enhance Immunogenicity"
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003266#s5
"""

from __future__ import print_function, division, absolute_import
from os.path import join

import pandas as pd

from .static_data import DATA_DIR
from .common import bad_amino_acids
from .features import make_ngram_dataset_from_args

def load_s1(
        human=True,
        peptide_length=None,
        hla_type=None,
        exclude_hla_type=None):
    """
    Immunogenic and non-immunogenic pMHCs that were found in the IEDB,
    Vaccinia, Arena and Coxiella data sets.
    """
    s1_path = join(DATA_DIR, 'calis_s1.csv')
    df = pd.read_csv(s1_path)
    mask = ~df.Peptide.str.contains(bad_amino_acids)
    if human:
        mask &= df.Species == 'Homo'
    if peptide_length:
        mask &= df.Peptide.str.len() == peptide_length
    if hla_type:
        mask &= df.MHC.str.contains(hla_type, na=False).astype('bool')
    if exclude_hla_type:
        contains_hla_type = df.MHC.str.contains(exclude_hla_type, na=False)
        mask &= ~(contains_hla_type.astype('bool'))
    return df[mask]

def load_s1_values(*args, **kwargs):
    group_by_allele = kwargs.pop('group_by_allele', False)
    df = load_s1(*args, **kwargs)
    pos = (df.Immunogenicity == 'immunogenic')
    if group_by_allele:
        groups = pos.groupby([df.Peptide, df.MHC])
    else:
        groups = pos.groupby(df.Peptide)
    return groups.mean()

def load_s1_classes(*args, **kwargs):
    values = load_s1_values(*args, **kwargs)
    imm_mask = values > 0
    non_mask = ~imm_mask
    imm = set(values.index[imm_mask])
    non = set(values.index[non_mask])
    return imm, non

def load_s1_ngrams(*args, **kwargs):
    return make_ngram_dataset_from_args(load_s1_classes, *args, **kwargs)

def load_s2(
        human=True,
        hla_type=None,
        exclude_hla_type=None):
    """
    Non-redundant murine and human Dengue epitopes and non-epitopes
    """
    s2_path = join(DATA_DIR, 'calis_s2.csv')
    df = pd.read_csv(s2_path)
    mask = ~df.peptide.str.contains(bad_amino_acids)
    if human:
        mask &= df.host == 'Homo'
    if hla_type:
        mask &= df.HLA.str.contains(hla_type, na=False).astype('bool')
    if exclude_hla_type:
        mask &= ~df.HLA.str.contains(exclude_hla_type, na=False).astype('bool')
    return df[mask]

def load_s2_values(*args, **kwargs):
    group_by_allele = kwargs.pop('group_by_allele', False)
    df = load_s2(*args, **kwargs)
    pos = (df['epitope status'] == 'epitope')
    if group_by_allele:
        groups = pos.groupby([df.peptide, df.HLA])
    else:
        groups = pos.groupby(df.peptide)
    return groups.mean()

def load_s2_classes(*args, **kwargs):
    values = load_s2_values(*args, **kwargs)
    imm_mask = values > 0
    non_mask = ~imm_mask
    imm = set(values.index[imm_mask])
    non = set(values.index[non_mask])
    return imm, non


def load_s2_ngrams(*args, **kwargs):
    return make_ngram_dataset_from_args(load_s2_classes, *args, **kwargs)
