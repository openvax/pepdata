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
T-cell epitopes from Dana-Farber Repository for Machine Learning in Immunology
http://bio.dfci.harvard.edu/DFRMLI/HTML/TCellEpitopes.php
"""

from __future__ import print_function, division, absolute_import
import pandas as pd

from .common import bad_amino_acids, fetch_file

def load_tumor(
        peptide_length=None,
        hla_type=None,
        source_protein=None,
        nrows=None):
    """
    Tumor antigens.
    This data set is a list of 718 T cell epitopes, 8-31 amino acids in length,
    derived from human tumor antigens. They are restricted by multiple HLA class
    I and class II alleles. The epitopes were experimentally validated based on
    immune recognition of HLA molecules and T cell stimulation. It could be used
    as a validation dataset for computational models to predict T cell epitopes.
    """
    path = fetch_file(
        filename="tumor_epitopes.csv",
        download_url="http://bio.dfci.harvard.edu/DFRMLI/datasets/tumor_epitopes.htm",
        subdir="pepdata")
    df = pd.read_csv(path, skipinitialspace=True, nrows=nrows)
    mask = ~df.Peptide.str.contains(bad_amino_acids)
    if peptide_length:
        mask &= df.Peptide.str.len() == peptide_length
    if hla_type:
        mask &= df.Allele.str.contains(hla_type)
    if source_protein:
        mask &= df.Protein.str.contains(source_protein)
    return df[mask]

def load_tumor_set(*args, **kwargs):
    df = load_tumor(*args, **kwargs)
    return set(df.Peptide)


def load_virus(
        peptide_length=None,
        hla_type=None,
        source_protein=None,
        nrows=None):
    """
    Virus antigens.
    This data set is a list of 44 HLA-A2 restricted T cell epitopes,
    9 amino acids in length, derived from human medically important viruses.
    The epitopes were experimentally characterized based on immune recognition
    of HLA molecule and T cell stimulation. It could be used as a validation
    dataset for computational models to predict T cell epitopes.
    """
    path = fetch_file(
        filename="virus_epitopes_A2.csv",
        download_url="http://bio.dfci.harvard.edu/DFRMLI/datasets/virus_epitopes_A2.htm",
        subdir="pepdata")
    df = pd.read_csv(path, skipinitialspace=True, nrows=nrows)
    mask = ~df.Epitope.str.contains(bad_amino_acids)
    if peptide_length:
        mask &= df.Epitope.str.len() == peptide_length
    if hla_type:
        mask &= df.Allele.str.contains(hla_type)
    if source_protein:
        mask &= df['Virus protein'].str.contains(source_protein)
    return df[mask]

def load_virus_set(*args, **kwargs):
    df = load_virus(*args, **kwargs)
    return set(df.Epitope)


def load_cef(
        peptide_length=None,
        hla_type=None,
        source_protein=None,
        nrows=None):
    """
    This dataset is a list of 32 T cell epitopes, 8-12 amino acids in length,
    with sequences derived from the human Cytomegalovirus, Epstein-Barr Virus
    and Influenza Virus.
    """
    path = fetch_file(
        filename="CEF.csv",
        download_url="http://bio.dfci.harvard.edu/DFRMLI/datasets/CEF.htm")
    df = pd.read_csv(path, skipinitialspace=True, nrows=nrows)
    mask = ~df.Peptide.str.contains(bad_amino_acids)
    if peptide_length:
        mask &= df.Peptide.str.len() == peptide_length
    if hla_type:
        mask &= df.Allele.str.contains(hla_type)
    if source_protein:
        mask &= df.Protein.str.contains(source_protein)
    return df[mask]

def load_cef_set(*args, **kwargs):
    df = load_cef(*args, **kwargs)
    return set(df.Peptide)
