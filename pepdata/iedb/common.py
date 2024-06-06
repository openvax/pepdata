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

from __future__ import annotations

import pandas as pd

import datacache

cache = datacache.Cache("pepdata")



bad_amino_acids = 'U|X|J|B|Z'

def find_column(df : pd.DataFrame, group_candidates : list[str], column_candidates : list[str]) -> pd.Series | None: 
    """
    Try to find a column that contains a combination of the two candidate lists. 

    Motivation: format for MHC ligand CSV used to have:
        epitope_column_key = ("Epitope", "Description")
        mhc_allele_column_key = ("MHC", "Allele Name")
        mhc_class_column_key = ("MHC", "MHC allele class")
        mhc_assay_column_key = ("Assay", "Method/Technique")

    Now it's:
        epitope_column_key = ("Epitope", "Name")
        mhc_allele_column_key = ("MHC Restriction", "Name")
        mhc_class_column_key = ("MHC Restriction", "Class")
        mhc_assay_column_key = ("Assay", "Method")

    ...who knows what it will be next!
    """
    group_candidates = [s.lower() for s in group_candidates]
    column_candidates = [s.lower() for s in column_candidates]

    possible_matches = []
    for a in group_candidates:
        for b in column_candidates:
            for pair in df.columns:
                assert type(pair) is tuple and len(pair) == 2
                group, col = pair
                if a in group.lower() and b in col.lower():
                    possible_matches.append(pair)
    
    if len(possible_matches) == 0:
        return None
    # get the shortest matches
    


MHC_GROUP_CANDIDATES : list[str] = ["MHC", "MHC Restriction"]
EPITOPE_GROUP_CANDIDATES : list[str] = ["Epitope"] 
ASSAY_GROUP_CANDIDATES : list[str] = ["Assay"]

def get_mhc_allele_column(
        df : pd.DataFrame, 
        group_candidates : list[str] = MHC_GROUP_CANDIDATES,
        column_candidates : list[str] = ["Allele", "Allele name", "Name"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_mhc_class_column(
        df : pd.DataFrame, 
        group_candidates : list[str] = MHC_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Class", "MHC allele class"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_mhc_assay_column(
        df : pd.Series, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] =["method"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_epitope_name_column(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["name"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_epitope_type_column(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Object Type", "Type"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)

def get_epitope_modifications_column(
         df : pd.DataFrame,
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES,
        column_candidates : list[str]  = ["Modified Residue(s)"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_epitope_IRI_column(
        df : pd.DataFrame,
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Epitope IRI"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)


def get_epitope_source_molecule_column(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates=["Source Molecule"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)

def get_epitope_source_molecule_iri_column(
        df : pd.DataFrame, 
        group_candidates : listr[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Source Molecule IRI"]) -> pd.Series | None:
    return find_column(df, group_candidates, column_candidates)

Assay Columns:

'Method', 'Response measured',
       'Units', 'IRI', 'Qualitative Measurement', 'Measurement Inequality',
       'Quantitative measurement', 'Number of Subjects Tested',
       'Number of Subjects Responded', 'Response Frequency (%)', 'PDB ID',
       'Comments'],
"IEDB IRI","Object Type",Name,"Modified Residue(s)",Modifications,"Starting Position","Ending Position",IRI,Synonyms,"Source Molecule","Source Molecule IRI"
    