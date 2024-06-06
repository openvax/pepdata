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

def find(df : pd.DataFrame, group_candidates : list[str], column_candidates : list[str]) -> pd.Series | None: 
    """
    Try to find a column that contains a combination of the two candidate lists. 

    Motivation: format for MHC ligand CSV used to have:
        epitope_key = ("Epitope", "Description")
        mhc_allele_key = ("MHC", "Allele Name")
        mhc_class_key = ("MHC", "MHC allele class")
        mhc_assay_key = ("Assay", "Method/Technique")

    Now it's:
        epitope_key = ("Epitope", "Name")
        mhc_allele_key = ("MHC Restriction", "Name")
        mhc_class_key = ("MHC Restriction", "Class")
        mhc_assay_key = ("Assay", "Method")

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
    return find(df, group_candidates, column_candidates)


def get_mhc_class(
        df : pd.DataFrame, 
        group_candidates : list[str] = MHC_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Class", "MHC allele class"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_mhc_assay(
        df : pd.Series, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] =["method"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_epitope_name(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["name"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_epitope_type(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Object Type", "Type"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

def get_epitope_modifications(
         df : pd.DataFrame,
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES,
        column_candidates : list[str]  = ["Modified Residue(s)"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_epitope_IRI(
        df : pd.DataFrame,
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] =["Epitope IRI"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_epitope_source_molecule(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates=["Source Molecule"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

def get_epitope_source_molecule_iri(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Source Molecule IRI"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_epitope_source_molecule_iri(
        df : pd.DataFrame, 
        group_candidates : list[str] = EPITOPE_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Source Molecule IRI"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)



def get_assay_method(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Method", "Method/Technique"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

def get_assay_response_measured(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Response measured"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_assay_units(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Units"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)


def get_assay_qualitative(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Qualitative Measurement"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

def get_assay_num_tested(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Number of Subjects Tested"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

def get_assay_num_responded(
        df : pd.DataFrame, 
        group_candidates : list[str] = ASSAY_GROUP_CANDIDATES, 
        column_candidates : list[str] = ["Number of Subjects Responded"]) -> pd.Series | None:
    return find(df, group_candidates, column_candidates)

