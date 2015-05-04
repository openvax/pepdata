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
Construct DataFrames which join multiple IEDB datasets
"""

from __future__ import print_function, division, absolute_import

import pandas as pd

from ..common import memoize
from .mhc import load_groups as load_mhc_values
from .tcell import load_groups as load_tcell_values

@memoize
def load_tcell_vs_mhc(
        mhc_class=None,  # 1, 2, or None for neither
        hla=None,
        exclude_hla=None,
        peptide_length=None,
        min_count=0,
        mhc_assay_method=None,
        mhc_assay_group=None,
        tcell_assay_method=None,
        tcell_assay_group=None,
        nrows=None,
        group_by_allele=False):
    """
    Percentage positive results for both T-cell response assays
    and MHC binding assays (keyed by epitopes for which we have data
    for both)
    """
    mhc = load_mhc_values(
            mhc_class=mhc_class,
            hla=hla,
            exclude_hla=exclude_hla,
            peptide_length=peptide_length,
            assay_method=mhc_assay_method,
            assay_group=mhc_assay_group,
            nrows=nrows,
            min_count=min_count,
            group_by_allele=group_by_allele)

    tcell = load_tcell_values(
                mhc_class=mhc_class,
                hla=hla,
                exclude_hla=exclude_hla,
                assay_method=tcell_assay_method,
                assay_group=tcell_assay_group,
                peptide_length=peptide_length,
                nrows=nrows,
                min_count=min_count,
                group_by_allele=group_by_allele)

    df_combined = pd.DataFrame({'mhc': mhc.value, 'tcell': tcell.value})
    both = ~(df_combined.mhc.isnull() | df_combined.tcell.isnull())
    return df_combined[both]
