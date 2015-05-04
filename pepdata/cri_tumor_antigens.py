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
Cancer Research Institute's
Table of Tumor Antigens Resulting from Mutations
http://cancerimmunity.org/peptide/mutations/
"""

from __future__ import print_function, division, absolute_import
from os.path import join

import pandas as pd

from .static_data import DATA_DIR

def load_dataframe(
        mhc_class=None,
        hla_type=None):
    path = join(DATA_DIR, 'cri_mutations.csv')
    df = pd.read_csv(path)
    mhc2 = df.HLA.str.startswith('D')
    if mhc_class == 1:
        return df.ix[~mhc2]
    elif mhc_class == 2:
        return df.ix[mhc2]
    elif hla_type:
        return df.ix[df.HLA == hla_type]
    else:
        return df

def load_peptides(*args, **kwargs):
    df = load_dataframe(*args, **kwargs)
    return set(df.Peptide)
