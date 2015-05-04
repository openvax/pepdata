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
Data from 'HLA-binding properties of tumor neoepitopes in humans'
by Edward F Fritsch, Mohini Rajasagi, Patrick A Ott, et al.
Cancer Immunology Research
"""

from __future__ import print_function, division, absolute_import
from os.path import join
import pandas as pd

from .static_data import DATA_DIR

def load_dataframe(
        hla_type=None,
        exclude_hla_type=None):
    path = join(DATA_DIR, 'fritsch2014_neoepitopes.csv')
    df = pd.read_csv(path, skipinitialspace=True)
    hla = df['HLA Allele']
    if hla_type:
        df = df[hla.str.contains(hla_type, na=False).astype('bool')]
    if exclude_hla_type:
        df = df[~(hla.str.contains(exclude_hla_type, na=True).astype('bool'))]
    return df
