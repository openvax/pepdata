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
Peptide response data from:
"Consistent CTL targeting of immunodominant regions in HIV across different ethnicities." by Frahm et al.
http://www.hiv.lanl.gov/content/immunology/hlatem/study1/index.html
"""

import numpy as np
import pandas as pd

from download import fetch_data

def load_dataframe(min_count = None, max_count = None):
    url = \
        "http://www.hiv.lanl.gov/content/immunology/hlatem/study1/peptides.html"
    df = pd.read_csv(fetch_data('frahm.csv', url))
    peptides = df['Sequence']
    # header parsing of the table messes up, so
    # reactions gets labeled as nan and
    # entropy  gets labeled as 'nan.1'

    reactions = df['Unnamed: 7']
    entropy = df['nan.1']
    df = pd.DataFrame({
        'Peptide' : peptides,
        'Reactions' : reactions,
        'Entropy' : entropy
    })
    if min_count is not None:
        df = df.ix[df.Reactions >= min_count]
    if max_count is not None:
        df = df.ix[df.Reactions <= max_count]
    return df

def load_set(*args, **kwargs):
    df = load_dataframe(*args, **kwargs)
    return set(df.Peptide)


