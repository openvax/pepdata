# Copyright (c) 2014-2016. Mount Sinai School of Medicine
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


from __future__ import print_function, division, absolute_import

import numpy as np

def transform_peptide(peptide, property_dict):
    return np.array([property_dict[amino_acid] for amino_acid in peptide])

def transform_peptides(peptides, property_dict):
    return np.array([
        [property_dict[aa] for aa in peptide]
        for peptide in peptides])

