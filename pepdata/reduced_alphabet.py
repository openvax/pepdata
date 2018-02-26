# Copyright (c) 2014-2018. Mount Sinai School of Medicine
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
Amino acid groupings from
'Reduced amino acid alphabets improve the sensitivity...' by
Peterson, Kondev, et al.
http://www.rpgroup.caltech.edu/publications/Peterson2008.pdf
"""
from __future__ import print_function, division, absolute_import

def dict_from_list(groups):
    aa_to_group = {}
    for i, group in enumerate(groups):
        for c in group:
            aa_to_group[c] = group[0]
    return aa_to_group

gbmr4 = dict_from_list(["ADKERNTSQ", "YFLIVMCWH", "G", "P"])

sdm12 = dict_from_list([
    "A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"
])

hsdm17 = dict_from_list([
    "A", "D", "KE", "R", "N", "T", "S", "Q", "Y",
    "F", "LIV", "M", "C", "W", "H", "G", "P"
])

"""
Other alphabets from
http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html
"""

# hydrophilic vs. hydrophobic
hp2 = dict_from_list(["AGTSNQDEHRKP", "CMFILVWY"])

murphy10 = dict_from_list([
    "LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"
])

alex6 = dict_from_list(["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"])

aromatic2 = dict_from_list(["FHWY", "ADKERNTSQLIVMCGP"])

hp_vs_aromatic = dict_from_list(["H", "CMILV", "FWY", "ADKERNTSQGP"])
