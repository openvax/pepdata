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
Amino acid groupings from
'Reduced amino acid alphabets improve the sensitivity...' by
Peterson, Kondev, et al.
http://www.rpgroup.caltech.edu/publications/Peterson2008.pdf
"""
from __future__ import print_function, division, absolute_import

def dict_from_list(groups):
    result = {}
    for i, group in enumerate(groups):
        for c in group:
            result[c.upper()] = i
            result[c.lower()] = i
    return result

gbmr4 = dict_from_list(["ADKERNTSQ", "YFLIVMCWH", "G", "P"])

sdm12 = dict_from_list(
  ["A", "D", "KER", "N", "TSQ", "YF", "LIVM", "C", "W", "H", "G", "P"]
)

hsdm17 = dict_from_list(
  ["A", "D", "KE", "R", "N", "T", "S", "Q", "Y",
   "F", "LIV", "M", "C", "W", "H", "G", "P"])

"""
Other alphabets from
http://bio.math-inf.uni-greifswald.de/viscose/html/alphabets.html
"""

# hydrophilic vs. hydrophobic
hp2 = dict_from_list(["AGTSNQDEHRKP", "CMFILVWY"])

murphy10 = dict_from_list(
  ["LVIM", "C", "A", "G", "ST", "P", "FYW", "EDNQ", "KR", "H"]
)

alex6 = dict_from_list(["C", "G", "P", "FYW", "AVILM", "STNQRHKDE"])

aromatic2 = dict_from_list(["FHWY", "ADKERNTSQLIVMCGP"])

hp_vs_aromatic = dict_from_list(["H", "CMILV", "FWY", "ADKERNTSQGP"])

class AlphabetTransformer(object):
    def __init__(self, reduced_alphabet):
        assert isinstance(reduced_alphabet, dict), \
            "Expected dictionary, got %s" % (reduced_alphabet,)
        self.reduced_alphabet = reduced_alphabet

    def __call__(self, s):
        return self.transform(s)

    def transform(self, s):
        d = self.reduced_alphabet
        return ''.join([chr(48 + d[char]) for char in s])

    def __getstate__(self):
        return {'reduced_alphabet': self.reduced_alphabet}

def make_alphabet_transformer(reduced_alphabet):
    if isinstance(reduced_alphabet, (str, unicode)):
        reduced_alphabet = globals()[reduced_alphabet]
    return AlphabetTransformer(reduced_alphabet)
