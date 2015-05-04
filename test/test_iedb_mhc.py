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

from __future__ import print_function, division, absolute_import

from pepdata import iedb, reduced_alphabet

def test_mhc_hla_a2():
    """
    IEDB MHC: Test that HLA restriction actually decreases number of results and
    that regular expression patterns are being used correctly
    """
    df_all = iedb.mhc.load_dataframe(nrows=1000)
    df_a2_1 = iedb.mhc.load_dataframe(hla='HLA-A2', nrows=1000)
    df_a2_2 = iedb.mhc.load_dataframe(hla='HLA-A\*02', nrows=1000)
    df_a2_combined = iedb.mhc.load_dataframe(hla='HLA-A2|HLA-A\*02', nrows=1000)
    assert len(df_a2_1) < len(df_all)
    assert len(df_a2_2) < len(df_all)
    assert len(df_a2_combined) <= len(df_a2_1) + len(df_a2_2), \
        "Expected %d <= %d + %d" % \
        (len(df_a2_combined), len(df_a2_1), len(df_a2_2))

def test_mhc_reduced_alphabet():
    pos, neg = iedb.mhc.load_classes(nrows=100)
    pos2, neg2 = iedb.mhc.load_classes(
        nrows=100,
        reduced_alphabet=reduced_alphabet.hp2)
    assert len(pos) + len(neg) > len(pos2) + len(neg2)
