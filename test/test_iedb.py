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

from epitopes import iedb, reduced_alphabet

def test_tcell_hla_a24():
    """
    IEDB T-cell:
    Test that HLA restriction actually decreases number of results and
    that regular expression patterns are being used correctly
    """
    df_all = iedb.load_tcell(nrows=1000)
    df_a24_1 = iedb.load_tcell(hla_type='HLA-A24', nrows=1000)
    df_a24_2 = iedb.load_tcell(hla_type='HLA-A\*24', nrows=1000)
    df_a24_combined = \
        iedb.load_tcell(hla_type = 'HLA-A24|HLA-A\*24', nrows=1000)
    assert len(df_a24_1) < len(df_all)
    assert len(df_a24_2) < len(df_all)
    assert len(df_a24_combined) <= \
        len(df_a24_1) + len(df_a24_2), \
        "Expected %d <= %d + %d" % \
        (len(df_a24_combined), len(df_a24_1), len(df_a24_2))

def test_tcell_reduced_alphabet():
    """
    IEBD T-cell:
    Changing to a binary amino acid alphabet should reduce the number of
    samples since some distinct 20-letter strings collide as 2-letter strings
    """
    imm, non = iedb.load_tcell_classes(nrows = 100)
    imm2, non2 = \
        iedb.load_tcell_classes(
            nrows = 100,
            reduced_alphabet = reduced_alphabet.hp2)
    assert len(imm) + len(non) > len(imm2) + len(non2)
def test_mhc_hla_a2():
    """
    IEDB MHC: Test that HLA restriction actually decreases number of results and
    that regular expression patterns are being used correctly
    """
    df_all = iedb.load_mhc()
    df_a2_1 = iedb.load_mhc(hla_type='HLA-A2', nrows=1000)
    df_a2_2 = iedb.load_mhc(hla_type='HLA-A\*02', nrows=1000)
    df_a2_combined = iedb.load_mhc(hla_type = 'HLA-A2|HLA-A\*02', nrows=1000)
    assert len(df_a2_1) < len(df_all)
    assert len(df_a2_2) < len(df_all)
    assert len(df_a2_combined) <= len(df_a2_1) + len(df_a2_2), \
        "Expected %d <= %d + %d" % \
        (len(df_a2_combined), len(df_a2_1), len(df_a2_2))

def test_mhc_reduced_alphabet():
    pos, neg = iedb.load_mhc_classes(nrows = 100)
    pos2, neg2 = iedb.load_mhc_classes(
        nrows = 100,
        reduced_alphabet = reduced_alphabet.hp2)
    assert len(pos) + len(neg) > len(pos2) + len(neg2)
