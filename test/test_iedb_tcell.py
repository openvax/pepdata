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

from pepdata import iedb, reduced_alphabet

def test_tcell_hla_restrict_a24():
    """
    IEDB T-cell: Test that HLA restriction actually decreases
    number of results and that regular expression patterns
    are being used correctly
    """
    df_all = iedb.tcell.load_dataframe(nrows=1000)
    df_a24_1 = iedb.tcell.load_dataframe(hla='HLA-A24', nrows=1000)
    df_a24_2 = iedb.tcell.load_dataframe(hla='HLA-A\*24', nrows=1000)
    df_a24_combined = \
        iedb.tcell.load_dataframe(hla='HLA-A24|HLA-A\*24', nrows=1000)
    assert len(df_a24_1) < len(df_all)
    assert len(df_a24_2) < len(df_all)
    assert len(df_a24_combined) <= \
        len(df_a24_1) + len(df_a24_2), \
        "Expected %d <= %d + %d" % \
        (len(df_a24_combined), len(df_a24_1), len(df_a24_2))

def test_tcell_hla_exclude_a0201():
    """
    Test that excluding HLA allele A*02:01
    actually returns a DataFrame not containing
    that allele
    """
    df_all = iedb.tcell.load_dataframe(nrows=1000)
    assert (df_all['MHC Allele Name'] == "HLA-A*02:01").any()

    df_exclude = iedb.tcell.load_dataframe(
        nrows=1000,
        exclude_hla="HLA-A\*02:01")

    n_A0201_entries = (df_exclude['MHC Allele Name'] == "HLA-A*02:01").sum()
    assert n_A0201_entries == 0, \
        ("Not supposed to contain HLA-A*02:01, "
         " but found %d rows of that allele") % n_A0201_entries

def test_tcell_reduced_alphabet():
    """
    IEBD T-cell: Changing to a binary amino acid alphabet should reduce
    the number of samples since some distinct 20-letter strings collide
    as 2-letter strings
    """
    imm, non = iedb.tcell.load_classes(nrows=100)
    imm2, non2 = \
        iedb.tcell.load_classes(
            nrows=100,
            reduced_alphabet=reduced_alphabet.hp2)
    assert len(imm) + len(non) > len(imm2) + len(non2)
