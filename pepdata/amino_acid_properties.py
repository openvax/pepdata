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

from .amino_acid import letter_to_index

"""
Quantify amino acids by their physical/chemical properties
"""


class AminoAcidPropertyTable(object):
    @classmethod
    def _aa_dict_to_positional_list(value_dict):
        value_list = [None] * 20
        for letter, value in value_dict.items():
            idx = letter_to_index(letter)
            assert idx >= 0
            assert idx < 20
            value_list[idx] = value
        assert all(elt is not None for elt in value_list), \
            "Missing amino acids in:\n%s" % value_dict.keys()
        return value_list


    @classmethod
    def _parse_table(cls, table_string):
        value_dict = {}
        for line in table_string.splitlines():
            line = line.strip()
            if not line:
                continue
            fields = line.split(" ")
            fields = [f for f in fields if len(f.strip()) > 0]
            assert len(fields) >= 2
            value, letter = fields[:2]
            assert letter not in value_dict, "Repeated amino acid " + line
            value_dict[letter] = float(value)
        return value_dict

    def __init__(self, table):
        self.table = table
        self.value_dict = self._parse_table(self.table)
        self.value_list = self._aa_dict_to_positional_list(self.value_dict)

    def __call__(self, letter):
        return self.value_dict[letter]

    def __getitem__(self, letter):
        return self.value_dict[letter]

    def transform_string(self, letters):
        return np.array([self.value_dict[amino_acid] for amino_acid in letters])

    def transform_strings(self, strings):
        d = self.value_dict
        return np.array([[d[x] for x in strings] for x in strings])

"""
Amino acids property tables copied from CRASP website
"""

hydropathy = AminoAcidPropertyTable("""
1.80000 A ALA
-4.5000 R ARG
-3.5000 N ASN
-3.5000 D ASP
2.50000 C CYS
-3.5000 Q GLN
-3.5000 E GLU
-0.4000 G GLY
-3.2000 H HIS
4.50000 I ILE
3.80000 L LEU
-3.9000 K LYS
1.90000 M MET
2.80000 F PHE
-1.6000 P PRO
-0.8000 S SER
-0.7000 T THR
-0.9000 W TRP
-1.3000 Y TYR
4.20000 V VAL
""")

volume = AminoAcidPropertyTable("""
91.5000 A ALA
202.0000 R ARG
135.2000 N ASN
124.5000 D ASP
118.0000 C CYS
161.1000 Q GLN
155.1000 E GLU
66.40000 G GLY
167.3000 H HIS
168.8000 I ILE
167.9000 L LEU
171.3000 K LYS
170.8000 M MET
203.4000 F PHE
129.3000 P PRO
99.10000 S SER
122.1000 T THR
237.6000 W TRP
203.6000 Y TYR
141.7000 V VAL
""")

polarity = AminoAcidPropertyTable("""
0.0000 A ALA
52.000 R ARG
3.3800 N ASN
40.700 D ASP
1.4800 C CYS
3.5300 Q GLN
49.910 E GLU
0.0000 G GLY
51.600 H HIS
0.1500 I ILE
0.4500 L LEU
49.500 K LYS
1.4300 M MET
0.3500 F PHE
1.5800 P PRO
1.6700 S SER
1.6600 T THR
2.1000 W TRP
1.6100 Y TYR
0.1300 V VAL
""")

pK_side_chain = AminoAcidPropertyTable("""
0.0000 A ALA
12.480 R ARG
0.0000 N ASN
3.6500 D ASP
8.1800 C CYS
0.0000 Q GLN
4.2500 E GLU
0.0000 G GLY
6.0000 H HIS
0.0000 I ILE
0.0000 L LEU
10.530 K LYS
0.0000 M MET
0.0000 F PHE
0.0000 P PRO
0.0000 S SER
0.0000 T THR
0.0000 W TRP
10.700 Y TYR
0.0000 V VAL
""")

prct_exposed_residues = AminoAcidPropertyTable("""
15.0000 A ALA
67.0000 R ARG
49.0000 N ASN
50.0000 D ASP
5.00000 C CYS
56.0000 Q GLN
55.0000 E GLU
10.0000 G GLY
34.0000 H HIS
13.0000 I ILE
16.0000 L LEU
85.0000 K LYS
20.0000 M MET
10.0000 F PHE
45.0000 P PRO
32.0000 S SER
32.0000 T THR
17.0000 W TRP
41.0000 Y TYR
14.0000 V VAL
""")

hydrophilicity = AminoAcidPropertyTable("""
-0.5000 A ALA
3.00000 R ARG
0.20000 N ASN
3.00000 D ASP
-1.0000 C CYS
0.20000 Q GLN
3.00000 E GLU
0.00000 G GLY
-0.5000 H HIS
-1.8000 I ILE
-1.8000 L LEU
3.00000 K LYS
-1.3000 M MET
-2.5000 F PHE
0.00000 P PRO
0.30000 S SER
-0.4000 T THR
-3.4000 W TRP
-2.3000 Y TYR
-1.5000 V VAL
""")

accessible_surface_area = AminoAcidPropertyTable("""
27.8000 A ALA
94.7000 R ARG
60.1000 N ASN
60.6000 D ASP
15.5000 C CYS
68.7000 Q GLN
68.2000 E GLU
24.5000 G GLY
50.7000 H HIS
22.8000 I ILE
27.6000 L LEU
103.000 K LYS
33.5000 M MET
25.5000 F PHE
51.5000 P PRO
42.0000 S SER
45.0000 T THR
34.7000 W TRP
55.2000 Y TYR
23.7000 V VAL
""")

local_flexibility = AminoAcidPropertyTable("""
705.42000 A ALA
1484.2800 R ARG
513.46010 N ASN
34.960000 D ASP
2412.5601 C CYS
1087.8300 Q GLN
1158.6600 E GLU
33.180000 G GLY
1637.1300 H HIS
5979.3701 I ILE
4985.7300 L LEU
699.69000 K LYS
4491.6602 M MET
5203.8599 F PHE
431.96000 P PRO
174.76000 S SER
601.88000 T THR
6374.0698 W TRP
4291.1001 Y TYR
4474.4199 V VAL
""")

accessible_surface_area_folded = AminoAcidPropertyTable("""
31.5000 A ALA
93.8000 R ARG
62.2000 N ASN
60.9000 D ASP
13.9000 C CYS
74.0000 Q GLN
72.3000 E GLU
25.2000 G GLY
46.7000 H HIS
23.0000 I ILE
29.0000 L LEU
110.300 K LYS
30.5000 M MET
28.7000 F PHE
53.7000 P PRO
44.2000 S SER
46.0000 T THR
41.7000 W TRP
59.1000 Y TYR
23.5000 V VAL
""")

refractivity = AminoAcidPropertyTable("""
4.34000 A ALA
26.6600 R ARG
13.2800 N ASN
12.0000 D ASP
35.7700 C CYS
17.5600 Q GLN
17.2600 E GLU
0.00000 G GLY
21.8100 H HIS
19.0600 I ILE
18.7800 L LEU
21.2900 K LYS
21.6400 M MET
29.4000 F PHE
10.9300 P PRO
6.35000 S SER
11.0100 T THR
42.5300 W TRP
31.5300 Y TYR
13.9200 V VAL
""")

###
# Amino acid residue masses copied from:
# "Table 1. Amino acid residues sorted by name."
# http://www1.bioinfor.com/peaks/downloads/masstable.html

mass = dict(
    A=71.08,
    R=156.2,
    N=114.1,
    D=115.1,
    C=103.1,
    E=129.1,
    Q=128.1,
    G=57.05,
    H=137.1,
    I=113.2,
    L=113.2,
    K=128.2,
    M=131.2,
    F=147.2,
    P=97.12,
    S=87.08,
    T=101.1,
    W=186.2,
    Y=163.2,
    V=99.13)

# Chou-Fasman of structural properties from
# http://prowl.rockefeller.edu/aainfo/chou.htm
chou_fasman_table = """
Alanine        142     83       66      0.06    0.076   0.035   0.058
Arginine        98     93       95      0.070   0.106   0.099   0.085
Aspartic Acid  101     54      146      0.147   0.110   0.179   0.081
Asparagine      67     89      156      0.161   0.083   0.191   0.091
Cysteine        70    119      119      0.149   0.050   0.117   0.128
Glutamic Acid  151    037       74      0.056   0.060   0.077   0.064
Glutamine      111    110       98      0.074   0.098   0.037   0.098
Glycine         57     75      156      0.102   0.085   0.190   0.152
Histidine      100     87       95      0.140   0.047   0.093   0.054
Isoleucine     108    160       47      0.043   0.034   0.013   0.056
Leucine        121    130       59      0.061   0.025   0.036   0.070
Lysine         114     74      101      0.055   0.115   0.072   0.095
Methionine     145    105       60      0.068   0.082   0.014   0.055
Phenylalanine  113    138       60      0.059   0.041   0.065   0.065
Proline         57     55      152      0.102   0.301   0.034   0.068
Serine          77     75      143      0.120   0.139   0.125   0.106
Threonine       83    119       96      0.086   0.108   0.065   0.079
Tryptophan     108    137       96      0.077   0.013   0.064   0.167
Tyrosine        69    147      114      0.082   0.065   0.114   0.125
Valine         106    170       50      0.062   0.048   0.028   0.053
"""


def parse_chou_fasman(table):
    alpha_helix_score_dict = {}
    beta_sheet_score_dict = {}
    turn_score_dict = {}

    for line in table.split("\n"):
        fields = [field for field in line.split(" ") if len(field.strip()) > 0]
        if len(fields) == 0:
            continue

        if fields[1] == 'Acid':
            name = fields[0] + " " + fields[1]
            fields = fields[1:]
        else:
            name = fields[0]

        assert name in long_amino_acid_names, "Invalid amino acid name %s" % name
        idx = long_amino_acid_names.index(name)
        letter = amino_acid_letters[idx]
        alpha = int(fields[1])
        beta = int(fields[2])
        turn = int(fields[3])
        alpha_helix_score_dict[letter] = alpha
        beta_sheet_score_dict[letter] = beta
        turn_score_dict[letter] = turn

    assert len(alpha_helix_score_dict) == 20
    assert len(beta_sheet_score_dict) == 20
    assert len(turn_score_dict) == 20
    return alpha_helix_score_dict, beta_sheet_score_dict, turn_score_dict

alpha_helix_score, beta_sheet_score, turn_score = \
    parse_chou_fasman(chou_fasman_table)



###
# Values copied from:
# "Solvent accessibility of AA in known protein structures"
# http://prowl.rockefeller.edu/aainfo/access.htm
###
"""
Solvent accessibility of AA in known protein structures

Figure 1.

S   0.70    0.20    0.10
T   0.71    0.16    0.13
A   0.48    0.35    0.17
G   0.51    0.36    0.13
P   0.78    0.13    0.09
C   0.32    0.54    0.14
D   0.81    0.09    0.10
E   0.93    0.04    0.03
Q   0.81    0.10    0.09
N   0.82    0.10    0.08
L   0.41    0.49    0.10
I   0.39    0.47    0.14
V   0.40    0.50    0.10
M   0.44    0.20    0.36
F   0.42    0.42    0.16
Y   0.67    0.20    0.13
W   0.49    0.44    0.07
K   0.93    0.02    0.05
R   0.84    0.05    0.11
H   0.66    0.19    0.15
"""

solvent_exposed_area = dict(
    S=0.70,
    T=0.71,
    A=0.48,
    G=0.51,
    P=0.78,
    C=0.32,
    D=0.81,
    E=0.93,
    Q=0.81,
    N=0.82,
    L=0.41,
    I=0.39,
    V=0.40,
    M=0.44,
    F=0.42,
    Y=0.67,
    W=0.49,
    K=0.93,
    R=0.84,
    H=0.66,
)
