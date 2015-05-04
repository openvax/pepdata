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
Quantify amino acids by their physical/chemical properties
"""

from __future__ import print_function, division, absolute_import
from os.path import join

import numpy as np

from .static_data import MATRIX_DIR


long_amino_acid_names = [
    "Alanine",
    "Arginine",
    "Asparagine",
    "Aspartic Acid",
    "Cysteine",
    "Glutamic Acid",
    "Glutamine",
    "Glycine",
    "Histidine",
    "Isoleucine",
    "Leucine",
    "Lysine",
    "Methionine",
    "Phenylalanine",
    "Proline",
    "Serine",
    "Threonine",
    "Tryptophan",
    "Tyrosine",
    "Valine",
]

short_amino_acid_names = [
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Glu",
    "Gln",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val"
]

amino_acid_letters = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V"
]


amino_acid_letter_pairs = [
  "%s%s" % (x, y) for y in amino_acid_letters for x in amino_acid_letters
]

amino_acid_pair_positions = {
    pair: i for (i, pair) in enumerate(amino_acid_letter_pairs)
}

def index_to_long_name(idx):
    return long_amino_acid_names[idx]

def index_to_short_name(idx):
    return short_amino_acid_names[idx]

def index_to_letter(idx):
    return amino_acid_letters[idx]

def letter_to_index(x):
    """
    Convert from an amino acid's letter code to its position index
    """
    # assert len(x) == 1
    x = x.upper()
    assert x in amino_acid_letters, x
    return amino_acid_letters.index(x)


def peptide_to_indices(xs):
    return [letter_to_index(x) for x in xs if x != 'X' and x != 'U']

def letter_to_short_name(x):
    return short_amino_acid_names[letter_to_index(x)]

def long_name_to_letter(name):
    assert name in long_amino_acid_names, "%s not found" % name
    idx = long_amino_acid_names.index(x)
    return amino_acid_letters[idx]


def peptide_toshort_amino_acid_names(xs):
    return [letter_to_short_name(x) for x in xs]

def short_name_to_index(x):
    return short_amino_acid_names.index(x.capitalize())

def short_name_to_letter(x):
    return amino_acid_letters[short_name_to_index(x)]

def parse_table(table_string):
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

def aa_value_dict_to_positional_list(value_dict):
    value_list = [None] * 20
    for letter, value in value_dict.items():
        idx = letter_to_index(letter)
        assert idx >= 0
        assert idx < 20
        value_list[idx] = value
    assert all(elt is not None for elt in value_list), "Missing amino acids in:\n%s" % value_dict.keys()
    return value_list

def get_idx(x):
    if isinstance(x, int):
        return x
    assert isinstance(x, str) and len(x) in (1, 3), "Unexpected %s" % x
    x = x.upper()
    if len(x) == 3:
        return short_name_to_index(x)
    else:
        return letter_to_index(x)

class SequenceTransformer(object):
    def __init__(self, table):
        self.table = table
        self.value_dict = parse_table(self.table)
        self.value_list = aa_value_dict_to_positional_list(self.value_dict)

    def __call__(self, letter):
        return self.value_dict[letter]

    def __getitem__(self, letter):
        return self.value_dict[letter]

    def transform_string(self, letters):
        return np.array([self.value_dict[amino_acid] for amino_acid in letters])

    def transform_strings(self, strings):
        d = self.value_dict
        return np.array([[d[x] for x in strings] for x in strings])

def transformation_from_table(table):
    return SequenceTransformer(table)

"""
Amino acids property tables copied from CRASP website
"""

hydropathy = transformation_from_table("""
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

volume = transformation_from_table("""
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

polarity = transformation_from_table("""
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

pK_side_chain = transformation_from_table("""
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

prct_exposed_residues = transformation_from_table("""
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

hydrophilicity = transformation_from_table("""
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

accessible_surface_area = transformation_from_table("""
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

local_flexibility = transformation_from_table("""
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

accessible_surface_area_folded = transformation_from_table("""
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

refractivity = transformation_from_table("""
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


def parse_interaction_table(table):
    table = table.strip()
    while "  " in table:
        table = table.replace("  ", " ")

    lines = [l.strip() for l in table.split("\n")]
    lines = [l for l in lines if len(l) > 0 and not l.startswith("#")]
    assert len(lines) == 20, "Malformed amino acid interaction table"
    d = {}
    for i, line in enumerate(lines):
        coeff_strings = line.split(" ")
        assert len(coeff_strings) == 20, \
          "Malformed row in amino acid interaction table"
        x = amino_acid_letters[i]
        d[x] = {}
        for j, coeff_str in enumerate(coeff_strings):
            value = float(coeff_str)
            y = amino_acid_letters[j]
            d[x][y] = value
    return d

def transpose_interaction_dict(d):
    transposed = {}
    for x in amino_acid_letters:
        transposed[x] = {}
        for y in amino_acid_letters:
            transposed[x][y] = d[y][x]
    return transposed


with open(join(MATRIX_DIR, 'strand_vs_coil.txt'), 'r') as f:
    strand_vs_coil_table = f.read()
    strand_vs_coil = parse_interaction_table(strand_vs_coil_table)
    coil_vs_strand = transpose_interaction_dict(strand_vs_coil)

with open(join(MATRIX_DIR, 'helix_vs_strand.txt'), 'r') as f:
    helix_vs_strand_table = f.read()
    helix_vs_strand = parse_interaction_table(helix_vs_strand_table)
    strand_vs_helix = transpose_interaction_dict(helix_vs_strand)

with open(join(MATRIX_DIR, 'helix_vs_coil.txt'), 'r') as f:
    helix_vs_coil_table = f.read()
    helix_vs_coil = parse_interaction_table(helix_vs_coil_table)
    coil_vs_helix = transpose_interaction_dict(helix_vs_coil)


def parse_blosum_table(table, coeff_type=int, key_type='row'):
    """
    Parse a table of pairwise amino acid coefficient (e.g. BLOSUM50)
    """

    lines = table.split("\n")
    # drop comments
    lines = [line for line in lines if not line.startswith("#")]
    # drop CR endline characters
    lines = [line.replace("\r", "") for line in lines]
    # skip empty lines
    lines = [line for line in lines if line]

    labels = lines[0].split()

    assert len(labels) >= 20, \
        "Expected 20+ amino acids but first line '%s' has %d fields" % (
          lines[0],
          len(labels)
        )
    coeffs = {}
    for line in lines[1:]:

        fields = line.split()
        assert len(fields) >= 21, \
            "Expected AA and 20+ coefficients but '%s' has %d fields" % (
                line, len(fields))
        x = fields[0]
        for i, coeff_str in enumerate(fields[1:]):
            y = labels[i]
            coeff = coeff_type(coeff_str)
            if key_type == 'pair':
                coeffs[(x, y)] = coeff
            elif key_type == 'pair_string':
                coeffs[x + y] = coeff
            else:
                assert key_type == 'row', "Unknown key type: %s" % key_type
                if x not in coeffs:
                    coeffs[x] = {}
                coeffs[x][y] = coeff
    return coeffs


with open(join(MATRIX_DIR, 'BLOSUM30'), 'r') as f:
    blosum30 = parse_blosum_table(f.read())

with open(join(MATRIX_DIR, 'BLOSUM50'), 'r') as f:
    blosum50 = parse_blosum_table(f.read())

with open(join(MATRIX_DIR, 'BLOSUM62'), 'r') as f:
    blosum62 = parse_blosum_table(f.read())

