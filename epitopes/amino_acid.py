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

import numpy as np

"""
Quantify amino acids by their physical/chemical properties
"""

_long_names = [
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

_short_names = [
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

_letters = [
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

def index_to_long_name(idx):
  return _long_names[idx]

def index_to_short_name(idx):
  return _short_names[idx]

def index_to_letter(idx):
  return _letters[idx]

def letter_to_index(x):
  """
  Convert from an amino acid's letter code to its position index
  """
  # assert len(x) == 1
  x = x.upper()
  assert x in _letters, x
  return _letters.index(x)


def peptide_to_indices(xs):
  return [letter_to_index(x) for x in xs if x != 'X' and x != 'U'  ]

def letter_to_short_name(x):
  return _short_names[letter_to_index(x)]

def long_name_to_letter(name):
  assert name in _long_names, "%s not found" % name
  idx = _long_names.index(x)
  return _letters[idx]


def peptide_to_short_names(xs):
  return [letter_to_short_name(x) for x in xs]


def short_name_to_index(x):
  return _short_names.index(x.capitalize())

def short_name_to_letter(x):
  return _letters[short_name_to_index(x)]

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
    assert letter not in value_dict,  "Repeated amino acid " + line
    value_dict[letter] = float(value)
  return value_dict

def aa_value_dict_to_positional_list(value_dict):
  value_list = [None] * 20
  for letter, value in value_dict.iteritems():
    idx = letter_to_index(letter)
    assert idx >= 0
    assert idx < 20
    value_list[idx] = value
  assert all(elt is not None for elt in value_list), "Missing amino acids in:\n%s" % value_dict.keys()
  return value_list

def get_idx(x):
  if isinstance(x, int):
    return x
  assert isinstance(x, str) and len(x) in (1,3),  "Unexpected %s"  % x
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

  def transform_string(self, letters):
    return np.array([self.value_dict[amino_acid] for amino_acid in letters])

  def transform_strings(self, strings):
    d = self.value_dict
    return np.array([[d[x] for x in s] for x in strings])

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



#Chou-Fasman of structural properties from
#http://prowl.rockefeller.edu/aainfo/chou.htm
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

        assert name in _long_names, "Invalid amino acid name %s" % name
        idx = _long_names.index(name)
        letter = _letters[idx]
        alpha = int(fields[1])
        beta = int(fields[2])
        turn = int(fields[3])
        alpha_helix_score_dict[letter] = alpha
        beta_sheet_score_dict[letter] = beta
        turn_score_dict[letter]= turn

    assert len(alpha_helix_score_dict) == 20
    assert len(beta_sheet_score_dict) == 20
    assert len(turn_score_dict) == 20
    return alpha_helix_score_dict, beta_sheet_score_dict, turn_score_dict

alpha_helix_score_dict, beta_sheet_score_dict, turn_score_dict = \
    parse_chou_fasman(chou_fasman_table)


_interaction_letters = "ARNDCQEGHILKMFPSTWYV"
def parse_interaction_table(table):
  table = table.strip()
  while "  " in table:
    table = table.replace("  ", " ")

  lines = [l.strip() for l in table.split("\n")]
  lines = [l for l in lines if len(l) > 0]
  assert len(lines) == 20, "Malformed amino acid interaction table"
  d = {}
  for i, line in enumerate(lines):
    coeff_strings = line.split(" ")
    assert len(coeff_strings) == 20, \
      "Malformed row in amino acid interaction table"
    x = _interaction_letters[i]
    d[x] = {}
    for j, coeff_str in enumerate(coeff_strings):
      value = float(coeff_str)
      y = _interaction_letters[j]
      d[x][y] = value
  return d
def transpose_interaction_dict(d):
  transposed = {}
  for x in _interaction_letters:
    transposed[x] = {}
    for y in _interaction_letters:
      transposed[x][y] = d[y][x]
  return transposed


#H ZHAC000105
#D Environment-dependent residue contact energies (rows = strand, cols = coil)
#R PMID:10706611
#A Zhang, C. and Kim, S.H.
#T Environment-dependent residue contact energies for proteins
#J Proc. Natl. Acad. Sci. USA 97, 2550-2555 (2000)
#M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
strand_vs_coil_table = \
"""
  -0.57  0.47 0.30 0.62 -1.60 0.45 0.61 -0.24  0.07 -1.64 -1.63 0.62 -1.03 -1.55 -0.11 -0.10 -0.34 -1.44 -0.39 -1.55
   0.23  0.79 0.76 0.39 -0.41 0.92 0.76  0.52  0.51 -0.30  0.13 1.58  0.88 -0.07  0.60  0.65  0.37  0.14  0.32  0.17
  -0.28  0.74 0.57 0.87 -0.68 0.52 1.00 -0.07  0.32 -0.31 -0.08 0.87  0.29 -0.17  0.57  0.11  0.19  0.04  0.24 -0.23
   0.15 -0.25 0.46 0.69 -0.46 0.41 1.34  0.56 -0.51 -0.23  0.27 0.59  0.60 -0.38  1.02  0.08  0.05 -0.48  0.02  0.34
  -1.19 -0.46 0.21 0.51 -3.30 0.26 0.20 -1.03 -0.72 -1.55 -1.71 0.27 -1.24 -1.70 -0.50 -0.55 -0.97 -0.67 -1.26 -1.62
   0.63  1.18 0.92 1.37 -0.30 0.93 1.27  0.56  0.91 -0.28 -0.11 0.98  0.15 -0.30  0.64  0.88  0.68 -0.44  0.66  0.15
   0.97  0.89 1.37 1.89  0.30 1.25 2.34  0.98  0.58  0.20  0.50 0.67  1.23  0.58  1.26  0.95  1.06  0.04  0.87  0.48
  -0.64  0.12 0.27 0.31 -1.37 0.38 0.98 -0.40 -0.12 -1.58 -1.40 0.78 -0.46 -1.38 -0.21  0.05 -0.26 -1.41 -0.61 -1.13
  -0.02  0.75 0.68 0.14 -0.58 0.73 0.84  0.41 -0.64 -0.75  0.03 1.46 -0.16 -0.49  0.52  0.31 -0.11 -1.00 -0.58  0.03
  -0.94 -0.14 0.31 0.26 -1.70 0.07 0.46 -0.37 -0.50 -1.88 -1.79 0.84 -0.99 -1.82 -0.47 -0.05 -0.54 -1.65 -1.09 -1.64
  -0.76  0.32 0.43 0.25 -1.63 0.22 0.68 -0.17 -0.40 -1.84 -1.70 0.47 -1.06 -1.76 -0.39  0.09 -0.42 -1.81 -1.15 -1.64
   1.02  1.99 1.18 0.59  0.08 1.10 0.60  0.61  0.95  0.24  0.34 2.69  0.97 -0.03  1.23  1.07  0.83  0.00  0.26  0.36
  -0.16  0.83 0.47 0.92 -1.63 0.36 0.71 -0.20  0.90 -1.00 -1.12 1.55 -0.31 -1.35 -0.01  0.34  0.20 -1.70 -0.60 -0.79
  -0.70  0.03 0.63 0.15 -1.26 0.29 0.35 -0.11 -0.36 -1.73 -1.55 0.71 -0.97 -1.55 -0.28 -0.09 -0.32 -1.23 -0.91 -1.30
   0.17  0.50 0.60 0.67 -1.31 0.50 0.94  0.02 -0.45 -1.26 -0.91 1.08  0.83 -0.87  0.63  0.31  0.26 -0.50 -0.55 -0.79
  -0.06  0.99 0.73 0.86 -0.89 0.85 0.67  0.08  0.06 -0.22 -0.29 0.94 -0.08 -0.41  0.67  0.33  0.13 -1.01  0.13 -0.24
   0.26  0.93 0.70 0.87 -0.78 0.58 1.20  0.12  0.52 -0.30 -0.24 1.11  0.01 -0.08  0.65  0.47  0.41 -0.31  0.12 -0.32
  -0.03 -0.11 0.27 0.66 -1.50 0.65 0.50 -0.12 -0.32 -1.13 -1.01 0.52 -1.08 -1.04 -0.32 -0.03 -0.10 -0.67 -0.73 -0.64
  -0.44  0.20 0.20 0.20 -1.26 0.16 0.10 -0.21 -0.52 -1.26 -1.30 0.60 -0.76 -1.17 -0.42  0.05 -0.27 -1.20 -0.75 -0.84
  -0.83  0.20 0.48 0.62 -1.44 0.17 0.73 -0.12 -0.26 -1.64 -1.59 0.52 -0.70 -1.55 -0.28  0.12 -0.17 -1.16 -0.85 -1.42
  """
strand_vs_coil_dict = parse_interaction_table(strand_vs_coil_table)
coil_vs_strand_dict = transpose_interaction_dict(strand_vs_coil_dict)
#H ZHAC000102
#D Environment-dependent residue contact energies (rows = helix, cols = strand)
#R PMID:10706611
#A Zhang, C. and Kim, S.H.
#T Environment-dependent residue contact energies for proteins
#J Proc. Natl. Acad. Sci. USA 97, 2550-2555 (2000)
#M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
helix_vs_strand_table = \
"""
  -0.94  1.26  0.55 0.76 -1.54  1.14 1.57 -0.78  0.44 -1.59 -1.64  1.91 -0.90 -1.49  0.28  0.20 -0.04  -0.92 -0.75 -1.45
   0.56  1.79  2.31 0.79 -0.67  2.54 0.72  1.09  0.94 -0.01  0.01  3.68  0.89 -0.05  1.37  0.83  1.35   0.00  0.33  0.44
   0.59  2.21  1.82 0.77 -0.90  0.46 3.06 -0.16  0.63 -0.33  0.20  2.43  0.99  0.63  0.54  0.24  0.63   0.11 -0.19  0.23
   0.66  0.76  0.76 1.19 -0.21  1.66 2.22  0.29  0.57  0.59  0.79  1.13  1.41  0.49  1.70  1.03  1.19   1.85  0.18  0.86
  -1.75  0.78 -1.00 0.32 -3.64  0.48 0.87 -1.67 -0.62 -2.77 -2.32  0.19 -1.22 -2.67 -1.62 -0.83 -1.14  -0.52 -1.94 -2.35
   0.33  2.15  1.22 1.26  1.37  1.17 2.56  0.92  1.02  0.11  0.00  2.58  0.79 -0.26  0.53  1.19  1.11   0.21  0.39  0.15
   0.82  1.05  2.18 2.11  0.01  2.42 2.58  1.15  0.97  0.20  0.31  1.31  1.25  0.12  2.00  1.09  1.13   0.58  0.31  0.39
  -0.40  0.95  0.03 0.14 -1.00  0.34 0.99 -1.32  0.13 -1.40 -1.36  1.58 -0.90 -1.41  0.82 -0.27  0.21  -0.59 -1.27 -1.09
  -0.75  2.19  0.13 0.68 -1.37  1.98 1.13  0.01  1.52 -0.83 -0.58  2.26 -0.82 -1.01  0.53 -0.17  0.02 -49.00 -0.61 -0.56
  -1.99  0.25 -0.20 1.00 -2.44 -0.12 0.88 -1.54 -0.05 -2.64 -2.33  0.75 -1.85 -2.46 -1.06 -0.59 -0.65  -1.82 -1.88 -2.45
  -2.02  0.34 -0.04 0.13 -2.29  0.24 0.73 -1.27 -0.46 -2.53 -2.44  0.67 -1.80 -2.28 -1.29 -0.40 -0.34  -1.76 -1.66 -2.26
   0.60  3.11  2.23 1.06  0.50  1.80 1.65  0.82  1.25  0.10  0.34  3.51  0.98 -0.21  1.15  2.09  1.30  -0.14  0.28  0.13
  -1.54 -0.06 -0.63 1.76 -2.51  0.14 0.72 -1.74  0.07 -2.27 -2.22  1.27 -1.77 -1.87  0.34 -0.02 -0.21  -0.93 -1.54 -1.81
  -2.12  0.33 -0.70 0.17 -2.30 -0.59 0.26 -1.60 -0.88 -2.53 -2.44 -0.42 -1.83 -2.68 -1.40 -0.82 -0.61  -1.63 -1.83 -2.25
   0.63  2.43 -0.19 1.31 -1.63  1.46 1.91  0.08  1.11 -0.20  0.47  1.94 -0.34  0.15  0.57  0.00  1.15   0.06  0.26 -0.06
  -0.41  0.88  1.02 1.04 -0.21  1.27 0.94  0.04  0.75 -0.48 -0.67  2.28  0.45 -0.92  0.75  0.50  0.96   0.22 -0.19 -0.54
  -0.32  1.48  0.35 0.43 -1.44  0.38 1.36 -0.38  0.20 -1.14 -1.00  1.38 -0.35 -0.97 -0.05 -0.16  0.29  -0.53 -0.76 -0.73
  -1.85  0.45 -0.03 0.80 -1.64 -0.23 0.11 -0.95  0.67 -1.58 -2.13  0.61 -1.75 -1.59 -1.07 -0.34 -0.40  -1.29 -1.27 -1.79
  -0.88 -0.20 -0.29 0.14 -1.31  0.09 0.71 -0.56 -0.57 -1.66 -1.38  1.40 -1.60 -1.97 -0.73 -0.32 -0.37  -1.40 -0.96 -1.38
  -1.74  0.85  0.24 0.72 -2.25  0.45 0.81 -1.29 -0.24 -2.46 -2.38  0.37 -1.21 -2.16 -1.00 -0.10 -0.57  -1.34 -1.52 -2.31
"""

helix_vs_strand_dict = parse_interaction_table(helix_vs_strand_table)
strand_vs_helix_dict = transpose_interaction_dict(helix_vs_strand_dict)

#H ZHAC000103
#D Environment-dependent residue contact energies (rows = helix, cols = coil)
#R PMID:10706611
#A Zhang, C. and Kim, S.H.
#T Environment-dependent residue contact energies for proteins
#J Proc. Natl. Acad. Sci. USA 97, 2550-2555 (2000)
#M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
helix_vs_coil_table = \
"""
   0.12  1.17  0.84 0.90 -0.81  1.16 1.44  0.10  0.69 -0.81 -0.78 1.16 -0.22 -0.67  0.61  0.47  0.36 -0.72 -0.37 -0.43
   0.98  1.65  1.16 0.60 -0.21  1.26 1.12  1.09  1.16 -0.04 -0.09 2.37  0.47 -0.04  1.22  1.05  0.92 -0.09  0.06  0.32
   0.69  1.16  1.16 1.22 -0.06  1.23 1.45  0.96  0.88  0.26  0.12 1.48  0.32  0.03  1.14  0.73  0.62  0.62  0.53  0.23
   0.90  0.40  1.06 1.45  0.58  1.88 2.18  1.13  0.69  0.43  0.65 0.95  0.75  0.33  1.41  0.39  0.54 -0.10  0.12  0.77
  -0.83  0.10  0.40 0.12 -2.65 -0.24 0.96 -0.26 -0.26 -1.61 -1.77 0.80 -1.02 -1.47 -0.31 -0.31 -0.49 -1.30 -0.98 -1.62
   1.13  1.10  1.28 1.37  0.14  1.62 1.84  1.29  1.31  0.05 -0.05 1.50  0.41  0.20  1.14  0.86  0.62  0.45  0.31  0.48
   1.33  0.91  1.33 1.60  0.31  1.60 1.93  1.62  1.01  0.33  0.38 1.12  0.82  0.55  1.54  0.78  0.54  0.23  0.52  0.86
  -0.22  0.72  0.27 0.47 -0.95  0.42 1.39 -0.23  0.40 -0.48 -0.81 1.04 -0.62 -0.36  0.41  0.23 -0.04 -0.71  0.08 -0.35
   0.47  0.81  0.95 0.51 -1.56  0.90 0.89  0.86  0.20 -0.43 -0.48 1.31 -0.63 -0.41  0.56  0.40  0.28 -0.20 -0.22 -0.21
  -0.58  0.17  0.61 0.46 -1.17  0.24 0.80  0.04 -0.16 -1.64 -1.66 0.87 -0.89 -1.56 -0.27  0.02 -0.32 -1.40 -1.13 -1.36
  -0.44  0.20  0.50 0.71 -1.56  0.11 0.82  0.28 -0.15 -1.67 -1.62 0.72 -0.96 -1.55  0.02  0.19 -0.09 -1.46 -0.95 -1.32
   1.07  2.48  1.75 0.98  0.42  1.68 1.04  1.31  1.39  0.41  0.29 2.95  0.98  0.27  1.63  1.51  1.48  0.32  0.60  0.64
  -0.22  0.65  0.76 0.88 -0.95  0.68 1.92  0.27  0.31 -1.32 -1.04 1.02 -0.57 -1.60  0.07  0.47  0.04 -1.29 -0.85 -0.82
  -0.33 -0.06  0.42 0.42 -1.90  0.25 0.64  0.12 -0.01 -1.64 -1.50 0.58 -1.36 -1.77 -0.30  0.02  0.04 -1.41 -1.36 -1.34
   0.78  1.30  1.31 1.27 -0.04  1.44 1.71  0.69  0.84  0.05  0.15 1.68  0.38  0.27  1.05  1.19  0.83 -0.24  0.23  0.12
   0.46  1.07  1.04 0.73 -0.31  1.47 1.23  0.57  0.58 -0.11 -0.24 1.37  0.08 -0.34  0.76  0.51  0.48 -0.04  0.47  0.18
   0.50  0.90  0.75 0.91 -0.26  1.03 1.25  0.55  0.55 -0.20 -0.26 1.42  0.50 -0.22  0.88  0.69  0.56  0.41  0.11 -0.15
  -0.41 -0.06 -0.19 0.32 -0.79 -0.14 0.58  0.07 -0.62 -1.58 -1.16 0.18 -1.03 -1.33 -0.56  0.15 -0.19 -1.83 -0.67 -0.92
  -0.22 -0.07  0.52 0.46 -0.87  0.38 0.59  0.40 -0.17 -1.29 -1.15 0.83 -0.98 -1.16 -0.16  0.34 -0.12 -0.79 -0.77 -0.78
  -0.51  0.49  0.48 0.67 -1.40  0.66 0.63 -0.06  0.28 -1.25 -1.50 1.14 -0.93 -1.36 -0.04  0.10 -0.01 -1.11 -0.82 -1.14
"""
helix_vs_coil_dict = parse_interaction_table(helix_vs_coil_table)
coil_vs_helix_dict = transpose_interaction_dict(helix_vs_coil_dict)