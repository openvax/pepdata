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

from .amino_acid_alphabet import letter_to_index

"""
Quantify amino acids by their physical/chemical properties
"""


def aa_dict_to_positional_list(aa_property_dict):
    value_list = [None] * 20
    for letter, value in aa_property_dict.items():
        idx = letter_to_index(letter)
        assert idx >= 0
        assert idx < 20
        value_list[idx] = value
    assert all(elt is not None for elt in value_list), \
        "Missing amino acids in:\n%s" % aa_property_dict.keys()
    return value_list

def parse_property_table(table_string):
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


"""
Amino acids property tables copied from CRASP website
"""

hydropathy = parse_property_table("""
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

volume = parse_property_table("""
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

polarity = parse_property_table("""
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

pK_side_chain = parse_property_table("""
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

prct_exposed_residues = parse_property_table("""
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

hydrophilicity = parse_property_table("""
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

accessible_surface_area = parse_property_table("""
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

local_flexibility = parse_property_table("""
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

accessible_surface_area_folded = parse_property_table("""
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

refractivity = parse_property_table("""
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


mass = parse_property_table("""
70.079 A ALA
156.188 R ARG
114.104 N ASN
115.089 D ASP
103.144 C CYS
128.131 Q GLN
129.116 E GLU
57.052 G GLY
137.142 H HIS
113.160 I ILE
113.160 L LEU
128.174 K LYS
131.198 M MET
147.177 F PHE
97.177 P PRO
87.078 S SER
101.105 T THR
186.213 W TRP
163.170 Y TYR
99.133 V VAL
""")

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
