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

from os.path import join

from .amino_acid_alphabet import canonical_amino_acid_letters, dict_to_amino_acid_matrix
from .static_data import MATRIX_DIR


def parse_interaction_table(table, amino_acid_order="ARNDCQEGHILKMFPSTWYV"):
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
        x = amino_acid_order[i]
        d[x] = {}
        for j, coeff_str in enumerate(coeff_strings):
            value = float(coeff_str)
            y = amino_acid_order[j]
            d[x][y] = value
    return d

def transpose_interaction_dict(d):
    transposed = {}
    for x in canonical_amino_acid_letters:
        transposed[x] = {}
        for y in canonical_amino_acid_letters:
            transposed[x][y] = d[y][x]
    return transposed


with open(join(MATRIX_DIR, 'strand_vs_coil.txt'), 'r') as f:
    # Strand vs. Coil
    strand_vs_coil_dict = parse_interaction_table(f.read())
    strand_vs_coil_array = dict_to_amino_acid_matrix(strand_vs_coil_dict)

    # Coil vs. Strand
    coil_vs_strand_dict = transpose_interaction_dict(strand_vs_coil_dict)
    coil_vs_strand_array = dict_to_amino_acid_matrix(coil_vs_strand_dict)

with open(join(MATRIX_DIR, 'helix_vs_strand.txt'), 'r') as f:
    # Helix vs. Strand
    helix_vs_strand_dict = parse_interaction_table(f.read())
    helix_vs_strand_array = dict_to_amino_acid_matrix(helix_vs_strand_dict)

    # Strand vs. Helix
    strand_vs_helix_dict = transpose_interaction_dict(helix_vs_strand_dict)
    strand_vs_helix_array = dict_to_amino_acid_matrix(strand_vs_helix_dict)

with open(join(MATRIX_DIR, 'helix_vs_coil.txt'), 'r') as f:
    # Helix vs. Coil
    helix_vs_coil_dict = parse_interaction_table(f.read())
    helix_vs_coil_array = dict_to_amino_acid_matrix(helix_vs_coil_dict)

    # Coil vs. Helix
    coil_vs_helix_dict = transpose_interaction_dict(helix_vs_coil_dict)
    coil_vs_helix_array = dict_to_amino_acid_matrix(coil_vs_helix_dict)
