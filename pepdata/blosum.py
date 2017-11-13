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

from os.path import join

from .static_data import MATRIX_DIR

from .amino_acid_alphabet import dict_to_amino_acid_matrix

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

    if len(labels) < 20:
        raise ValueError(
            "Expected 20+ amino acids but first line '%s' has %d fields" % (
                lines[0],
                len(labels)))
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
    blosum30_dict = parse_blosum_table(f.read())
    blosum30_matrix = dict_to_amino_acid_matrix(blosum30_dict)

with open(join(MATRIX_DIR, 'BLOSUM50'), 'r') as f:
    blosum50_dict = parse_blosum_table(f.read())
    blosum50_matrix = dict_to_amino_acid_matrix(blosum50_dict)

with open(join(MATRIX_DIR, 'BLOSUM62'), 'r') as f:
    blosum62_dict = parse_blosum_table(f.read())
    blosum62_matrix = dict_to_amino_acid_matrix(blosum62_dict)

