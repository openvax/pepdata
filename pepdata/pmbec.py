# Copyright (c) 2014-2016. Mount Sinai School of Medicine
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

from os.path import join

from .static_data import MATRIX_DIR

from .amino_acid_alphabet import dict_to_amino_acid_matrix

def read_pmbec_coefficients(
        key_type='row',
        verbose=True,
        filename=join(MATRIX_DIR, 'pmbec.mat')):
    """
    Parameters
    ------------

    filename : str
        Location of PMBEC coefficient matrix

    key_type : str
        'row' : every key is a single amino acid,
           which maps to a dictionary for that row
        'pair' : every key is a tuple of amino acids
        'pair_string' : every key is a string of two amino acid characters

    verbose : bool
        Print rows of matrix as we read them
    """
    d = {}
    if key_type == 'row':
        def add_pair(row_letter, col_letter, value):
            if row_letter not in d:
                d[row_letter] = {}
            d[row_letter][col_letter] = value
    elif key_type == 'pair':
        def add_pair(row_letter, col_letter, value):
            d[(row_letter, col_letter)] = value

    else:
        assert key_type == 'pair_string', \
            "Invalid dictionary key type: %s" % key_type

        def add_pair(row_letter, col_letter, value):
            d["%s%s" % (row_letter, col_letter)] = value

    with open(filename, 'r') as f:
        lines = [line for line in f.read().split('\n') if len(line) > 0]
        header = lines[0]
        if verbose:
            print(header)
        residues = [
            x for x in header.split()
            if len(x) == 1 and x != ' ' and x != '\t'
        ]
        assert len(residues) == 20
        if verbose:
            print(residues)
        for line in lines[1:]:
            cols = [
                x
                for x in line.split(' ')
                if len(x) > 0 and x != ' ' and x != '\t'
            ]
            assert len(cols) == 21, "Expected 20 values + letter, got %s" % cols
            row_letter = cols[0]
            for i, col in enumerate(cols[1:]):
                col_letter = residues[i]
                assert col_letter != ' ' and col_letter != '\t'
                value = float(col)
                add_pair(row_letter, col_letter, value)
    return d

# dictionary of PMBEC coefficient accessed like pmbec_dict["V"]["R"]
pmbec_dict = read_pmbec_coefficients(key_type="row")
pmbec_matrix = dict_to_amino_acid_matrix(pmbec_dict)
