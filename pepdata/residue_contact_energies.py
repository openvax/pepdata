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

import warnings
from os.path import join

from .amino_acid_alphabet import canonical_amino_acid_letters, dict_to_amino_acid_matrix
from .pairwise_table import parse_pairwise_table, read_pairwise_table
from .static_data import MATRIX_DIR

# Row and column order shared by the three shipped contact energy tables
CONTACT_TABLE_AMINO_ACID_ORDER = "ARNDCQEGHILKMFPSTWYV"


def parse_interaction_table(table, amino_acid_order=CONTACT_TABLE_AMINO_ACID_ORDER):
    """
    Parse a headerless 20x20 table of residue contact energies.

    Deprecated since 1.3.8: use
    `pepdata.pairwise_table.parse_pairwise_table` with its `column_order`
    argument. This wrapper will be removed in a future release.

    Parameters
    ----------
    table : str
        Contents of a table with one row of coefficients per residue.

    amino_acid_order : str
        Residue letters naming both the rows and the columns.

    Returns
    -------
    dict of str to dict of str to float
        Coefficients accessed like `d[row_letter][column_letter]`.

    Raises
    ------
    ValueError
        If the table does not have one row per residue or a row has the
        wrong number of coefficients.
    """
    warnings.warn(
        "parse_interaction_table is deprecated; use "
        "pepdata.pairwise_table.parse_pairwise_table instead",
        DeprecationWarning,
        stacklevel=2)
    return parse_pairwise_table(table, column_order=amino_acid_order)


def transpose_interaction_dict(d):
    transposed = {}
    for x in canonical_amino_acid_letters:
        transposed[x] = {}
        for y in canonical_amino_acid_letters:
            transposed[x][y] = d[y][x]
    return transposed


# Strand vs. Coil
strand_vs_coil_dict = read_pairwise_table(
    join(MATRIX_DIR, 'strand_vs_coil.txt'),
    column_order=CONTACT_TABLE_AMINO_ACID_ORDER)
strand_vs_coil_array = dict_to_amino_acid_matrix(strand_vs_coil_dict)

# Coil vs. Strand
coil_vs_strand_dict = transpose_interaction_dict(strand_vs_coil_dict)
coil_vs_strand_array = dict_to_amino_acid_matrix(coil_vs_strand_dict)

# Helix vs. Strand
helix_vs_strand_dict = read_pairwise_table(
    join(MATRIX_DIR, 'helix_vs_strand.txt'),
    column_order=CONTACT_TABLE_AMINO_ACID_ORDER)
helix_vs_strand_array = dict_to_amino_acid_matrix(helix_vs_strand_dict)

# Strand vs. Helix
strand_vs_helix_dict = transpose_interaction_dict(helix_vs_strand_dict)
strand_vs_helix_array = dict_to_amino_acid_matrix(strand_vs_helix_dict)

# Helix vs. Coil
helix_vs_coil_dict = read_pairwise_table(
    join(MATRIX_DIR, 'helix_vs_coil.txt'),
    column_order=CONTACT_TABLE_AMINO_ACID_ORDER)
helix_vs_coil_array = dict_to_amino_acid_matrix(helix_vs_coil_dict)

# Coil vs. Helix
coil_vs_helix_dict = transpose_interaction_dict(helix_vs_coil_dict)
coil_vs_helix_array = dict_to_amino_acid_matrix(coil_vs_helix_dict)
