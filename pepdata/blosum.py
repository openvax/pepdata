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

from .static_data import MATRIX_DIR

from .amino_acid_alphabet import dict_to_amino_acid_matrix
from .pairwise_table import (
    parse_pairwise_table,
    read_pairwise_table,
    rekey_pairwise_table,
)


def parse_blosum_table(table, coeff_type=int, key_type='row'):
    """
    Parse a table of pairwise amino acid coefficient (e.g. BLOSUM50)

    Deprecated since 1.3.8: use
    `pepdata.pairwise_table.parse_pairwise_table` or
    `pepdata.pairwise_table.read_pairwise_table`. This wrapper will be
    removed in a future release.

    Parameters
    ----------
    table : str
        Contents of a table with a header line of column labels.

    coeff_type : callable, default int
        Converts each coefficient string to a number.

    key_type : {'row', 'pair', 'pair_string'}
        Layout of the returned dictionary; see
        `pepdata.pairwise_table.rekey_pairwise_table`.

    Returns
    -------
    dict
        Coefficients in the requested layout.

    Raises
    ------
    ValueError
        If the table is malformed, has fewer than 20 column labels, or
        `key_type` is unknown.
    """
    warnings.warn(
        "parse_blosum_table is deprecated; use "
        "pepdata.pairwise_table.parse_pairwise_table instead",
        DeprecationWarning,
        stacklevel=2)
    coeffs = parse_pairwise_table(table, coeff_type=coeff_type)
    n_labels = len(next(iter(coeffs.values())))
    if n_labels < 20:
        raise ValueError(
            "Expected 20+ amino acids but header has %d fields" % n_labels)
    return rekey_pairwise_table(coeffs, key_type=key_type)


blosum30_dict = read_pairwise_table(join(MATRIX_DIR, 'BLOSUM30'), coeff_type=int)
blosum30_matrix = dict_to_amino_acid_matrix(blosum30_dict)

blosum50_dict = read_pairwise_table(join(MATRIX_DIR, 'BLOSUM50'), coeff_type=int)
blosum50_matrix = dict_to_amino_acid_matrix(blosum50_dict)

blosum62_dict = read_pairwise_table(join(MATRIX_DIR, 'BLOSUM62'), coeff_type=int)
blosum62_matrix = dict_to_amino_acid_matrix(blosum62_dict)
