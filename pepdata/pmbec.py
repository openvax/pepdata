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

import warnings
from os.path import join

from .static_data import MATRIX_DIR

from .amino_acid_alphabet import dict_to_amino_acid_matrix
from .pairwise_table import read_pairwise_table, rekey_pairwise_table

PMBEC_PATH = join(MATRIX_DIR, 'pmbec.mat')


def read_pmbec_coefficients(
        key_type='row',
        verbose=False,
        filename=PMBEC_PATH):
    """
    Read the PMBEC coefficient matrix.

    Deprecated since 1.3.8: use
    `pepdata.pairwise_table.read_pairwise_table`. This wrapper will be
    removed in a future release.

    Parameters
    ----------
    key_type : str
        'row' : every key is a single amino acid,
           which maps to a dictionary for that row
        'pair' : every key is a tuple of amino acids
        'pair_string' : every key is a string of two amino acid characters

    verbose : bool
        Ignored. Earlier releases printed the header and residue list to
        stdout by default; nothing is printed any more.

    filename : str
        Location of PMBEC coefficient matrix

    Returns
    -------
    dict
        Coefficients in the requested layout.

    Raises
    ------
    ValueError
        If the table is malformed, does not have exactly 20 columns, or
        `key_type` is unknown.
    """
    warnings.warn(
        "read_pmbec_coefficients is deprecated; use "
        "pepdata.pairwise_table.read_pairwise_table instead",
        DeprecationWarning,
        stacklevel=2)
    coeffs = read_pairwise_table(filename)
    n_labels = len(next(iter(coeffs.values())))
    if n_labels != 20:
        raise ValueError(
            "Expected 20 amino acids but header has %d fields" % n_labels)
    return rekey_pairwise_table(coeffs, key_type=key_type)


# dictionary of PMBEC coefficient accessed like pmbec_dict["V"]["R"]
pmbec_dict = read_pairwise_table(PMBEC_PATH)
pmbec_matrix = dict_to_amino_acid_matrix(pmbec_dict)
