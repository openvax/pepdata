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

import numpy as np

from .amino_acid import AminoAcid

canonical_amino_acids = [
    AminoAcid("Alanine", "Ala", "A"),
    AminoAcid("Arginine", "Arg", "R"),
    AminoAcid("Asparagine","Asn", "N"),
    AminoAcid("Aspartic Acid", "Asp", "D"),
    AminoAcid("Cysteine", "Cys", "C"),
    AminoAcid("Glutamic Acid", "Glu", "E"),
    AminoAcid("Glutamine", "Gln", "Q"),
    AminoAcid("Glycine", "Gly", "G"),
    AminoAcid("Histidine", "His", "H"),
    AminoAcid("Isoleucine",  "Ile", "I"),
    AminoAcid("Leucine", "Leu", "L"),
    AminoAcid("Lysine", "Lys", "K"),
    AminoAcid("Methionine",  "Met", "M"),
    AminoAcid("Phenylalanine", "Phe", "F"),
    AminoAcid("Proline", "Pro", "P"),
    AminoAcid("Serine", "Ser", "S"),
    AminoAcid("Threonine", "Thr", "T"),
    AminoAcid("Tryptophan", "Trp", "W"),
    AminoAcid("Tyrosine", "Tyr", "Y"),
    AminoAcid("Valine", "Val", "V")
]

canonical_amino_acid_letters = [aa.letter for aa in canonical_amino_acids]

###
# Post-translation modifications commonly detected by mass-spec
###

# TODO: figure out three letter codes for modified AAs

modified_amino_acids = [
    AminoAcid("Phospho-Serine", "Sep", "s"),
    AminoAcid("Phospho-Threonine", "???", "t"),
    AminoAcid("Phospho-Tyrosine", "???", "y"),
    AminoAcid("Cystine", "???", "c"),
    AminoAcid("Methionine sulfoxide", "???", "m"),
    AminoAcid("Pyroglutamate", "???", "q"),
    AminoAcid("Pyroglutamic acid", "???", "n"),
]

###
# Amino acid tokens which represent multiple canonical amino acids
###
wildcard_amino_acids = [
    AminoAcid("Unknown", "Xaa", "X", contains=set(canonical_amino_acid_letters)),
    AminoAcid("Asparagine-or-Aspartic-Acid", "Asx",  "B", contains={"D", "N"}),
    AminoAcid("Glutamine-or-Glutamic-Acid", "Glx", "Z", contains={"E", "Q"}),
    AminoAcid("Leucine-or-Isoleucine", "Xle", "J", contains={"I", "L"})
]

###
# Canonical amino acids + wilcard tokens
###

canonical_amino_acids_with_unknown = canonical_amino_acids + wildcard_amino_acids


###
# Rare amino acids which aren't considered part of the core 20 "canonical"
###

rare_amino_acids = [
    AminoAcid("Selenocysteine", "Sec", "U"),
    AminoAcid("Pyrrolysine", "Pyl", "O"),
]

###
# Extended amino acids + wildcard tokens
###

extended_amino_acids = canonical_amino_acids + rare_amino_acids + wildcard_amino_acids
extended_amino_acid_letters = [
    aa.letter for aa in extended_amino_acids
]
extended_amino_acids_with_unknown_names = [
    aa.full_name for aa in extended_amino_acids
]


amino_acid_letter_indices = {
    c: i for (i, c) in
    enumerate(extended_amino_acid_letters)
}


amino_acid_letter_pairs = [
    "%s%s" % (x, y)
    for y in extended_amino_acids
    for x in extended_amino_acids
]


amino_acid_name_indices = {
    aa_name: i for (i, aa_name)
    in enumerate(extended_amino_acids_with_unknown_names)
}

amino_acid_pair_positions = {
    pair: i for (i, pair) in enumerate(amino_acid_letter_pairs)
}

def index_to_full_name(idx):
    return extended_amino_acids[idx].full_name

def index_to_short_name(idx):
    return extended_amino_acids[idx].short_name

def index_to_letter(idx):
    return extended_amino_acids[idx]

def letter_to_index(x):
    """
    Convert from an amino acid's letter code to its position index
    """
    assert x in amino_acid_letter_indices, "Unknown amino acid: %s" % x
    return amino_acid_letter_indices[x]

def peptide_to_indices(xs):
    return [amino_acid_letter_indices[x] for x in xs]

def letter_to_short_name(x):
    return index_to_short_name(letter_to_index(x))

def peptide_to_short_amino_acid_names(xs):
    return [amino_acid_letter_indices[x] for x in xs]

def dict_to_amino_acid_matrix(d, alphabet=canonical_amino_acids):
    n_aa = len(d)
    result_matrix = np.zeros((n_aa, n_aa), dtype="float32")
    for i, aa_row in enumerate(alphabet):
        d_row = d[aa_row.letter]
        for j, aa_col in enumerate(alphabet):
            value = d_row[aa_col.letter]
            result_matrix[i, j] = value
    return result_matrix

