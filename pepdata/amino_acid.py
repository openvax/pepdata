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

class AminoAcid(object):
    def __init__(self, full_name, short_name, letter, contains=None):
        self.letter = letter
        self.full_name = full_name
        self.short_name = short_name
        if not contains:
            contains = [letter]
        self.contains = contains

    def __str__(self):
        return (
            ("AminoAcid(full_name='%s', short_name='%s', letter='%s', "
             "contains=%s)") % (
            self.letter, self.full_name, self.short_name, self.contains))

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return other.__class__ is AminoAcid and self.letter == other.letter

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
# Amino acid tokens which represent multiple canonical amino acids
###
wildcard_amino_acids = [
    AminoAcid("Unknown", "Xaa", "X", contains=set(canonical_amino_acid_letters)),
    AminoAcid("Asparagine or aspartic acid", "Asx",  "B", contains={"D", "N"}),
    AminoAcid("Glutamine or glutamic acid", "Glx", "Z", contains={"E", "Q"}),
    AminoAcid("Leucine or Isoleucine", "Xle", "J", contains={"I", "L"})
]

###
# Canonical amino acids + wilcard tokens
###
canonical_amino_acids_with_unknown = canonical_amino_acids + wildcard_amino_acids


###
# Rare amino acids which aren't considered part of the core 20 "canonical"
###

extended_amino_acids = canonical_amino_acids + [
    AminoAcid("Selenocysteine", "Sec", "U"),
    AminoAcid("Pyrrolysine", "Pyl", "O"),
]

extended_amino_acid_letters = [aa.letter for aa in extended_amino_acids]

###
# Extended amino acids + wildcard tokens
###
extended_amino_acids_with_unknown = extended_amino_acids + wildcard_amino_acids
extended_amino_acids_with_unknown_letters = [
    aa.letter for aa in extended_amino_acids_with_unknown
]


amino_acid_letter_indices = {
    c: i for (i, c) in
    enumerate(extended_amino_acids_with_unknown_letters)
}

amino_acid_letter_pairs = [
    "%s%s" % (x, y)
    for y in extended_amino_acids_with_unknown_letters
    for x in extended_amino_acids_with_unknown_letters
]

amino_acid_pair_positions = {
    pair: i for (i, pair) in enumerate(amino_acid_letter_pairs)
}

def index_to_full_name(idx):
    return extended_amino_acids_with_unknown[idx].full_name

def index_to_short_name(idx):
    return extended_amino_acids_with_unknown[idx].short_name

def index_to_letter(idx):
    return extended_amino_acids_with_unknown_letters[idx]

def letter_to_index(x):
    """
    Convert from an amino acid's letter code to its position index
    """
    # assert len(x) == 1
    x = x.upper()
    assert x in amino_acid_letter_indices, "Unknown amino acid: %s" % x
    return amino_acid_letter_indices[x]

def peptide_to_indices(xs):
    xs = xs.upper()
    return [amino_acid_letter_indices[x] for x in xs]

def letter_to_short_name(x):
    return index_to_short_name(letter_to_index(x))

def peptide_to_short_amino_acid_names(xs):
    xs = xs.upper()
    return [amino_acid_letter_indices[x] for x in xs]

def aa_value_dict_to_positional_list(value_dict):
    value_list = [None] * 20
    for letter, value in value_dict.items():
        idx = letter_to_index(letter)
        assert idx >= 0
        assert idx < 20
        value_list[idx] = value
    assert all(elt is not None for elt in value_list), \
        "Missing amino acids in:\n%s" % value_dict.keys()
    return value_list
