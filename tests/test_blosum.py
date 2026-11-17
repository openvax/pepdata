from pathlib import Path

import pytest

from pepdata.amino_acid_alphabet import canonical_amino_acid_letters
from pepdata.blosum import (
    blosum30_dict,
    blosum30_matrix,
    blosum50_dict,
    blosum50_matrix,
    blosum62_dict,
    blosum62_matrix,
)

BLOSUM_TABLES = {
    30: (blosum30_dict, blosum30_matrix),
    50: (blosum50_dict, blosum50_matrix),
    62: (blosum62_dict, blosum62_matrix),
}

# Published self-substitution scores for alanine, tryptophan and cysteine.
DIAGONAL_SCORES = {
    30: {"A": 4, "W": 20, "C": 17},
    50: {"A": 5, "W": 15, "C": 13},
    62: {"A": 4, "W": 11, "C": 9},
}


@pytest.mark.parametrize("number", sorted(BLOSUM_TABLES))
def test_blosum_covers_canonical_residues(number):
    table, _ = BLOSUM_TABLES[number]
    for row in canonical_amino_acid_letters:
        assert set(canonical_amino_acid_letters) <= set(table[row])


@pytest.mark.parametrize("number", sorted(BLOSUM_TABLES))
def test_blosum_diagonal_scores(number):
    table, _ = BLOSUM_TABLES[number]
    for residue, score in DIAGONAL_SCORES[number].items():
        assert table[residue][residue] == score
    for residue in canonical_amino_acid_letters:
        row = [table[residue][other] for other in canonical_amino_acid_letters]
        assert table[residue][residue] == max(row)


@pytest.mark.parametrize("number", sorted(BLOSUM_TABLES))
def test_blosum_is_symmetric(number):
    table, _ = BLOSUM_TABLES[number]
    for x in canonical_amino_acid_letters:
        for y in canonical_amino_acid_letters:
            assert table[x][y] == table[y][x]


def test_blosum50_matches_ncbi_canonical_scores():
    reference = Path(__file__).parent / "data/ncbi_BLOSUM50.txt"
    rows = [line.split() for line in reference.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    header = rows[0]
    expected = {row[0]: dict(zip(header, map(int, row[1:]))) for row in rows[1:]}
    for x in canonical_amino_acid_letters:
        for y in canonical_amino_acid_letters:
            assert blosum50_dict[x][y] == expected[x][y], (x, y)


@pytest.mark.parametrize("number", sorted(BLOSUM_TABLES))
def test_blosum_matrix_matches_dict_in_canonical_order(number):
    table, matrix = BLOSUM_TABLES[number]
    for i, x in enumerate(canonical_amino_acid_letters):
        for j, y in enumerate(canonical_amino_acid_letters):
            assert matrix[i, j] == table[x][y]


@pytest.mark.parametrize("number", sorted(BLOSUM_TABLES))
def test_blosum_matrix_has_only_the_requested_alphabet(number):
    _, matrix = BLOSUM_TABLES[number]
    assert matrix.shape == (len(canonical_amino_acid_letters),) * 2


def test_dict_to_matrix_respects_a_subset_and_its_order():
    from pepdata.amino_acid import AminoAcid
    from pepdata.amino_acid_alphabet import dict_to_amino_acid_matrix

    alphabet = [AminoAcid("Cysteine", "Cys", "C"), AminoAcid("Alanine", "Ala", "A")]
    matrix = dict_to_amino_acid_matrix(blosum62_dict, alphabet=alphabet)
    assert matrix.tolist() == [[9, 0], [0, 4]]
