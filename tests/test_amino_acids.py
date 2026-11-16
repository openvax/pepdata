from pepdata.amino_acid_alphabet import (
    canonical_amino_acids,
    canonical_amino_acid_letters,
    extended_amino_acids,
    extended_amino_acid_letters,
)

def test_canonical_amino_acids():
    assert len(canonical_amino_acids) == 20

def test_canonical_amino_acids_letters():
    assert len(canonical_amino_acid_letters) == 20
    assert "X" not in canonical_amino_acid_letters
    expected_letters = [aa.letter for aa in canonical_amino_acids]
    assert expected_letters == canonical_amino_acid_letters

def test_extended_amino_acids():
    assert len(extended_amino_acids) > 20

def test_extended_amino_acids_letters():
    assert len(extended_amino_acid_letters) > 20
    assert "X" in extended_amino_acid_letters
    assert "J" in extended_amino_acid_letters
    expected_letters = [aa.letter for aa in extended_amino_acids]
    assert expected_letters == extended_amino_acid_letters


def test_amino_acid_representation_labels_match_values():
    alanine = canonical_amino_acids[0]
    expected = "AminoAcid(full_name='Alanine', short_name='Ala', letter='A', contains=['A'])"
    assert str(alanine) == expected
    assert repr(alanine) == expected


def test_alphabet_helpers_return_letters_and_names():
    from pepdata.amino_acid_alphabet import (
        index_to_letter, letter_to_index, peptide_to_short_amino_acid_names,
    )
    for letter in extended_amino_acid_letters:
        assert index_to_letter(letter_to_index(letter)) == letter
    assert peptide_to_short_amino_acid_names("ACX") == ["Ala", "Cys", "Xaa"]
    assert peptide_to_short_amino_acid_names("") == []


def test_pair_positions_use_two_letter_keys_in_existing_order():
    from pepdata.amino_acid_alphabet import (
        amino_acid_letter_pairs, amino_acid_pair_positions,
    )
    letters = extended_amino_acid_letters
    assert len(amino_acid_pair_positions) == len(letters) ** 2
    for j, y in enumerate(letters):
        for i, x in enumerate(letters):
            position = j * len(letters) + i
            assert amino_acid_pair_positions[x + y] == position
            assert amino_acid_letter_pairs[position] == x + y
