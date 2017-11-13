from nose.tools import eq_
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
    eq_(expected_letters, canonical_amino_acid_letters)

def test_extended_amino_acids():
    assert len(extended_amino_acids) > 20

def test_extended_amino_acids_letters():
    assert len(extended_amino_acid_letters) > 20
    assert "X" in extended_amino_acid_letters
    assert "J" in extended_amino_acid_letters
    expected_letters = [aa.letter for aa in extended_amino_acids]
    eq_(expected_letters, extended_amino_acid_letters)
