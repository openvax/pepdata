from pepdata.amino_acid_alphabet import canonical_amino_acid_letters
from pepdata.reduced_alphabet import (
    ALPHABETS,
    alex6,
    aromatic2,
    gbmr4,
    hp2,
    hp_vs_aromatic,
    hsdm17,
    murphy8,
    murphy10,
    murphy15,
    sdm12,
)


def test_all_reduced_alphabets_cover_canonical_amino_acids():
    expected_counts = {
        "gbmr4": 4,
        "sdm12": 12,
        "hsdm17": 17,
        "hp2": 2,
        "murphy8": 8,
        "murphy10": 10,
        "murphy15": 15,
        "alex6": 6,
        "aromatic2": 2,
        "hp_vs_aromatic": 4,
    }
    assert list(ALPHABETS) == list(expected_counts)
    canonical = set(canonical_amino_acid_letters)
    for name, expected_groups in expected_counts.items():
        mapping = ALPHABETS[name]
        assert set(mapping) == canonical
        assert len(set(mapping.values())) == expected_groups


def test_named_exports_match_alphabets_registry():
    assert ALPHABETS["gbmr4"] is gbmr4
    assert ALPHABETS["sdm12"] is sdm12
    assert ALPHABETS["hsdm17"] is hsdm17
    assert ALPHABETS["hp2"] is hp2
    assert ALPHABETS["murphy8"] is murphy8
    assert ALPHABETS["murphy10"] is murphy10
    assert ALPHABETS["murphy15"] is murphy15
    assert ALPHABETS["alex6"] is alex6
    assert ALPHABETS["aromatic2"] is aromatic2
    assert ALPHABETS["hp_vs_aromatic"] is hp_vs_aromatic
