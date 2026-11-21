from pepdata.amino_acid_alphabet import canonical_amino_acid_letters
from pepdata.chou_fasman import (
    alpha_helix_score,
    beta_sheet_score,
    turn_score,
)


def test_chou_fasman_scores_keyed_by_one_letter_codes():
    expected_letters = set(canonical_amino_acid_letters)
    for table in (alpha_helix_score, beta_sheet_score, turn_score):
        assert set(table.keys()) == expected_letters
        assert all(isinstance(k, str) and len(k) == 1 for k in table)
        assert all(isinstance(v, int) for v in table.values())


def test_chou_fasman_known_values():
    assert alpha_helix_score["A"] == 142
    assert beta_sheet_score["A"] == 83
    assert turn_score["A"] == 66

    assert alpha_helix_score["D"] == 101
    assert beta_sheet_score["D"] == 54
    assert turn_score["D"] == 146

    assert alpha_helix_score["E"] == 151
    assert beta_sheet_score["E"] == 37
    assert turn_score["E"] == 74

    assert alpha_helix_score["G"] == 57
    assert beta_sheet_score["G"] == 75
    assert turn_score["G"] == 156

    assert alpha_helix_score["V"] == 106
    assert beta_sheet_score["V"] == 170
    assert turn_score["V"] == 50
