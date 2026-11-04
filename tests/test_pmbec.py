import numpy as np

from pepdata.amino_acid_alphabet import canonical_amino_acid_letters
from pepdata.pmbec import (
    pmbec_dict,
    pmbec_matrix,
)


def test_pmbec_covers_canonical_residues():
    assert sorted(pmbec_dict) == sorted(canonical_amino_acid_letters)
    for row in pmbec_dict.values():
        assert sorted(row) == sorted(canonical_amino_acid_letters)


def test_pmbec_matrix_is_symmetric_and_matches_dict():
    assert pmbec_matrix.shape == (20, 20)
    assert np.allclose(pmbec_matrix, pmbec_matrix.T)
    for i, x in enumerate(canonical_amino_acid_letters):
        for j, y in enumerate(canonical_amino_acid_letters):
            assert np.isclose(pmbec_matrix[i, j], pmbec_dict[x][y])


def test_importing_pmbec_prints_nothing():
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "-c", "import pepdata.pmbec"],
        capture_output=True, text=True, check=True)
    assert result.stdout == ""
