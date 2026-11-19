import subprocess
import sys

import pytest

from pepdata.blosum import blosum62_dict, parse_blosum_table
from pepdata.pairwise_table import (
    parse_pairwise_table,
    read_pairwise_table,
    rekey_pairwise_table,
)
from pepdata.pmbec import PMBEC_PATH, pmbec_dict, read_pmbec_coefficients
from pepdata.residue_contact_energies import (
    helix_vs_coil_dict,
    parse_interaction_table,
)
from pepdata.static_data import MATRIX_DIR

HEADER_TABLE = """
# comment
   A  C
A  4 -1
C -1  9
"""

HEADERLESS_TABLE = """
#M rows = AC, cols = AC
 0.5 -1.25
 2.0  3.0
"""


def test_header_table_uses_first_line_as_column_labels():
    assert parse_pairwise_table(HEADER_TABLE) == {
        "A": {"A": 4.0, "C": -1.0},
        "C": {"A": -1.0, "C": 9.0},
    }


def test_coeff_type_controls_the_value_type():
    table = parse_pairwise_table(HEADER_TABLE, coeff_type=int)
    assert table["C"]["C"] == 9
    assert type(table["C"]["C"]) is int
    assert type(parse_pairwise_table(HEADER_TABLE)["C"]["C"]) is float


def test_column_order_labels_a_headerless_table():
    assert parse_pairwise_table(HEADERLESS_TABLE, column_order="AC") == {
        "A": {"A": 0.5, "C": -1.25},
        "C": {"A": 2.0, "C": 3.0},
    }


def test_carriage_returns_and_indented_comments_are_ignored():
    text = "  # indented comment\r\n   A  C \r\nA  4 -1 \r\n\r\nC -1  9\r\n"
    assert parse_pairwise_table(text, coeff_type=int) == {
        "A": {"A": 4, "C": -1},
        "C": {"A": -1, "C": 9},
    }


def test_rectangular_header_table_is_allowed():
    table = parse_pairwise_table(" A C\nA 1 2\n")
    assert table == {"A": {"A": 1.0, "C": 2.0}}


@pytest.mark.parametrize("text, column_order", [
    ("", None),
    ("# only a comment\n", None),
    (" A C\n", None),
    (" A C\nA 1\n", None),
    (" A C\nA 1 2 3\n", None),
    (" A C\nA 1 x\n", None),
    (" A A\nA 1 2\n", None),
    (" A C\nA 1 2\nA 3 4\n", None),
    ("1 2\n", "AC"),
    ("1 2\n3\n", "AC"),
    ("1 2\n3 4\n5 6\n", "AC"),
    ("1 2\n3 4\n", "AA"),
])
def test_malformed_tables_raise_value_error(text, column_order):
    with pytest.raises(ValueError):
        parse_pairwise_table(text, column_order=column_order)


def test_malformed_tables_are_rejected_with_assertions_disabled(tmp_path):
    path = tmp_path / "short_row.txt"
    path.write_text(" A C\nA 1\n")
    code = (
        "import sys\n"
        "from pepdata.pairwise_table import read_pairwise_table\n"
        "try:\n"
        "    read_pairwise_table(sys.argv[1])\n"
        "except ValueError:\n"
        "    sys.exit(7)\n"
    )
    result = subprocess.run([sys.executable, "-O", "-c", code, str(path)])
    assert result.returncode == 7


def test_read_pairwise_table_reads_a_file(tmp_path):
    path = tmp_path / "table.txt"
    path.write_text(HEADERLESS_TABLE)
    assert read_pairwise_table(path, column_order="AC")["C"]["A"] == 2.0
    assert read_pairwise_table(str(path), column_order=["A", "C"])["A"]["C"] == -1.25


def test_rekey_pairwise_table_layouts():
    table = parse_pairwise_table(HEADER_TABLE, coeff_type=int)
    assert rekey_pairwise_table(table) is table
    assert rekey_pairwise_table(table, "pair") == {
        ("A", "A"): 4, ("A", "C"): -1, ("C", "A"): -1, ("C", "C"): 9}
    assert rekey_pairwise_table(table, "pair_string") == {
        "AA": 4, "AC": -1, "CA": -1, "CC": 9}
    with pytest.raises(ValueError):
        rekey_pairwise_table(table, "column")


def read_matrix_text(name):
    with open(MATRIX_DIR + "/" + name) as f:
        return f.read()


def test_parse_blosum_table_is_a_deprecated_wrapper():
    text = read_matrix_text("BLOSUM62")
    with pytest.warns(DeprecationWarning, match="parse_pairwise_table"):
        assert parse_blosum_table(text) == blosum62_dict
    with pytest.warns(DeprecationWarning):
        pairs = parse_blosum_table(text, key_type="pair")
    with pytest.warns(DeprecationWarning):
        pair_strings = parse_blosum_table(text, coeff_type=float, key_type="pair_string")
    assert pairs[("W", "W")] == 11 and type(pairs[("W", "W")]) is int
    assert pair_strings["WA"] == -3.0 and type(pair_strings["WA"]) is float
    assert len(pairs) == len(pair_strings) == len(blosum62_dict) ** 2
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError):
        parse_blosum_table(HEADER_TABLE)
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError):
        parse_blosum_table(text, key_type="column")


def test_read_pmbec_coefficients_is_a_deprecated_silent_wrapper(capsys):
    with pytest.warns(DeprecationWarning, match="read_pairwise_table"):
        assert read_pmbec_coefficients() == pmbec_dict
    with pytest.warns(DeprecationWarning):
        assert read_pmbec_coefficients(verbose=True, filename=PMBEC_PATH) == pmbec_dict
    with pytest.warns(DeprecationWarning):
        pairs = read_pmbec_coefficients(key_type="pair")
    with pytest.warns(DeprecationWarning):
        pair_strings = read_pmbec_coefficients("pair_string")
    assert pairs[("V", "R")] == pair_strings["VR"] == pmbec_dict["V"]["R"]
    assert len(pairs) == len(pair_strings) == 400
    assert capsys.readouterr().out == ""


def test_parse_interaction_table_is_a_deprecated_wrapper():
    text = read_matrix_text("helix_vs_coil.txt")
    with pytest.warns(DeprecationWarning, match="parse_pairwise_table"):
        assert parse_interaction_table(text) == helix_vs_coil_dict
    with pytest.warns(DeprecationWarning):
        assert parse_interaction_table(HEADERLESS_TABLE, amino_acid_order="AC") == {
            "A": {"A": 0.5, "C": -1.25},
            "C": {"A": 2.0, "C": 3.0},
        }
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError):
        parse_interaction_table(HEADERLESS_TABLE)


def test_importing_table_modules_emits_no_deprecation_warning():
    code = "import pepdata.blosum, pepdata.pmbec, pepdata.residue_contact_energies"
    subprocess.run(
        [sys.executable, "-W", "error::DeprecationWarning", "-c", code], check=True)
