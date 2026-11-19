import subprocess
import sys

import pytest

from pepdata import amino_acid_properties
from pepdata.amino_acid_alphabet import canonical_amino_acid_letters
from pepdata.amino_acid_properties import (
    CRASP_PROPERTIES_PATH,
    crasp_properties,
    parse_property_table,
    read_crasp_properties,
)

# Module-level names and the CRASP table each one is read from.
NAMED_TABLES = {
    "hydropathy": "19. Hydropathy (Kyte and Doolittle 1982)",
    "volume": "7. Volume (Chothia 1984)",
    "polarity": "32. Polarity (Ponnuswamy et al. 1980)",
    "pK_side_chain": "35. pK of side chain (White et al., 1978)",
    "prct_exposed_residues": "27. Percentage of exposed residues (Janin et al. 1978)",
    "hydrophilicity": "18. Hydrophilicity (Hopp and Woods 1981)",
    "accessible_surface_area": "12. Average accessibility surface area (Janin et al. 1978)",
    "local_flexibility": "8. Local flexibility (Ragone et al. 1989)",
    "accessible_surface_area_folded":
        "16. Average accessible surface area in folded proteins (Rose et al. 1985)",
    "refractivity": "22. Refractivity (Jones 1975)",
    "mass": "33. Average isotopic mass (Biemann, 1990)",
}


def write_property_file(tmp_path, text):
    path = tmp_path / "properties.txt"
    path.write_text(text)
    return path


def property_block(header, letters=canonical_amino_acid_letters):
    lines = [header, "Some Journal (1999) 1 2-3"]
    lines += ["%d.5000 %s XXX" % (i, letter) for i, letter in enumerate(letters)]
    return "\n".join(lines) + "\n"


def test_shipped_file_has_36_numbered_tables_in_order():
    assert crasp_properties == read_crasp_properties(CRASP_PROPERTIES_PATH)
    numbers = [int(name.split(".")[0]) for name in crasp_properties]
    assert numbers == list(range(1, 37))
    for name, table in crasp_properties.items():
        assert list(table) == list("ARNDCQEGHILKMFPSTWYV"), name
        assert all(type(value) is float for value in table.values()), name


def test_header_without_space_after_number_is_normalized():
    assert "34. Isoelectric point (Zimmerman et al., 1968)" in crasp_properties
    isoelectric = crasp_properties["34. Isoelectric point (Zimmerman et al., 1968)"]
    assert isoelectric["K"] == 9.74


@pytest.mark.parametrize("name", sorted(NAMED_TABLES))
def test_named_tables_are_lookups_into_the_file(name):
    assert getattr(amino_acid_properties, name) == crasp_properties[NAMED_TABLES[name]]


def test_tables_not_previously_exported_are_available():
    table = crasp_properties[
        "36. Energy of transfer from water to ethanol kcal/mol (Nozaki & Tanford, 1971)"]
    assert table["W"] == -3.0


def test_read_crasp_properties_skips_blocks_without_a_numbered_header(tmp_path):
    text = (
        "http://example.org/source\n\n"
        + property_block("1. First (Author 1999)") + "\n\n"
        + property_block("2.Second (Author 2000)"))
    properties = read_crasp_properties(write_property_file(tmp_path, text))
    assert list(properties) == ["1. First (Author 1999)", "2. Second (Author 2000)"]
    assert properties["1. First (Author 1999)"]["R"] == 1.5


@pytest.mark.parametrize("text", [
    property_block("1. Same") + "\n" + property_block("1. Same"),
    property_block("1. Missing valine", canonical_amino_acid_letters[:-1]),
    property_block("1. Extra residue", canonical_amino_acid_letters + ["X"]),
    property_block("1. Repeated residue", canonical_amino_acid_letters + ["A"]),
    property_block("1. Not a number").replace("3.5000", "abc"),
    property_block("1. Missing letter").replace("3.5000 D XXX", "3.5000"),
])
def test_malformed_property_files_raise_value_error(tmp_path, text):
    with pytest.raises(ValueError):
        read_crasp_properties(write_property_file(tmp_path, text))


def test_parse_property_table():
    assert parse_property_table("\n1.80000 A ALA\n-4.5000  R\n\n") == {"A": 1.8, "R": -4.5}
    for text in ["1.8", "1.8 A\n2.0 A", "x A"]:
        with pytest.raises(ValueError):
            parse_property_table(text)


def test_malformed_tables_are_rejected_with_assertions_disabled():
    code = (
        "import sys\n"
        "from pepdata.amino_acid_properties import parse_property_table\n"
        "try:\n"
        "    parse_property_table('1.8 A\\n2.0 A')\n"
        "except ValueError:\n"
        "    sys.exit(7)\n"
    )
    assert subprocess.run([sys.executable, "-O", "-c", code]).returncode == 7
