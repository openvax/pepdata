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

import re
from os.path import join

from .amino_acid_alphabet import canonical_amino_acid_letters, letter_to_index
from .static_data import MATRIX_DIR

# Copy of the CRASP amino-acid property page,
# http://www.sgi.sscc.ru/mgs/programs/crasp/texts/AA_Properties.htm
CRASP_PROPERTIES_PATH = join(MATRIX_DIR, "amino_acid_properties.txt")


def aa_dict_to_positional_list(aa_property_dict):
    value_list = [None] * 20
    for letter, value in aa_property_dict.items():
        idx = letter_to_index(letter)
        assert idx >= 0
        assert idx < 20
        value_list[idx] = value
    assert all(elt is not None for elt in value_list), \
        "Missing amino acids in:\n%s" % aa_property_dict.keys()
    return value_list


def parse_property_table(table_string):
    """
    Parse lines like `1.80000 A ALA` into a residue property dictionary.

    Parameters
    ----------
    table_string : str
        One line per residue: the value, the one-letter code and, optionally,
        further fields that are ignored. Blank lines are skipped.

    Returns
    -------
    dict of str to float
        Property value for each one-letter residue code, in line order.

    Raises
    ------
    ValueError
        If a line has fewer than two fields, a value is not a number, or a
        residue is repeated.
    """
    value_dict = {}
    for line in table_string.splitlines():
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 2:
            raise ValueError("Expected a value and a residue in '%s'" % line)
        value, letter = fields[:2]
        if letter in value_dict:
            raise ValueError("Repeated amino acid in '%s'" % line)
        value_dict[letter] = float(value)
    return value_dict


def read_crasp_properties(path=CRASP_PROPERTIES_PATH):
    """
    Read every table of the CRASP amino-acid property file.

    The file is a sequence of blocks separated by blank lines. Each block is
    a numbered header such as `19. Hydropathy (Kyte and Doolittle 1982)`, a
    line with the literature citation, and one `value letter code` line per
    canonical residue. Blocks without a numbered header, such as the source
    URL on the first line, are skipped.

    Parameters
    ----------
    path : str or path-like
        Location of the property file; defaults to the copy shipped with
        pepdata.

    Returns
    -------
    dict of str to dict of str to float
        Residue property dictionaries keyed by their header, written as
        `"<number>. <title>"` with a single space after the period.

    Raises
    ------
    ValueError
        If a header is repeated, a table is malformed, or a table does not
        cover exactly the 20 canonical residues.
    """
    with open(path, "r") as f:
        blocks = re.split(r"\n\s*\n", f.read())
    properties = {}
    for block in blocks:
        lines = block.strip().splitlines()
        header = re.fullmatch(r"(\d+)\.\s*(.+)", lines[0].strip()) if lines else None
        if header is None:
            continue
        name = "%s. %s" % (header.group(1), header.group(2))
        if name in properties:
            raise ValueError("Repeated property table '%s'" % name)
        table = parse_property_table("\n".join(lines[2:]))
        if sorted(table) != sorted(canonical_amino_acid_letters):
            raise ValueError(
                "Property table '%s' does not cover the 20 canonical "
                "amino acids: %s" % (name, "".join(table)))
        properties[name] = table
    return properties


# All 36 tables of the CRASP property file, keyed by their numbered header
crasp_properties = read_crasp_properties()

# The names below are the CRASP tables this module has always exported.
# Each comment gives the literature citation as written in the CRASP file.

# Kyte and Doolittle, J. Mol. Biol. (1982) 157 105-132
hydropathy = crasp_properties["19. Hydropathy (Kyte and Doolittle 1982)"]

# Chothia, Annu. Rev. Biochem. (1984) 53 537-572
volume = crasp_properties["7. Volume (Chothia 1984)"]

# Ponnuswamy et al., Biochim.Biophys.Acta(1980) 623 301-316
polarity = crasp_properties["32. Polarity (Ponnuswamy et al. 1980)"]

# White A,Handler P,Smith EL,Hill RL,Lehman IR;Principles of Biochemistry;1978
pK_side_chain = crasp_properties["35. pK of side chain (White et al., 1978)"]

# Janin et al., J. Mol. Biol. (1978) 125 357-386
prct_exposed_residues = crasp_properties[
    "27. Percentage of exposed residues (Janin et al. 1978)"]

# Hopp and Woods, Proc. Natl. Acad. Sci. USA (1981) 78 3824-3828
hydrophilicity = crasp_properties["18. Hydrophilicity (Hopp and Woods 1981)"]

# Janin et al., J. Mol. Biol. (1978) 125 357-386
accessible_surface_area = crasp_properties[
    "12. Average accessibility surface area (Janin et al. 1978)"]

# Ragone et al., Prot. Engineering(1989) 2 497-504
local_flexibility = crasp_properties["8. Local flexibility (Ragone et al. 1989)"]

# Rose et al., Science (1975) 229 834-838
accessible_surface_area_folded = crasp_properties[
    "16. Average accessible surface area in folded proteins (Rose et al. 1985)"]

# Jones, J.theor.Biol.(1975) 50 167-183
refractivity = crasp_properties["22. Refractivity (Jones 1975)"]

# Biemann, Meth.Enzimol. (1990) v.193, p.888
mass = crasp_properties["33. Average isotopic mass (Biemann, 1990)"]

###
# Values copied from:
# "Solvent accessibility of AA in known protein structures"
# http://prowl.rockefeller.edu/aainfo/access.htm
###
"""
Solvent accessibility of AA in known protein structures

Figure 1.

S   0.70    0.20    0.10
T   0.71    0.16    0.13
A   0.48    0.35    0.17
G   0.51    0.36    0.13
P   0.78    0.13    0.09
C   0.32    0.54    0.14
D   0.81    0.09    0.10
E   0.93    0.04    0.03
Q   0.81    0.10    0.09
N   0.82    0.10    0.08
L   0.41    0.49    0.10
I   0.39    0.47    0.14
V   0.40    0.50    0.10
M   0.44    0.20    0.36
F   0.42    0.42    0.16
Y   0.67    0.20    0.13
W   0.49    0.44    0.07
K   0.93    0.02    0.05
R   0.84    0.05    0.11
H   0.66    0.19    0.15
"""

solvent_exposed_area = dict(
    S=0.70,
    T=0.71,
    A=0.48,
    G=0.51,
    P=0.78,
    C=0.32,
    D=0.81,
    E=0.93,
    Q=0.81,
    N=0.82,
    L=0.41,
    I=0.39,
    V=0.40,
    M=0.44,
    F=0.42,
    Y=0.67,
    W=0.49,
    K=0.93,
    R=0.84,
    H=0.66,
)
