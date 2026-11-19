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
Read whitespace-delimited tables of amino-acid pair coefficients.

One parser serves the BLOSUM substitution matrices, the PMBEC matrix and the
residue contact energy tables that ship in `pepdata/matrices`.
"""


def parse_pairwise_table(text, column_order=None, coeff_type=float):
    """
    Parse the text of a table of pairwise amino-acid coefficients.

    Blank lines and lines starting with `#` are ignored.

    Parameters
    ----------
    text : str
        Contents of the table.

    column_order : sequence of str, optional
        Residue labels of a table that has no header line and no row labels,
        such as the residue contact energy tables. The same labels name the
        rows and the columns, so the table must be square. When omitted, the
        first line holds the column labels and every later line starts with
        its row label.

    coeff_type : callable, default float
        Converts each coefficient string to a number, e.g. `int` for
        BLOSUM scores.

    Returns
    -------
    dict of str to dict of str to number
        Coefficients accessed like `table[row_label][column_label]`.

    Raises
    ------
    ValueError
        If the table is empty, a label is repeated, a row has the wrong
        number of fields, or a coefficient is not a number.
    """
    rows = [
        line.split()
        for line in text.splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    if column_order is None:
        if not rows:
            raise ValueError("Pairwise table has no header line")
        column_labels = rows[0]
        labeled_rows = [(fields[0], fields[1:]) for fields in rows[1:]]
    else:
        column_labels = list(column_order)
        if len(rows) != len(column_labels):
            raise ValueError(
                "Expected %d rows in pairwise table but found %d" % (
                    len(column_labels), len(rows)))
        labeled_rows = list(zip(column_labels, rows))

    if len(set(column_labels)) != len(column_labels):
        raise ValueError(
            "Repeated column label in pairwise table: %s" % (column_labels,))
    if not labeled_rows:
        raise ValueError("Pairwise table has no rows")

    table = {}
    for row_label, coeff_strings in labeled_rows:
        if row_label in table:
            raise ValueError(
                "Repeated row label '%s' in pairwise table" % row_label)
        if len(coeff_strings) != len(column_labels):
            raise ValueError(
                "Expected %d coefficients in row '%s' but found %d" % (
                    len(column_labels), row_label, len(coeff_strings)))
        table[row_label] = {
            column_label: coeff_type(coeff_string)
            for column_label, coeff_string in zip(column_labels, coeff_strings)
        }
    return table


def read_pairwise_table(path, column_order=None, coeff_type=float):
    """
    Read a file holding a table of pairwise amino-acid coefficients.

    Parameters
    ----------
    path : str or path-like
        Location of the table.

    column_order : sequence of str, optional
        Residue labels of a table without a header line or row labels; see
        `parse_pairwise_table`.

    coeff_type : callable, default float
        Converts each coefficient string to a number.

    Returns
    -------
    dict of str to dict of str to number
        Coefficients accessed like `table[row_label][column_label]`.

    Raises
    ------
    ValueError
        If the table is malformed; see `parse_pairwise_table`.
    """
    with open(path, "r") as f:
        return parse_pairwise_table(
            f.read(), column_order=column_order, coeff_type=coeff_type)


def rekey_pairwise_table(table, key_type="row"):
    """
    Convert a nested pairwise table to one of the legacy key layouts.

    Exists to support the deprecated `key_type` argument of
    `pepdata.blosum.parse_blosum_table` and
    `pepdata.pmbec.read_pmbec_coefficients`.

    Parameters
    ----------
    table : dict of str to dict of str to number
        Coefficients accessed like `table[row_label][column_label]`.

    key_type : {'row', 'pair', 'pair_string'}
        `'row'` returns the nested layout unchanged, `'pair'` keys each
        coefficient by a `(row_label, column_label)` tuple and
        `'pair_string'` by the two labels joined into one string.

    Returns
    -------
    dict
        The coefficients in the requested layout.

    Raises
    ------
    ValueError
        If `key_type` is not one of the three layouts.
    """
    if key_type == "row":
        return table
    if key_type == "pair":
        return {
            (row_label, column_label): coeff
            for row_label, row in table.items()
            for column_label, coeff in row.items()
        }
    if key_type == "pair_string":
        return {
            row_label + column_label: coeff
            for row_label, row in table.items()
            for column_label, coeff in row.items()
        }
    raise ValueError("Unknown key type: %s" % key_type)
