# Residue and pairwise tables

This guide shows how to read the property tables and pairwise matrices that
ship with pepdata, and how to choose a residue alphabet. All of it works
offline.

## Read a single-residue property

```python
from pepdata.amino_acid_properties import hydropathy, volume

print(hydropathy["A"])
print(volume["A"])
```

```text
1.8
91.5
```

Each table is a dictionary. Keep the scale and source context with the values
when using them as features; the numbers from different properties are not
interchangeable measurements.

## Read a pairwise coefficient

```python
from pepdata.blosum import blosum62_dict

print(blosum62_dict["A"]["A"])
```

```text
4
```

The first and second keys select the two residue letters. The exported NumPy
matrices are 20 by 20 in `canonical_amino_acid_letters` order; extra symbols in
the source dictionaries do not add empty rows or columns. BLOSUM tables and
the other pairwise matrices describe different quantities; consult the
shipped matrix and its original source before interpreting coefficients.

## Read a pairwise table from a file

`read_pairwise_table` is the parser behind every shipped pairwise matrix. A
table with a header line of column labels and a label at the start of each row
needs only its path:

```python
from os.path import join

from pepdata.pairwise_table import read_pairwise_table
from pepdata.static_data import MATRIX_DIR

blosum62 = read_pairwise_table(join(MATRIX_DIR, "BLOSUM62"), coeff_type=int)
print(blosum62["W"]["W"])
```

```text
11
```

The residue contact energy tables have neither a header nor row labels, so
pass the residue order as `column_order`:

```python
contact = read_pairwise_table(
    join(MATRIX_DIR, "helix_vs_coil.txt"),
    column_order="ARNDCQEGHILKMFPSTWYV")
print(contact["A"]["R"])
```

```text
1.17
```

Malformed tables raise `ValueError`. The older `parse_blosum_table`,
`read_pmbec_coefficients` and `parse_interaction_table` still work but emit a
`DeprecationWarning` and will be removed in a future release.

## Use alphabets explicitly

`canonical_amino_acid_letters` describes the standard residue alphabet.
Extended and reduced alphabets are separate interfaces: do not assume every
ambiguous letter has a value in every property table. Select and document an
encoding policy before processing ambiguous sequences.
