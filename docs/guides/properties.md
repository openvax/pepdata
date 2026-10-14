# Residue and pairwise tables

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

The first and second keys select the two residue letters. BLOSUM tables and
the other pairwise matrices describe different quantities; consult the
shipped matrix and its original source before interpreting coefficients.

## Use alphabets explicitly

`canonical_amino_acid_letters` describes the standard residue alphabet.
Extended and reduced alphabets are separate interfaces: do not assume every
ambiguous letter has a value in every property table. Select and document an
encoding policy before processing ambiguous sequences.
