# Peptide features and IEDB tables

## N-gram features

`PeptideVectorizer` fits a vocabulary of residue n-grams, then transforms new
sequences using that vocabulary. `max_ngram` controls the largest n-gram;
`normalize_row` selects count normalization. Fit before calling `transform()`.
These are sequence features, not binding predictions.

## Load assay data

These templates require network access on the first use and can download large
IEDB exports:

```python
from pepdata.iedb import mhc, tcell

binding_assays = mhc.load_dataframe()
tcell_assays = tcell.load_dataframe()
```

By default, both loaders keep only peptides made of the 20 standard amino
acids (`only_standard_amino_acids=True`). The
[loader reference](../reference.md#assay-tables) lists the other filtering
arguments, such as `mhc_class`, `hla` and `human_only`.

The MHC-binding and T-cell tables represent different assays. A positive
binding measurement is not automatically a positive T-cell response. Preserve
assay, species and allele metadata when selecting records; the
[IEDB](https://www.iedb.org/) is the upstream data source.

Each module exposes `download(force=False)` and `local_path(auto_download=True)`.
Call `download()` to fetch the export before loading a table, for example when
preparing an environment. `force=True` requests a refresh; a cached export is a
snapshot and does not by itself identify a pinned IEDB release.
