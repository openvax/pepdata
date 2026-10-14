# Peptide features and IEDB tables

## N-gram features

`PeptideVectorizer` fits a vocabulary of residue n-grams, then transforms new
sequences using that vocabulary. `max_ngram` controls the largest n-gram;
`normalize_row` selects count normalization. Fit before calling `transform()`.
These are sequence features, not binding predictions.

## Load assay data

These templates require network access on the first use and can download large
IEDB exports. Review the loader's actual filters before interpreting a subset:

```python
from pepdata.iedb import mhc, tcell

binding_assays = mhc.load_dataframe()
tcell_assays = tcell.load_dataframe()
```

The MHC-binding and T-cell tables represent different assays. A positive
binding measurement is not automatically a positive T-cell response. Preserve
assay, species and allele metadata when selecting records; the
[IEDB](https://www.iedb.org/) is the upstream data source.

Each module exposes `download(force=False)` and `local_path(auto_download=True)`.
Use those explicit acquisition interfaces when preparing an environment before
loading a table. `force=True` requests a refresh; a cached export is a snapshot
and does not by itself identify a pinned IEDB release.

See the [loader reference](../reference.md#assay-tables) for filtering arguments.
