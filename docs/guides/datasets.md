# Peptide features and IEDB tables

This guide covers turning peptides into n-gram count vectors and loading the
IEDB's MHC-binding and T-cell assay tables, which need a download.

## N-gram features

```python
from pepdata import PeptideVectorizer

vectorizer = PeptideVectorizer(max_ngram=2, normalize_row=False)
counts = vectorizer.fit_transform(["SIINFEKL", "SIINFEKV"])
print(counts.shape)
print(vectorizer.transform(["SIINFEKA"]).sum())
```

```text
(2, 16)
13.0
```

`PeptideVectorizer` fits a vocabulary of residue n-grams, then transforms new
sequences using that vocabulary. `max_ngram` controls the largest n-gram, so
the two 8-mers above give 16 distinct 1- and 2-grams. `transform()` counts only
n-grams seen during fitting: `SIINFEKA` has 15, but `A` and `KA` are new.
Fit before calling `transform()`. These are sequence features, not binding
predictions.

`normalize_row=True`, the default, divides each row by its total so the
features are frequencies. On current scikit-learn releases (checked with 1.9)
that path raises `TypeError: np.matrix is not supported`; pass
`normalize_row=False` and normalize the counts yourself until it is fixed.

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
