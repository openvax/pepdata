# PepData

Look up shipped amino-acid properties and pairwise matrices, encode peptide
sequences, or load peptide assay tables from IEDB. Property lookups work offline;
IEDB acquisition requires a download and has separate assay filters.

## Install

```sh
python -m pip install pepdata
```

Python 3.9 or later is required.

## Look up a peptide's residue properties

```python
from pepdata.amino_acid_properties import hydropathy

peptide = "ACD"
print([hydropathy[residue] for residue in peptide])
```

```text
[1.8, 2.5, -3.5]
```

The values come from the shipped residue-property table. The dictionary maps
single-letter residues to numeric values; it does not predict MHC binding or
clinical immunogenicity. Use uppercase canonical residue letters. A missing
letter raises `KeyError`.

## Find the interface you need

| Task | Module or object |
| --- | --- |
| Single-residue properties | `pepdata.amino_acid_properties` |
| BLOSUM substitution matrices | `pepdata.blosum` |
| PMBEC and contact tables | `pepdata.pmbec`, `pepdata.residue_contact_energies` |
| Residue alphabets | `pepdata.amino_acid_alphabet`, `pepdata.reduced_alphabet` |
| N-gram features | `pepdata.PeptideVectorizer` |
| MHC-binding assay tables | `pepdata.iedb.mhc` |
| T-cell assay tables | `pepdata.iedb.tcell` |

The property dictionaries are in `amino_acid_properties`, rather than
`amino_acid`. There is no `StringTransformer` interface in this package.

See [residue and pairwise tables](docs/guides/properties.md), then
[peptide features and IEDB tables](docs/guides/datasets.md). The
[API reference](docs/reference.md) follows those same tasks.

## Documentation development

Install `requirements-docs.txt`, run `./docs.sh`, then `python scripts/check_docs_examples.py`.
