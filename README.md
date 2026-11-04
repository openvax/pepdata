[![Tests](https://github.com/openvax/pepdata/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/openvax/pepdata/actions/workflows/tests.yml)
[![Documentation](https://github.com/openvax/pepdata/actions/workflows/documentation.yml/badge.svg?branch=master)](https://openvax.github.io/pepdata/)
[![PyPI](https://img.shields.io/pypi/v/pepdata.svg)](https://pypi.org/project/pepdata/)

# pepdata

pepdata ships amino-acid property tables and pairwise residue matrices
(BLOSUM, PMBEC, contact energies), turns peptide sequences into n-gram feature
vectors, and loads MHC-binding and T-cell assay tables from the
[IEDB](https://www.iedb.org/). It is the [OpenVax](https://github.com/openvax)
library for peptide-level reference data. It does not predict MHC binding; for
that, see [mhctools](https://github.com/openvax/mhctools) and
[mhcflurry](https://github.com/openvax/mhcflurry).

## Install

```sh
python -m pip install pepdata
```

pepdata needs Python 3.11 or later. Property tables and matrices ship with the
package and work offline. The IEDB tables are downloaded on first use and
cached by [datacache](https://github.com/openvax/datacache).

## Quick start

```python
from pepdata.amino_acid_properties import hydropathy

peptide = "ACD"
print([hydropathy[residue] for residue in peptide])
```

```text
[1.8, 2.5, -3.5]
```

Each property table is a dictionary from single-letter residue to value. Use
uppercase canonical residue letters; a missing letter raises `KeyError`.

## Documentation

The [documentation](https://openvax.github.io/pepdata/) has two guides and an
API reference:

1. [Residue and pairwise tables](https://openvax.github.io/pepdata/guides/properties/):
   property dictionaries, BLOSUM and other matrices, alphabets.
2. [Peptide features and IEDB tables](https://openvax.github.io/pepdata/guides/datasets/):
   n-gram vectors and the MHC-binding and T-cell assay loaders.
3. [API reference](https://openvax.github.io/pepdata/reference/).

| Task | Module or object |
| --- | --- |
| Single-residue properties | `pepdata.amino_acid_properties` |
| BLOSUM substitution matrices | `pepdata.blosum` |
| PMBEC and contact tables | `pepdata.pmbec`, `pepdata.residue_contact_energies` |
| Residue alphabets | `pepdata.amino_acid_alphabet`, `pepdata.reduced_alphabet` |
| N-gram features | `pepdata.PeptideVectorizer` |
| MHC-binding assay tables | `pepdata.iedb.mhc` |
| T-cell assay tables | `pepdata.iedb.tcell` |

## Development

```sh
python -m pip install -e . pytest pytest-cov ruff
./lint.sh   # ruff check
./test.sh   # pytest with coverage
python -m pip install -r requirements-docs.txt
./docs.sh   # mkdocs build --strict, then a rendered-link check
python scripts/check_docs_examples.py   # run docs examples, compare output
```

## License

Apache 2.0; see [LICENSE.txt](https://github.com/openvax/pepdata/blob/master/LICENSE.txt).
