<a href="https://travis-ci.org/openvax/pepdata">
    <img src="https://travis-ci.org/openvax/pepdata.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/openvax/pepdata?branch=master">
    <img src="https://coveralls.io/repos/openvax/pepdata/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/pepdata/">
    <img src="https://img.shields.io/pypi/v/pepdata.svg?maxAge=1000" alt="PyPI" />
</a>

PepData
=======

Currently an amalgam of IEDB interfaces and amino acid properties. This package
will probably be eventually split and the IEDB portions placed into something
named `pyiedb`.

**Amino Acid Properties**

The `amino_acid` module contains a variety of physical/chemical properties for both single amino residues and interactions between pairs of residues.

Single residue feature tables are parsed into `StringTransformer` objects, which can be treated as dictionaries or will vectorize a string when you call their method `transform_string`.

Examples of single residue features:
- `hydropathy`
- `volume`
- `polarity`
- `pK_side_chain`
- `prct_exposed_residues`
- `hydrophilicity`
- `accessible_surface_area`
- `refractivity`
- `local_flexibility`
- `accessible_surface_area_folded`
- `alpha_helix_score` (Chou-Fasman)
- `beta_sheet_score` (Chou-Fasman)
- `turn_score` (Chou-Fasman)

Pairwise interaction tables are parsed into nested dictionaries, so that the interaction between amino acids `x` and `y` can be determined from `d[x][y]`.

Pairwise interaction dictionaries:
- `strand_vs_coil` (and its transpose `coil_vs_strand`)
- `helix_vs_strand` (and its transpose `strand_vs_helix`)
- `helix_vs_coil` (and its transpose `coil_vs_helix`)
- `blosum30`
- `blosum50`
- `blosum62`

There is also a function to parse the coefficients of the [PMBEC similarity matrix](http://www.biomedcentral.com/1471-2105/10/394), though this currently lives in the separate `pmbec` module.

