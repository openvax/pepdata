# PepData 1.3.9

Load the single-residue property tables from the shipped, cited
`pepdata/matrices/amino_acid_properties.txt` instead of from re-typed copies
in `pepdata/amino_acid_properties.py` (#66). The new public
`read_crasp_properties` parses the file and `crasp_properties` exposes all 36
of its tables keyed by their numbered header, e.g.
`crasp_properties["34. Isoelectric point (Zimmerman et al., 1968)"]`; 25 of
them were shipped before but could not be loaded.

The 11 existing names (`hydropathy`, `volume`, `polarity`, `pK_side_chain`,
`prct_exposed_residues`, `hydrophilicity`, `accessible_surface_area`,
`local_flexibility`, `accessible_surface_area_folded`, `refractivity`, `mass`)
are now lookups into that file, each with its literature citation. No value
changes: they are compared with a snapshot captured from 1.3.7.

`parse_property_table` raises `ValueError` rather than failing an `assert` on
malformed input.
