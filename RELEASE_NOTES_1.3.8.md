# PepData 1.3.8

Read every shipped amino-acid pair table (BLOSUM30/50/62, PMBEC and the three
residue contact energy tables) with one parser, the new public
`pepdata.pairwise_table.read_pairwise_table` and `parse_pairwise_table` (#67).
Malformed tables now raise `ValueError` instead of failing an `assert`, so
they are also rejected under `python -O`.

`parse_blosum_table`, `read_pmbec_coefficients` and `parse_interaction_table`
remain as wrappers that emit a `DeprecationWarning` and will be removed in a
future release. Their `key_type` layouts still work.
`read_pmbec_coefficients` no longer prints to stdout; its `verbose` argument
is accepted and ignored.

No table value changes: every public dictionary and matrix is compared with a
snapshot captured from 1.3.7.

The deprecated wrappers read every shipped file exactly as before, but they
are stricter than 1.3.7 on inputs no shipped file has. They now raise
`ValueError` for a header with a multi-character token, for repeated row or
column labels (the last one used to win), for a header-only table (formerly
`{}`), and for an `amino_acid_order` longer than 20 letters or with a repeated
letter. An unknown `key_type` and empty input also raise `ValueError` rather
than `AssertionError` or `IndexError`. Tab-delimited tables and indented `#`
comments are now accepted. `rekey_pairwise_table` exists only to support the
deprecated `key_type` argument and will be removed with the wrappers.
