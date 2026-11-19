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
