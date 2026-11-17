# PepData 1.3.6

Fix all-zero transforms after fitting already-reduced peptide sequences (#69).
Transforms apply the alphabet reduction once, then use the same case handling
as training. Existing lowercase vocabularies and their feature-column order
are preserved, including legacy fitted vectorizers. Newly fitted
already-reduced vectorizers also support standard pickle roundtrips.
