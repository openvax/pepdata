# PepData 1.3.7

Correct the shipped BLOSUM50 W-to-T substitution score from -4 to the NCBI
reference value -3 (#59). It now agrees with T-to-W and restores symmetry
over the canonical amino-acid alphabet. Applications using this entry will
receive the corrected numerical score. All 400 canonical entries are checked
against an independently downloaded NCBI reference fixture.
