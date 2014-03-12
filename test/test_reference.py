import pandas as pd
from epitopes import reference

def test_reference_peptide_counts():
    counts = reference.load_peptide_counts(peptide_lengths = [8], nrows = 20)
    assert counts is not None
    assert isinstance(counts, pd.DataFrame)
    assert (counts.Peptide.str.len() == 8).all()

def test_reference_peptide_set():
    peptides = reference.load_peptide_set(peptide_lengths = [8], nrows = 20)
    assert peptides is not None
    assert isinstance(peptides, set)
    assert all(len(p) == 8 for p in peptides)