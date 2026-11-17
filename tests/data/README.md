# Synthetic IEDB-format fixtures

These six-row CSVs are artificial input for allele filtering tests. They use
the two-row IEDB-style column header; sequences and assay labels are invented.
They are not IEDB observations, binding measurements or T-cell responses.

The rows cover broad and specific A2/A24 labels, another human allele and a
nonhuman allele, allowing the existing regex and exclusion tests to run
without downloading the full assay tables.

`alleles_test_fixture.xml` is a synthetic four-entry catalogue for the allele
parser tests, covering human and mouse class I and II labels. It retains the
parser's historical `Organsim` field spelling. It is not an IEDB download and
does not validate the current live catalogue or its availability.

## BLOSUM50 reference

`ncbi_BLOSUM50.txt` is the NCBI BLAST reference matrix downloaded from
https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM50 on 2026-11-17.
SHA-256: `a152e30b8afcab42281ff8604bc51358be98253155fde62ae20707805a0d3d67`.
Tests compare all 400 canonical-residue entries against this independent
reference; ambiguity symbols are outside this regression.
