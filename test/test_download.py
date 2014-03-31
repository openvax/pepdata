from epitopes.download import fetch_fasta_dict, fetch_fasta_db, fetch_data

FASTA_FILENAME = 'Homo_sapiens.GRCh37.75.dna_rm.chromosome.MT.fa'
DB_FILENAME = 'Homo_sapiens.GRCh37.75.dna_rm.chromosome.MT.db'
URL = \
'ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna_rm.chromosome.MT.fa.gz'

def test_fetch_data():
    path = fetch_data(FASTA_FILENAME, URL)
    assert path.endswith(FASTA_FILENAME)

    # if we change the subdir then data should end up in
    # something like /Users/me/Library/Caches/epitopes_test/
    other_path = fetch_data(FASTA_FILENAME, URL, subdir = "epitopes_test")
    assert other_path.endswith(FASTA_FILENAME)
    assert other_path != path, other_path

def test_download_fasta_dict():
    d = fetch_fasta_dict(FASTA_FILENAME, URL)
    assert hasattr(d, 'keys'), d
    assert hasattr(d, 'values'), d
    assert len(d) > 0

def test_downloadd_fasta_db():
    db = fetch_fasta_db(DB_FILENAME, "DNA", FASTA_FILENAME, URL)
    assert db is not None