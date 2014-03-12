from download import fetch_data

def load_reference_peptides():
    """
    List of all reference peptides encoded in a reference human exome
    """
    base_url = "ftp://ftp.ensembl.org"
    ftp_dir = "/pub/release-75/fasta/homo_sapiens/pep/"
    fasta_name =  "Homo_sapiens.GRCh37.75.pep.all.fa"
    gz_name =  fasta_name + ".gz"
    full_url = base_url + ftp_dir + gz_name
    local_path = fetch_data(fasta_name, full_url)
    with open(local_path, 'r') as f:
        print f.read()