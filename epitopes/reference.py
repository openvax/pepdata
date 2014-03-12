from collections import Counter

import Bio.SeqIO
import pandas as pd
from progressbar import ProgressBar

from download import fetch_data

def load_peptide_counts(peptide_lengths = [8,9,10,11]):
    """
    List of all reference peptides encoded in a reference human exome
    """
    base_url = "ftp://ftp.ensembl.org"
    ftp_dir = "/pub/release-75/fasta/homo_sapiens/pep/"
    fasta_name =  "Homo_sapiens.GRCh37.75.pep.all.fa"
    gz_name =  fasta_name + ".gz"
    full_url = base_url + ftp_dir + gz_name
    local_path = fetch_data(fasta_name, full_url)
    epitope_counts = Counter()
    with open(local_path, 'r') as f:
        seqs = [str(record.seq) for record in Bio.SeqIO.parse(f, "fasta")]
        print "Generating substrings"
        pbar = ProgressBar(maxval = len(seqs)).start()
        for seq_num, seq in enumerate(seqs):
            for size in peptide_lengths:
                for i in range(len(seq) - size + 1):
                    epitope_counts[seq[i:i+size]] += 1
            pbar.update(seq_num+1)
        pbar.finish()
    return pd.DataFrame(
        epitopes.most_common(),
        columns=["Peptide", "Count"])

def load_peptides(peptide_lengths = [8,9,10,11]):
    df = load_peptide_counts(peptide_lengths)
    return set(df.Peptide)
