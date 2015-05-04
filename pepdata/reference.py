# Copyright (c) 2014. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from __future__ import print_function, division, absolute_import
import cPickle
from gzip import GzipFile
import re

from datacache import fetch_and_transform
import Bio.SeqIO
import pandas as pd
from progressbar import ProgressBar

from .common import int_or_seq, dataframe_from_counts, cache

BASE_URL = "ftp://ftp.ensembl.org"
FTP_DIR = "/pub/release-75/fasta/homo_sapiens/pep/"
FASTA_FILENAME = "Homo_sapiens.GRCh37.75.pep.all.fa"
GZ_FILENAME = FASTA_FILENAME + ".gz"
FULL_URL = BASE_URL + FTP_DIR + GZ_FILENAME

def load_dataframe():
    """
    Loads the protein products of the reference genome
    in a dataframe with columns:
        - protein : amino acid string
        - protein_id
        - gene_id
        - transcript_id
    """
    filename = cache.fetch(filename=FASTA_FILENAME, url=FULL_URL)
    sequences = []
    protein_ids = []
    gene_ids = []
    transcript_ids = []
    gene_pattern = re.compile('gene:(ENSG[0-9]*)')
    transcript_pattern = re.compile('transcript:(ENST[0-9]*)')

    with open(filename, 'r') as f:
        for record in Bio.SeqIO.parse(f, 'fasta'):
            protein_ids.append(record.id)
            sequences.append(str(record.seq))
            desc = record.description
            gene_matches = gene_pattern.findall(desc)
            if gene_matches:
                gene_id = gene_matches[0]
            else:
                gene_id = None
            gene_ids.append(gene_id)
            transcript_matches = transcript_pattern.findall(desc)
            if transcript_matches:
                transcript_id = transcript_matches[0]
            else:
                transcript_id = None
            transcript_ids.append(transcript_id)
    df = pd.DataFrame({
        'protein': sequences,
        'gene_id': gene_ids,
        'protein_id': protein_ids,
        'transcript_id': transcript_ids})
    return df


def _generate_counts(src_filename, peptide_lengths, nrows):
    epitope_counts = {}
    with open(src_filename, 'r') as f:
        seqs = [str(record.seq) for record in Bio.SeqIO.parse(f, "fasta")]
        print("Generating substrings of length %s" % (peptide_lengths,))
        pbar = ProgressBar(maxval=len(seqs)).start()
        for seq_num, seq in enumerate(seqs):
            seq_len = len(seq)
            if nrows and seq_num > nrows:
                break
            for size in peptide_lengths:
                for i in xrange(seq_len - size + 1):
                    epitope = seq[i:i + size]
                    if epitope in epitope_counts:
                        epitope_counts[epitope] += 1
                    else:
                        epitope_counts[epitope] = 1
            pbar.update(seq_num + 1)
        pbar.finish()
    return dataframe_from_counts(epitope_counts)

def _generate_set(src_filename, peptide_lengths, nrows):
    peptides = set([])
    with open(src_filename, 'r') as f:
        seqs = [str(record.seq) for record in Bio.SeqIO.parse(f, "fasta")]
        print("Generating substrings of length %s" % (peptide_lengths,))
        pbar = ProgressBar(maxval=len(seqs)).start()
        for seq_num, seq in enumerate(seqs):
            if nrows and seq_num > nrows:
                break
            for size in peptide_lengths:
                for i in xrange(len(seq) - size + 1):
                    peptides.add(seq[i:i + size])
            pbar.update(seq_num + 1)
        pbar.finish()
    return peptides


def load_peptide_counts(peptide_length=[8, 9, 10, 11], nrows=None):
    """
    List of all reference peptides encoded in a reference human exome
    """
    peptide_lengths = int_or_seq(peptide_length)
    lens = "_".join(str(n) for n in peptide_lengths)
    cache_filename = \
        "reference_peptide_counts_" + lens + "_nrows_" + str(nrows) + ".csv"

    def save_counts(src_path, dst_path):
        counts = _generate_counts(src_path, peptide_lengths, nrows)
        print("Saving %s" % dst_path)
        counts.to_csv(dst_path)
        return counts

    return fetch_and_transform(
        transformed_filename=cache_filename,
        transformer=save_counts,
        loader=pd.read_csv,
        source_filename=FASTA_FILENAME,
        source_url=FULL_URL,
        subdir=cache.subdir)

def load_peptide_set(peptide_length=[8, 9, 10, 11], nrows=None):
    peptide_lengths = int_or_seq(peptide_length)
    lens = "_".join(str(n) for n in peptide_lengths)
    cache_filename = \
        "reference_peptide_set_" + lens + "_nrows_" + str(nrows) + ".pickle.gz"

    def save_set(src_path, dst_path):
        string_set = _generate_set(src_path, peptide_lengths, nrows)
        with GzipFile(dst_path, 'w') as out_file:
            out_file.write(cPickle.dumps(string_set))
        return string_set

    def load_set(path):
        result = None
        with GzipFile(path, 'r') as in_file:
            result = cPickle.loads(in_file.read())
        return result

    return fetch_and_transform(
        transformed_filename=cache_filename,
        transformer=save_set,
        loader=load_set,
        source_filename=FASTA_FILENAME,
        source_url=FULL_URL,
        subdir=cache.subdir)
