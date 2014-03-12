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

from gzip import GzipFile
from cPickle import dumps, dump, loads, load

import Bio.SeqIO
import pandas as pd
from progressbar import ProgressBar

from download import fetch_data, fetch_and_transform_data

def _most_common(d):
    return list(sorted(d.iteritems(), key=lambda (x,y): y))

def _generate_counts(src_filename, peptide_lengths, nrows):
    epitope_counts = {}
    get_count = epitope_counts.get
    with open(src_filename, 'r') as f:
        seqs = [str(record.seq) for record in Bio.SeqIO.parse(f, "fasta")]
        print "Generating substrings of length %s" % (peptide_lengths,)
        pbar = ProgressBar(maxval = len(seqs)).start()
        for seq_num, seq in enumerate(seqs):
            if nrows and seq_num > nrows:
                break
            for size in peptide_lengths:
                for i in xrange(len(seq) - size + 1):
                    epitope = seq[i:i+size]
                    epitope_counts[epitope] = get_count(epitope, 0) + 1
            pbar.update(seq_num+1)
        pbar.finish()
    return pd.DataFrame(
        _most_common(epitope_counts),
        columns=["Peptide", "Count"])

def _generate_set(src_filename, peptide_lengths, nrows):
    peptides = set([])
    with open(src_filename, 'r') as f:
        seqs = [str(record.seq) for record in Bio.SeqIO.parse(f, "fasta")]
        print "Generating substrings of length %s" % (peptide_lengths,)
        pbar = ProgressBar(maxval = len(seqs)).start()
        for seq_num, seq in enumerate(seqs):
            if nrows and seq_num > nrows:
                break
            for size in peptide_lengths:
                for i in xrange(len(seq) - size + 1):
                    peptides.add(seq[i:i+size])
            pbar.update(seq_num+1)
        pbar.finish()
    return peptides

BASE_URL = "ftp://ftp.ensembl.org"
FTP_DIR = "/pub/release-75/fasta/homo_sapiens/pep/"
FASTA_FILENAME = "Homo_sapiens.GRCh37.75.pep.all.fa"
GZ_FILENAME = FASTA_FILENAME + ".gz"
FULL_URL = BASE_URL + FTP_DIR + GZ_FILENAME

def load_peptide_counts(peptide_lengths = [8,9,10,11], nrows = None):
    """
    List of all reference peptides encoded in a reference human exome
    """
    lens = "_".join(str(n) for n in peptide_lengths)
    cache_filename = \
        "reference_peptide_counts_" + lens + "_nrows_" + str(nrows) + ".csv"

    def save_counts(src_path, dst_path):
        counts = _generate_counts(src_path, peptide_lengths, nrows)
        print "Saving %s" % dst_path
        counts.to_csv(dst_path)
        return counts

    return fetch_and_transform_data(
        transformed_filename = cache_filename,
        transformer = save_counts,
        loader = pd.read_csv,
        source_filename = FASTA_FILENAME,
        source_url = FULL_URL)

def load_peptide_set(peptide_lengths = [8,9,10,11], nrows = None):
    lens = "_".join(str(n) for n in peptide_lengths)
    cache_filename = \
        "reference_peptide_set_" + lens + "_nrows_" + str(nrows) + ".pickle.gz"

    def save_set(src_path, dst_path):
        string_set = _generate_set(src_path, peptide_lengths, nrows)
        with GzipFile(dst_path, 'w') as out_file:
            out_file.write(dumps(string_set))
        return string_set

    def load_set(path):
        result = None
        with GzipFile(path, 'r') as in_file:
            result = loads(in_file.read())
        return result

    return fetch_and_transform_data(
        transformed_filename = cache_filename,
        transformer = save_set,
        loader = load_set,
        source_filename = FASTA_FILENAME,
        source_url = FULL_URL)

