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

import re

import pandas as pd
import Bio.SeqIO
from tcga_sources import maf_urls, REFSEQ_PROTEIN_URL
from download import fetch_data

def open_maf(filename):
    """
    Load a TCGA MAF file into a Pandas DataFrame
    """
    with open(filename) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1
    return pd.read_csv(
        filename,
        skiprows=lines_to_skip,
        sep="\t",
        low_memory=False)

def load_maf_files(sources_dict):
    """
    Given a dictionary mapping cancer types to download urls,
    get all the source MAFs, load them as DataFrames, and then
    concatenate into a single DataFrame
    """
    data_frames = []
    for cancer_type, maf_url in sources_dict.iteritems():
        maf_filename = cancer_type + ".maf"
        path = fetch_data(maf_filename, maf_url)
        df = open_maf(path)
        df['Cancer Type'] = cancer_type
        data_frames.append(df)
    return pd.concat(data_frames)

def build_refseq_id_to_protein(refseq_path):
    """
    Given the path to a local FASTA file containing
    RefSeq ID's and their protein transcripts,
    build a dictionary from IDs to transcripts
    """
    result = {}
    with open(refseq_path, 'r') as f:
        for record in Bio.SeqIO.parse(f, "fasta"):
            protein = str(record.seq)
            try:
                name = record.id.split("|")[3]
                # TODO: handle multiple entries more intelligently than this.
                if name.startswith("NP_"):
                    before_dot = name.split('.')[0]
                    result[before_dot] = protein
            except IndexError:
                pass
    return result

SINGLE_AMINO_ACID_SUBSTITUTION = re.compile("p.([A-Z])([0-9]+)([A-Z])")

def load_tcga():
    combined_df = load_maf_files(maf_urls)
    filtered = combined_df[["Refseq_prot_Id", "Protein_Change"]].dropna()
    refseq_path = fetch_data('refseq_protein.faa', REFSEQ_PROTEIN_URL)
    refseq_id_to_protein = build_refseq_id_to_protein(refseq_path)
    return filtered
