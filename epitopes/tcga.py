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
import numpy as np
import Bio.SeqIO
from progressbar import ProgressBar

from common import int_or_seq, dataframe_from_counts
from tcga_sources import TCGA_SOURCES, REFSEQ_PROTEIN_URL
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

def _load_maf_files(sources_dict, cancer_type = None):
    """
    Given a dictionary mapping cancer types to download urls,
    get all the source MAFs, load them as DataFrames, and then
    concatenate into a single DataFrame
    """
    data_frames = []
    if cancer_type is None:
        cancer_types = sources_dict.keys()
    elif isinstance(cancer_type, str):
        cancer_types = [cancer_type]
    else:
        assert isinstance(cancer_type, list), \
            "Cancer type must be None, str, or list but got %s" % cancer_type
        cancer_types = cancer_type

    for key in cancer_types:
        assert key in sources_dict, "Unknown cancer type %s" % key
        maf_url = sources_dict[key]
        maf_filename = key + ".maf"
        path = fetch_data(maf_filename, maf_url)
        df = open_maf(path)
        df['Cancer Type'] = key
        data_frames.append(df)
    return pd.concat(data_frames, ignore_index = True)

def load_dataframe(cancer_type = None):
    return _load_maf_files(TCGA_SOURCES, cancer_type)

def _build_refseq_id_to_protein(refseq_path, predicted_proteins = False):
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
                if name.startswith("NP_") or name.startswith("YP_"):
                    before_dot = name.split('.')[0]
                    assert before_dot not in result, \
                        "Unexpected refseq ID repeat %s" % before_dot
                    result[before_dot] = protein
                    result[name] = protein
                elif name.startswith("XP_"):
                    if predicted_proteins:
                        before_dot = name.split('.')[0]
                        assert before_dot not in result, \
                            "Unexpected refseq ID repeat %s" % before_dot
                        result[before_dot] = protein
                        result[name] = protein
                else:
                    print "Unexpected refseq ID", name
            except IndexError:
                pass
    return result

# patterns for how protein changes are encoded i.e. p.Q206E
SINGLE_AMINO_ACID_SUBSTITUTION = "p.([A-Z])([0-9]+)([A-Z])"
DELETION = "p.([A-Z])([0-9]+)del"

def load_peptide_counts(
        peptide_length = [8,9,10,11],
        cancer_type = None,
        verbose = False):
    """
    Call given functions for mutation type
        subst_fn(protein, position, mutant_aa)
    """
    peptide_lengths = int_or_seq(peptide_length)

    combined_df = load_dataframe(cancer_type = cancer_type)
    filtered = combined_df[["Refseq_prot_Id", "Protein_Change"]].dropna()

    # discard indices of NA entries
    filtered.index = np.arange(len(filtered))
    refseq_path = fetch_data('refseq_protein.faa', REFSEQ_PROTEIN_URL)
    refseq_ids_to_protein = _build_refseq_id_to_protein(refseq_path)
    refseq_ids = filtered.Refseq_prot_Id

    n_failed = 0
    peptide_counts = {}

    subst_matches = \
        filtered.Protein_Change.str.extract(SINGLE_AMINO_ACID_SUBSTITUTION)

    # drop non-matching rows
    subst_matches = subst_matches.dropna()

    #
    # generate substrings from all single amino acid substitutions
    #
    subst_count = len(subst_matches)
    print "Simple substitutions (%d)" % subst_count
    pbar = ProgressBar(maxval = np.array(subst_matches.index).max() ).start()
    for match_num, wildtype, str_position, mutation in \
            subst_matches.itertuples():
        if wildtype == mutation:
            # silent mutation, skipping
            continue
        refseq_id = refseq_ids[match_num]
        protein = refseq_ids_to_protein.get(refseq_id)
        if protein is None:
            if verbose:
                print "Couldn't find refseq ID %s" % refseq_id
            n_failed += 1
            continue
        amino_acid_pos = int(str_position) - 1
        if len(protein) <= amino_acid_pos:
            if verbose:
                print "Protein %s too short, needed position %s but len %d" % \
                    (refseq_id, str_position, len(protein))
            n_failed += 1
            continue
        old_aa = protein[amino_acid_pos]
        if old_aa != wildtype:
            if verbose:
                print "Expected %s but got %s at position %s in %s (%s%s%s)" % \
                    (wildtype,
                    old_aa,
                    str_position,
                    refseq_id,
                    wildtype,
                    str_position,
                    mutation)
            n_failed += 1
            continue

        m = len(protein)

        for n in peptide_lengths:
            for i in xrange(n):
                start = amino_acid_pos - i
                stop = start + n
                if start >= 0 and stop <= m:
                    substr = \
                        protein[start:start+i] + \
                        mutation + \
                        protein[start+i+1:stop]
                    if substr in peptide_counts:
                        peptide_counts[substr] += 1
                    else:
                        peptide_counts[substr] = 1
        pbar.update(match_num)
    pbar.finish()


    del_matches = \
        filtered.Protein_Change.str.extract(DELETION)

    # drop non-matching rows
    del_matches = del_matches.dropna()

    #
    # generate substrings from all single amino acid substitutions
    #
    del_count = len(del_matches)
    print "Simple deletions (%d)" % del_count
    pbar = ProgressBar(maxval = np.array(del_matches.index).max()).start()
    for match_num, wildtype, str_position in \
            del_matches.itertuples():
        refseq_id = refseq_ids[match_num]
        protein = refseq_ids_to_protein.get(refseq_id)
        if protein is None:
            if verbose:
                print "Couldn't find refseq ID %s" % refseq_id
            n_failed += 1
            continue
        amino_acid_pos = int(str_position) - 1
        if len(protein) <= amino_acid_pos:
            if verbose:
                print "Protein %s too short, needed position %s but len %d" % \
                    (refseq_id, str_position, len(protein))
            n_failed += 1
            continue
        old_aa = protein[amino_acid_pos]
        if old_aa != wildtype:
            if verbose:
                print \
                    "Expected %s but got %s at position %s in %s (%s%sdel)" % \
                    (wildtype,
                    old_aa,
                    str_position,
                    refseq_id,
                    wildtype,
                    str_position,
                    )
            n_failed += 1
            continue

        m = len(protein)

        for n in peptide_lengths:
            for i in xrange(n):
                start = amino_acid_pos - i
                stop = start + n + 1
                if start >= 0 and stop <= m:
                    substr = \
                        protein[start:start+i] + \
                        protein[start+i+1:stop]
                    if substr in peptide_counts:
                        peptide_counts[substr] += 1
                    else:
                        peptide_counts[substr] = 1
        pbar.update(match_num)
    pbar.finish()
    print "Failed to parse %d / %d mutations" % (n_failed, len(filtered))
    return dataframe_from_counts(peptide_counts)

def load_peptide_set(
        peptide_length = [8,9,10,11],
        cancer_type = None):
    counts = load_peptide_counts(
        peptide_length = peptide_length,
        cancer_type = cancer_type)
    return set(counts.Peptide)