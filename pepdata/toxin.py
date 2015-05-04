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
from os.path import join

import numpy as np

from .static_data import DATA_DIR

TOXIN_TABLE_FILENAME = join(DATA_DIR, 'Toxin_Protein_Table.txt')
TOXIN_LIST_FILENAME = join(DATA_DIR, 'toxins.txt')

def read_toxins_from_table(filename=TOXIN_TABLE_FILENAME):
    f = open(filename)
    result = [line.split('\t')[29] for line in f]
    f.close()
    return result

def gen_toxin_list(input_filename=TOXIN_TABLE_FILENAME,
                   output_filename=TOXIN_LIST_FILENAME):
    proteins = read_toxins_from_table(input_filename)
    with open(output_filename, 'w') as f:
        for protein in proteins:
            f.write(protein)
            f.write('\n')
    return proteins

def read_toxin_list(filename=TOXIN_LIST_FILENAME):
    f = open(filename)
    result = f.read().splitlines()
    f.close()
    return result


def count_toxin_substring_matches(peptides, toxins=None, length=3):
    if toxins is None:
        toxins = read_toxin_list()
    total_hits = 0
    total_comparisons = 0
    toxin_count = 0
    for t in toxins:
        toxin_substrings = set([])
        for i in xrange(0, len(t) - length):
            toxin_substrings.add(t[i:i + length])
        for p in peptides:
            for i in xrange(0, len(p) - length):
                substr = p[i:i + length]
                if substr in toxin_substrings:
                    print("Peptide %s's substring %s matches!" % (p, substr))
                    total_hits += 1
                    incr_toxin_count = True
                total_comparisons += 1
        toxin_count += incr_toxin_count
    print("Total hits", total_hits)
    print("Total comparisons", total_comparisons)
    print("Total toxins w/ hits", toxin_count, "/", len(toxins))

def toxin_features(peptides, toxins=None, length=3, reverse=False):
    """
    Unordered binary features indicating whether a substring of a peptide
    occurs *anywhere* in the sequence of some known toxin
    """
    if toxins is None:
        toxins = read_toxin_list()
    substr_sets = []
    for t in toxins:
        toxin_substrings = set([])
        for i in xrange(0, len(t) - length):
            substr = t[i:i + length]
            toxin_substrings.add(substr)
            if reverse:
                toxin_substrings.add(substr[::-1])
        substr_sets.append(toxin_substrings)
    X = np.zeros((len(peptides), len(toxins))).astype('int')
    for i, p in enumerate(peptides):
        for pos in xrange(0, 9 - length):
            substr = p[pos:pos + length]
            for j, substr_set in enumerate(substr_sets):
                if substr in substr_set:
                    X[i, j] = 1
    return X

def positional_toxin_features(peptides, toxins=None, length=3, reverse=True):
    """
    For each substring of the peptide, count how many toxins this occurs in
    """
    if toxins is None:
        toxins = read_toxin_list()
    substr_sets = []
    for t in toxins:
        toxin_substrings = set([])
        for i in xrange(0, len(t) - length):
            substr = t[i:i + length]
            toxin_substrings.add(substr)
            if reverse:
                toxin_substrings.add(substr[::-1])
        substr_sets.append(toxin_substrings)

    X = np.zeros((len(peptides), 9 - length)).astype('int')
    for i, p in enumerate(peptides):
        for pos in xrange(0, 9 - length):
            substr = p[pos:pos + length]
            for substr_set in substr_sets:
                if substr in substr_set:
                    X[i, pos] += 1
    return X
