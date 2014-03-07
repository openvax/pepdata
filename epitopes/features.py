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


import numpy as np
from amino_acid import peptide_to_indices

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize


def make_ngram_dataset(
        imm, non = [],
        max_ngram = 1,
        normalize_row = True,
        reduced_alphabet = None,
        subsample_bigger_class = False,
        return_transformer = False,
        verbose = True):

    if reduced_alphabet is None:
        preprocessor = None
    else:
        preprocessor = \
            lambda s: ''.join([chr(48 + reduced_alphabet[char]) for char in s])

    c = CountVectorizer(analyzer='char',
                        ngram_range=(1,max_ngram),
                        dtype=np.float,
                        preprocessor = preprocessor)

    total = list(imm) + list(non)
    # returns a sparse matrix
    X = c.fit_transform(total).todense()

    n_imm = len(imm)
    n_non = len(non)

    if subsample_bigger_class:
        X_true = X[:n_imm]
        X_false = X[n_imm:]
        n_min = min(n_imm, n_non)
        idx = np.arange(n_imm)
        np.random.shuffle(idx)
        X_true = X_true[idx[:n_min]]
        idx = np.arange(n_non)
        np.random.shuffle(idx)
        X_false = X_false[idx[:n_min]]
        X = np.vstack([X_true, X_false])
        Y = np.ones(2*n_min, dtype='bool')
        Y[n_min:] = False
        if verbose:
            print "Rebalancing %d, %d -> %d" % (n_imm, n_non, n_min)
    else:
        Y = np.ones(len(total), dtype='bool')
        Y[n_imm:] = 0

    if reduced_alphabet and verbose:
        print "Alphabet", c.get_feature_names()
    if normalize_row:
        X = normalize(X, norm='l1')
    if verbose:
        print "Dataset size", X.shape
    if return_transformer:
        def transform(test_strings):
            X_test = c.transform(test_strings).todense()
            if normalize_row:
                X_test = normalize(X_test, norm='l1')
            return X_test
        return X, Y, transform
    else:
        return X, Y



import amino_acid

amino_acid_features = [
    amino_acid.hydropathy,
    amino_acid.volume,
    amino_acid.pK_side_chain,
    amino_acid.polarity,
    amino_acid.prct_exposed_residues,
    amino_acid.hydrophilicity,
    amino_acid.accessible_surface_area,
    amino_acid.local_flexibility,
    amino_acid.accessible_surface_area_folded,
    amino_acid.refractivity
]

def transform_rows(
        X,
        fns = amino_acid_features,
        positions = None,
        mean = False,
        pairwise_ratios = False):
    X2 = []
    for x in X:
        row = []
        for fn in fns:
            row_entries = [fn(xi) for xi in x]
            if positions:
                row_entries = \
                    [xi
                     for i, xi in enumerate(row_entries)
                     if i in positions]
            if mean:
                row.append(np.mean(row_entries))
            else:
                row.extend(row_entries)
            if pairwise_ratios:
                for i, y in enumerate(row_entries):
                    for j, z in enumerate(row_entries):
                        if i < j:
                            if z == 0:
                                ratio = 0
                            else:
                                ratio = y / z
                            row.append(ratio)
        X2.append(row)
    X2 = np.array(X2)
    return X2

import toxin
def toxin_features(imm, non = [], substring_length=3, positional=False):
    if positional:
        imm_toxin_features = \
            toxin.positional_toxin_features(imm, length=substring_length)
        non_toxin_features = \
            toxin.positional_toxin_features(non, length=substring_length)
    else:
        imm_toxin_features = toxin.toxin_features(imm, length=substring_length)
        non_toxin_features = toxin.toxin_features(non, length=substring_length)

    X = np.vstack([imm_toxin_features, non_toxin_features])
    Y = np.ones(len(X), dtype='bool')
    Y[len(imm):] = 0
    return X, Y
