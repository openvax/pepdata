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

import numpy as np

from .amino_acid import peptide_to_indices
from .peptide_vectorizer import PeptideVectorizer
from . import (toxin, amino_acid)

def make_ngram_dataset(
        imm,
        non=[],
        max_ngram=1,
        normalize_row=True,
        subsample_bigger_class=False,
        reduced_alphabet=None,
        training_already_reduced=False,
        return_transformer=False,
        verbose=True):
    """
    Given a set of positive amino acid string examples and a set
    of negative examples, transform them into an array of n-gram counts
    and a vector of target outputs.

    Parameters
    -----------
    max_ngram : int
        Degree of ngrams

    normalize_row : bool
        Return frequencies (default) or counts

    subsample_bigger_class : bool
        Default False

    reduced_alphabet : dictionary, optional
        Remap amino acid characters into some alternative alphabet

    training_already_reduced : bool
        Only apply reduced alphabet to testing set (training data is already
        transformed)

    return_transformer : bool
        In addition to the data X, Y, also return a feature transformer
        (default False)

    verbose : bool
        Print debug messages (default True)
    """

    imm = list(imm)
    non = list(non)
    n_imm = len(imm)
    n_non = len(non)

    # if class sizes are unequal
    # may want to subsample the larger
    # class so that training set becomes balanced
    if subsample_bigger_class:
        n_min = min(n_imm, n_non)
        idx = np.arange(n_imm)
        np.random.shuffle(idx)
        imm = [imm[i] for i in idx[:n_min]]
        idx = np.arange(n_non)
        np.random.shuffle(idx)
        non = [non[i] for i in idx[:n_min]]
        n_imm = n_min
        n_non = n_min

    combined = imm + non

    V = PeptideVectorizer(
        max_ngram=max_ngram,
        normalize_row=normalize_row,
        reduced_alphabet=reduced_alphabet,
        training_already_reduced=training_already_reduced)

    X = V.fit_transform(combined)
    Y = np.ones(len(combined), dtype='bool')
    Y[n_imm:] = 0

    if verbose:
        print("Dataset size", X.shape)

    if return_transformer:
        return X, Y, V
    else:
        return X, Y

def make_unlabeled_ngram_dataset(strings, **kwargs):
    # reuse the vectorizer for labeled dataset but give
    # an empty list for the negative class
    # then throw away the label array Y
    kwargs['subsample_bigger_class'] = False
    result = make_ngram_dataset(strings, [], **kwargs)
    # expect to get back either just the labeled data (X, Y)(
    # or (X, Y, vectorizer)
    if len(result) == 2:
        return result[0]
    else:
        assert len(result) == 3
        return result[0], result[2]

def _split_ngram_kwargs(kwargs):
    """
    Functions which call both a data loader and then vectorizer
    need to split their kwargs into two dictionaries
    """
    d = {
        'max_ngram': kwargs.pop('max_ngram', 1),
        'normalize_row': kwargs.pop('normalize_row', True),
        'subsample_bigger_class': kwargs.pop('subsample_bigger_class', False),
        'return_transformer': kwargs.pop('return_transformer', False),
        'reduced_alphabet': kwargs.get('reduced_alphabet'),
        'training_already_reduced':
            kwargs.pop('training_already_reduced', False),
        'verbose': kwargs.get('verbose')
    }
    return kwargs, d

def make_ngram_dataset_from_args(loader, *args, **kwargs):
    """
    Parameters
    -----------
    max_ngram : int
        Degree of ngrams

    normalize_row : bool
        Return frequencies (default) or counts

    subsample_bigger_class : bool
        Default False

    return_transformer : bool
        Default False

    reduced_alphabet : dictionary, optional
        Remap amino acid characters into some alternative alphabet

    training_already_reduced : bool
        Only apply reduced alphabet to testing set
        (training data is already transformed)

    verbose : bool
        Default True

    loader : fn
        Takes all the remaining arguments and returns two sets
        of positive and negative samples
    """
    loader_kwargs, kwargs = _split_ngram_kwargs(kwargs)
    pos, neg = loader(*args, **loader_kwargs)
    return make_ngram_dataset(pos, neg, **kwargs)

def make_unlabeled_ngram_dataset_from_args(loader, *args, **kwargs):
    loader_kwargs, kwargs = _split_ngram_kwargs(kwargs)
    strings = loader(*args, **loader_kwargs)
    return make_unlabeled_ngram_dataset(strings, **kwargs)


def array_from_kmers(strings):
    result = [peptide_to_indices(s) for s in strings]
    return np.array(result)


def make_kmer_dataset(imm_strings, non_strings, verbose=True):
    X_imm = array_from_kmers(imm_strings)
    X_non = array_from_kmers(non_strings)
    X = np.vstack([X_imm, X_non])
    Y = np.ones(len(X), dtype='bool')
    Y[len(X_imm):] = 0
    if verbose:
        print("[make_kmer_dataset] X.shape = %s" % (X.shape,))
    return X, Y

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
        fns=amino_acid_features,
        positions=None,
        mean=False,
        pairwise_ratios=False):
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

def toxin_features(imm, non=[], substring_length=3, positional=False):
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
