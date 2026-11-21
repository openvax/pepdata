# Copyright (c) 2014-2016. Mount Sinai School of Medicine
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


from copy import copy

import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

def make_count_vectorizer(reduced_alphabet, max_ngram):
    if reduced_alphabet is None:
        preprocessor = None
    else:
        def preprocessor(s):
            return "".join([reduced_alphabet[si] for si in s])

    return CountVectorizer(
        analyzer='char',
        ngram_range=(1, max_ngram),
        dtype=np.float64,
        preprocessor=preprocessor)

class PeptideVectorizer(object):
    """
    Make n-gram frequency vectors from peptide sequences
    """
    def __init__(
            self,
            max_ngram=1,
            normalize_row=True,
            reduced_alphabet=None,
            training_already_reduced=False):
        self.reduced_alphabet = reduced_alphabet
        self.max_ngram = max_ngram
        self.normalize_row = normalize_row
        self.training_already_reduced = training_already_reduced
        self.count_vectorizer = None

    def __getstate__(self):
        return {
            'reduced_alphabet': self.reduced_alphabet,
            'count_vectorizer': self.count_vectorizer,
            'training_already_reduced': self.training_already_reduced,
            'normalize_row': self.normalize_row,
            'max_ngram': self.max_ngram,
        }

    def fit_transform(self, amino_acid_strings):
        """
        Learn the n-gram vocabulary of the given sequences and vectorize them.

        Parameters
        ----------
        amino_acid_strings : list of str
            Peptide or protein sequences.

        Returns
        -------
        numpy.ndarray
            Array of shape (len(amino_acid_strings), number of distinct
            n-grams) with n-gram counts, or per-row frequencies (each row
            divided by its total) when normalize_row is True.
        """
        # Already-reduced training must use the same tokenization later used
        # for reduced transform input. Keep its existing lowercase vocabulary
        # and feature-column order, including in saved vectorizers.
        alphabet = None if self.training_already_reduced else self.reduced_alphabet
        self.count_vectorizer = make_count_vectorizer(alphabet, self.max_ngram)
        X = self.count_vectorizer.fit_transform(amino_acid_strings).toarray()

        if self.normalize_row:
            X = normalize(X, norm='l1')
        return X

    def fit(self, amino_acid_strings):
        self.fit_transform(amino_acid_strings)

    def transform(self, amino_acid_strings):
        """
        Vectorize sequences with the vocabulary learned by fit.

        n-grams that weren't seen during fitting are not counted.
        Input sequences are unreduced, even when training_already_reduced was
        used during fitting.

        Parameters
        ----------
        amino_acid_strings : list of str
            Peptide or protein sequences.

        Returns
        -------
        numpy.ndarray
            Array of shape (len(amino_acid_strings), size of the fitted
            vocabulary), normalized per row as in fit_transform.
        """
        assert self.count_vectorizer, "Must call 'fit' before 'transform'"
        vectorizer = self.count_vectorizer
        if self.training_already_reduced:
            if self.reduced_alphabet is not None:
                amino_acid_strings = [
                    "".join(self.reduced_alphabet[residue] for residue in sequence)
                    for sequence in amino_acid_strings
                ]
            # Older saved models attached the raw-sequence preprocessor to a
            # lowercase vocabulary fitted by a different CountVectorizer.
            # Preserve that vocabulary without reducing twice or mutating it.
            if vectorizer.preprocessor is not None:
                vectorizer = copy(vectorizer)
                vectorizer.preprocessor = None
        X = vectorizer.transform(amino_acid_strings).toarray()
        if self.normalize_row:
            X = normalize(X, norm='l1')
        return X
