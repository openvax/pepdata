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
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

from .reduced_alphabet import make_alphabet_transformer

def make_count_vectorizer(reduced_alphabet, max_ngram):
    if reduced_alphabet is None:
        preprocessor = None
    else:
        preprocessor = make_alphabet_transformer(reduced_alphabet)

    return CountVectorizer(
        analyzer='char',
        ngram_range=(1, max_ngram),
        dtype=np.float,
        preprocessor=preprocessor)

class PeptideVectorizer(object):
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
        self.count_vectorizer = \
            make_count_vectorizer(self.reduced_alphabet, self.max_ngram)

        if self.training_already_reduced:
            c = make_count_vectorizer(None, self.max_ngram)
            X = c.fit_transform(amino_acid_strings).todense()
            self.count_vectorizer.vocabulary_ = c.vocabulary_
        else:
            c = self.count_vectorizer
            X = c.fit_transform(amino_acid_strings).todense()

        if self.normalize_row:
            X = normalize(X, norm='l1')
        return X

    def fit(self, amino_acid_strings):
        self.fit_transform(amino_acid_strings)

    def transform(self, amino_acid_strings):
        assert self.count_vectorizer, "Must call 'fit' before 'transform'"
        X = self.count_vectorizer.transform(amino_acid_strings).todense()
        if self.normalize_row:
            X = normalize(X, norm='l1')
        return X
