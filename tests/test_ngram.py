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

import pickle

import numpy as np
import pytest

from pepdata import PeptideVectorizer
from pepdata.reduced_alphabet import hp2

# isoforms of two different proteins a, b

a1 = (
    "MSPHPTALLGLVLCLAQTIHTQEEDLPRPSISAEPGTVIPLGSHVTFVCRGPVGVQTFRLERESRSTYND"
    "TEDVSQASPSESEARFRIDSVSEGNAGPYRCIYYKPPKWSEQSDYLELLVKETSGGPDSPDTEPGSSAGPT"
    "QRPSDNSHNEHAPASQGLKAEHLYILIGVSVVFLFCLLLLVLFCLHRQNQIKQGPPRSKDEEQKPQQRPDL"
    "AVDVLERTADKATVNGLPEKDRETDTSALAAGSSQEVTYAQLDHWALTQRTARAVSPQSTKPMAESITYAA"
    "VARH"
)

a2 = (
    "MSLMVVSMACVGFFLLQGAWPHEGVHRKPSLLAHPGPLVKSEETVILQCWSDVRFEHFLLHREGKYKDTLH"
    "LIGEHHDGVSKANFSIGPMMQDLAGTYRCYGSVTHSPYQLSAPSDPLDIVITGLYEKPSLSAQPGPTVLAG"
    "ESVTLSCSSRSSYDMYHLSREGEAHERRFSAGPKVNGTFQADFPLGPATHGGTYRCFGSFRDSPYEWSNSS"
    "DPLLVSVTGNPSNSWPSPTEPSSKTGNPRHLHVLIGTSVVKIPFTILLFFLLHRWCSNKKNAAVMDQEPAG"
    "NRTVNSEDSDEQDHQEVSYA"
)

a3 = (
    "MSLMVVSMACVGFFLLEGPWPHVGGQDKPFLSAWPGTVVSEGQHVTLQCRSRLGFNEFSLSKEDGMPVPEL"
    "YNRIFRNSFLMGPVTPAHAGTYRCCSSHPHSPTGWSAPSNPVVIMVTGVHRKPSLLAHPGPLVKSEETVIL"
    "QCWSDVRFEHFLLHREGKYKDTLHLIGEHHDGVSKANFSIGPMMQDLAGTYRCYGSVTHSPYQLSAPSDPL"
    "DIVITGLYEKPSLSAQPGPTVLAGESVTLSCSSRSSYDMYHLSREGEAHERRFSAGPKVNGTFQADFPLGP"
    "ATHGGTYRCFGSFRDSPYEWSNSSDPLLVSVTAFLSVKSSGHKYIY"
)

A = [a1, a2, a3]

b1 = (
    "MPKGRAGSLPTTSIGWRFQLWFLGLTCPERHLARRLKNNSFYPFVQQEPNVFVLEYYLDTLWKGMLLFII"
    "SVVLVSFSSLREVQKQETWVFLVYGVGVGLWLVISSLPRRRLVLNHTRGVYHFSIQGRTVCQGPLHLVYV"
    "RLALSSDAHGRCFFHLVLGGHRLEPLVLVQLSEHYEQMEYLGRYIARKLNINYFDYLATSYRHVVRHWPP"
    "PGAGTVMGKSPMGHKPSSSQSSLEV"
)

b2 = (
    "MPKGRAGSLPTTSIGWRFQLWFLGLTCPERHLARRLKNNSFYPFVQQEPNVFVLEYYLDTLWKGMLLFII"
    "SVVLVSFSSLREVQKQETWVFLVYGVGVGLWLVISSLPRRRLVLNHTRGVYHFSIQGRTVCQGPLHLVYV"
    "RLALSSDAHGRCFFHLVLGGHRLEPLVLVQLSEHYEQMEYLGRYIARKLNINYFDYLATSYRHVVRHWPPP"
    "GAGTVMGKSPMGHKPSSSQSSLEV"
)

B = [b1, b2]


@pytest.mark.parametrize("normalize_row", [True, False])
def test_fit_transform_returns_one_array_row_per_sequence(normalize_row):
    vectorizer = PeptideVectorizer(max_ngram=2, normalize_row=normalize_row)
    X = vectorizer.fit_transform(A + B)
    assert type(X) is np.ndarray
    assert X.shape == (len(A + B), len(vectorizer.count_vectorizer.vocabulary_))
    Y = vectorizer.transform(A + B)
    assert type(Y) is np.ndarray
    np.testing.assert_allclose(Y, X)


def test_default_normalize_row_gives_row_frequencies():
    # The default path raised "np.matrix is not supported" on current
    # scikit-learn (#62).
    vectorizer = PeptideVectorizer()
    X = vectorizer.fit_transform(["SIINFEKL", "GILGFVFTL"])
    assert X.shape == (2, 10)
    np.testing.assert_allclose(X.sum(axis=1), [1.0, 1.0])
    column_i = vectorizer.count_vectorizer.vocabulary_["i"]
    np.testing.assert_allclose(X[:, column_i], [2 / 8, 1 / 9])


def test_normalize_row_false_gives_counts():
    vectorizer = PeptideVectorizer(normalize_row=False)
    X = vectorizer.fit_transform(["SIINFEKL", "GILGFVFTL"])
    np.testing.assert_array_equal(X.sum(axis=1), [8, 9])


@pytest.mark.parametrize("normalize_row, expected_total", [(False, 13.0), (True, 1.0)])
def test_transform_counts_only_fitted_ngrams(normalize_row, expected_total):
    # SIINFEKA has 15 1- and 2-grams; "A" and "KA" were not seen during fitting.
    vectorizer = PeptideVectorizer(max_ngram=2, normalize_row=normalize_row)
    vectorizer.fit(["SIINFEKL", "SIINFEKV"])
    X = vectorizer.transform(["SIINFEKA"])
    assert X.shape == (1, 16)
    assert X.sum() == pytest.approx(expected_total)


def test_fit_transform_on_already_reduced_training_sequences():
    reduced = ["".join(hp2[residue] for residue in sequence) for sequence in A + B]
    vectorizer = PeptideVectorizer(reduced_alphabet=hp2, training_already_reduced=True)
    X = vectorizer.fit_transform(reduced)
    assert type(X) is np.ndarray
    assert X.shape == (len(A + B), 2)
    np.testing.assert_allclose(X.sum(axis=1), np.ones(len(A + B)))


def test_pickled_vectorizer_transforms_identically():
    vectorizer = PeptideVectorizer(max_ngram=2)
    vectorizer.fit(A)
    restored = pickle.loads(pickle.dumps(vectorizer))
    np.testing.assert_allclose(restored.transform(B), vectorizer.transform(B))
