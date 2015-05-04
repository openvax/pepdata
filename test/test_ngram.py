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

from pepdata import features

# isoforms of two different proteins a, b

a1 = \
    "MSPHPTALLGLVLCLAQTIHTQEEDLPRPSISAEPGTVIPLGSHVTFVCRGPVGVQTFRLERESRSTYNDTEDVSQASPSESEARFRIDSVSEGNAGPYRCIYYKPPKWSEQSDYLELLVKETSGGPDSPDTEPGSSAGPTQRPSDNSHNEHAPASQGLKAEHLYILIGVSVVFLFCLLLLVLFCLHRQNQIKQGPPRSKDEEQKPQQRPDLAVDVLERTADKATVNGLPEKDRETDTSALAAGSSQEVTYAQLDHWALTQRTARAVSPQSTKPMAESITYAAVARH"

a2 = \
    "MSLMVVSMACVGFFLLQGAWPHEGVHRKPSLLAHPGPLVKSEETVILQCWSDVRFEHFLLHREGKYKDTLHLIGEHHDGVSKANFSIGPMMQDLAGTYRCYGSVTHSPYQLSAPSDPLDIVITGLYEKPSLSAQPGPTVLAGESVTLSCSSRSSYDMYHLSREGEAHERRFSAGPKVNGTFQADFPLGPATHGGTYRCFGSFRDSPYEWSNSSDPLLVSVTGNPSNSWPSPTEPSSKTGNPRHLHVLIGTSVVKIPFTILLFFLLHRWCSNKKNAAVMDQEPAGNRTVNSEDSDEQDHQEVSYA"

a3 = 'MSLMVVSMACVGFFLLEGPWPHVGGQDKPFLSAWPGTVVSEGQHVTLQCRSRLGFNEFSLSKEDGMPVPELYNRIFRNSFLMGPVTPAHAGTYRCCSSHPHSPTGWSAPSNPVVIMVTGVHRKPSLLAHPGPLVKSEETVILQCWSDVRFEHFLLHREGKYKDTLHLIGEHHDGVSKANFSIGPMMQDLAGTYRCYGSVTHSPYQLSAPSDPLDIVITGLYEKPSLSAQPGPTVLAGESVTLSCSSRSSYDMYHLSREGEAHERRFSAGPKVNGTFQADFPLGPATHGGTYRCFGSFRDSPYEWSNSSDPLLVSVTAFLSVKSSGHKYIY'

A = [a1, a2, a3]

b1 = \
    'MPKGRAGSLPTTSIGWRFQLWFLGLTCPERHLARRLKNNSFYPFVQQEPNVFVLEYYLDTLWKGMLLFIISVVLVSFSSLREVQKQETWVFLVYGVGVGLWLVISSLPRRRLVLNHTRGVYHFSIQGRTVCQGPLHLVYVRLALSSDAHGRCFFHLVLGGHRLEPLVLVQLSEHYEQMEYLGRYIARKLNINYFDYLATSYRHVVRHWPPPGAGTVMGKSPMGHKPSSSQSSLEV'

b2 = \
    'MPKGRAGSLPTTSIGWRFQLWFLGLTCPERHLARRLKNNSFYPFVQQEPNVFVLEYYLDTLWKGMLLFIISVVLVSFSSLREVQKQETWVFLVYGVGVGLWLVISSLPRRRLVLNHTRGVYHFSIQGRTVCQGPLHLVYVRLALSSDAHGRCFFHLVLGGHRLEPLVLVQLSEHYEQMEYLGRYIARKLNINYFDYLATSYRHVVRHWPPPGAGTVMGKSPMGHKPSSSQSSLEV'

B = [b1, b2]

def test_labeled_unigrams():
    """
    Unigram amino acid counts should result in length 20 vectors
    """
    n = len(A) + len(B)
    X, Y = features.make_ngram_dataset(A, B, max_ngram=1)
    assert len(X) == n, "Expected %d vectors but got %d" % (n, len(X))
    assert len(X.shape) == 2
    assert len(Y.shape) == 1
    assert Y.sum() == len(A), "Wrong number of positive labels: %d" % Y.sum()
    assert X.shape[1] == 20, \
        "Expected feature vectors to be of length 20, got %d" % X.shape[1]

def test_unlabeled_unigrams():
    """
    Unigram amino acid counts should result in length 20 vectors
    """
    n = len(A)
    X = features.make_unlabeled_ngram_dataset(A, max_ngram=1)
    assert len(X) == n, "Expected %d vectors but got %d" % (n, len(X))
    assert len(X.shape) == 2
    assert X.shape[1] == 20, \
        "Expected feature vectors to be of length 20, got %d" % X.shape[1]

def test_pickle_ngram_vectorizer():
    X, V = features.make_unlabeled_ngram_dataset(
        A,
        max_ngram=1,
        return_transformer=True)
    s = cPickle.dumps(V)
    V2 = cPickle.loads(s)
    X2 = V2.transform(A)
    assert (X == X2).all()
