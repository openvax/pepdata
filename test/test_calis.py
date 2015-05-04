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

from pepdata import calis

def test_calis_s1():
    """
    Calis S1: Make sure humans are a subset of all species
    """
    df_all = calis.load_s1(human=False)
    df_human = calis.load_s1(human=True)
    assert len(df_all) > len(df_human)

def test_calis_s2():
    """
    Calis S2: Make sure humans are a subset of all species
    """
    df_all = calis.load_s2(human=False)
    df_human = calis.load_s2(human=True)
    assert len(df_all) > len(df_human)

def test_calis_s1_ngrams():
    """
    Calis S1: Number of samples in ngram dataset should be same as strings
    """
    X, Y = calis.load_s1_ngrams()
    imm, non = calis.load_s1_classes()
    assert len(X) == len(Y)
    assert len(X) == len(imm) + len(non)

def test_calis_s2_ngrams():
    """
    Calis S2: Number of samples in ngram dataset should be same as strings
    """
    X, Y = calis.load_s2_ngrams()
    imm, non = calis.load_s2_classes()
    assert len(X) == len(Y)
    assert len(X) == len(imm) + len(non)
