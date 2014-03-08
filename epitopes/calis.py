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


"""
Datasets from Calis, et al. "Properties of MHC Class I Presented Peptides That Enhance Immunogenicity"
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003266#s5
"""
from os.path import join

import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

from base import DATA_DIR

S1_FILE = join(DATA_DIR, 'calis_s1.csv')
S2_FILE = join(DATA_DIR, 'calis_s2.csv')

def load_s1():
    """
    Immunogenic and non-immunogenic pMHCs that were found in the IEDB,
    Vaccinia, Arena and Coxiella data sets.
    """
    return pd.read_csv(S1_FILE)

def load_s1_classes():
    data = load_s1()
    imm = data[ (data.Immunogenicity == 'immunogenic') & (data.Species == 'Mus')].Peptide
    non = data[ (data.Immunogenicity == 'non-immunogenic') & (data.Species == 'Mus')].Peptide
    return imm, non

def load_s1_ngrams():
    imm, non = load_s1_classes()

def load_s2():
    """
    Non-redundant murine and human Dengue epitopes and non-epitopes
    """
    data = pd.read_csv(S2_FILE)

def load_s2_classes():
    data = load_s2()
    human = data.host == 'Homo'
    pos = data['epitope status'] == 'epitope'
    neg = data['epitope status'] == 'non-epitope'
    imm = data[human & pos].peptide
    non = data[human & neg].peptide
    return imm, non

def load_s2_ngrams():
    imm, non = load_s2_classes()

def load_counts(normalize_row = False, redundant = False, max_ngram = 1):
    c = CountVectorizer(analyzer='char', ngram_range=(1, max_ngram), dtype=np.float)
    if redudant:
        imm, non = load_s1()
    else:
        imm, non = load_s2()
    imm, non = load_csv()
    total = list(imm) + list(non)
    X = c.fit_transform(total)
    if normalize_row:
        X = normalize(X, norm='l1')
    Y = np.ones(len(total), dtype='bool')
    Y[len(imm):] = 0
    return X, Y
