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

def load_s1(imm_file = S1_FILE):
    """
    Immunogenic and non-immunogenic pMHCs that were found in the IEDB, Vaccinia, Arena and Coxiella data sets
    """
    data = pd.read_csv(imm_file)
    imm = data[ (data.Immunogenicity == 'immunogenic') & (data.Species == 'Mus')].Peptide
    non = data[ (data.Immunogenicity == 'non-immunogenic') & (data.Species == 'Mus')].Peptide
    return imm, non

def load_s2(imm_file = S2_FILE):
    """
    Non-redundant murine and human Dengue epitopes and non-epitopes
    """
    data = pd.read_csv(imm_file)
    imm = data[ (data['epitope status'] == 'epitope') & (data.host == 'Homo')].peptide
    non = data[ (data['epitope status'] == 'non-epitope') & (data.host == 'Homo')].peptide
    return imm, non

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
