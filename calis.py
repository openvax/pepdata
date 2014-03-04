"""
Datasets from Calis, et al. "Properties of MHC Class I Presented Peptides That Enhance Immunogenicity"
http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1003266#s5



"""

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
import pandas as pd
import numpy as np

def load_s1(imm_file = 'calis_s1.csv'):
    """
    Immunogenic and non-immunogenic pMHCs that were found in the IEDB, Vaccinia, Arena and Coxiella data sets
    """
    data = pd.read_csv(imm_file)
    imm = data[ (data.Immunogenicity == 'immunogenic') & (data.Species == 'Mus')].Peptide
    non = data[ (data.Immunogenicity == 'non-immunogenic') & (data.Species == 'Mus')].Peptide
    return imm, non

def load_s2(imm_file = 'calis_s2.csv'):
    """
    Non-redundant murine and human Dengue epitopes and non-epitopes
    """
    data = pd.read_csv(imm_file)
    imm = data[ (data['epitope status'] == 'epitope') & (data.host == 'Homo')].peptide
    non = data[ (data['epitope status'] == 'non-epitope') & (data.host == 'Homo')].peptide
    return imm, non

def load_counts(normalizeRow = False, redundant = False):
    c = CountVectorizer(analyzer='char', ngram_range=(1,2), dtype=np.float)
    if redudant:
        imm, non = load_s1()
    else:
        imm, non = load_s2()
    imm, non = load_csv()
    total = list(imm) + list(non)
    X = c.fit_transform(total)
    if normalizeRow:
        X = normalize(X, norm='l1')
    Y = np.ones(len(total), dtype='bool')
    Y[len(imm):] = 0
    return X, Y
