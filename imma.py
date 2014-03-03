from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
import pandas as pd
import numpy as np

def load_csv(imm_file = 'IMMA2_imm.txt', non_file = 'IMMA2_non.txt'):
  imm = [line.strip() for line in open(imm_file)]
  non = [line.strip() for line in open(non_file)]
  return imm, non

def load_dataset(normalize = False):
  c = CountVectorizer(analyzer='char', ngram_range=(1,2), dtype=np.float)
  imm, non = load_csv()
  total = list(imm) + list(non)
  X = c.fit_transform(total)
  if normalize:
  	X = normalize(X, norm='l1')
  Y = np.ones(len(total), dtype='bool')
  Y[len(imm):] = 0
  return X, Y
