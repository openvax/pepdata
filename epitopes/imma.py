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


from os.path  import join

import pandas as pd
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize

from base import DATA_DIR

IMMA2_IMM_FILE = join(DATA_DIR, 'IMMA2_imm.txt')
IMMA2_NON_FILE = join(DATA_DIR, 'IMMA2_non.txt')

def load_imma2(imm_file = IMMA2_IMM_FILE, non_file = IMMA2_NON_FILE):
    with open(join(data_dir, imm_file)) as f:
        imm = [line.strip() for line in f]
    with open(join(data_dir, non_file)) as f:
        non = [line.strip() for line in f]
    return imm, non

def load_ngrams(normalize = False, max_ngram = 1):
    c = CountVectorizer(analyzer='char', ngram_range=(1,max_ngram), dtype=np.float)
    imm, non = load_csv()
    total = list(imm) + list(non)
    X = c.fit_transform(total)
    if normalize:
        X = normalize(X, norm='l1')
    Y = np.ones(len(total), dtype='bool')
    Y[len(imm):] = 0
    return X, Y
