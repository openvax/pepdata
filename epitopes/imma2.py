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

from base import DATA_DIR
from features import make_ngram_dataset_from_args, make_alphabet_transformer

IMMA2_IMM_FILE = join(DATA_DIR, 'IMMA2_imm.txt')
IMMA2_NON_FILE = join(DATA_DIR, 'IMMA2_non.txt')

def load_imma2(reduced_alphabet = None):
    if reduced_alphabet:
        transformer = make_alphabet_transformer(reduced_alphabet)
    else:
        transformer = lambda x: x
    with open(IMMA2_IMM_FILE) as f:
        imm = [transformer(line.strip()) for line in f]
    with open(IMMA2_NON_FILE) as f:
        non = [transformer(line.strip()) for line in f]
    return set(imm), set(non)


def load_imma2_ngrams(*args, **kwargs):
    """
    Load IMMA2 dataset and transform into n-gram vector representation

    Parameters
    -----------
    max_ngram : int
        Degree of ngrams

    normalize_row : bool
        Return frequencies (default) or counts

    subsample_bigger_class : bool
        Default False

    return_transformer : bool
        Default False

    verbose : bool
        Default True
    """
    return make_ngram_dataset_from_args(load_imma2, *args, **kwargs)



