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

from __future__ import print_function, division, absolute_import
from os.path import join

from .features import make_ngram_dataset_from_args
from .reduced_alphabet import make_alphabet_transformer
from .static_data import DATA_DIR

def load_imm_list(reduced_alphabet=None):
    path = join(DATA_DIR, 'IMMA2_imm.txt')
    if reduced_alphabet:
        transformer = make_alphabet_transformer(reduced_alphabet)
    else:
        def transformer(x):
            return x

    with open(path) as f:
        imm = [transformer(line.strip()) for line in f]
    return imm

def load_non_list(reduced_alphabet=None):
    path = join(DATA_DIR, 'IMMA2_non.txt')
    if reduced_alphabet:
        transformer = make_alphabet_transformer(reduced_alphabet)
    else:
        def transformer(x):
            return x

    with open(path) as f:
        non = [transformer(line.strip()) for line in f]
    return non

def load_classes(reduced_alphabet=None):
    imm = load_imm_list(reduced_alphabet)
    non = load_non_list(reduced_alphabet)
    return set(imm), set(non)


def load_ngrams(*args, **kwargs):
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
    return make_ngram_dataset_from_args(load_classes, *args, **kwargs)
