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

from pepdata import imma2

def test_imma2_no_overlap():
    """
    IMMA2: Make sure the immunogenic and non-immunogenic sets have no common
    strings
    """
    imm, non = imma2.load_classes()
    assert len(imm.intersection(non)) == 0

def test_imma2_ngrams_same_size():
    """
    IMMA2: The number of samples in the n-gram dataset should be the same
    as the original sets
    """
    X, Y = imma2.load_ngrams()
    assert len(X) == len(Y)
    imm, non = imma2.load_classes()
    assert len(X) == len(imm) + len(non)
