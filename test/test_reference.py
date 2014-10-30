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

import pandas as pd
from pepdata import reference

def test_reference_peptide_counts():
    counts = reference.load_peptide_counts(peptide_length = 8, nrows = 20)
    assert counts is not None
    assert isinstance(counts, pd.DataFrame)
    assert (counts.Peptide.str.len() == 8).all()

def test_reference_peptide_set():
    peptides = reference.load_peptide_set(peptide_length = 8, nrows = 20)
    assert peptides is not None
    assert isinstance(peptides, set)
    assert all(len(p) == 8 for p in peptides)
