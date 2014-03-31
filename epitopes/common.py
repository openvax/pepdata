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
import logging

bad_amino_acids = 'U|X|J|B|Z'

def int_or_seq(x):
    if isinstance(x, int):
        return [x]
    else:
        return list(x)

def dataframe_from_counts(counts):
    invert =  {
        'Peptide': counts.keys(),
        'Count' : counts.values()
    }

    df = pd.DataFrame(invert)
    return df.sort("Count", ascending=False)

def split_classes(
        values,
        noisy_labels = 'majority',
        unique_sequences = True,
        verbose = False):
    """
    Given a dataframe mapping epitope strings to percentages in [0,1],
    split them into two sets (positive and negative examples.

    Values for `noisy_labels`:
    - majority = epitope is Positive if majority of its results are Positive
    - negative = epitope is Negative if any result is Negative
    - positive = epitope is Positive if any result is Positive
    - drop = remove any epitopes with contradictory results
    - keep = leave contradictory results in both positive and negative sets
    """
    if noisy_labels == 'majority':
        pos_mask = values >= 0.5
    elif noisy_labels == 'positive':
        pos_mask = values > 0
    else:
        pos_mask = values == 1.0

    neg_mask = ~pos_mask
    pos = pos_mask.index[pos_mask]
    neg = neg_mask.index[neg_mask]

    pos_set = set(pos)
    neg_set = set(neg)

    if verbose:
        print "# positive sequences", len(pos)
        print "# negative sequences", len(neg)

    noisy_set = pos_set.intersection(neg_set)

    if verbose:
        print "# unique positive", len(pos_set)
        print "# unique negative", len(neg_set)
        print "# overlap %d (%0.4f)" % (len(noisy_set), \
          float(len(noisy_set)) / len(pos_set))

    if noisy_labels != 'majority':
        if (noisy_labels == 'drop') or (noisy_labels == 'negative'):
            pos_set = pos_set.difference(noisy_set)
        if (noisy_labels == 'drop') or (noisy_labels == 'positive'):
            neg_set = neg_set.difference(noisy_set)
    if unique_sequences:
        return pos_set, neg_set
    else:
        pos = [epitope for epitope in pos if epitope not in pos_set]
        neg = [epitope for epitope in neg if epitope not in neg_set]
        return pos, neg
