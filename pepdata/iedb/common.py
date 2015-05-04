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
import pandas as pd

def group_peptides(
        peptides,
        pos_mask,
        mhc_alleles=None,
        min_count=1):
    """Takes three series of equal length and returns a DataFrame.

    Parameters
    ----------
    peptides : Series
        Amino acid sequences, will serve as unique key for entries in the
        result DataFrame

    pos_mask : Series
        Boolean series indicating whether assay result was positive.

    mhc_alleles : Series, optional
        MHC allele names, if present then results are pairs of
        (peptide, mhc_allele) entries, otherwise just use peptide sequences

    min_count : int, optional
        Minimum count per group to be included in the result
    """

    if mhc_alleles is not None:
        groups = pos_mask.groupby([peptides, mhc_alleles])
    else:
        groups = pos_mask.groupby(peptides)

    values = groups.mean()
    counts = groups.count()
    pos_counts = groups.sum()
    neg_counts = counts - pos_counts

    if min_count > 1:
        mask = counts >= min_count
        values = values[mask]
        counts = counts[mask]
        pos_counts = pos_counts[mask]
        neg_counts = neg_counts[mask]

    result = pd.DataFrame(
        data={
            'value': values,
            'count': counts,
            'pos': pos_counts,
            'neg': neg_counts,
        },
        index=values.index)
    return result


def split_classes(
        values,
        noisy_labels='majority',
        unique_sequences=True,
        verbose=False):
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
        print("# positive sequences", len(pos))
        print("# negative sequences", len(neg))

    noisy_set = pos_set.intersection(neg_set)

    if verbose:
        print("# unique positive", len(pos_set))
        print("# unique negative", len(neg_set))
        print("# overlap %d (%0.4f)" % (len(noisy_set),
          float(len(noisy_set)) / len(pos_set)))

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
