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

def group_peptides(peptides, mhc_alleles, pos_mask, group_by_allele, min_count):

    if group_by_allele:
        groups = pos_mask.groupby([peptides, mhc_alleles])
    else:
        groups = pos_mask.groupby(peptides)

    values = groups.mean()
    counts = groups.count()
    pos_counts = groups.sum()
    neg_counts = counts - pos_counts

    if min_count:
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
