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

from .common import cache
from .tcga_sources import TCGA_SOURCES

def open_maf(filename):
    """
    Load a TCGA MAF file into a Pandas DataFrame
    """
    with open(filename) as fd:
        lines_to_skip = 0
        while next(fd).startswith('#'):
            lines_to_skip += 1
    return pd.read_csv(
        filename,
        skiprows=lines_to_skip,
        sep="\t",
        low_memory=False)

def _load_maf_files(sources_dict, cancer_type=None):
    """
    Given a dictionary mapping cancer types to download urls,
    get all the source MAFs, load them as DataFrames, and then
    concatenate into a single DataFrame
    """
    data_frames = []
    if cancer_type is None:
        cancer_types = sources_dict.keys()
    elif isinstance(cancer_type, str):
        cancer_types = [cancer_type]
    else:
        assert isinstance(cancer_type, list), \
            "Cancer type must be None, str, or list but got %s" % cancer_type
        cancer_types = cancer_type

    for key in cancer_types:
        assert key in sources_dict, "Unknown cancer type %s" % key
        maf_url = sources_dict[key]
        maf_filename = key + ".maf"
        path = cache.fetch(
            url=maf_url,
            filename=maf_filename)
        df = open_maf(path)
        df['Cancer Type'] = key
        data_frames.append(df)
    return pd.concat(data_frames, ignore_index=True)

def load_dataframe(cancer_type=None):
    return _load_maf_files(TCGA_SOURCES, cancer_type)
