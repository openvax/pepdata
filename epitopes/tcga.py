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
from tcga_sources import TCGA_SOURCES
from download import fetch_data

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

def load_maf_files(sources_dict):
    """
    Given a dictionary mapping cancer types to download urls,
    get all the source MAFs, load them as DataFrames, and then
    concatenate into a single DataFrame
    """
    data_frames = []
    for cancer_type, maf_url in sources_dict.iteritems():
        maf_filename = cancer_type + ".maf"
        path = fetch_data(maf_filename, maf_url)
        df = open_maf(path)
        df['Cancer Type'] = cancer_type
        data_frames.append(df)
    return pd.concat(data_frames)

def load_tcga():
    combined_df = load_maf_files(TCGA_SOURCES)
    return combined_df
