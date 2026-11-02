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

import pytest

from pepdata import iedb

# Bound at import, before conftest's autouse fixture replaces it with a fixture path.
original_mhc_local_path = iedb.mhc.local_path

def test_mhc_hla_a2():
    """
    IEDB MHC: Test that HLA restriction actually decreases number of results and
    that regular expression patterns are being used correctly
    """
    df_all = iedb.mhc.load_dataframe(nrows=1000)
    df_a2_1 = iedb.mhc.load_dataframe(hla='HLA-A2', nrows=1000)
    df_a2_2 = iedb.mhc.load_dataframe(hla=r'HLA-A\*02', nrows=1000)
    df_a2_combined = iedb.mhc.load_dataframe(hla=r'HLA-A2|HLA-A\*02', nrows=1000)
    assert len(df_a2_1) < len(df_all)
    assert len(df_a2_2) < len(df_all)
    assert len(df_a2_combined) <= len(df_a2_1) + len(df_a2_2), \
        "Expected %d <= %d + %d" % \
        (len(df_a2_combined), len(df_a2_1), len(df_a2_2))


def test_mhc_load_dataframe_prints_nothing(capsys):
    # nrows=999 avoids the memoized results of the other tests.
    df = iedb.mhc.load_dataframe(nrows=999)
    assert len(df) > 0
    assert capsys.readouterr().out == ""


def test_mhc_missing_file_error_names_the_download_function(monkeypatch, tmp_path):
    monkeypatch.setattr(
        iedb.mhc.cache, "local_path", lambda **kwargs: str(tmp_path / "missing.csv"))
    with pytest.raises(ValueError, match=r"pepdata\.iedb\.mhc\.download\(\)"):
        original_mhc_local_path(auto_download=False)
