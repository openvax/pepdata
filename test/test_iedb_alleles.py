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

from nose.tools import eq_

from pepdata import iedb

def test_iedb_human_class1_allele():
    allele_dict = iedb.alleles.load_alleles_dict()
    allele = allele_dict["HLA-C*07:02"]
    eq_(allele.mhc_class, "I")
    eq_(allele.locus, "C")

def test_iedb_human_class2_allele():
    allele_dict = iedb.alleles.load_alleles_dict()
    allele = allele_dict["HLA-DRA*01:01/DRB1*04:04"]
    eq_(allele.mhc_class, "II")
    eq_(allele.locus, "DR")


def test_iedb_mouse_class1_allele():
    allele_dict = iedb.alleles.load_alleles_dict()
    allele = allele_dict["H-2-Ds"]
    eq_(allele.mhc_class, "I")
    eq_(allele.locus, "D")

def test_iedb_mouse_class2_allele():
    allele_dict = iedb.alleles.load_alleles_dict()
    allele = allele_dict["H-2-IAq"]
    eq_(allele.mhc_class, "II")
    eq_(allele.locus, "IA")
