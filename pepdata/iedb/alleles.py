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
from collections import namedtuple
import os
import xml

from ..common import memoize, cache

ALLELE_XML_FILENAME = "MhcAlleleNames.xml"
ALLELE_XML_URL = "http://www.iedb.org/doc/MhcAlleleNameList.zip"
ALLELE_XML_DECOMPRESS = True

def local_path(force_download=False):
    """Downloads allele database from IEDB, returns local path to XML file."""
    return cache.fetch(
        filename=ALLELE_XML_FILENAME,
        url=ALLELE_XML_URL,
        decompress=ALLELE_XML_DECOMPRESS,
        force=force_download)

def delete():
    """Deletes local XML file"""
    path = cache.local_path(
        filename=ALLELE_XML_FILENAME,
        url=ALLELE_XML_URL,
        decompress=ALLELE_XML_DECOMPRESS)
    os.remove(path)

Allele = namedtuple("Allele", [
    "name",
    "mhc_class",
    "locus",
    "organism",
    "synonyms"
])

@memoize
def load_alleles():
    """Parses the IEDB MhcAlleleName XML file and returns a list of Allele
    namedtuple objects containing information about that each allele's HLA
    class and source organism.
    """
    result = []
    path = local_path()
    etree = xml.etree.ElementTree.parse(path)
    for allele in etree.iterfind("MhcAlleleName"):
        name_element = allele.find("DisplayedRestriction")
        mhc_class_element = allele.find("Class")
        # need at least a name and an HLA class
        if name_element is None or mhc_class_element is None:
            continue
        name = name_element.text

        synonyms = set([])
        for synonym_element in allele.iterfind("Synonyms"):
            for synonym in synonym_element.text.split(","):
                synonyms.add(synonym.strip())
        mhc_class = mhc_class_element.text
        organism_element = allele.find("Organsim")
        if organism_element is None:
            organism = None
        else:
            organism = organism_element.text

        locus_element = allele.find("Locus")

        if locus_element is None:
            locus = None
        else:
            locus = locus_element.text

        allele_object = Allele(
            name=name,
            mhc_class=mhc_class,
            locus=locus,
            organism=organism,
            synonyms=synonyms)
        result.append(allele_object)
    return result

@memoize
def load_alleles_dict():
    """Create a dictionary mapping each unique allele name to a namedtuple
    containing information about that alleles class, locus, species, &c.
    """
    alleles = load_alleles()
    result = {}
    for allele in alleles:
        for name in {allele.name}.union(allele.synonyms):
            result[name] = allele
    return result
