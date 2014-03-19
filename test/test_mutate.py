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

from epitopes import mutate
from Bio.Seq import Seq

def test_snp_mutation():
    seq = Seq("AACCTT")

    mutated = mutate.mutate(seq, 1, 'A', 'G')
    assert(mutated[0] == seq[0])
    assert(mutated[1] == 'G')
    assert(mutated[2] == seq[2])

def test_insertion_mutation():
    seq = Seq("AACT")

    mutated = mutate.mutate(seq, 1, 'A', 'AG')
    assert(len(mutated) == len(seq) + 1)
    assert(mutated[0] == seq[0])
    assert(mutated[1] == 'A')
    assert(mutated[2] == 'G')
    assert(mutated[3] == 'C')

def test_deletion_mutation():
    seq = Seq("AACT")
    mutated = mutate.mutate(seq, 1, 'ACT', 'T')
    assert(len(mutated) == 2)
    assert(mutated[0] == 'A')
    assert(mutated[1] == 'T')

def test_get_mutation_region_snp():
    seq = Seq("EDLTVKIGDFGLATEKSRWSGSHQFEQLS")
    region = mutate.get_mutation_region(seq, 10, variant_length=1, min_padding = 2)
    print(str(seq))
    assert(len(region) == 5)

def test_get_mutation_region_snp_end_codon():
    seq = Seq("EDLTVKIGDFGLA*TEKSRWSGSHQFEQLS")
    region = mutate.get_mutation_region(seq, 10, variant_length=1, min_padding = 2, max_length=10)
    print str(region)
    assert(len(region) == 5)

def test_get_mutation_region_indel():
    seq = Seq("EDLTVKIGDFGLATEKSRWSGSHQFEQLS")
    region = mutate.get_mutation_region(seq, 3, variant_length = 4, max_length=20, min_padding=2)
    print(str(region))
    assert(region[0] == 'D')
    assert(len(region) == 20)

def test_get_mutation_region_indel_end_codon():
    seq = Seq("EDLTVKIGDFGLA*TEKSRWSGSHQFEQLS")
    region = mutate.get_mutation_region(seq, 3, variant_length = 4, max_length=20, min_padding=2)
    print str(region)
    assert(region[0] == 'D')
    assert(len(region) == 12)

def test_mutate_protein_from_transcript_snp():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")
    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AATGCTATTCGTAGT").translate()

    region = mutate.mutate_protein_from_transcript(seq, 1, 'C', 'A', min_padding = 8)
    print(str(region))
    assert(region[0] == 'N')
    assert(str(mutated_seq) == str(region))
    assert(len(region) == 5)

def test_mutate_protein_from_transcript_indel():
    seq = Seq("ACTGCTATTCGTAGT")
    prot_seq = Seq("TAIRS")

    assert(str(seq.translate()) == str(prot_seq))

    mutated_seq = Seq("AAATGCTATTCGTAGT").translate()
    print (str(mutated_seq))

    region = mutate.mutate_protein_from_transcript(seq, 1, 'C', 'AA', min_padding = 8)
    print(str(region))
    assert(region[0] == 'K')
    assert(region[1] == 'C')
    assert(region[2] == 'Y')
    assert(str(region) == str(mutated_seq[:-1]))



