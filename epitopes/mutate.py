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

import math
from collections import namedtuple

from Bio.Seq import Seq



def mutate_split(sequence, position, ref, alt):
    """
    Return the prefix, mutated region, and suffix of
    applying a mutation to a sequence.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    ref : sequence or str
        What do we expect to find at the position?

    alt : sequence or str
        Alternate allele to insert
    """

    transcript_ref_base = sequence[position:position+len(ref)]
    if ref != '.':
        assert str(transcript_ref_base) == ref, \
            "Transcript ref base %s at position %d != given reference %s" % \
            (transcript_ref_base, position, ref)
    prefix = sequence[:position]
    suffix = sequence[position+len(ref):]
    return prefix, alt, suffix

def mutate(sequence, position, ref, alt):
    """
    Mutate a sequence by substituting given `alt` at instead of `ref` at the
    given `position`.

    Parameters
    ----------
    sequence : sequence
        String of amino acids or DNA bases

    position : int
        Position in the sequence, starts from 0

    ref : sequence or str
        What do we expect to find at the position?

    alt : sequence or str
        Alternate allele to insert
    """
    prefix, alt, suffix = mutate_split(sequence, position, ref, alt)
    return prefix + alt + suffix

def is_frameshift(variant_length):
    return variant_length > 1 and variant_length % 3 != 0

MutationRegion = \
    namedtuple("MutationRegion",
    (
        "seq",   # amino acid sequence of region around mutation
        "start"  # where in the original protein did we start?
        "stop"   # where in the original protein did we end?
        "mutation_start", # where in the region does the mutation start?
        "mutation_stop",  # where in the region does the mutation end?
        "frameshift" # was the mutation a frameshift?
    )

def get_mutation_region(
        seq,
        position,
        variant_length = 1,
        max_length = None,
        min_padding = None,
        with_mutation_coordinates = False):
    """
    Get surrounding region of a sequence around a specific position.

    Parameters
    ----------
    seq : BioPython sequence
        Full amino acid string from which we're extracting a substring

    position : int
        Position around variant

    variant_length : int
        Length of the mutation

    max_length : int, optional
        Maximum length peptide to return

    min_padding : int, optional
        Minimum amount of residues before and after variant affected residues

    return_info : bool (default = False)
        Return start/stop for region, start/stop of mutation in region, and
        boolean flag for frameshifts
    """
    n = len(seq)
    if max_length is None:
        max_length = n
    if min_padding is None:
        min_padding = max_length

    end_pos = min(position + min_padding + 1, len(seq))
    start_pos = max(0, position - min_padding)
    num_residue_affected = int(math.ceil(variant_length / 3.0))

    # move start pos to not include stop codons
    prefix_stop_codon = str(seq[:position]).rfind("*") + 1
    start_pos = max(start_pos, prefix_stop_codon)

    # Choose end position to fit max length
    if is_frameshift(variant_length):
        end_pos = start_pos + max_length
        num_residue_affected = end_pos - position
        frameshift = True
    else:
        frameshift = False

    end_codon = str(seq[position:]).find("*")
    if end_codon > 0:
        end_pos = min(end_pos, position + end_codon)

    seq_region = seq[start_pos : end_pos]

    mutation_end_pos = \
        min(len(seq_region), mutation_start_pos + num_residue_affected)
    mutation_start_pos = position - start_pos
    return MutationRegion(
        seq_region, start_pos, end_pos,
        mutation_start_pos, mutation_end_pos,
        frameshift)



Mutation = namedtuple("Mutation",
    (
        "seq",  # mutated sequence
        "original",  # original sequence
        "pos",



def mutate_protein_from_transcript(
        transcript_seq, position, ref, alt,
            min_padding = None,
            max_length = None):
    """
    Mutate a sequence by inserting the allele into the genomic transcript
    and translate to protein sequence

    Parameters
    ----------
    transcript_seq :  sequence
        Protein sequence we're going to mutate

    position : int
        Position in `transcript_seq`, starts from 0

    ref : sequence or str
        What's already supposed to be at the given position

    alt : sequence or str
        Alternate substring to insert

    """
    # turn any character sequence into a BioPython sequence
    transcript_seq = Seq(str(transcript_seq))

    transcript_ref_base = transcript_seq[position:position+len(ref)]
    if ref != '.':
        assert str(transcript_ref_base) == ref, \
            "Transcript reference base %s at position %d != reference %s" % \
            (transcript_ref_base, position, ref)

    mutated_dna = mutate(transcript_seq, position, ref, alt)
    mutated_peptide = mutated_dna.translate()

    aa_position = int(position / 3)  # genomic position to codon position
    variant_length = max(len(ref), len(alt))
    if max_length:
        variant_length = min(variant_length, max_length)
   region = \
        get_mutation_region(
            mutated_peptide,
            aa_position,
            variant_length,
            min_padding = min_padding,
            with_mutation_coordinates = True)
    original_peptide = transcript_seq.translate()_
    aa_ref =
    return Mutat