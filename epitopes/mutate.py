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

from Bio.Seq import Seq

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
    transcript_ref_base = sequence[position:position+len(ref)]

    assert str(transcript_ref_base) == ref, \
        "Transcript reference base %s at position %d != given reference %s" % \
        (transcript_ref_base, position, ref)
    prefix = sequence[:position]
    suffix = sequence[position+len(ref):]
    return prefix + alt + suffix


def get_mutation_region(
        seq,
        position,
        variant_length = 1,
        max_length = 50,
        min_padding = 2):
    """
    Get surrounding region of a sequence around a specific position.

    Parameters
    ----------
    seq : BioPython sequence
        Full amino acid string from which we're extracting a substring

    position : int
        Position around variant

    variant_legnth : int
        Length of the mutation

    max_length : int
        Maximum length peptide to return

    min_padding : int
        Minimum amount of residues before and after variant affected residues
    """
    end_pos = min(position + min_padding + 1, len(seq))
    start_pos = max(0, position - min_padding)

    # Choose end position to fit max length
    if variant_length > 1:
        end_pos = start_pos + max_length

    end_codon = str(seq[position:]).find("*")
    if end_codon > 0:
        end_pos = min(end_pos, position + end_codon)
    return seq[start_pos : end_pos]


def mutate_protein_from_transcript(
        transcript_seq, position, ref, alt, min_padding = 8, max_length = None):
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
    assert str(transcript_ref_base) == ref, \
        "Transcript reference base %s at position %d != reference %s" % \
        (transcript_ref_base, position, ref)

    mutated_dna = mutate(transcript_seq, position, ref, alt)
    mutated_peptide = mutated_dna.translate()

    aa_position = int(position / 3)  # genomic position to codon position
    variant_length = max(len(ref), len(alt))
    if max_length:
        variant_length = min(variant_length, max_length)
    seq = get_mutation_region(
        mutated_peptide,
        aa_position,
        variant_length,
        min_padding = min_padding)
    return seq
