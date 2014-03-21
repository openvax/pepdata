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


def mutate(sequence, position, ref, alt):
    """ mutate a sequence by inserting the allele

    Parameters
    ----------
    sequence : sequence or BioPython sequence
    position : int 
    allele : sequence or str, alternate allele to insert

    """
    transcript_ref_base = sequence[position:position+len(ref)]
    assert str(transcript_ref_base) == ref, \
        "Transcript reference base %s at position %d does not match reference %s" % (transcript_ref_base, position, ref)

    mutable_sequence = sequence.tomutable()
    mutable_sequence[position : position + len(ref)] = alt
    return mutable_sequence

def mutate_protein_from_transcript(transcript_seq, position, ref, alt, max_length = 50, min_padding = 2):
    """ mutate a sequence by inserting the allele into the genomic transcript
        and translate to protein sequence

    Parameters
    ----------
    sequence : sequence or BioPython sequence
    position : int 
    allele : sequence or str, alternate allele to insert

    """
    transcript_ref_base = transcript_seq[position:position+len(ref)]
    assert str(transcript_ref_base) == ref, \
        "Transcript reference base %s at position %d does not match reference %s" % (transcript_ref_base, position, ref)

    mutated_sequence = mutate(transcript_seq, position, ref, alt).toseq().translate()
    aa_position = int(position / 3)  # genomic position to codon position
    variant_length = max(len(ref), len(alt))
    seq = get_mutation_region(mutated_sequence.tomutable(), aa_position, variant_length, max_length = max_length, min_padding = min_padding)

    return seq

def get_mutation_region(seq, position, variant_length = 1, max_length = 50, min_padding = 2):
    """ get surronding region of a sequence around a specific position.

    Parameters
    ----------
    sequence : sequence or BioPython sequence
    position : int, position around variant
    variant_legnth : int, length of the mutation 
    max_length : int, maximum length peptide to return
    min_padding : int, minimum amount of residues before and after variant affected residues
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
