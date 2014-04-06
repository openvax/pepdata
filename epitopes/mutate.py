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


def annotate(
        original_protein,
        mutated_protein,
        aa_position,
        n_deleted,
        n_inserted,
        frameshift):

    ref_start = aa_position
    ref_stop = max(aa_position + n_deleted, 1)
    aa_ref = original_protein[ref_start : ref_stop]

    mut_start = aa_position
    mut_stop = max(aa_position + n_inserted, 1)
    aa_mut = mutated_protein[mut_start : mut_stop]

    if frameshift:
        return "%s%dfs" % (aa_ref[0], aa_position+1)
    elif n_inserted == 0:
        return "%s%ddel" % (aa_ref, aa_position+1)
    else:
        return "%s%d%s" % \
            (aa_ref, aa_position+1, aa_mut)

Mutation = \
    namedtuple(
        "Mutation",
        (
            "seq",    # mutated region sequence
            "start",  # start position of region in the protein
            "stop",   # stop position of region in the protein
            "mutation_start", # where in the region is the first mutated AA?
            "n_wildtype_removed",  # how many wildtype residues removed?
            "n_mutant_inserted",  # how many new residues in the seq?
            "annot", # mutation annotation i.e. "V600E"
        ))


def mutate_protein_from_transcript(
        transcript_seq,
        position,
        dna_ref,
        dna_alt,
        min_padded_length  = 0,
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

    dna_ref : sequence or str
        What's already supposed to be at the given position

    dna_alt : sequence or str
        Alternate substring to insert

    min_padded_length : int, optional
        Try to extend the result to this length, unless stop codons
        get in the way.

    max_length : int, optional
        Maximum length of the returned string

    """
    # turn any character sequence into a BioPython sequence
    transcript_seq = Seq(str(transcript_seq))

    transcript_ref_base = transcript_seq[position:position+len(dna_ref)]

    if ref != '.':
        assert str(transcript_ref_base) == dna_ref, \
            "Transcript reference base %s at position %d != reference %s" % \
            (transcript_ref_base, position, dna_ref)

    original_protein = transcript_seq.translate()
    original_len = len(original_protein)

    mutated_dna = mutate(transcript_seq, position, dna_ref, dna_alt)
    mutated_protein = mutated_dna.translate()
    mutated_len = len(mutated_protein)

    if max_length is None:
        max_length = mutated_len

    aa_position = int(position / 3)  # genomic position to codon position

    n_dna_ref = len(dna_ref)
    n_dna_alt = len(dna_alt)

    # is this a frameshift mutation?
    if abs(n_dna_ref - n_dna_alt) % 3 != 0:
        # frameshifts 'delete' the rest of the origin protein string
        # and 'insert' the frameshifted residues.
        # These numbers, since they encompass the entire protein,
        # will typically be larger than n_wildtype_deleted/n_mutant_inserted
        # on the region struct below, since that restricts the residues counts
        # to those within a particular region
        n_aa_deleted = original_len - aa_position
        n_aa_inserted = mutated_len - aa_position
        frameshift = True
    else:
        n_aa_deleted = int(math.ceil(n_dna_ref / 3.0))
        n_aa_inserted = int(math.ceil(n_dna_alt / 3.0))
        frameshift = False


    # move start pos to not include stop codons
    prefix_stop_codon = str(mutated_protein[:aa_position]).rfind("*") + 1

    suffix_stop_codon = str(mutated_protein[aa_position:]).find("*")

    # try to come up with padded length
    end_pos = min(aa_position + padding_right + 1, mutated_len)
    start_pos = max(0, aa_position - padding_left)

    start_pos = max(start_pos, prefix_stop_codon)


    if end_codon >= 0:
        end_pos = min(end_pos, aa_position + end_codon)

    seq_region = mutated_protein[start_pos : end_pos]
    mutation_start_pos_in_region = aa_position - start_pos

    annot = \
        _annotate(
            original_protein,
            mutated_protein,
            aa_position,
            n_aa_deleted,
            n_aa_inserted,
            frameshift)

    return Mutation(
        seq = seq_region,
        start = start_pos,
        stop = end_pos,
        mutation_start = aa_position,
        n_wildtype_removed = n_aa_deleted,
        n_mutant_inserted = n_aa_inserted,
        annot = annot)
