from typing import Union
import os
import modules.nucleic_acids_functions as na
import modules.fastq_filters as ff
import modules.amino_acids_functions as aa
NUCLEOTIDES = {'U', 'A', 'g', 't', 'G', 'T', 'a', 'c', 'C', 'u'}
AMINO_ACIDS = {'M', 'O', 'v', 'D', 'f', 'N', 'c', 'A', 'R', 'W', 'I', 'm', 'L', 's', 'H', 'q', 'w', 'V', 'n', 'i',
               'g', 'F', 'S', 'e', 'l', 'U', 'P', 'Q', 'K', 'Y', 'u', 'y', 'd', 'h', 'k', 'r', 't', 'G', 'o', 'E',
               'p', 'T', 'C', 'a'}


def dna_rna_tools(*args: str) -> str | list:
    """
        Performs functions for working with poly- and oligonucleotide sequences.

        Parameters:
            The function must accept a number of nucleotide sequences (str) as input,
            the last  variable must be the function (str) you want to execute.
            The nucleotide sequence can consist of both uppercase and lowercase letters.
        Input example:
            dna_rna_tools('GCGGT','auccuc','GCTatGc','complement')
        Function:
            complement: makes complement sequence to seq
            reverse: reverses sequence, (from 5'-3' to 3' -5' there or back)
            reverse_complement: makes complement sequence and reverse it
            transcribe: make transcript of DNA sequence. in the case of RNA sequence, no changes will be returned
        Returns:
            If one sequence is supplied, a string with the result is returned.
            If several are submitted, a list of strings is returned.
            Depending on the function performed, the following returns will occur:
            complement: (str) or (list) complement sequence to seq
            reverse: (str) or (list) reversed sequence to seq
            reverse_complement: (str) or (list) complement and reversed sequence
            transcribe (str) or (list): transcript of DNA sequence.
                                        in the case of RNA sequence, no changes will be returned,
            If seq contains both U(u) and T(t) the result for this str will be
            "U and T in one seq" regardless of function
            If seq contains symbols differ from standardized oligonucleotide notation, the result will be
            'unexpected symbols in sequence' regardless of function
            If action is not in the function the message 'unexpected action' will occur
    """
    *seqs, action = args
    result = list()
    for seq in seqs:
        if not set(seq).issubset(NUCLEOTIDES):
            result.append('unexpected symbols in sequence')
        elif (('U' in seq) or ('u' in seq)) and (('T' in seq) or ('t' in seq)):
            result.append('T and U in one seq')
        else:
            if action == 'reverse':
                result.append(na.reverse(seq))
            elif action == 'complement':
                result.append(na.complement(seq))
            elif action == 'transcribe':
                result.append(na.transcribe(seq))
            elif action == 'reverse_complement':
                result.append(na.reverse_complement(seq))
            else:
                raise ValueError("unexpected action")
    if len(result) == 1:
        result = result[0]
    return result


def amino_acid_tools(*args: str) -> list | int | float | str :
    """
    Performs functions for working with protein sequences.

        Parameters:
            The function must accept an unlimited number of protein sequences (str) as input,
            the last  variable must be the function (str) you want to execute.
            The amino acid sequence can consist of both uppercase and lowercase letters.
        Input example:
            amino_acid_tools('PLPKVEL','VDviRIkLQ','PPDFGKT','folding')
        Function:
            molecular_weight: calculates molecular weight of the amino acid chain
            three_letter_code: converts single letter translations to three letter translations
            length: counts the number of amino acids in the given sequence
            folding: counts the number of amino acids characteristic separately for alpha helixes and beta sheets,
                    and gives out what will be the structure of the protein more
            seq_charge: evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7)

        Returns:
            If one sequence is supplied, a string with the result is returned.
            If several are submitted, a list of strings is returned.
            Depending on the function performed, the following returns will occur:
                molecular_weight (int) or (list): amino acid sequence molecular weight number or list of numbers
                three_letter_code (str) or (list): translated sequence from one-letter in three-letter code
                length (int) or (list): integer number of amino acid residues
                folding (str) or (list): 'alpha_helix', if there are more alpha helices
                                        'beta_sheet', if there are more beta sheets
                                        'equally', if the probability of alpha spirals and beta sheets are the same
                seq_charge(str) or (list): "positive", "negative" or "neutral"
            If seq contains symbols differ from IUPAC 1letter code for 22 proteinogenic amino acids,
            the result for this sequence will be 'unexpected symbols in sequence' regardless of function
            If action is not in the function the message 'unexpected action' will occur
    """
    *seqs, function = args
    d_of_functions = {'molecular_weight': aa.molecular_weight,
                      'three_letter_code': aa.three_letter_code,
                      'length': aa.length,
                      'folding': aa.folding,
                      'seq_charge': aa.seq_charge}
    answer = []
    if function not in d_of_functions.keys():
        raise ValueError("unexpected action")
    for sequence in seqs:
        if not set(sequence).issubset(AMINO_ACIDS):
            answer.append('unexpected symbols in sequence')
        else:
            answer.append(d_of_functions[function](sequence))
    if len(answer) == 1:
        return answer[0]
    else:
        return answer


def fastq_filtration(input_fastq, gc_bounds=(0, 100), length_bounds=(0, 2 ** 32), quality_treshold=0, output_fastq=''):
    """
    This function provides you the opportunity to filter the FASTQ file to select sequences
    that fit  requirements on 5 parameters: input and output(optional) files, length, GC composition,
    and quality of the reed
    :parameters
        input_fastq: path to fastq file
        gc_bounds: (tuple) interval for the of acceptable GC content, in %
        length_bounds: (tuple) interval for the of acceptable sequense length in number of nucleotide
        quality_treshold: (float) threshold value for average quality per nucleotide (phred33 scale)
        output_fastq = name of output file, ./fastq_filtrator_resuls/output_fastq, if it is not defined,
        it will be the same of the input file

    """
    if not os.path.isdir('fastq_filtrator_resuls'):
        os.mkdir('fastq_filtrator_resuls')
    if output_fastq == '':
        output_fastq = os.path.join('fastq_filtrator_resuls', os.path.basename(input_fastq))
    else:
        output_fastq = os.path.join('fastq_filtrator_resuls', output_fastq + ".fasta")
    seqs = ff.convert_fastq_to_dict(input_fastq)
    if type(gc_bounds) == float or type(gc_bounds) == int:
        gc_bounds_both_side = (0, gc_bounds)
    else:
        gc_bounds_both_side = gc_bounds
    if type(length_bounds) == int:
        length_bounds_both_side = (0, length_bounds)
    else:
        length_bounds_both_side = length_bounds
    good_seqs = ff.filter_length(ff.filter_quality(ff.filter_gc(seqs, gc_bounds_both_side), quality_treshold),
                              length_bounds_both_side)

    ff.write_dict_file_to_fastq(good_seqs, output_fastq)
    return good_seqs


