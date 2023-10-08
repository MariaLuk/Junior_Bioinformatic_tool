DICT_COMPLEMENT_DNA = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', 'a': 't', 'c': 'g', 't': 'a', 'g': 'c'}
DICT_COMPLEMENT_RNA = {'A': 'U', 'C': 'G', 'U': 'A', 'G': 'C', 'a': 'u', 'c': 'g', 'u': 'a', 'g': 'c'}
DICT_TRANSCRIPTION = {'A': 'A', 'C': 'C', 'T': 'U', 'G': 'G', 'a': 'a', 'c': 'c', 't': 'u', 'g': 'g'}


def complement(seq: str) -> str:
    """
    Function return complement sequence to sec
        Parameters:
            seq (str): oligo- and polynucleotide  sequence
        Returns:
            (str): oligo- and polynucleotide  sequence
    """
    if ('U' in seq) or ('u' in seq):
        complementary_rna = seq.maketrans(DICT_COMPLEMENT_RNA)
        modified_seq = seq.translate(complementary_rna)
    else:
        complementary_dna = seq.maketrans(DICT_COMPLEMENT_DNA)
        modified_seq = seq.translate(complementary_dna)
    return modified_seq


def transcribe(seq: str) -> str:
    """
    Function return transcript of DNA sequence, if the case of RNA sequence, no changes will be returned
        Parameters:
            seq (str): oligo- and polynucleotide  sequence
        Returns:
            (str): oligo- and polynucleotide  sequence
    """
    if ('U' in seq) or ('u' in seq):
        modified_seq = seq
    else:
        transcribe_dna = seq.maketrans(DICT_TRANSCRIPTION)
        modified_seq = seq.translate(transcribe_dna)
    return modified_seq


def reverse(seq: str) -> str:
    """
    Function return reversed sequence, (from 5'-3' to 3' -5' there or back)
        Parameters:
            seq (str): oligo- and polynucleotide  sequence
        Returns:
            (str): oligo- and polynucleotide  sequence
    """
    return seq[::-1]


def reverse_complement(seq: str) -> str:
    """
    Function makes complement sequence and reverse it, (from 5'-3' to 3' -5' there or back)
        Parameters:
            seq (str): oligo- and polynucleotide  sequence
        Returns:
            (str): oligo- and polynucleotide  sequence
    """
    return complement(seq)[::-1]
