def gc_content(seq: str) -> float:
    """
    Supporting function
    Function counts GC-content in sequence, and returns result in %
        Parameters:
            seq (str): oligo- and polynucleotide  sequence
        Returns:
            (float): % of GC-content
    """
    n = 0
    for nucl in seq:
        if nucl == 'c' or nucl == 'g' or nucl == 'C' or nucl == 'G':
            n += 1
    return 100 * n / len(seq)


def filter_gc(seqs: dict, gc_bounds_both_side=(0, 100)) -> dict:
    """
    This function selects sequences with the GC content of your interest
    :parameters:
        seqs: dictionary of FASTQ sequences {name: (sequence, quality)}
        gc_bound: interval for the of acceptable GC content, in %
    :return:(dict) new dictionary consists of selected sequences
    """
    d_new = {}
    for key, params in seqs.items():
        if gc_bounds_both_side[1] >= gc_content(params[0]) >= gc_bounds_both_side[0]:
            d_new[key] = (params[0], params[1])
    return d_new


def filter_length(seqs: dict, length_bounds_both_side=(0, 2 ** 32)) -> dict:
    """
    This function selects sequences with the length of your interest
    :parameters:
        seqs: dictionary of FASTQ sequences {name: (sequence, quality)}
        length_bound: interval for the of acceptable sequense length in number of nucleotide
    :return:(dict) new dictionary consists of selected sequences
    """
    d_new = {}
    for key, params in seqs.items():
        if length_bounds_both_side[1] >= len(params[0]) >= length_bounds_both_side[0]:
            d_new[key] = (params[0], params[1])
    return d_new


def filter_quality(seqs: dict, quality_treshold=0) -> dict:
    """
    This function selects  FASTQ sequences with appropriate average nucleotide read quality
    parameters:
        seqs: dictionary of FASTQ sequences {name: (sequence, quality)}
        quality_treshold: threshold value for average quality per nucleotide (phred33 scale)
    :return:(dict) new dictionary consists of selected sequences
    """
    d_new = {}
    for key, params in seqs.items():
        quality_sum = 0
        for n in params[1]:
            quality_sum += ord(n) - 33
        if quality_sum / len(params[1]) >= quality_treshold:
            d_new[key] = (params[0], params[1])
    return d_new


def convert_fastq_to_dict(input_path) -> dict:
    """
    This function converts fastq file to dictionary with name as a key, and sequence, comment and quality as values
    :parameters:
        input_path: path for data
    :return:(dict) dictionary {name: (sequence, comment, quality}}
    """
    with open(input_path) as input_file:
        seqs = {}
        while True:
            key = input_file.readline().strip()
            if key == '':
                break
            seq = input_file.readline().strip()
            com = input_file.readline().strip()
            qual = input_file.readline().strip()
            seqs[key] = (seq, com, qual)
    return seqs


def write_dict_file_to_fastq(seqs, output_fastq):
    """
    This function forms output file from the  fastq dictionary {name: (sequence, comment, quality}}
     with name as a key, and sequence, comment and quality as values
    :parameters:
        seqs: dictionary of FASTQ sequences {name: (sequence, quality)}
        output_file: filename for output data
    :return:(dict) output_file.fastq
    """
    with open(output_fastq, 'w') as output_file:
        for key, params in seqs.items():
            output_file.write(key + '\n')
            output_file.write(params[0] + '\n')
            output_file.write(params[1] + '\n')
            output_file.write(params[2] + '\n')
