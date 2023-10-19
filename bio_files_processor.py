import os


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = ''):
    """
    This function covert multiple string FASTA-file to one-long-string format
    :param input_fasta: the pass to the FASTA file
    :param output_fasta: name for output FASTA-file, if not defined it will be saved to the directory 
    'Result'/input_fasta
    Example input: convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta', 'output2'): 
    """
    if not os.path.isdir('multiple_to_online_results'):
        os.mkdir('multiple_to_online_results')
    if output_fasta == '':
        output_fasta = os.path.join('multiple_to_online_results', os.path.basename(input_fasta))
    else:
        output_fasta = os.path.join('multiple_to_online_results', output_fasta + ".fasta")
    with open(input_fasta) as input_file, open(output_fasta, 'w') as output_file:
        current = []
        output_file.write(input_file.readline())
        while True:
            line = input_file.readline()
            current.append(line.strip())
            if line.startswith('>'):
                output_file.write(''.join(current) + '\n')
                output_file.write(line)
                current = []
                break
    
        for line in input_file:
            if line.startswith('>'):
                output_file.write(''.join(current) + '\n')
                output_file.write(line)
                current = []
            else:
                current.append(line.strip())
        output_file.write(''.join(current) + '\n')
        

def change_fasta_start_pos(input_fasta, n: int, output_fasta=''):
    """
    This function moves the position to a new starting nucleotide.
    :param input_fasta: data for processing (1-string single fasta)
    :param n: position of new start nucleotide (started with 0)
    :param output_fasta: name for output fasta file
    :return:
    """
    if not os.path.isdir('shifted_fasta_results'):
        os.mkdir('shifted_fasta_results')
    if output_fasta == '':
        output_fasta = os.path.join('shifted_fasta_results', os.path.basename(input_fasta))
    else:
        output_fasta = os.path.join('shifted_fasta_results', output_fasta + ".fasta")
    with open(input_fasta) as input_file:
        name = input_file.readline().strip()
        string = input_file.readline().strip()
    s1 = ''
    s2 = ''
    if n <= 0:
        n = len(string) + n
    for i in range(0, len(string)):
        if i < n:
            s2 = s2 + string[i]
        else:
            s1 = s1 + string[i]
    new_fasta = s1 +s2
    with open(output_fasta, 'w') as output_file:
        output_file.write(name +'\n')
        output_file.write(new_fasta[0] + '\n')
