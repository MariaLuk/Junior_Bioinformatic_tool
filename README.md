# My junior biochemical tool
This package for junior bioinformaticians is designed to work with nucleotide and amino acid sequences and FASTQ files, and consists of three parts. Each module `amino_acid_tools`, `dna_rna_tools` and `fastq_filtration` presents its own possibilities for sequence processing.  
 
# `amino_acids_tools`
This tool is designed  to work with amino acid sequences consisting of _22 proteinogenic amino acid_ residues (including pyrrolizine and selenocysteine) recorded in a standard one-letter format. It is not intended to process sequences with post-translational and other amino acid modifications

## Usage  
You call the `amino_acid_tools` function, which takes as input an arbitrary number of arguments with amino-acid sequences (str), as well as the name of the procedure to be executed (it is always the last argument, str). After that the command performs the specified action on all the given sequences. If one sequence is submitted, a string with the result is returned. If several sequences are submitted, a list of strings is returned.  
Input sequences can contain both uppercase and lowercase letters, but the last argument with the function name must correspond to the listed functions.


## Options  
The following options for aminoacid sequence processing are available at the moment:

- **molecular_weight**: calculate the molecular weight of the amino acid chain in Da, according to the average amino acid residues molecular masses rounded to 1 or 2 decimal places.  
- **three_letter_code**: converts standard single letter translations to three letter translations  
- **show_length**: count the overall number of amino acids in the given  
- **sequence folding**: count the number of amino acids characteristic separately for alpha helixes and beta sheets,and give out what will be the structure of the protein more. This function has been tested on proteins such as 2M3X, 6DT4 (PDB ID) and MHC, CRP. The obtained results corresponded to reality.  
- **seq_charge**: evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7), according to the pKa of amino acid side chains, lysine, pyrrolizine and arginine contribute +1, while asparagine and glutamic amino acids contribute -1. The total charge of a protein is evaluated as positive, negative, or neutral as the sum of these contributions

### Remarks
If sequense contains symbols differ from IUPAC 1-letter code for 22 proteinogenic amino acids,
the result *for this sequence* will be `'unexpected symbols in sequence'` regardless of function
If action is not in the function the message `'unexpected action'` will occur
  
## Examples  
Below is an example of processing an amino acid sequence for different input data
```Python
amino_acid_tools('DNA', 'molecular_weight')
```
Output: 300.28

```Python
amino_acid_tools('DNA', 'russco', 'LOLkEk', 'azaz', 'three_letter_code')
```
Output: ['AspAsnAla', 'ArgSecSerSerCysPyl', 'LeuPylLeuLysGluLys', 'unexpected symbols in sequence']

```Python
amino_acid_tools('DNA', 'russco', 'LOLkEk', 'azaz', 'letter_code')
```
ValueError: Unexpected action

### Using the function for molecular weight calculation

```Python  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'molecular_weight')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'molecular_weight'  
Output: '[1228.66, 1447.8400000000001, 1224.6399999999999]'

### Using the function to convert one-letter translations to three-letter translations

```Python  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'three_letter_code')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'three_letter_code'  
Output: '['GluGlyValIleMetSerGluLeuLysLeuLys', 'ProLeuProLysValGluLeuProProAspPheValAsp', 'AspValIleGlyIleSerIleLeuGlyLysGluVal']'

### Using the function to counts the number of amino acids in the given sequence

```Python  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'show_length')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'show_length'  
Output: '[11, 13, 12]'

### Using the function to determine the predominant secondary structure

```Python  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'folding')  
```  
Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'folding'  
Output: '['alfa_helix', 'equally', 'equally']'

### Using the function to estimate relative charge

```Python  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'seq_charge')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'seq_charge'  
Output: '['neutral', 'negative', 'negative']'

# `DNA_RNA_tools`

This program is designed to work with nucleic acid sequences. We expect RNA and DNA chains in the standard NNNNNNNN form, without phosphate groups, in uppercase or lowercase letters

### Usage
You call the `dna_rna_tools` function, which takes as input anumber of arguments with amino-acid sequences (str), as well as the name of the procedure to be executed (it is always the last argument, str).

The nucleotide sequence can consist of both uppercase and lowercase letters.
Input example:
```Python
dna_rna_tools('GCGGT','auccuc','GCTatGc','complement')
```

### Functions:
Folowing options are avaliable now:

 - **complement:** makes complement sequence to seq
 - **reverse**: reverses sequence, (from 5'-3' to 3' -5' there or back)
 - **reverse_complement**: makes complement sequence and reverse it
 - **transcribe**: make transcript of DNA sequence. in the case of RNA
   sequence, no changes will be returned

*Returns:*
If one sequence is supplied, a string with the result is returned.
If several are submitted, a list of strings is returned.

#### Remarks
If seq contains both U(u) and T(t) the result *for this sequence* will be `"U and T in one seq"` regardless of function

If seq contains symbols differ from standardized oligonucleotide notation, the result *for this sequence* will be `'unexpected symbols in sequence'` regardless of function
If action is not in the function the message `'unexpected action'` will occur

### Examples:
```Python
dna_rna_tools('ATTC', 'CGcGc', 'AZA', 'atuc', 'Au', 'transcribe')
```
Output: ['AUUC', 'CGcGc', 'unexpected symbols in sequence', 'T and U in one seq', 'Au']


```Python
dna_rna_tools('ATTC', 'CGcGc', 'AZA', 'atuc', 'Au', 'revese')
```
ValueError: Unexpected action

```Python
dna_rna_tools('UuCG','complement')
```
Output: 'AaGC'

# `FASTQ_filtration tool`

This function provides you the opportunity to filter the FASTQ list to select sequences according to requirements on three parameters: length, GC composition, and quality of the reed

### Usage
It is required to input the list of sequences in dictionary format, as well as values for filtering parameters: interval by GC-composition, interval by sequence length, threshold value of average quality of reed
Input example:
```Python
fastq_filtration(seqs, gc_bounds=(20, 40), length_bounds=(0, 2 *8), quality_treshold=10)
```

####  Filtering parameters

 - seqs: dictionary of FASTQ sequences *{name: (sequence, quality)}*
 - **gc_bounds:**  interval for the of acceptable GC content, in %, *Default
   value = (0,100)*
 - **length_bounds**:  interval for the of acceptable sequense length in
   number of nucleotide, *Default value = (2,2**32)*
 - **quality_treshold**:  threshold value for average quality per nucleotide
   (phred33 scale), *Default value = 0*

 
### Result:
New dictionary consists of selected sequences after 3-step filtration
#### Remarks
After running without specifying the filtering parameters, all sequences will be be selected as appropriate
You also can specify only the upper limit for GC-content and length filtering

### Examples

```Python
d = {'a': ('atcaaa', '@@@@@@'), 'b': ('gcc', '@@!'), 'c': ('ga', '!!'), 'd': ('ga', '@@')}
```
```Python
fastq_filtration(d) 
```
Output: {'a': ('atcaaa', '@@@@@@'), 'b': ('gcc', '@@!'), 'c': ('ga', '!!'), 'd': ('ga', '@@')}
No filtration with default values

```Python
fastq_filtration(d, 40, (0, 3), 10)
```
Output: {}

```Python
fastq_filtration(d, 50, (0, 7), 10)
```
Output: {'a': ('atcaaa', '@@@@@@'), 'd': ('ga', '@@')}

```Python
fastq_filtration(d, 50, (0, 7), 0)
```
Output: {'a': ('atcaaa', '@@@@@@'), 'c': ('ga', '!!'), 'd': ('ga', '@@')}


##### Contacts
Maria Lukina
maria.v.luk@gmail.com



