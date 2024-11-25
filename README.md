# fastqfil
bionf utility for filtering sequences

## fastqfil.py contains 2 functions:
1. `filter_fastq(input_fastq, output_fastq)`  
  Accepts path of fastq-file "*.fastq" and name/path of output-file, filters data according to the specified parameters of the percentage of GC nucleotides, the length of the sequences and the average quality of reading. Creates a filtered output "*.fastq"-file inside of directory "filtered"
3. `run_dna_rna_tools(*nucleotide_chains_and_procedure)`  
  Receives the nucleotide chains and the procedure as last argument. Сhecks the correctness of the input. Returns the changed circuits
## bio_files_processor.py contains 1 function:
1. `convert_multiline_fasta_to_oneline(input_fasta, output_fasta)`  
   Accepts 2 arguments as input (input_fasta and output_fasta). Reads the fasta file submitted for input, in which the sequence (DNA/RNA/protein/ ... ) can be divided into several lines, after which it saves it to a new fasta file in which each sequence fits into one line. Output_fasta is an optional argument

### How to use:
You can import the needed functions from fastqfil, using import:  
for example, `from fastqfil.fastqfil import filter_fastq`  
or `from fastqfil.bio_files_processor import convert_multiline_fasta_to_oneline`

#### Authors: ItValeria (github)
