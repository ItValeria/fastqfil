# fastqfil
bionf utility for filtering sequences

Utility fastqfil contains 2 functions:
1. filter_fastq()
  Accepts dictionaries with the sequence name, sequence and its quality of reading,
filters according to the specified parameters of the percentage of GC nucleotides,
the length of the sequences and the average quality of reading. Returns a filtered dictionary
2. run_dna_rna_tools()
  Receives the nucleotide chains and the procedure
Сhecks the correctness of the input
Returns the changed circuits
