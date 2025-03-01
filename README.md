# fastqfil
Bioinformatics utility for filtering and processing biological sequences.

## Contents
This package contains functions for:
- Filtering FastQ files based on GC content, sequence length, and quality score.
- Performing DNA/RNA sequence manipulations.
- Converting multi-line FASTA files into a single-line format.

---

## Files and Functions

### fastqfil.py  

1. filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 232), quality_threshold=0)
   - Description: Filters sequences from a FastQ file based on specified GC content, sequence length, and average quality score.  
   - Parameters:  
     - input_fastq *(str)* – Path to the input FastQ file.  
     - output_fastq *(str)* – Name/path of the output filtered FastQ file.  
     - gc_bounds *(tuple or int, default=(0,100))* – Allowed range of GC content in %. If an integer is passed
     - length_bounds *(tuple or int, default=(0, 2**32))* – Allowed range of sequence lengths. If an integer is passed
     - quality_threshold *(float, default=0)* – Minimum average quality score required for sequences.  
   - Output: A filtered FastQ file inside the filtered directory.
2.  Classes

#### InvalidBiologicalSequence
- Exception for invalid biological sequence alphabet.

#### BiologicalSequence (Abstract class)
- Base class for biological sequences.
- Methods: len(), getitem(), str(), repr(), abstract method _check_alphabet().

#### NucleicAcidSequence (Abstract class)
- Base class for nucleic acid sequences (DNA/RNA).
- Methods: reverse(), complement(), reverse_complement().

#### DNASequence
- Class for working with DNA sequences.
- Methods: _check_alphabet(), transcribe().

#### RNASequence
- Class for working with RNA sequences.
- Methods: _check_alphabet().

#### AminoAcidSequence
- Class for working with amino acid sequences.
- Methods: _check_alphabet(), amino_acid_composition().
---

### bio_files_processor.py
This script contains one function:

1. convert_multiline_fasta_to_oneline(input_fasta, output_fasta)  
   - Description: Converts a FASTA file where sequences are split across multiple lines into a new FASTA file where each sequence is written on a single line.  
   - Parameters:  
     - input_fasta *(str)* – Path to the input FASTA file.  
     - output_fasta *(str, optional)* – Path to save the reformatted FASTA file. If not provided, the function will overwrite the input file.  

---

### Example_of_use.ipynb
This file contains examples of using fastqfil.py

---

### input.fastq
An example fastq-sequence for running filter_fastq

---

## Installation
To use this package, install libraries from requirements.txt

#### Authors: ItValeria (github)
