from additional_scripts.filtration_fastq import filt_fastq
import os
from additional_scripts.procedures_with_dna_and_rna import (
    check_of_procedure,
    complement,
    reverse,
    reverse_complement,
    transcribe,
    type_of_chain,
)

"""Accepts dictionaries with the sequence name,
sequence and its quality of reading,
filters according to the specified parameters
of the percentage of GC nucleotides,
the length of the sequences and the average quality
of reading. Returns a filtered dictionary"""


def filter_fastq(
    input_fastq,
    output_fastq,
    gc_bounds=(0, 100),
    length_bounds=(0, 2**32),
    quality_threshold=0,
):
    if not os.path.isdir("filtered"):
        os.mkdir("filtered")
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
    with open(input_fastq, "r") as read_fastq, open(
        os.path.join("filtered", output_fastq), "w"
    ) as write_fastq:
        name, sequence, commentary, quality = "", "", "", ""
        count_of_line = 0
        # Reading 4 lines of the same sequence
        for line in read_fastq:
            if count_of_line == 0:
                name = line
                count_of_line += 1
            elif count_of_line == 1:
                sequence = line
                count_of_line += 1
            elif count_of_line == 2:
                commentary = line
                count_of_line += 1
            elif count_of_line == 3:
                quality = line
                count_of_line = 0
                # Checking the length, GC content and reading quality of nucleotide chains
                if filt_fastq(
                    name, sequence, quality, gc_bounds, length_bounds, quality_threshold
                ):
                    # Writing the verified sequences to an external file
                    write_fastq.write(name)
                    write_fastq.write(sequence)
                    write_fastq.write(commentary)
                    write_fastq.write(quality)
                name, sequence, commentary, quality = "", "", "", ""


# Receives the nucleotide chains and the procedure
# Checks the correctness of the input
# Returns the changed circuits
def run_dna_rna_tools(*nucleotide_chains_and_procedure):
    procedures = {
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }
    nucleotide_chains = []
    procedure_name = nucleotide_chains_and_procedure[-1]
    if not check_of_procedure(procedure_name):
        return "ERROR: The requested procedure does not exist"
    for chain in nucleotide_chains_and_procedure[:-1]:
        if type_of_chain(chain) != "error":
            nucleotide_chains += [chain]
        else:
            return f'ERROR: "{chain}" is not RNA or DNA'
    if len(nucleotide_chains) == 1:
        return procedures[procedure_name](nucleotide_chains[0])
    else:
        result_chains = []
        for chain in nucleotide_chains:
            result_chains += [procedures[procedure_name](chain)]
        return result_chains
