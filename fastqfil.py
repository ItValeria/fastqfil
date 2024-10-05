from additional_scripts.check_of_procedure import check_of_procedure
from additional_scripts.type_of_chain import type_of_chain
from additional_scripts.transcribe import transcribe
from additional_scripts.reverse import reverse
from additional_scripts.complement import complement
from additional_scripts.reverse_complement import reverse_complement

'''Accepts dictionaries with the sequence name,
sequence and its quality of reading,
filters according to the specified parameters
of the percentage of GC nucleotides,
the length of the sequences and the average quality 
of reading. Returns a filtered dictionary'''
def filter_fastq(seqs, gc_bounds = (0, 100), length_bounds = (0, 2**32), quality_threshold = 0):
    filtered_fastq = dict()
    if type(gc_bounds) == int:
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) == int:
        length_bounds = (0, length_bounds)
    for name, sequal in seqs.items():
        sequence, quality = sequal[0], sequal[1]
        if length_bounds[0] <= len(sequence) <= length_bounds[1]:
            if gc_bounds[0] <= (sequence.count("C") + sequence.count("G")) / len(sequence) * 100 <= gc_bounds[1]:
                median_quality = 0
                for symbol in quality:
                    median_quality += ord(symbol) - 33
                median_quality = median_quality/len(sequence)
                if median_quality > quality_threshold:
                    filtered_fastq[name] = (sequence, quality)
    return filtered_fastq

#Receives the nucleotide chains and the procedure
#Checks the correctness of the input
#Returns the changed circuits
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
