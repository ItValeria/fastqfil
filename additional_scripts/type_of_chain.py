# Determine the type of chain - DNA, RNA or erroneous
def type_of_chain(chain):
    dna_rna_nucl = ["A", "G", "C", "a", "g", "c"]
    dna_nucl = ["T", "t"]
    rna_nucl = ["U", "u"]
    type_of_chain = "unknown"
    for nucleotide in chain:
        if nucleotide in dna_rna_nucl:
            type_of_chain = type_of_chain
        elif nucleotide in rna_nucl and type_of_chain != "dna":
            type_of_chain = "rna"
        elif nucleotide in dna_nucl and type_of_chain != "rna":
            type_of_chain = "dna"
        else:
            return "error"
    return type_of_chain
