# Verify the existence of the requested procedure
def check_of_procedure(procedure_name):
    return procedure_name in [
        "transcribe",
        "reverse",
        "complement",
        "reverse_complement",
    ]


# To return a complementary sequence
def complement(chain):
    complementary_chain = ""
    for nucleotide in chain:
        if nucleotide == "A":
            complementary_chain += "T"
        elif nucleotide == "a":
            complementary_chain += "t"
        elif nucleotide == "U":
            complementary_chain += "A"
        elif nucleotide == "u":
            complementary_chain += "a"
        elif nucleotide == "T":
            complementary_chain += "A"
        elif nucleotide == "t":
            complementary_chain += "a"
        elif nucleotide == "C":
            complementary_chain += "G"
        elif nucleotide == "c":
            complementary_chain += "g"
        elif nucleotide == "G":
            complementary_chain += "C"
        elif nucleotide == "g":
            complementary_chain += "c"
    return complementary_chain


# Return the expanded sequence
def reverse(chain):
    return chain[::-1]


# Return the reverse complementary sequence
def reverse_complement(chain):
    return complement(reverse(chain))


# Return the transcribed sequence
def transcribe(chain):
    transcribed_chain = ""
    for nucleotide in chain:
        if nucleotide == "T":
            transcribed_chain += "U"
        elif nucleotide == "t":
            transcribed_chain += "u"
        else:
            transcribed_chain += nucleotide
    return transcribed_chain


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
