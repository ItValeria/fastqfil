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
