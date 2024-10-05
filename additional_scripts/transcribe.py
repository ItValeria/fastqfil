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
