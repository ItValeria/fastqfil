#Creates a file with combined sequences
def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = "oneline_fasta.fasta"):
    with open(input_fasta, "r") as read_fasta, open(output_fasta, "w") as write_fasta:
        full_chain = ""
        for line in read_fasta:
            if line[0] == ">":
                full_chain += "\n"
                full_chain += line
            else:
                full_chain += line[:-2]
        write_fasta.write(full_chain[2:])