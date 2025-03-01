import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


#### Ex 1


class InvalidBiologicalSequence(ValueError):
    pass
 
class BiologicalSequence(ABC):
    """Abstract class for biological sequences"""

    def __init__(self, sequence: str):
        if not self._check_alphabet(sequence):
            raise InvalidBiologicalSequence(f"Wrong alphapet of biological sequence: {sequence}")
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        return self.sequence[index]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

    @abstractmethod
    def _check_alphabet(self, sequence: str) -> bool:
        """Cheks symbols of biological sequence alphabet"""
        pass


class NucleicAcidSequence(BiologicalSequence):
    """Class for nucleic acid sequences"""

    complements = {}
    
    def reverse(self):
        """Returns reversed sequence"""
        return self.__class__(self.sequence[::-1])

    def complement(self):
        """Returns complementary sequence"""
        return self.__class__("".join(self.complements[n] for n in self.sequence))

    def reverse_complement(self):
        """Returns reversed and complementary sequence"""
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """Class for DNA sequences"""

    complements = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def _check_alphabet(self, sequence: str) -> bool:
        return all(n in "ATGC" for n in sequence)

    def transcribe(self):
        """Transcribes DNA to RNA"""
        return RNASequence(self.sequence.replace("T", "U"))

    
class RNASequence(NucleicAcidSequence):
    """Class for RNA sequences"""

    complements = {"A": "U", "U": "A", "G": "C", "C": "G"}

    def _check_alphabet(self, sequence: str) -> bool:
        return all(n in "AUGC" for n in sequence)

    
class AminoAcidSequence(BiologicalSequence):
    """Class for amino acid sequences"""

    valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    def _check_alphabet(self, sequence: str) -> bool:
        return all(aa in self.valid_amino_acids for aa in sequence)

    def amino_acid_composition(self):
        """Returns dictionary with count of each amino acid"""
        return {aa: self.sequence.count(aa) for aa in set(self.sequence)}
    
    
#### Ex 2


def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):
    """Filtrates FastQ-file by length, GC-count and quality"""
    
    if not os.path.isdir("filtered"):
        os.mkdir("filtered")
    if isinstance(gc_bounds, int):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)
   
    with open(input_fastq, "r") as read_fastq, open(os.path.join("filtered", output_fastq), "w") as write_fastq:
       
        # Reading fastq-file
        for record in SeqIO.parse(read_fastq, "fastq"):
            sequence = record.seq
            quality = record.letter_annotations["phred_quality"]
           
            # checks length
            if not (length_bounds[0] <= len(sequence) <= length_bounds[1]):
                continue

            # checks GC-count
            gc_content = gc_fraction(sequence)
            if not (gc_bounds[0] <= gc_content*100 <= gc_bounds[1]):
                continue

            # checks quality
            avg_quality = sum(quality) / len(quality)
            if avg_quality < quality_threshold:
                continue

            # Records output file
            SeqIO.write(record, write_fastq, "fastq")