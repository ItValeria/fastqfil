# Checking the length, GC content and reading quality of nucleotide chains
def filt_fastq(name, sequence, quality, gc_bounds, length_bounds, quality_threshold):
    if length_bounds[0] <= len(sequence) <= length_bounds[1]:
        if (
            gc_bounds[0]
            <= (sequence.count("C") + sequence.count("G")) / len(sequence) * 100
            <= gc_bounds[1]
        ):
            median_quality = 0
            for symbol in quality:
                median_quality += ord(symbol) - 33
            median_quality = median_quality / len(sequence)
            if median_quality > quality_threshold:
                return True
