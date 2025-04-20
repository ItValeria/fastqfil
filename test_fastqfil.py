import pytest
import os
import logging
from fastqfil import filter_fastq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import shutil


class TestFastqFilter:
    """Tests for filter_fastq"""

    @pytest.fixture
    def create_test_fastq(self, tmp_path):
        """Creates test fastq file"""
        test_file = tmp_path / "test.fastq"
        records = [
            SeqRecord(
                Seq("ACGT"),
                id="good",
                description="",
                letter_annotations={"phred_quality": [40]*4},
            ),
            SeqRecord(
                Seq("TTTT"),
                id="low_gc",
                description="",
                letter_annotations={"phred_quality": [30]*4},
            ),
            SeqRecord(
                Seq("GGGG"),
                id="low_qual",
                description="",
                letter_annotations={"phred_quality": [10]*4},
            ),
        ]
        SeqIO.write(records, test_file, "fastq")
        return test_file

    def test_output_file_is_created(self, create_test_fastq, tmp_path):
        """Checks creating output file"""
        output = os.path.join(tmp_path, "output.fastq")
        filter_fastq(create_test_fastq, output)
        assert os.path.exists(os.path.join("filtered", output))

    def test_all_records_pass_with_defaults(self, create_test_fastq, tmp_path):
        """Makes output without filtration"""
        output = os.path.join(tmp_path, "default_output.fastq")
        filter_fastq(create_test_fastq, output)
        records = list(SeqIO.parse(os.path.join("filtered", output), "fastq"))
        assert len(records) == 3

    def test_gc_filtering(self, create_test_fastq, tmp_path):
        """Checks filtration by GC count"""
        output = os.path.join(tmp_path, "output.fastq")
        filter_fastq(create_test_fastq, output, gc_bounds=(40, 60))

        records = list(SeqIO.parse(os.path.join("filtered", output), "fastq"))
        assert len(records) == 1
        assert records[0].id == "good"

    def test_quality_filtering(self, create_test_fastq, tmp_path):
        """Checks filtration by quality_threshold"""
        output = os.path.join(tmp_path, "output.fastq")
        filter_fastq(create_test_fastq, output, quality_threshold=20)

        records = list(SeqIO.parse(os.path.join("filtered", output), "fastq"))
        assert len(records) == 2

    def test_length_filtering(self, create_test_fastq, tmp_path):
        """Checks filtration by record lenght"""
        output = os.path.join(tmp_path, "output.fastq")
        filter_fastq(create_test_fastq, output, length_bounds=(5, 100))

        records = list(SeqIO.parse(os.path.join("filtered", output), "fastq"))
        assert len(records) == 0

    def test_nonexistent_file_error(self):
        """Checks error of nonexistent file"""
        with pytest.raises(FileNotFoundError):
            filter_fastq("nonexistent.fastq", "output.fastq")

    def test_filtered_dir_created(self, create_test_fastq, tmp_path):
        """Checks creating of directory filtered"""
        output = os.path.join(tmp_path, "output.fastq")
        if os.path.exists("filtered"):
            shutil.rmtree("filtered")
        filter_fastq(create_test_fastq, output)
        assert os.path.isdir("filtered")

    def test_no_records_pass(self, create_test_fastq, tmp_path):
        """Checks empty output with hard filtering parameters"""
        output = os.path.join(tmp_path, "output.fastq")
        filter_fastq(
            create_test_fastq,
            output,
            gc_bounds=(100, 100),
            quality_threshold=50,
        )

        records = list(SeqIO.parse(os.path.join("filtered", output), "fastq"))
        assert len(records) == 0
        

## Command for test: 
# pytest test_fastqfil.py -v

