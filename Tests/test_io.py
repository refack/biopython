# Copyright 2024 by The Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for the new Bio.io module."""

import unittest
from Bio import io


class TestIoApi(unittest.TestCase):
    def test_fasta_parsing(self):
        """Test parsing a FASTA file with the new API."""
        path = "Tests/Fasta/f002"
        records = list(io.parse(path, "fasta"))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].id, "gi|1348912|gb|G26680|G26680")
        self.assertEqual(len(records[0].seq), 633)

    def test_fastq_parsing(self):
        """Test parsing a FASTQ file with the new API."""
        path = "Tests/Quality/example.fastq"
        records = list(io.parse(path, "fastq"))
        self.assertEqual(len(records), 3)
        self.assertEqual(records[0].id, "EAS54_6_R1_2_1_413_324")
        self.assertEqual(len(records[0].seq), 25)
        self.assertEqual(len(records[0].letter_annotations["phred_quality"]), 25)

    def test_fastq_sanger_parsing(self):
        """Test parsing a fastq-sanger file."""
        path = "Tests/Quality/sanger_faked.fastq"
        record = list(io.parse(path, "fastq-sanger"))[0]
        self.assertEqual(record.id, "Test")
        self.assertEqual(len(record.seq), 41)
        self.assertEqual(
            record.letter_annotations["phred_quality"],
            [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        )

    def test_fastq_solexa_parsing(self):
        """Test parsing a fastq-solexa file."""
        path = "Tests/Quality/solexa_faked.fastq"
        record = list(io.parse(path, "fastq-solexa"))[0]
        self.assertEqual(record.id, "slxa_0001_1_0001_01")
        self.assertEqual(len(record.seq), 46)
        self.assertEqual(
            record.letter_annotations["solexa_quality"],
            [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5]
        )

    def test_fastq_illumina_parsing(self):
        """Test parsing a fastq-illumina file."""
        path = "Tests/Quality/illumina_faked.fastq"
        record = list(io.parse(path, "fastq-illumina"))[0]
        self.assertEqual(record.id, "Test")
        self.assertEqual(len(record.seq), 41)
        self.assertEqual(
            record.letter_annotations["phred_quality"],
            [40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        )


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
