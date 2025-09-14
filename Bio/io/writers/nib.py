# Copyright 2019 by Michiel de Hoon.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.io support for writing the UCSC nib file format."""

import binascii
import struct
import sys

from ..Interfaces import SequenceWriter


class NibWriter(SequenceWriter):
    """Nib file writer."""

    modes = "b"

    def write_header(self):
        """Write the file header."""
        handle = self.handle
        byteorder = sys.byteorder
        if byteorder == "little":  # little-endian
            signature = "3a3de96b"
        elif byteorder == "big":  # big-endian
            signature = "6be93d3a"
        else:
            raise RuntimeError(f"unexpected system byte order {byteorder}")
        handle.write(bytes.fromhex(signature))

    def write_record(self, record):
        """Write a single record to the output file."""
        handle = self.handle
        sequence = record.seq
        nucleotides = bytes(sequence)
        length = len(sequence)
        handle.write(struct.pack("i", length))
        table = bytes.maketrans(b"TCAGNtcagn", b"0123489abc")
        padding = length % 2
        suffix = padding * b"T"
        nucleotides += suffix
        if not set(nucleotides).issubset(b"ACGTNacgtn"):
            raise ValueError("Sequence should contain A,C,G,T,N,a,c,g,t,n only")
        indices = nucleotides.translate(table)
        handle.write(binascii.unhexlify(indices))

    def write_records(self, records):
        """Write records to the output file, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        self.write_header()
        count = 0
        for record in records:
            if count == 1:
                raise ValueError("More than one sequence found")
            self.write_record(record)
            count += 1
        if count != 1:
            raise ValueError("Must have one sequence")
        return count
