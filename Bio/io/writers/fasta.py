# Copyright 2006-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for writing FASTA format files as SeqRecord
# objects.

from Bio.SeqRecord import SeqRecord
from ..Interfaces import _clean, _get_seq_string, SequenceWriter


class FastaWriter(SequenceWriter):
    """FASTA file writer."""

    modes = "t"

    def __init__(self, target, wrap=60, record2title=None):
        """Create a Fasta writer."""
        super().__init__(target)
        if wrap:
            if wrap < 1:
                raise ValueError
        self.wrap = wrap
        self.record2title = record2title

    @classmethod
    def to_string(cls, record):
        """Turn a SeqRecord into a FASTA formatted string, and return it."""
        id = _clean(record.id)
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id:
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id} {description}"
        else:
            title = id
        assert "\n" not in title
        assert "\r" not in title
        lines = [f">{title}\n"]

        data = _get_seq_string(record)  # Catches sequence being None
        assert "\n" not in data
        assert "\r" not in data
        for i in range(0, len(data), 60):
            lines.append(data[i : i + 60] + "\n")

        return "".join(lines)

    def write_record(self, record):
        """Write a single Fasta record to the file."""
        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            id = self.clean(record.id)
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == id:
                # The description includes the id at the start
                title = description
            elif description:
                title = f"{id} {description}"
            else:
                title = id

        assert "\n" not in title
        assert "\r" not in title
        self.handle.write(f">{title}\n")

        data = _get_seq_string(record)  # Catches sequence being None

        assert "\n" not in data
        assert "\r" not in data

        if self.wrap:
            for i in range(0, len(data), self.wrap):
                self.handle.write(data[i : i + self.wrap] + "\n")
        else:
            self.handle.write(data + "\n")


class FastaTwoLineWriter(FastaWriter):
    """Class to write 2-line per record Fasta format files."""

    def __init__(self, handle, record2title=None):
        """Create a 2-line per record Fasta writer."""
        super().__init__(handle, wrap=None, record2title=record2title)

    @classmethod
    def to_string(cls, record):
        """Return a string in FASTA format with the sequence as one line."""
        id = _clean(record.id)
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id:
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id} {description}"
        else:
            title = id
        assert "\n" not in title
        assert "\r" not in title

        data = _get_seq_string(record)  # Catches sequence being None
        assert "\n" not in data
        assert "\r" not in data

        return f">{title}\n{data}\n"
