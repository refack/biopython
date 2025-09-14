# Copyright 2009-2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.io support for writing FASTQ and QUAL file formats."""

import warnings

from Bio.SeqRecord import SeqRecord
from ..Interfaces import _clean, _get_seq_string, SequenceWriter
from ..qual import _get_phred_quality
from ..qual import _get_sanger_quality_str
from ..qual import _get_illumina_quality_str
from ..qual import _get_solexa_quality_str


class FastqPhredWriter(SequenceWriter):
    """Class to write standard FASTQ format files (using PHRED quality scores)."""

    modes = "t"

    @classmethod
    def to_string(cls, record):
        """Turn a SeqRecord into a Sanger FASTQ formatted string, and return it."""
        seq_str = _get_seq_string(record)
        qualities_str = _get_sanger_quality_str(record)
        if len(qualities_str) != len(seq_str):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq_str), len(qualities_str))
            )
        id_ = _clean(record.id) if record.id else ""
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id_:
            title = description
        elif description:
            title = f"{id_} {description}"
        else:
            title = id_
        return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"

    def write_record(self, record: SeqRecord) -> None:
        """Write a single FASTQ record to the file."""
        self.handle.write(self.to_string(record))


class QualPhredWriter(SequenceWriter):
    """Class to write QUAL format files (using PHRED quality scores)."""

    modes = "t"

    def __init__(self, handle, wrap=60, record2title=None):
        """Create a QUAL writer."""
        super().__init__(handle)
        self.wrap = None
        if wrap:
            if wrap < 1:
                raise ValueError
            self.wrap = wrap
        self.record2title = record2title

    @classmethod
    def to_string(cls, record: SeqRecord) -> str:
        """Turn a SeqRecord into a QUAL formatted string."""
        id_ = _clean(record.id) if record.id else ""
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id_:
            title = description
        elif description:
            title = f"{id_} {description}"
        else:
            title = id_
        lines = [f">{title}\n"]

        qualities = _get_phred_quality(record)
        try:
            qualities_strs = [("%i" % round(q, 0)) for q in qualities]
        except TypeError:
            if None in qualities:
                raise TypeError("A quality value of None was found") from None
            else:
                raise

        while qualities_strs:
            line = qualities_strs.pop(0)
            while qualities_strs and len(line) + 1 + len(qualities_strs[0]) < 60:
                line += " " + qualities_strs.pop(0)
            lines.append(line + "\n")
        return "".join(lines)

    def write_record(self, record: SeqRecord) -> None:
        """Write a single QUAL record to the file."""
        self.handle.write(self.to_string(record))


class FastqSolexaWriter(SequenceWriter):
    """Write old style Solexa/Illumina FASTQ format files."""

    modes = "t"

    @classmethod
    def to_string(cls, record: SeqRecord) -> str:
        """Turn a SeqRecord into a Solexa FASTQ formatted string."""
        seq_str = _get_seq_string(record)
        qualities_str = _get_solexa_quality_str(record)
        if len(qualities_str) != len(seq_str):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq_str), len(qualities_str))
            )
        id_ = _clean(record.id) if record.id else ""
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id_:
            title = description
        elif description:
            title = f"{id_} {description}"
        else:
            title = id_
        return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"

    def write_record(self, record: SeqRecord) -> None:
        """Write a single FASTQ record to the file."""
        self.handle.write(self.to_string(record))


class FastqIlluminaWriter(SequenceWriter):
    """Write Illumina 1.3+ FASTQ format files."""

    modes = "t"

    @classmethod
    def to_string(cls, record: SeqRecord) -> str:
        """Turn a SeqRecord into an Illumina FASTQ formatted string."""
        seq_str = _get_seq_string(record)
        qualities_str = _get_illumina_quality_str(record)
        if len(qualities_str) != len(seq_str):
            raise ValueError(
                "Record %s has sequence length %i but %i quality scores"
                % (record.id, len(seq_str), len(qualities_str))
            )
        id_ = _clean(record.id) if record.id else ""
        description = _clean(record.description)
        if description and description.split(None, 1)[0] == id_:
            title = description
        elif description:
            title = f"{id_} {description}"
        else:
            title = id_
        return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"

    def write_record(self, record: SeqRecord) -> None:
        """Write a single FASTQ record to the file."""
        self.handle.write(self.to_string(record))
