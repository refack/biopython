# Copyright 2006-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now removed module Bio.SeqIO.FASTA
"""Bio.SeqIO support for the "fasta" (aka FastA or Pearson) file format.

You are expected to use this module via the Bio.SeqIO functions.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import BiopythonDeprecationWarning


from .Interfaces import _clean
from .Interfaces import _get_seq_string
from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter

import warnings


def SimpleFastaParser(handle):
    """Iterate over Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    For each record a tuple of two strings is returned, the FASTA title
    line (without the leading '>' character), and the sequence (with any
    whitespace removed). The title line is not divided up into an
    identifier (the first word) and comment or description.

    >>> with open("Fasta/dups.fasta") as handle:
    ...     for values in SimpleFastaParser(handle):
    ...         print(values)
    ...
    ('alpha', 'ACGTA')
    ('beta', 'CGTC')
    ('gamma', 'CCGCC')
    ('alpha (again - this is a duplicate entry to test the indexing code)', 'ACGTA')
    ('delta', 'CGCGC')

    """
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return

    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")


def FastaTwoLineParser(handle):
    """Iterate over no-wrapping Fasta records as string tuples.

    Arguments:
     - handle - input stream opened in text mode

    Functionally the same as SimpleFastaParser but with a strict
    interpretation of the FASTA format as exactly two lines per
    record, the greater-than-sign identifier with description,
    and the sequence with no line wrapping.

    Any line wrapping will raise an exception, as will excess blank
    lines (other than the special case of a zero-length sequence
    as the second line of a record).

    Examples
    --------
    This file uses two lines per FASTA record:

    >>> with open("Fasta/aster_no_wrap.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    gi|3298468|dbj|BAA31520.1| SAMIPF = GGH...

    This equivalent file uses line wrapping:

    >>> with open("Fasta/aster.pro") as handle:
    ...     for title, seq in FastaTwoLineParser(handle):
    ...         print("%s = %s..." % (title, seq[:3]))
    ...
    Traceback (most recent call last):
       ...
    ValueError: Expected FASTA record starting with '>' character. Perhaps this file is using FASTA line wrapping? Got: 'MTFGLVYTVYATAIDPKKGSLGTIAPIAIGFIVGANI'

    """
    idx = -1  # for empty file
    for idx, line in enumerate(handle):
        if idx % 2 == 0:  # title line
            if line[0] != ">":
                raise ValueError(
                    "Expected FASTA record starting with '>' character. "
                    "Perhaps this file is using FASTA line wrapping? "
                    f"Got: '{line}'"
                )
            title = line[1:].rstrip()
        else:  # sequence line
            if line[0] == ">":
                raise ValueError(
                    "Two '>' FASTA lines in a row. Missing sequence line "
                    "if this is strict two-line-per-record FASTA format. "
                    f"Have '>{title}' and '{line}'"
                )
            yield title, line.strip()

    if idx == -1:
        pass  # empty file
    elif idx % 2 == 0:  # on a title line
        raise ValueError(
            "Missing sequence line at end of file if this is strict "
            f"two-line-per-record FASTA format. Have title line '{line}'"
        )
    else:
        assert line[0] != ">", "line[0] == '>' ; this should be impossible!"






if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)


class FastaPyfastxIterator(SequenceIterator):
    """Parser for Fasta files using pyfastx."""

    modes = "t"

    def __init__(self, source, alphabet=None):
        """Iterate over Fasta records as SeqRecord objects."""
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        # Don't call super().__init__ as it opens the file.
        # pyfastx will handle the file opening.
        try:
            import pyfastx
        except ImportError:
            raise ImportError(
                "Please install pyfastx to use the pyfastx-based FASTA parser."
            ) from None
        self._iterator = self._create_iterator(source)

    def _create_iterator(self, source):
        import pyfastx
        for name, seq, comment in pyfastx.Fastx(source, comment=True):
            if comment:
                description = f"{name} {comment}"
            else:
                description = name
            yield SeqRecord(Seq(seq), id=name, description=description)

    def __next__(self):
        return next(self._iterator)
