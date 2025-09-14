# Copyright 2024 by The Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

import warnings
from Bio import BiopythonDeprecationWarning
from ..Interfaces import SequenceWriter, _get_seq_string, _clean


class TabWriter(SequenceWriter):
    """Class to write simple tab separated format files.

    Each line consists of "id(tab)sequence" only.

    Any description, name or other annotation is not recorded.

    This class is not intended to be used directly. Instead, please use
    the top level ``Bio.SeqIO.write()`` function with ``format="tab"``.
    """

    modes = "t"

    @classmethod
    def to_string(cls, record):
        """Return record as tab separated (id(tab)seq) string."""
        title = _clean(record.id)
        seq = _get_seq_string(record)  # Catches sequence being None
        assert "\t" not in title
        assert "\n" not in title
        assert "\r" not in title
        assert "\t" not in seq
        assert "\n" not in seq
        assert "\r" not in seq
        return f"{title}\t{seq}\n"

    def write_record(self, record):
        """Write a single tab line to the file."""
        self.handle.write(self.to_string(record))


def as_tab(record):
    """Return record as tab separated (id(tab)seq) string."""
    warnings.warn(
        """\
TabIO.as_tab is deprecated.

Instead of

TabIO.as_tab(record)

please use

format(record, "tab")
""",
        DeprecationWarning,
    )
    return TabWriter.to_string(record)


if __name__ == "__main__":
    pass