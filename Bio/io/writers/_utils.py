# Copyright 2024 by The Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Internal miscellaneous functions for sequence input/output.
"""


class UndefinedSequenceError(ValueError):
    """Sequence content is not defined."""

    pass


def _get_seq_string(record):
    """Get the sequence as a string (PRIVATE)."""
    if record.seq is None:
        raise UndefinedSequenceError("Sequence content is not defined.")
    return str(record.seq)
