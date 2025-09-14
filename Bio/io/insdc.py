# Copyright 2007-2016 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the "genbank" and "embl" file formats.

You are expected to use this module via the Bio.SeqIO functions.
Note that internally this module calls Bio.GenBank to do the actual
parsing of GenBank, EMBL and IMGT files.

See Also:
International Nucleotide Sequence Database Collaboration
http://www.insdc.org/

GenBank
http://www.ncbi.nlm.nih.gov/Genbank/

EMBL Nucleotide Sequence Database
http://www.ebi.ac.uk/embl/

DDBJ (DNA Data Bank of Japan)
http://www.ddbj.nig.ac.jp/

IMGT (use a variant of EMBL format with longer feature indents)
http://imgt.cines.fr/download/LIGM-DB/userman_doc.html
http://imgt.cines.fr/download/LIGM-DB/ftable_doc.html
http://www.ebi.ac.uk/imgt/hla/docs/manual.html

"""

import warnings
from datetime import datetime, date as datetime_date
from string import ascii_letters
from string import digits

from Bio import BiopythonWarning
from Bio import SeqFeature
from Bio import io
from Bio.GenBank.Scanner import _ImgtScanner
from Bio.GenBank.Scanner import EmblScanner
from Bio.GenBank.Scanner import GenBankScanner
from Bio.Seq import UndefinedSequenceError

from .Interfaces import _get_seq_string
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter

# Set containing all characters allowed in feature qualifier keys. See
# https://www.insdc.org/submitting-standards/feature-table/#3.1
_allowed_table_component_name_chars = set(ascii_letters + digits + "_-'*")

# NOTE
# ====
# The "brains" for parsing GenBank, EMBL and IMGT files (and any
# other flat file variants from the INSDC in future) is in
# Bio.GenBank.Scanner (plus the _FeatureConsumer in Bio.GenBank)
# However, all the writing code is in this file.


class GenBankIterator(SequenceIterator):
    """Parser for GenBank files."""

    modes = "t"

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the GenBank file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("GenBank/cor6_6.gb", "gb"):
        ...     print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        Equivalently,

        >>> with open("GenBank/cor6_6.gb") as handle:
        ...     for record in GenBankIterator(handle):
        ...         print(record.id)
        ...
        X55053.1
        X62281.1
        M81224.1
        AJ237582.1
        L31939.1
        AF297471.1

        """
        super().__init__(source, fmt="GenBank")
        self.records = GenBankScanner(debug=0).parse_records(self.stream)

    def __next__(self):
        """Return the next SeqRecord."""
        return next(self.records)


class EmblIterator(SequenceIterator):
    """Parser for EMBL files."""

    modes = "t"

    def __init__(self, source):
        """Break up an EMBL file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.

        This gets called internally by Bio.SeqIO for the EMBL file format:

        >>> from Bio import SeqIO
        >>> for record in SeqIO.parse("EMBL/epo_prt_selection.embl", "embl"):
        ...     print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        Equivalently,

        >>> with open("EMBL/epo_prt_selection.embl") as handle:
        ...     for record in EmblIterator(handle):
        ...         print(record.id)
        ...
        A00022.1
        A00028.1
        A00031.1
        A00034.1
        A00060.1
        A00071.1
        A00072.1
        A00078.1
        CQ797900.1

        """
        super().__init__(source, fmt="EMBL")
        self.records = EmblScanner(debug=0).parse_records(self.stream)

    def __next__(self):
        """Return the next SeqRecord."""
        return next(self.records)


class ImgtIterator(SequenceIterator):
    """Parser for IMGT files."""

    modes = "t"

    def __init__(self, source):
        """Break up an IMGT file into SeqRecord objects.

        Argument source is a file-like object opened in text mode or a path to a file.
        Every section from the LOCUS line to the terminating // becomes
        a single SeqRecord with associated annotation and features.

        Note that for genomes or chromosomes, there is typically only
        one record.
        """
        super().__init__(source, fmt="IMGT")
        self.records = _ImgtScanner(debug=0).parse_records(self.stream)

    def __next__(self):
        """Return the next SeqRecord."""
        return next(self.records)


class GenBankCdsFeatureIterator(SequenceIterator):
    """Parser for GenBank files, creating a SeqRecord for each CDS feature."""

    modes = "t"

    def __init__(self, source):
        """Break up a Genbank file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        super().__init__(source, fmt="GenBank")
        self.records = GenBankScanner(debug=0).parse_cds_features(self.stream)

    def __next__(self):
        """Return the next SeqRecord."""
        return next(self.records)


class EmblCdsFeatureIterator(SequenceIterator):
    """Parser for EMBL files, creating a SeqRecord for each CDS feature."""

    modes = "t"

    def __init__(self, source):
        """Break up a EMBL file into SeqRecord objects for each CDS feature.

        Argument source is a file-like object opened in text mode or a path to a file.

        Every section from the LOCUS line to the terminating // can contain
        many CDS features.  These are returned as with the stated amino acid
        translation sequence (if given).
        """
        super().__init__(source, fmt="EMBL")
        self.records = EmblScanner(debug=0).parse_cds_features(self.stream)

    def __next__(self):
        """Return the next SeqRecord."""
        return next(self.records)


def _insdc_feature_position_string(pos, offset=0):
    """Build a GenBank/EMBL position string (PRIVATE).

    Use offset=1 to add one to convert a start position from python counting.
    """
    if isinstance(pos, SeqFeature.ExactPosition):
        return "%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.WithinPosition):
        # TODO - avoid private variables
        return "(%i.%i)" % (
            pos._left + offset,
            pos._right + offset,
        )
    elif isinstance(pos, SeqFeature.BetweenPosition):
        # TODO - avoid private variables
        return "(%i^%i)" % (
            pos._left + offset,
            pos._right + offset,
        )
    elif isinstance(pos, SeqFeature.BeforePosition):
        return "<%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.AfterPosition):
        return ">%i" % (pos + offset)
    elif isinstance(pos, SeqFeature.OneOfPosition):
        return "one-of(%s)" % ",".join(
            _insdc_feature_position_string(p, offset) for p in pos.position_choices
        )
    elif isinstance(pos, SeqFeature.Position):
        raise NotImplementedError("Please report this as a bug in Biopython.")
    else:
        raise ValueError("Expected a SeqFeature position object.")


def _insdc_location_string_ignoring_strand_and_subfeatures(location, rec_length):
    if location.ref:
        ref = f"{location.ref}:"
    else:
        ref = ""
    assert not location.ref_db
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start == location.end
    ):
        # Special case, for 12:12 return 12^13
        # (a zero length slice, meaning the point between two letters)
        if location.end == rec_length:
            # Very special case, for a between position at the end of a
            # sequence (used on some circular genomes, Bug 3098) we have
            # N:N so return N^1
            return "%s%i^1" % (ref, rec_length)
        else:
            return "%s%i^%i" % (ref, location.end, location.end + 1)
    if (
        isinstance(location.start, SeqFeature.ExactPosition)
        and isinstance(location.end, SeqFeature.ExactPosition)
        and location.start + 1 == location.end
    ):
        # Special case, for 11:12 return 12 rather than 12..12
        # (a length one slice, meaning a single letter)
        return "%s%i" % (ref, location.end)
    elif isinstance(location.start, SeqFeature.UnknownPosition) or isinstance(
        location.end, SeqFeature.UnknownPosition
    ):
        # Special case for features from SwissProt/UniProt files
        if isinstance(location.start, SeqFeature.UnknownPosition) and isinstance(
            location.end, SeqFeature.UnknownPosition
        ):
            # warnings.warn("Feature with unknown location", BiopythonWarning)
            # return "?"
            raise ValueError("Feature with unknown location")
        elif isinstance(location.start, SeqFeature.UnknownPosition):
            # Treat the unknown start position as a BeforePosition
            return "%s<%i..%s" % (
                ref,
                location.end,
                _insdc_feature_position_string(location.end),
            )
        else:
            # Treat the unknown end position as an AfterPosition
            return "%s%s..>%i" % (
                ref,
                _insdc_feature_position_string(location.start, +1),
                location.start + 1,
            )
    else:
        # Typical case, e.g. 12..15 gets mapped to 11:15
        return (
            ref
            + _insdc_feature_position_string(location.start, +1)
            + ".."
            + _insdc_feature_position_string(location.end)
        )


def _insdc_location_string(location, rec_length):
    """Build a GenBank/EMBL location from a (Compound) SimpleLocation (PRIVATE).

    There is a choice of how to show joins on the reverse complement strand,
    GenBank used "complement(join(1,10),(20,100))" while EMBL used to use
    "join(complement(20,100),complement(1,10))" instead (but appears to have
    now adopted the GenBank convention). Notice that the order of the entries
    is reversed! This function therefore uses the first form. In this situation
    we expect the CompoundLocation and its parts to all be marked as
    strand == -1, and to be in the order 19:100 then 0:10.
    """
    try:
        parts = location.parts
        # CompoundLocation
        if location.strand == -1:
            # Special case, put complement outside the join/order/... and reverse order
            return "complement(%s(%s))" % (
                location.operator,
                ",".join(
                    _insdc_location_string_ignoring_strand_and_subfeatures(
                        p, rec_length
                    )
                    for p in parts[::-1]
                ),
            )
        else:
            return "%s(%s)" % (
                location.operator,
                ",".join(_insdc_location_string(p, rec_length) for p in parts),
            )
    except AttributeError:
        # SimpleLocation
        loc = _insdc_location_string_ignoring_strand_and_subfeatures(
            location, rec_length
        )
        if location.strand == -1:
            return f"complement({loc})"
        else:
            return loc



    from Bio._utils import run_doctest

    run_doctest(verbose=0)
