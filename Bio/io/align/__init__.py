# Copyright 2024 by The Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Alignment file format parsers and writers."""


from .clustal import ClustalIterator, ClustalWriter
from .emboss import EmbossIterator
from .fasta_m10 import FastaM10Iterator
from .maf import MafIterator, MafWriter
from .mauve import MauveIterator, MauveWriter
from .msf import MsfIterator
from .nexus import NexusIterator, NexusWriter
from .phylip import (
    PhylipIterator,
    PhylipWriter,
    RelaxedPhylipIterator,
    RelaxedPhylipWriter,
    SequentialPhylipIterator,
    SequentialPhylipWriter,
)
from .stockholm import StockholmIterator, StockholmWriter

_FormatToIterator = {
    "clustal": ClustalIterator,
    "emboss": EmbossIterator,
    "fasta-m10": FastaM10Iterator,
    "maf": MafIterator,
    "mauve": MauveIterator,
    "msf": MsfIterator,
    "nexus": NexusIterator,
    "phylip": PhylipIterator,
    "phylip-relaxed": RelaxedPhylipIterator,
    "phylip-sequential": SequentialPhylipIterator,
    "stockholm": StockholmIterator,
}

_FormatToWriter = {
    "clustal": ClustalWriter,
    "maf": MafWriter,
    "mauve": MauveWriter,
    "nexus": NexusWriter,
    "phylip": PhylipWriter,
    "phylip-relaxed": RelaxedPhylipWriter,
    "phylip-sequential": SequentialPhylipWriter,
    "stockholm": StockholmWriter,
}


def parse(handle, format, seq_count=None):
    """Parse an alignment file and return an iterator of MultipleSeqAlignment objects."""
    iterator = _FormatToIterator[format](handle, seq_count=seq_count)
    return iterator


def write(alignments, handle, format):
    """Write a set of alignments to a file."""
    writer = _FormatToWriter[format](handle)
    count = writer.write_file(alignments)
    return count


__all__ = ["parse", "write"]
