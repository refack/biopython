# Copyright 2006-2024 by Peter Cock.  All rights reserved.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from abc import ABC, abstractmethod
from collections.abc import Callable
from collections.abc import Iterable
from typing import Union

from Bio.SeqRecord import SeqRecord

# Import the parsers
from . import fasta
from . import fastq
from . import abi
from . import ace
from . import gck
from . import gfa
from . import ig
from . import insdc
from . import nib
from . import pdb
from . import phd
from . import pir
from . import seqxml
from . import sff
from . import snapgene
from . import swiss
from . import tab
from . import twobit
from . import uniprot
from . import xdna
from .writers import fasta as fasta_writer
from .writers import fastq as fastq_writer
from .writers import phd as phd_writer
from .writers import pir as pir_writer
from .writers import seqxml as seqxml_writer
from .writers import sff as sff_writer
from .writers import tab as tab_writer
from .writers import xdna as xdna_writer
from .writers import insdc as insdc_writer
from .writers import nib as nib_writer

# Registry of parsers
_FormatToIterator = {
    "ace": ace.AceIterator,
    "abi": abi.AbiIterator,
    "uniprot-xml": uniprot.UniprotIterator,
    "xdna": xdna.XdnaIterator,
    "tab": tab.TabIterator,
    "twobit": twobit.TwoBitIterator,
    "swiss": swiss.SwissIterator,
    "snapgene": snapgene.SnapGeneIterator,
    "phd": phd.PhdIterator,
    "pir": pir.PirIterator,
    "seqxml": seqxml.SeqXmlIterator,
    "sff": sff.SffIterator,
    "sff-trim": sff._SffTrimIterator,
    "gck": gck.GckIterator,
    "gfa1": gfa.Gfa1Iterator,
    "ig": ig.IgIterator,
    "nib": nib.NibIterator,
    "cif-seqres": pdb.CifSeqresIterator,
    "cif-atom": pdb.CifAtomIterator,
    "pdb-atom": pdb.PdbAtomIterator,
    "pdb-seqres": pdb.PdbSeqresIterator,
    "embl": insdc.EmblIterator,
    "embl-cds": insdc.EmblCdsFeatureIterator,
    "gb": insdc.GenBankIterator,
    "genbank": insdc.GenBankIterator,
    "genbank-cds": insdc.GenBankCdsFeatureIterator,
    "imgt": insdc.ImgtIterator,
    "gfa2": gfa.Gfa2Iterator,
    "abi-trim": abi._AbiTrimIterator,
    "fasta": fasta.FastaPyfastxIterator,
    "fastq": fastq.FastqSangerPyfastxIterator,
    "fastq-sanger": fastq.FastqSangerPyfastxIterator,
    "fastq-solexa": fastq.FastqSolexaPyfastxIterator,
    "fastq-illumina": fastq.FastqIlluminaPyfastxIterator,
    "qual": fastq.QualPhredIterator,
}

# Registry of writers
_FormatToWriter = {
    "phd": phd_writer.PhdWriter,
    "pir": pir_writer.PirWriter,
    "seqxml": seqxml_writer.SeqXmlWriter,
    "sff": sff_writer.SffWriter,
    "tab": tab_writer.TabWriter,
    "xdna": xdna_writer.XdnaWriter,
    "gb": insdc_writer.GenBankWriter,
    "genbank": insdc_writer.GenBankWriter,
    "embl": insdc_writer.EmblWriter,
    "imgt": insdc_writer.ImgtWriter,
    "fasta": fasta_writer.FastaWriter,
    "fasta-2line": fasta_writer.FastaTwoLineWriter,
    "nib": nib_writer.NibWriter,
    "phd": phd.PhdWriter,
    "pir": pir.PirWriter,
    "seqxml": seqxml.SeqXmlWriter,
    "sff": sff.SffWriter,
    "tab": tab.TabWriter,
    "xdna": xdna.XdnaWriter,
    "gb": insdc.GenBankWriter,
    "genbank": insdc.GenBankWriter,
    "embl": insdc.EmblWriter,
    "imgt": insdc.ImgtWriter,
    "fastq": fastq_writer.FastqPhredWriter,
    "fastq-sanger": fastq_writer.FastqPhredWriter,
    "fastq-solexa": fastq_writer.FastqSolexaWriter,
    "fastq-illumina": fastq_writer.FastqIlluminaWriter,
    "qual": fastq_writer.QualPhredWriter,
}


from Bio.File import as_handle

def write(
    records: Iterable[SeqRecord] | SeqRecord,
    handle,
    format: str,
) -> int:
    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if not format.islower():
        raise ValueError(f"Format string '{format}' should be lower case")

    if isinstance(handle, SeqRecord):
        raise TypeError("Check arguments, handle should NOT be a SeqRecord")
    if isinstance(handle, list):
        # e.g. list of SeqRecord objects
        raise TypeError("Check arguments, handle should NOT be a list")

    if isinstance(records, SeqRecord):
        # This raised an exception in older versions of Biopython
        records = [records]

    writer_class = _FormatToWriter.get(format)
    if writer_class is not None:
        # The writer class must be instantiated with the handle
        with as_handle(handle, "w") as f:
            writer = writer_class(f)
            count = writer.write_file(records)
        if not isinstance(count, int):
            raise RuntimeError(
                "Internal error - the underlying %s writer "
                "should have returned the record count, not %r" % (format, count)
            )
        return count

    if format in _FormatToIterator:
        raise ValueError(f"Reading format '{format}' is supported, but not writing")

    raise ValueError(f"Unknown format '{format}'")


def parse(handle, format, alphabet=None):
    # Try and give helpful error messages:
    if not isinstance(format, str):
        raise TypeError("Need a string for the file format (lower case)")
    if not format:
        raise ValueError("Format required (lower case string)")
    if not format.islower():
        raise ValueError(f"Format string '{format}' should be lower case")
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    iterator_generator = _FormatToIterator.get(format)
    if iterator_generator:
        return iterator_generator(handle)

    raise ValueError(f"Unknown format '{format}'")


def read(handle, format, alphabet=None):
    with parse(handle, format, alphabet) as records:
        try:
            record = next(records)
        except StopIteration:
            raise ValueError("No records found in handle") from None
        try:
            next(records)
            raise ValueError("More than one record found in handle")
        except StopIteration:
            pass
    return record


def to_dict(records, key_function=None):
    # This is to avoid a lambda function:

    def _default_key_function(rec):
        return rec.id

    if key_function is None:
        key_function = _default_key_function

    d = {}
    for record in records:
        key = key_function(record)
        if key in d:
            raise ValueError(f"Duplicate key '{key}'")
        d[key] = record
    return d


from . import align


def parse_alignment(handle, format, seq_count=None):
    """Parse an alignment file and return an iterator of MultipleSeqAlignment objects."""
    return align.parse(handle, format, seq_count=seq_count)


def write_alignment(alignments, handle, format):
    """Write a set of alignments to a file."""
    return align.write(alignments, handle, format)


__all__ = (
    "parse",
    "read",
    "to_dict",
    "write",
    "parse_alignment",
    "write_alignment",
)


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest()
