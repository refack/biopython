# Copyright 2009-2020 by Peter Cock.  All rights reserved.
# Copyright 2020 by Michael R. Crusoe
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.SeqIO support for the FASTQ and QUAL file formats.

Note that you are expected to use this code via the Bio.SeqIO interface, as
shown below.

The FASTQ file format is used frequently at the Wellcome Trust Sanger Institute
to bundle a FASTA sequence and its PHRED quality data (integers between 0 and
90).  Rather than using a single FASTQ file, often paired FASTA and QUAL files
are used containing the sequence and the quality information separately.

The PHRED software reads DNA sequencing trace files, calls bases, and
assigns a non-negative quality value to each called base using a logged
transformation of the error probability, Q = -10 log10( Pe ), for example::

    Pe = 1.0,         Q =  0
    Pe = 0.1,         Q = 10
    Pe = 0.01,        Q = 20
    ...
    Pe = 0.00000001,  Q = 80
    Pe = 0.000000001, Q = 90

In typical raw sequence reads, the PHRED quality valuea will be from 0 to 40.
In the QUAL format these quality values are held as space separated text in
a FASTA like file format.  In the FASTQ format, each quality values is encoded
with a single ASCI character using chr(Q+33), meaning zero maps to the
character "!" and for example 80 maps to "q".  For the Sanger FASTQ standard
the allowed range of PHRED scores is 0 to 93 inclusive. The sequences and
quality are then stored in pairs in a FASTA like format.

Unfortunately there is no official document describing the FASTQ file format,
and worse, several related but different variants exist. For more details,
please read this open access publication::

    The Sanger FASTQ file format for sequences with quality scores, and the
    Solexa/Illumina FASTQ variants.
    P.J.A.Cock (Biopython), C.J.Fields (BioPerl), N.Goto (BioRuby),
    M.L.Heuer (BioJava) and P.M. Rice (EMBOSS).
    Nucleic Acids Research 2010 38(6):1767-1771
    https://doi.org/10.1093/nar/gkp1137

The good news is that Roche 454 sequencers can output files in the QUAL format,
and sensibly they use PHREP style scores like Sanger.  Converting a pair of
FASTA and QUAL files into a Sanger style FASTQ file is easy. To extract QUAL
files from a Roche 454 SFF binary file, use the Roche off instrument command
line tool "sffinfo" with the -q or -qual argument.  You can extract a matching
FASTA file using the -s or -seq argument instead.

The bad news is that Solexa/Illumina did things differently - they have their
own scoring system AND their own incompatible versions of the FASTQ format.
Solexa/Illumina quality scores use Q = - 10 log10 ( Pe / (1-Pe) ), which can
be negative.  PHRED scores and Solexa scores are NOT interchangeable (but a
reasonable mapping can be achieved between them, and they are approximately
equal for higher quality reads).

Confusingly early Solexa pipelines produced a FASTQ like file but using their
own score mapping and an ASCII offset of 64. To make things worse, for the
Solexa/Illumina pipeline 1.3 onwards, they introduced a third variant of the
FASTQ file format, this time using PHRED scores (which is more consistent) but
with an ASCII offset of 64.

i.e. There are at least THREE different and INCOMPATIBLE variants of the FASTQ
file format: The original Sanger PHRED standard, and two from Solexa/Illumina.

The good news is that as of CASAVA version 1.8, Illumina sequencers will
produce FASTQ files using the standard Sanger encoding.

You are expected to use this module via the Bio.SeqIO functions, with the
following format names:

    - "qual" means simple quality files using PHRED scores (e.g. from Roche 454)
    - "fastq" means Sanger style FASTQ files using PHRED scores and an ASCII
      offset of 33 (e.g. from the NCBI Short Read Archive and Illumina 1.8+).
      These can potentially hold PHRED scores from 0 to 93.
    - "fastq-sanger" is an alias for "fastq".
    - "fastq-solexa" means old Solexa (and also very early Illumina) style FASTQ
      files, using Solexa scores with an ASCII offset 64. These can hold Solexa
      scores from -5 to 62.
    - "fastq-illumina" means newer Illumina 1.3 to 1.7 style FASTQ files, using
      PHRED scores but with an ASCII offset 64, allowing PHRED scores from 0
      to 62.

We could potentially add support for "qual-solexa" meaning QUAL files which
contain Solexa scores, but thus far there isn't any reason to use such files.

For example, consider the following short FASTQ file::

    @EAS54_6_R1_2_1_413_324
    CCCTTCTTGTCTTCAGCGTTTCTCC
    +
    ;;3;;;;;;;;;;;;7;;;;;;;88
    @EAS54_6_R1_2_1_540_792
    TTGGCAGGCCAAGGCCGATGGATCA
    +
    ;;;;;;;;;;;7;;;;;-;;;3;83
    @EAS54_6_R1_2_1_443_348
    GTTGCTTCTGGCGTGGGTGGGGGGG
    +
    ;;;;;;;;;;;9;7;;.7;393333

This contains three reads of length 25.  From the read length these were
probably originally from an early Solexa/Illumina sequencer but this file
follows the Sanger FASTQ convention (PHRED style qualities with an ASCII
offset of 33).  This means we can parse this file using Bio.SeqIO using
"fastq" as the format name:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Quality/example.fastq", "fastq"):
...     print("%s %s" % (record.id, record.seq))
EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

The qualities are held as a list of integers in each record's annotation:

>>> print(record)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Seq('GTTGCTTCTGGCGTGGGTGGGGGGG')
>>> print(record.letter_annotations["phred_quality"])
[26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

You can use the SeqRecord format method to show this in the QUAL format:

>>> print(record.format("qual"))
>EAS54_6_R1_2_1_443_348
26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
24 18 18 18 18
<BLANKLINE>

Or go back to the FASTQ format, use "fastq" (or "fastq-sanger"):

>>> print(record.format("fastq"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
;;;;;;;;;;;9;7;;.7;393333
<BLANKLINE>

Or, using the Illumina 1.3+ FASTQ encoding (PHRED values with an ASCII offset
of 64):

>>> print(record.format("fastq-illumina"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
ZZZZZZZZZZZXZVZZMVZRXRRRR
<BLANKLINE>

You can also get Biopython to convert the scores and show a Solexa style
FASTQ file:

>>> print(record.format("fastq-solexa"))
@EAS54_6_R1_2_1_443_348
GTTGCTTCTGGCGTGGGTGGGGGGG
+
ZZZZZZZZZZZXZVZZMVZRXRRRR
<BLANKLINE>

Notice that this is actually the same output as above using "fastq-illumina"
as the format! The reason for this is all these scores are high enough that
the PHRED and Solexa scores are almost equal. The differences become apparent
for poor quality reads. See the functions solexa_quality_from_phred and
phred_quality_from_solexa for more details.

If you wanted to trim your sequences (perhaps to remove low quality regions,
or to remove a primer sequence), try slicing the SeqRecord objects.  e.g.

>>> sub_rec = record[5:15]
>>> print(sub_rec)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Seq('TTCTGGCGTG')
>>> print(sub_rec.letter_annotations["phred_quality"])
[26, 26, 26, 26, 26, 26, 24, 26, 22, 26]
>>> print(sub_rec.format("fastq"))
@EAS54_6_R1_2_1_443_348
TTCTGGCGTG
+
;;;;;;9;7;
<BLANKLINE>

If you wanted to, you could read in this FASTQ file, and save it as a QUAL file:

>>> from Bio import SeqIO
>>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
>>> with open("Quality/temp.qual", "w") as out_handle:
...     SeqIO.write(record_iterator, out_handle, "qual")
3

You can of course read in a QUAL file, such as the one we just created:

>>> from Bio import SeqIO
>>> for record in SeqIO.parse("Quality/temp.qual", "qual"):
...     print("%s read of length %d" % (record.id, len(record.seq)))
EAS54_6_R1_2_1_413_324 read of length 25
EAS54_6_R1_2_1_540_792 read of length 25
EAS54_6_R1_2_1_443_348 read of length 25

Notice that QUAL files don't have a proper sequence present!  But the quality
information is there:

>>> print(record)
ID: EAS54_6_R1_2_1_443_348
Name: EAS54_6_R1_2_1_443_348
Description: EAS54_6_R1_2_1_443_348
Number of features: 0
Per letter annotation for: phred_quality
Undefined sequence of length 25
>>> print(record.letter_annotations["phred_quality"])
[26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

Just to keep things tidy, if you are following this example yourself, you can
delete this temporary file now:

>>> import os
>>> os.remove("Quality/temp.qual")

Sometimes you won't have a FASTQ file, but rather just a pair of FASTA and QUAL
files.  Because the Bio.SeqIO system is designed for reading single files, you
would have to read the two in separately and then combine the data.  However,
since this is such a common thing to want to do, there is a helper iterator
defined in this module that does this for you - PairedFastaQualIterator.

Alternatively, if you have enough RAM to hold all the records in memory at once,
then a simple dictionary approach would work:

>>> from Bio import SeqIO
>>> reads = SeqIO.to_dict(SeqIO.parse("Quality/example.fasta", "fasta"))
>>> for rec in SeqIO.parse("Quality/example.qual", "qual"):
...     reads[rec.id].letter_annotations["phred_quality"]=rec.letter_annotations["phred_quality"]

You can then access any record by its key, and get both the sequence and the
quality scores.

>>> print(reads["EAS54_6_R1_2_1_540_792"].format("fastq"))
@EAS54_6_R1_2_1_540_792
TTGGCAGGCCAAGGCCGATGGATCA
+
;;;;;;;;;;;7;;;;;-;;;3;83
<BLANKLINE>

It is important that you explicitly tell Bio.SeqIO which FASTQ variant you are
using ("fastq" or "fastq-sanger" for the Sanger standard using PHRED values,
"fastq-solexa" for the original Solexa/Illumina variant, or "fastq-illumina"
for the more recent variant), as this cannot be detected reliably
automatically.

To illustrate this problem, let's consider an artificial example:

>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> test = SeqRecord(Seq("NACGTACGTA"), id="Test", description="Made up!")
>>> print(test.format("fasta"))
>Test Made up!
NACGTACGTA
<BLANKLINE>
>>> print(test.format("fastq"))
Traceback (most recent call last):
 ...
ValueError: No suitable quality scores found in letter_annotations of SeqRecord (id=Test).

We created a sample SeqRecord, and can show it in FASTA format - but for QUAL
or FASTQ format we need to provide some quality scores. These are held as a
list of integers (one for each base) in the letter_annotations dictionary:

>>> test.letter_annotations["phred_quality"] = [0, 1, 2, 3, 4, 5, 10, 20, 30, 40]
>>> print(test.format("qual"))
>Test Made up!
0 1 2 3 4 5 10 20 30 40
<BLANKLINE>
>>> print(test.format("fastq"))
@Test Made up!
NACGTACGTA
+
!"#$%&+5?I
<BLANKLINE>

We can check this FASTQ encoding - the first PHRED quality was zero, and this
mapped to a exclamation mark, while the final score was 40 and this mapped to
the letter "I":

>>> ord('!') - 33
0
>>> ord('I') - 33
40
>>> [ord(letter)-33 for letter in '!"#$%&+5?I']
[0, 1, 2, 3, 4, 5, 10, 20, 30, 40]

Similarly, we could produce an Illumina 1.3 to 1.7 style FASTQ file using PHRED
scores with an offset of 64:

>>> print(test.format("fastq-illumina"))
@Test Made up!
NACGTACGTA
+
@ABCDEJT^h
<BLANKLINE>

And we can check this too - the first PHRED score was zero, and this mapped to
"@", while the final score was 40 and this mapped to "h":

>>> ord("@") - 64
0
>>> ord("h") - 64
40
>>> [ord(letter)-64 for letter in "@ABCDEJT^h"]
[0, 1, 2, 3, 4, 5, 10, 20, 30, 40]

Notice how different the standard Sanger FASTQ and the Illumina 1.3 to 1.7 style
FASTQ files look for the same data! Then we have the older Solexa/Illumina
format to consider which encodes Solexa scores instead of PHRED scores.

First let's see what Biopython says if we convert the PHRED scores into Solexa
scores (rounding to one decimal place):

>>> for q in [0, 1, 2, 3, 4, 5, 10, 20, 30, 40]:
...     print("PHRED %i maps to Solexa %0.1f" % (q, solexa_quality_from_phred(q)))
PHRED 0 maps to Solexa -5.0
PHRED 1 maps to Solexa -5.0
PHRED 2 maps to Solexa -2.3
PHRED 3 maps to Solexa -0.0
PHRED 4 maps to Solexa 1.8
PHRED 5 maps to Solexa 3.3
PHRED 10 maps to Solexa 9.5
PHRED 20 maps to Solexa 20.0
PHRED 30 maps to Solexa 30.0
PHRED 40 maps to Solexa 40.0

Now here is the record using the old Solexa style FASTQ file:

>>> print(test.format("fastq-solexa"))
@Test Made up!
NACGTACGTA
+
;;>@BCJT^h
<BLANKLINE>

Again, this is using an ASCII offset of 64, so we can check the Solexa scores:

>>> [ord(letter)-64 for letter in ";;>@BCJT^h"]
[-5, -5, -2, 0, 2, 3, 10, 20, 30, 40]

This explains why the last few letters of this FASTQ output matched that using
the Illumina 1.3 to 1.7 format - high quality PHRED scores and Solexa scores
are approximately equal.

"""

import warnings
from math import log
from abc import abstractmethod
from typing import Any
from collections.abc import Callable
from typing import IO
from collections.abc import Iterator
from collections.abc import Mapping
from typing import Optional
from collections.abc import Sequence
from typing import Union
from collections.abc import Iterable
import array
from dataclasses import dataclass

from Bio import BiopythonParserWarning
from Bio import BiopythonWarning
from Bio import BiopythonDeprecationWarning
from Bio import StreamModeError
from Bio.File import as_handle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .Interfaces import _clean
from .Interfaces import _get_seq_string
from .Interfaces import _TextIOSource
from .Interfaces import SequenceIterator
from .Interfaces import SequenceWriter
from .qual import _get_phred_quality
from .qual import _get_sanger_quality_str
from .qual import _get_illumina_quality_str
from .qual import _get_solexa_quality_str
from .qual import SANGER_SCORE_OFFSET
from .qual import SOLEXA_SCORE_OFFSET


# TODO - Default to nucleotide or even DNA?
def FastqGeneralIterator(source: _TextIOSource) -> Iterator[tuple[str, str, str]]:
    """Iterate over Fastq records as string tuples (not as SeqRecord objects).

    Arguments:
     - source - input stream opened in text mode, or a path to a file

    This code does not try to interpret the quality string numerically.  It
    just returns tuples of the title, sequence and quality as strings.  For
    the sequence and quality, any whitespace (such as new lines) is removed.

    Our SeqRecord based FASTQ iterators call this function internally, and then
    turn the strings into a SeqRecord objects, mapping the quality string into
    a list of numerical scores.  If you want to do a custom quality mapping,
    then you might consider calling this function directly.

    For parsing FASTQ files, the title string from the "@" line at the start
    of each record can optionally be omitted on the "+" lines.  If it is
    repeated, it must be identical.

    The sequence string and the quality string can optionally be split over
    multiple lines, although several sources discourage this.  In comparison,
    for the FASTA file format line breaks between 60 and 80 characters are
    the norm.

    **WARNING** - Because the "@" character can appear in the quality string,
    this can cause problems as this is also the marker for the start of
    a new sequence.  In fact, the "+" sign can also appear as well.  Some
    sources recommended having no line breaks in the  quality to avoid this,
    but even that is not enough, consider this example::

        @071113_EAS56_0053:1:1:998:236
        TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA
        +071113_EAS56_0053:1:1:998:236
        IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
        @071113_EAS56_0053:1:1:182:712
        ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG
        +
        @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
        @071113_EAS56_0053:1:1:153:10
        TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT
        +
        IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
        @071113_EAS56_0053:1:3:990:501
        TGGGAGGTTTTATGTGGA
        AAGCAGCAATGTACAAGA
        +
        IIIIIII.IIIIII1@44
        @-7.%<&+/$/%4(++(%

    This is four PHRED encoded FASTQ entries originally from an NCBI source
    (given the read length of 36, these are probably Solexa Illumina reads where
    the quality has been mapped onto the PHRED values).

    This example has been edited to illustrate some of the nasty things allowed
    in the FASTQ format.  Firstly, on the "+" lines most but not all of the
    (redundant) identifiers are omitted.  In real files it is likely that all or
    none of these extra identifiers will be present.

    Secondly, while the first three sequences have been shown without line
    breaks, the last has been split over multiple lines.  In real files any line
    breaks are likely to be consistent.

    Thirdly, some of the quality string lines start with an "@" character.  For
    the second record this is unavoidable.  However for the fourth sequence this
    only happens because its quality string is split over two lines.  A naive
    parser could wrongly treat any line starting with an "@" as the beginning of
    a new sequence!  This code copes with this possible ambiguity by keeping
    track of the length of the sequence which gives the expected length of the
    quality string.

    Using this tricky example file as input, this short bit of code demonstrates
    what this parsing function would return:

    >>> with open("Quality/tricky.fastq") as handle:
    ...     for (title, sequence, quality) in FastqGeneralIterator(handle):
    ...         print(title)
    ...         print("%s %s" % (sequence, quality))
    ...
    071113_EAS56_0053:1:1:998:236
    TTTCTTGCCCCCATAGACTGAGACCTTCCCTAAATA IIIIIIIIIIIIIIIIIIIIIIIIIIIIICII+III
    071113_EAS56_0053:1:1:182:712
    ACCCAGCTAATTTTTGTATTTTTGTTAGAGACAGTG @IIIIIIIIIIIIIIICDIIIII<%<6&-*).(*%+
    071113_EAS56_0053:1:1:153:10
    TGTTCTGAAGGAAGGTGTGCGTGCGTGTGTGTGTGT IIIIIIIIIIIICIIGIIIII>IAIIIE65I=II:6
    071113_EAS56_0053:1:3:990:501
    TGGGAGGTTTTATGTGGAAAGCAGCAATGTACAAGA IIIIIII.IIIIII1@44@-7.%<&+/$/%4(++(%

    Finally we note that some sources state that the quality string should
    start with "!" (which using the PHRED mapping means the first letter always
    has a quality score of zero).  This rather restrictive rule is not widely
    observed, so is therefore ignored here.  One plus point about this "!" rule
    is that (provided there are no line breaks in the quality sequence) it
    would prevent the above problem with the "@" character.
    """
    with as_handle(source) as handle:
        if handle.read(0) != "":
            raise StreamModeError("Fastq files must be opened in text mode") from None

        line = handle.readline()
        if line == "":
            return  # Premature end of file, or just empty?

        while True:
            if line[0] != "@":
                raise ValueError(
                    "Records in Fastq files should start with '@' character"
                )
            title_line = line[1:].rstrip()
            seq_string = ""
            # There will now be one or more sequence lines; keep going until we
            # find the "+" marking the quality line:
            for line in handle:
                if line[0] == "+":
                    break
                seq_string += line.rstrip()
            else:
                if seq_string:
                    raise ValueError("End of file without quality information.")
                else:
                    raise ValueError("Unexpected end of file")
            # The title here is optional, but if present must match!
            second_title = line[1:].rstrip()
            if second_title and second_title != title_line:
                raise ValueError("Sequence and quality captions differ.")
            # This is going to slow things down a little, but assuming
            # this isn't allowed we should try and catch it here:
            if " " in seq_string or "\t" in seq_string:
                raise ValueError("Whitespace is not allowed in the sequence.")
            seq_len = len(seq_string)

            # There will now be at least one line of quality data, followed by
            # another sequence, or EOF
            line = None
            quality_string = ""
            for line in handle:
                if line[0] == "@":
                    # This COULD be the start of a new sequence. However, it MAY just
                    # be a line of quality data which starts with a "@" character.  We
                    # should be able to check this by looking at the sequence length
                    # and the amount of quality data found so far.
                    if len(quality_string) >= seq_len:
                        # We expect it to be equal if this is the start of a new record.
                        # If the quality data is longer, we'll raise an error below.
                        break
                    # Continue - its just some (more) quality data.
                quality_string += line.rstrip()
            else:
                if line is None:
                    raise ValueError("Unexpected end of file")
                line = None

            if seq_len != len(quality_string):
                raise ValueError(
                    "Lengths of sequence and quality values differs for %s (%i and %i)."
                    % (title_line, seq_len, len(quality_string))
                )

            # Return the record and then continue...
            yield (title_line, seq_string, quality_string)

            if line is None:
                break


class FastqIteratorAbstractBaseClass(SequenceIterator[str]):
    """Abstract base class for FASTQ file parsers."""

    modes = "t"

    @property
    @abstractmethod
    def q_mapping(self):
        """Dictionary that maps letters in the quality string to quality values."""
        pass

    @property
    @abstractmethod
    def q_key(self):
        """Key name (string) of the quality values in record.letter_annotations."""
        pass

    def __init__(self, source):
        """Iterate over FASTQ records as SeqRecord objects.

        Arguments:
         - source - input stream opened in text mode, or a path to a file

        The quality values are stored in the `letter_annotations` dictionary
        attribute under the key `q_key`.
        """
        super().__init__(source, fmt="Fastq")
        self.line = None

    def __next__(self) -> SeqRecord:
        """Parse the file and generate SeqRecord objects."""

        line = self.line
        if line is None:
            line = self.stream.readline()
        if not line:
            raise StopIteration
        if line[0] != "@":
            raise ValueError("Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()
        seq_string = ""
        # There will now be one or more sequence lines; keep going until we
        # find the "+" marking the quality line:
        for line in self.stream:
            if line[0] == "+":
                break
            seq_string += line.rstrip()
        else:
            if seq_string:
                raise ValueError("End of file without quality information.")
            else:
                raise ValueError("Unexpected end of file")
        seq_len = len(seq_string)
        # The title here is optional, but if present must match!
        second_title = line[1:].rstrip()
        if second_title and second_title != title_line:
            raise ValueError("Sequence and quality captions differ.")

        # Note: str.isprintable is False for ASCII characters 0-32 and 127
        if not seq_string.isprintable() or " " in seq_string:  # type: ignore
            # first printable character
            raise ValueError("Whitespace is not allowed in the sequence.")

        # There will now be at least one line of quality data, followed by
        # another sequence, or EOF
        line = None
        quality_string = ""
        for line in self.stream:
            if line[0] == "@":
                # This COULD be the start of a new sequence. However, it MAY just
                # be a line of quality data which starts with a "@" character.  We
                # should be able to check this by looking at the sequence length
                # and the amount of quality data found so far.
                if len(quality_string) >= seq_len:
                    # We expect it to be equal if this is the start of a new record.
                    # If the quality data is longer, we'll raise an error below.
                    self.line = line
                    break
                # Continue - its just some (more) quality data.
            quality_string += line.rstrip()
        else:
            if line is None:
                raise ValueError("Unexpected end of file")
            self.line = None

        descr = title_line
        id = descr.split()[0]
        name = id

        if not quality_string.isascii():
            # Look for invalid non-ascii characters
            index = _find_index_where(quality_string, lambda c: not c.isascii())
            assert index >= 0, "Non-ascii char in qualities not found. Biopython bug?"

            details = "is not an ASCII character"
            raise InvalidCharError(quality_string, index, details)

        if len(quality_string) != seq_len:
            # This should happen after ascii check, because non-ascii characters will often trigger this
            raise ValueError(
                f"Lengths of sequence and quality values differs for {title_line} ({seq_len} and {len(quality_string)})."
            )

        byte_scores = quality_string.encode().translate(self.q_mapping)

        if INVALID_CHAR in byte_scores:
            # Look for invalid but still ascii characters
            invalid_index = byte_scores.find(INVALID_CHAR_CODE)

            details = "not in correct range (are you sure you're using the right QualityIO parser?)"
            raise InvalidCharError(quality_string, invalid_index, details)

        # Pass through (standard library) array to handle negative scores from old quality formats
        qualities = array.array("b", byte_scores).tolist()

        # SeqRecord._from_validated avoids length/type checking
        # .encode isn't strictly necessary (Seq init can handle a string), but it is faster to pre-encode
        record = SeqRecord._from_validated(
            Seq(seq_string.encode()),
            id=id,
            name=name,
            description=descr,
            letter_annotations={self.q_key: qualities},
        )
        return record




class QualPhredIterator(SequenceIterator):
    """Parser for QUAL files with PHRED quality scores but no sequence."""

    modes = "t"

    def __init__(
        self,
        source: _TextIOSource,
        alphabet: None = None,
    ) -> None:
        """For QUAL files which include PHRED quality scores, but no sequence.

        For example, consider this short QUAL file::

            >EAS54_6_R1_2_1_413_324
            26 26 18 26 26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26
            26 26 26 23 23
            >EAS54_6_R1_2_1_540_792
            26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26 26 12 26 26
            26 18 26 23 18
            >EAS54_6_R1_2_1_443_348
            26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
            24 18 18 18 18

        Using this module directly you might run:

        >>> with open("Quality/example.qual") as handle:
        ...     for record in QualPhredIterator(handle):
        ...         print("%s read of length %d" % (record.id, len(record.seq)))
        EAS54_6_R1_2_1_413_324 read of length 25
        EAS54_6_R1_2_1_540_792 read of length 25
        EAS54_6_R1_2_1_443_348 read of length 25

        Typically however, you would call this via Bio.SeqIO instead with "qual"
        as the format:

        >>> from Bio import SeqIO
        >>> with open("Quality/example.qual") as handle:
        ...     for record in SeqIO.parse(handle, "qual"):
        ...         print("%s read of length %d" % (record.id, len(record.seq)))
        EAS54_6_R1_2_1_413_324 read of length 25
        EAS54_6_R1_2_1_540_792 read of length 25
        EAS54_6_R1_2_1_443_348 read of length 25

        Only the sequence length is known, as the QUAL file does not contain
        the sequence string itself.

        The quality scores themselves are available as a list of integers
        in each record's per-letter-annotation:

        >>> print(record.letter_annotations["phred_quality"])
        [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

        You can still slice one of these SeqRecord objects:

        >>> sub_record = record[5:10]
        >>> print("%s %s" % (sub_record.id, sub_record.letter_annotations["phred_quality"]))
        EAS54_6_R1_2_1_443_348 [26, 26, 26, 26, 26]

        As of Biopython 1.59, this parser will accept files with negatives quality
        scores but will replace them with the lowest possible PHRED score of zero.
        This will trigger a warning, previously it raised a ValueError exception.
        """
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        super().__init__(source, fmt="QUAL")
        # Skip any text before the first record (e.g. blank lines, comments)
        for line in self.stream:
            if line[0] == ">":
                break
        else:
            line = None
        self._line = line

    def __next__(self) -> SeqRecord:
        """Parse the file and generate SeqRecord objects."""

        line = self._line
        if line is None:
            raise StopIteration
        while True:
            descr = line[1:].rstrip()
            id = descr.split()[0]
            name = id

            qualities: list[int] = []
            for line in self.stream:
                if line[0] == ">":
                    break
                qualities.extend(int(word) for word in line.split())
            else:
                line = None
            self._line = line

            if qualities and min(qualities) < 0:
                warnings.warn(
                    "Negative quality score %i found, substituting PHRED zero instead."
                    % min(qualities),
                    BiopythonParserWarning,
                )
                qualities = [max(0, q) for q in qualities]

            # Return the record and then continue...
            sequence = Seq(None, length=len(qualities))

            # Avoid unnecessary length/type checks
            record = SeqRecord._from_validated(
                sequence,
                id=id,
                name=name,
                description=descr,
                letter_annotations={"phred_quality": qualities},
            )
            return record


assert SANGER_SCORE_OFFSET == ord("!")


class FastqPhredWriter(SequenceWriter):
    """Class to write standard FASTQ format files (using PHRED quality scores).

    Although you can use this class directly, you are strongly encouraged
    to use the top level ``Bio.SeqIO.write()`` function instead via the format
    name "fastq" or the alias "fastq-sanger".

    For example, this code reads in a standard Sanger style FASTQ file
    (using PHRED scores) and re-saves it as another Sanger style FASTQ file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq")
    3

    You might want to do this if the original file included extra line breaks,
    which while valid may not be supported by all tools.  The output file from
    Biopython will have each sequence on a single line, and each quality
    string on a single line (which is considered desirable for maximum
    compatibility).

    In this next example, an old style Solexa/Illumina FASTQ file (using Solexa
    quality scores) is converted into a standard Sanger style FASTQ file using
    PHRED qualities:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/solexa_example.fastq", "fastq-solexa")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq")
    5

    This code is also called if you use the .format("fastq") method of a
    SeqRecord, or .format("fastq-sanger") if you prefer that alias.

    Note that Sanger FASTQ files have an upper limit of PHRED quality 93, which is
    encoded as ASCII 126, the tilde. If your quality scores are truncated to fit, a
    warning is issued.

    P.S. To avoid cluttering up your working directory, you can delete this
    temporary file now:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """

    modes = "t"

    @classmethod
    def to_string(cls, record):
        """Turn a SeqRecord into a Sanger FASTQ formatted string, and return it."""
        # TODO - Is an empty sequence allowed in FASTQ format?
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


def as_fastq(record: SeqRecord) -> str:
    """Turn a SeqRecord into a Sanger FASTQ formatted string, and return it."""
    warnings.warn(
        """\
QualityIO.as_fastq is deprecated.

Instead of

QualityIO.as_fastq(record)

please use

format(record, "fastq")
""",
        DeprecationWarning,
    )
    return FastqPhredWriter.to_string(record)


class QualPhredWriter(SequenceWriter):
    """Class to write QUAL format files (using PHRED quality scores).

    Although you can use this class directly, you are strongly encouraged
    to use the top level ``Bio.SeqIO.write()`` function instead.

    For example, this code reads in a FASTQ file and saves the quality scores
    into a QUAL file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/example.fastq", "fastq")
    >>> with open("Quality/temp.qual", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "qual")
    3

    This code is also called if you use the .format("qual") method of a
    SeqRecord.

    P.S. Don't forget to clean up the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.qual")
    """

    modes = "t"

    def __init__(
        self,
        handle: _TextIOSource,
        wrap: int = 60,
        record2title: Callable[[SeqRecord], str] | None = None,
    ) -> None:
        """Create a QUAL writer.

        Arguments:
         - handle - Handle to an output file, e.g. as returned
           by open(filename, "w")
         - wrap   - Optional line length used to wrap sequence lines.
           Defaults to wrapping the sequence at 60 characters. Use
           zero (or None) for no wrapping, giving a single long line
           for the sequence.
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default a
           combination of the record.id and record.description is
           used.  If the record.description starts with the record.id,
           then just the record.description is used.

        The record2title argument is present for consistency with the
        Bio.SeqIO.FastaIO writer class.
        """
        super().__init__(handle)
        # self.handle = handle
        self.wrap: int | None = None
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
            # This rounds to the nearest integer.
            # TODO - can we record a float in a qual file?
            qualities_strs = [("%i" % round(q, 0)) for q in qualities]
        except TypeError:
            if None in qualities:
                raise TypeError("A quality value of None was found") from None
            else:
                raise

        # Safe wrapping
        while qualities_strs:
            line = qualities_strs.pop(0)
            while qualities_strs and len(line) + 1 + len(qualities_strs[0]) < 60:
                line += " " + qualities_strs.pop(0)
            lines.append(line + "\n")
        return "".join(lines)

    def write_record(self, record: SeqRecord) -> None:
        """Write a single QUAL record to the file."""
        handle = self.handle
        wrap = self.wrap

        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            id_ = self.clean(record.id) if record.id else ""
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == id_:
                # The description includes the id at the start
                title = description
            elif description:
                title = f"{id_} {description}"
            else:
                title = id_
        handle.write(f">{title}\n")

        qualities = _get_phred_quality(record)
        try:
            # This rounds to the nearest integer.
            # TODO - can we record a float in a qual file?
            qualities_strs = [("%i" % round(q, 0)) for q in qualities]
        except TypeError:
            if None in qualities:
                raise TypeError("A quality value of None was found") from None
            else:
                raise

        if wrap is not None and wrap > 5:
            # Fast wrapping
            data = " ".join(qualities_strs)
            while True:
                if len(data) <= wrap:
                    self.handle.write(data + "\n")
                    break
                else:
                    # By construction there must be spaces in the first X chars
                    # (unless we have X digit or higher quality scores!)
                    i = data.rfind(" ", 0, wrap)
                    handle.write(data[:i] + "\n")
                    data = data[i + 1 :]
        elif wrap:
            # Safe wrapping
            while qualities_strs:
                line = qualities_strs.pop(0)
                while qualities_strs and len(line) + 1 + len(qualities_strs[0]) < wrap:
                    line += " " + qualities_strs.pop(0)
                handle.write(line + "\n")
        else:
            # No wrapping
            data = " ".join(qualities_strs)
            handle.write(data + "\n")


def as_qual(record: SeqRecord) -> str:
    """Turn a SeqRecord into a QUAL formatted string."""
    warnings.warn(
        """\
QualityIO.as_qual is deprecated.

Instead of

QualityIO.as_qual(record)

please use

format(record, "qual")
""",
        DeprecationWarning,
    )
    return QualPhredWriter.to_string(record)


class FastqSolexaWriter(SequenceWriter):
    r"""Write old style Solexa/Illumina FASTQ format files (with Solexa qualities).

    This outputs FASTQ files like those from the early Solexa/Illumina
    pipeline, using Solexa scores and an ASCII offset of 64. These are
    NOT compatible with the standard Sanger style PHRED FASTQ files.

    If your records contain a "solexa_quality" entry under letter_annotations,
    this is used, otherwise any "phred_quality" entry will be used after
    conversion using the solexa_quality_from_phred function. If neither style
    of quality scores are present, an exception is raised.

    Although you can use this class directly, you are strongly encouraged
    to use the ``as_fastq_solexa`` function, or top-level ``Bio.SeqIO.write()``
    function instead.  For example, this code reads in a FASTQ file and re-saves
    it as another FASTQ file:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("Quality/solexa_example.fastq", "fastq-solexa")
    >>> with open("Quality/temp.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "fastq-solexa")
    5

    You might want to do this if the original file included extra line breaks,
    which (while valid) may not be supported by all tools.  The output file
    from Biopython will have each sequence on a single line, and each quality
    string on a single line (which is considered desirable for maximum
    compatibility).

    This code is also called if you use the .format("fastq-solexa") method of
    a SeqRecord. For example,

    >>> record = SeqIO.read("Quality/sanger_faked.fastq", "fastq-sanger")
    >>> print(record.format("fastq-solexa"))
    @Test PHRED qualities from 40 to 0 inclusive
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN
    +
    hgfedcba`_^]\[ZYXWVUTSRQPONMLKJHGFECB@>;;
    <BLANKLINE>

    Note that Solexa FASTQ files have an upper limit of Solexa quality 62, which is
    encoded as ASCII 126, the tilde.  If your quality scores must be truncated to fit,
    a warning is issued.

    P.S. Don't forget to delete the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """

    modes = "t"

    @classmethod
    def to_string(cls, record: SeqRecord) -> str:
        """Turn a SeqRecord into a Solexa FASTQ formatted string.

        This is used internally by the SeqRecord's .format("fastq-solexa")
        method and by the SeqIO.write(..., ..., "fastq-solexa") function.
        """
        # TODO - Is an empty sequence allowed in FASTQ format?
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
            # The description includes the id at the start
            title = description
        elif description:
            title = f"{id_} {description}"
        else:
            title = id_
        return f"@{title}\n{seq_str}\n+\n{qualities_str}\n"

    def write_record(self, record: SeqRecord) -> None:
        """Write a single FASTQ record to the file."""
        self.handle.write(self.to_string(record))


def as_fastq_solexa(record: SeqRecord) -> str:
    """Turn a SeqRecord into a Solexa FASTQ formatted string."""
    warnings.warn(
        """\
QualityIO.as_fastq_solexa is deprecated.

Instead of

QualityIO.as_fastq_solexa(record)

please use

format(record, "fastq-solexa")
""",
        DeprecationWarning,
    )
    return FastqSolexaWriter.to_string(record)


class FastqIlluminaWriter(SequenceWriter):
    r"""Write Illumina 1.3+ FASTQ format files (with PHRED quality scores).

    This outputs FASTQ files like those from the Solexa/Illumina 1.3+ pipeline,
    using PHRED scores and an ASCII offset of 64. Note these files are NOT
    compatible with the standard Sanger style PHRED FASTQ files which use an
    ASCII offset of 32.

    Although you can use this class directly, you are strongly encouraged to
    use the ``as_fastq_illumina`` or top-level ``Bio.SeqIO.write()`` function
    with format name "fastq-illumina" instead. This code is also called if you
    use the .format("fastq-illumina") method of a SeqRecord. For example,

    >>> from Bio import SeqIO
    >>> record = SeqIO.read("Quality/sanger_faked.fastq", "fastq-sanger")
    >>> print(record.format("fastq-illumina"))
    @Test PHRED qualities from 40 to 0 inclusive
    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTN
    +
    hgfedcba`_^]\[ZYXWVUTSRQPONMLKJIHGFEDCBA@
    <BLANKLINE>

    Note that Illumina FASTQ files have an upper limit of PHRED quality 62, which is
    encoded as ASCII 126, the tilde. If your quality scores are truncated to fit, a
    warning is issued.
    """

    modes = "t"

    @classmethod
    def to_string(cls, record: SeqRecord) -> str:
        """Turn a SeqRecord into an Illumina FASTQ formatted string.

        This is used internally by the SeqRecord's .format("fastq-illumina")
        method and by the SeqIO.write(..., ..., "fastq-illumina") function.
        """
        # TODO - Is an empty sequence allowed in FASTQ format?
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


def as_fastq_illumina(record: SeqRecord) -> str:
    """Turn a SeqRecord into an Illumina FASTQ formatted string."""
    warnings.warn(
        """\
QualityIO.as_fastq_illumina is deprecated.

Instead of

QualityIO.as_fastq_illumina(record)

please use

format(record, "fastq-illumina")
""",
        DeprecationWarning,
    )
    return FastqIlluminaWriter.to_string(record)


def PairedFastaQualIterator(
    fasta_source: _TextIOSource,
    qual_source: _TextIOSource,
    alphabet: None = None,
) -> Iterator[SeqRecord]:
    """Iterate over matched FASTA and QUAL files as SeqRecord objects.

    For example, consider this short QUAL file with PHRED quality scores::

        >EAS54_6_R1_2_1_413_324
        26 26 18 26 26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26
        26 26 26 23 23
        >EAS54_6_R1_2_1_540_792
        26 26 26 26 26 26 26 26 26 26 26 22 26 26 26 26 26 12 26 26
        26 18 26 23 18
        >EAS54_6_R1_2_1_443_348
        26 26 26 26 26 26 26 26 26 26 26 24 26 22 26 26 13 22 26 18
        24 18 18 18 18

    And a matching FASTA file::

        >EAS54_6_R1_2_1_413_324
        CCCTTCTTGTCTTCAGCGTTTCTCC
        >EAS54_6_R1_2_1_540_792
        TTGGCAGGCCAAGGCCGATGGATCA
        >EAS54_6_R1_2_1_443_348
        GTTGCTTCTGGCGTGGGTGGGGGGG

    You can parse these separately using Bio.SeqIO with the "qual" and
    "fasta" formats, but then you'll get a group of SeqRecord objects with
    no sequence, and a matching group with the sequence but not the
    qualities.  Because it only deals with one input file handle, Bio.SeqIO
    can't be used to read the two files together - but this function can!
    For example,

    >>> with open("Quality/example.fasta") as f:
    ...     with open("Quality/example.qual") as q:
    ...         for record in PairedFastaQualIterator(f, q):
    ...             print("%s %s" % (record.id, record.seq))
    ...
    EAS54_6_R1_2_1_413_324 CCCTTCTTGTCTTCAGCGTTTCTCC
    EAS54_6_R1_2_1_540_792 TTGGCAGGCCAAGGCCGATGGATCA
    EAS54_6_R1_2_1_443_348 GTTGCTTCTGGCGTGGGTGGGGGGG

    As with the FASTQ or QUAL parsers, if you want to look at the qualities,
    they are in each record's per-letter-annotation dictionary as a simple
    list of integers:

    >>> print(record.letter_annotations["phred_quality"])
    [26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 24, 26, 22, 26, 26, 13, 22, 26, 18, 24, 18, 18, 18, 18]

    If you have access to data as a FASTQ format file, using that directly
    would be simpler and more straight forward.  Note that you can easily use
    this function to convert paired FASTA and QUAL files into FASTQ files:

    >>> from Bio import SeqIO
    >>> with open("Quality/example.fasta") as f:
    ...     with open("Quality/example.qual") as q:
    ...         SeqIO.write(PairedFastaQualIterator(f, q), "Quality/temp.fastq", "fastq")
    ...
    3

    And don't forget to clean up the temp file if you don't need it anymore:

    >>> import os
    >>> os.remove("Quality/temp.fastq")
    """
    if alphabet is not None:
        raise ValueError("The alphabet argument is no longer supported")

    from Bio.SeqIO.FastaIO import FastaIterator

    fasta_iter = FastaIterator(fasta_source)
    qual_iter = QualPhredIterator(qual_source)

    # Using zip wouldn't load everything into memory, but also would not catch
    # any extra records found in only one file.
    while True:
        try:
            f_rec = next(fasta_iter)
        except StopIteration:
            f_rec = None
        try:
            q_rec = next(qual_iter)
        except StopIteration:
            q_rec = None
        if f_rec is None and q_rec is None:
            # End of both files
            break
        if f_rec is None:
            raise ValueError("FASTA file has more entries than the QUAL file.")
        if q_rec is None:
            raise ValueError("QUAL file has more entries than the FASTA file.")
        if f_rec.id != q_rec.id:
            raise ValueError(
                f"FASTA and QUAL entries do not match ({f_rec.id} vs {q_rec.id})."
            )
        if len(f_rec) != len(q_rec.letter_annotations["phred_quality"]):
            raise ValueError(
                f"Sequence length and number of quality scores disagree for {f_rec.id}"
            )
        # Merge the data....
        f_rec.letter_annotations["phred_quality"] = q_rec.letter_annotations[
            "phred_quality"
        ]
        yield f_rec
    # Done


def _fastq_generic(
    in_file: _TextIOSource,
    out_file: _TextIOSource,
    mapping: Sequence[str] | Mapping[int, str | int | None],
) -> int:
    """FASTQ helper function where can't have data loss by truncation (PRIVATE)."""
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    with as_handle(out_file, "w") as out_handle:
        for title, seq, old_qual in FastqGeneralIterator(in_file):
            count += 1
            # map the qual...
            qual = old_qual.translate(mapping)
            if null in qual:
                invalid_index = qual.find(null)
                raise InvalidCharError(
                    old_qual,
                    invalid_index,
                    details="not in correct range (are you sure you're using the right QualityIO parser?)",
                )

            out_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")
    return count


def _fastq_generic2(
    in_file: _TextIOSource,
    out_file: _TextIOSource,
    mapping: Sequence[str] | Mapping[int, str | int | None],
    truncate_char: str,
    truncate_msg: str,
) -> int:
    """FASTQ helper function where there could be data loss by truncation (PRIVATE)."""
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    null = chr(0)
    with as_handle(out_file, "w") as out_handle:
        for title, seq, old_qual in FastqGeneralIterator(in_file):
            count += 1
            # map the qual...
            qual = old_qual.translate(mapping)
            if null in qual:
                invalid_index = qual.find(null)
                raise InvalidCharError(
                    old_qual,
                    invalid_index,
                    details="not in correct range (are you sure you're using the right QualityIO parser?)",
                )
            if truncate_char in qual:
                qual = qual.replace(truncate_char, chr(126))
                warnings.warn(truncate_msg, BiopythonWarning)
            out_handle.write(f"@{title}\n{seq}\n+\n{qual}\n")
    return count


def _fastq_sanger_convert_fastq_sanger(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Sanger FASTQ to Sanger FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(33)]
        + [chr(ascii) for ascii in range(33, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_solexa_convert_fastq_solexa(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Solexa FASTQ to Solexa FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(59)]
        + [chr(ascii) for ascii in range(59, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_illumina(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Illumina 1.3+ FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Useful for removing line wrapping and the redundant second identifier
    on the plus lines. Will check also check the quality string is valid.
    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(64)]
        + [chr(ascii) for ascii in range(64, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_sanger(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Illumina 1.3+ FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(64)]
        + [chr(33 + q) for q in range(62 + 1)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_sanger_convert_fastq_illumina(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Sanger FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Illumina 1.3+ FASTQ format)
    """
    # Map unexpected chars to null
    trunc_char = chr(1)
    mapping = "".join(
        [chr(0) for ascii in range(33)]
        + [chr(64 + q) for q in range(62 + 1)]
        + [trunc_char for ascii in range(96, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic2(
        in_file,
        out_file,
        mapping,
        trunc_char,
        "Data loss - max PHRED quality 62 in Illumina 1.3+ FASTQ",
    )


def _fastq_solexa_convert_fastq_sanger(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Solexa FASTQ to Sanger FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(59)]
        + [
            chr(33 + int(round(phred_quality_from_solexa(q))))
            for q in range(-5, 62 + 1)
        ]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_sanger_convert_fastq_solexa(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Sanger FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion. Will issue a warning if the scores had to be truncated at 62
    (maximum possible in the Solexa FASTQ format)
    """
    # Map unexpected chars to null
    trunc_char = chr(1)
    mapping = "".join(
        [chr(0) for ascii in range(33)]
        + [chr(64 + int(round(solexa_quality_from_phred(q)))) for q in range(62 + 1)]
        + [trunc_char for ascii in range(96, 127)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic2(
        in_file,
        out_file,
        mapping,
        trunc_char,
        "Data loss - max Solexa quality 62 in Solexa FASTQ",
    )


def _fastq_solexa_convert_fastq_illumina(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Solexa FASTQ to Illumina 1.3+ FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(59)]
        + [
            chr(64 + int(round(phred_quality_from_solexa(q))))
            for q in range(-5, 62 + 1)
        ]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_illumina_convert_fastq_solexa(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Illumina 1.3+ FASTQ to Solexa FASTQ conversion (PRIVATE).

    Avoids creating SeqRecord and Seq objects in order to speed up this
    conversion.
    """
    # Map unexpected chars to null
    mapping = "".join(
        [chr(0) for ascii in range(64)]
        + [chr(64 + int(round(solexa_quality_from_phred(q)))) for q in range(62 + 1)]
        + [chr(0) for ascii in range(127, 256)]
    )
    assert len(mapping) == 256
    return _fastq_generic(in_file, out_file, mapping)


def _fastq_convert_fasta(in_file: _TextIOSource, out_file: _TextIOSource) -> int:
    """Fast FASTQ to FASTA conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f">{title}\n")
            # Do line wrapping
            for i in range(0, len(seq), 60):
                out_handle.write(seq[i : i + 60] + "\n")
    return count


def _fastq_convert_tab(in_file: _TextIOSource, out_file: _TextIOSource) -> int:
    """Fast FASTQ to simple tabbed conversion (PRIVATE).

    Avoids dealing with the FASTQ quality encoding, and creating SeqRecord and
    Seq objects in order to speed up this conversion.

    NOTE - This does NOT check the characters used in the FASTQ quality string
    are valid!
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f"{title.split(None, 1)[0]}\t{seq}\n")
    return count


def _fastq_convert_qual(
    in_file: _TextIOSource,
    out_file: _TextIOSource,
    mapping: Mapping[str, str],
) -> int:
    """FASTQ helper function for QUAL output (PRIVATE).

    Mapping should be a dictionary mapping expected ASCII characters from the
    FASTQ quality string to PHRED quality scores (as strings).
    """
    # For real speed, don't even make SeqRecord and Seq objects!
    count = 0
    with as_handle(out_file, "w") as out_handle:
        for title, seq, qual in FastqGeneralIterator(in_file):
            count += 1
            out_handle.write(f">{title}\n")
            # map the qual... note even with Sanger encoding max 2 digits
            try:
                qualities_strs = [mapping[ascii_] for ascii_ in qual]
            except KeyError:
                invalid_index = _find_index_where(qual, lambda x: x not in mapping)
                assert invalid_index >= 0, "Invalid char not in mapping not found!"
                raise InvalidCharError(
                    qual,
                    invalid_index,
                    details="not in correct range (are you sure you're using the right QualityIO parser?)",
                ) from None
            data = " ".join(qualities_strs)
            while len(data) > 60:
                # Know quality scores are either 1 or 2 digits, so there
                # must be a space in any three consecutive characters.
                if data[60] == " ":
                    out_handle.write(data[:60] + "\n")
                    data = data[61:]
                elif data[59] == " ":
                    out_handle.write(data[:59] + "\n")
                    data = data[60:]
                else:
                    assert data[58] == " ", "Internal logic failure in wrapping"
                    out_handle.write(data[:58] + "\n")
                    data = data[59:]
            out_handle.write(data + "\n")
    return count


def _fastq_sanger_convert_qual(in_file: _TextIOSource, out_file: _TextIOSource) -> int:
    """Fast Sanger FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {chr(q + 33): str(q) for q in range(93 + 1)}
    return _fastq_convert_qual(in_file, out_file, mapping)


def _fastq_solexa_convert_qual(in_file: _TextIOSource, out_file: _TextIOSource) -> int:
    """Fast Solexa FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {
        chr(q + 64): str(int(round(phred_quality_from_solexa(q))))
        for q in range(-5, 62 + 1)
    }
    return _fastq_convert_qual(in_file, out_file, mapping)


def _fastq_illumina_convert_qual(
    in_file: _TextIOSource, out_file: _TextIOSource
) -> int:
    """Fast Illumina 1.3+ FASTQ to QUAL conversion (PRIVATE)."""
    mapping = {chr(q + 64): str(q) for q in range(62 + 1)}
    return _fastq_convert_qual(in_file, out_file, mapping)


@dataclass
class InvalidCharError(ValueError):
    """
    Custom error for strings that have a character that is invalid for whatever reason (eg: non-ascii, invalid range)

    Main attributes:
     - full_string    - the string which contains the invalid character (str)
     - index          - position of the invalid character in full_string (int)
     - details        - additional information to add to the error message. Like: 'not in correct range' (str)
     - r              - how many characters on each side of the invalid character to include in the error message (int)
    """

    full_string: str
    index: int
    details: str
    r: int = 3

    def __str__(self) -> str:
        char = self.full_string[self.index]

        surrounding_characters = self.full_string[
            max(self.index - self.r, 0) : self.index + self.r + 1
        ]
        left_complete = self.index - self.r < 0
        prefix = "" if left_complete else "..."
        right_complete = self.index + self.r + 1 >= len(self.full_string)
        suffix = "" if right_complete else "..."

        return f"Invalid character ({char}) or (0x{char.encode().hex()}) in quality string {self.details} with context: [{prefix}{surrounding_characters}{suffix}]"


def _find_index_where(iterable: Iterable, predicate: Callable[[Any], bool]) -> int:
    for i, x in enumerate(iterable):
        if predicate(x) is True:
            return i

    return -1


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)


class FastqSangerPyfastxIterator(SequenceIterator):
    """Parser for Sanger FASTQ files using pyfastx."""

    modes = "t"

    def __init__(self, source, alphabet=None):
        """Iterate over Fastq records as SeqRecord objects."""
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        try:
            import pyfastx
        except ImportError:
            raise ImportError(
                "Please install pyfastx to use the pyfastx-based FASTQ parser."
            ) from None
        self._iterator = self._create_iterator(source)

    def _create_iterator(self, source):
        import pyfastx
        for name, seq, qual, comment in pyfastx.Fastx(source, comment=True):
            if comment:
                description = f"{name} {comment}"
            else:
                description = name
            qualities = [ord(c) - 33 for c in qual]
            yield SeqRecord(
                Seq(seq),
                id=name,
                description=description,
                letter_annotations={"phred_quality": qualities},
            )

    def __next__(self):
        return next(self._iterator)


class FastqSolexaPyfastxIterator(SequenceIterator):
    """Parser for Solexa FASTQ files using pyfastx."""

    modes = "t"

    def __init__(self, source, alphabet=None):
        """Iterate over Fastq records as SeqRecord objects."""
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        try:
            import pyfastx
        except ImportError:
            raise ImportError(
                "Please install pyfastx to use the pyfastx-based FASTQ parser."
            ) from None
        self._iterator = self._create_iterator(source)

    def _create_iterator(self, source):
        import pyfastx
        for name, seq, qual, comment in pyfastx.Fastx(source, comment=True):
            if comment:
                description = f"{name} {comment}"
            else:
                description = name
            qualities = [ord(c) - 64 for c in qual]
            yield SeqRecord(
                Seq(seq),
                id=name,
                description=description,
                letter_annotations={"solexa_quality": qualities},
            )

    def __next__(self):
        return next(self._iterator)


class FastqIlluminaPyfastxIterator(SequenceIterator):
    """Parser for Illumina FASTQ files using pyfastx."""

    modes = "t"

    def __init__(self, source, alphabet=None):
        """Iterate over Fastq records as SeqRecord objects."""
        if alphabet is not None:
            raise ValueError("The alphabet argument is no longer supported")
        try:
            import pyfastx
        except ImportError:
            raise ImportError(
                "Please install pyfastx to use the pyfastx-based FASTQ parser."
            ) from None
        self._iterator = self._create_iterator(source)

    def _create_iterator(self, source):
        import pyfastx
        for name, seq, qual, comment in pyfastx.Fastx(source, comment=True):
            if comment:
                description = f"{name} {comment}"
            else:
                description = name
            qualities = [ord(c) - 64 for c in qual]
            yield SeqRecord(
                Seq(seq),
                id=name,
                description=description,
                letter_annotations={"phred_quality": qualities},
            )

    def __next__(self):
        return next(self._iterator)
