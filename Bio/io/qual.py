# Copyright 2009-2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Bio.io support for QUAL files and quality conversions."""

import warnings
from math import log

from Bio import BiopythonWarning
from Bio.SeqRecord import SeqRecord


SANGER_SCORE_OFFSET = 33
SOLEXA_SCORE_OFFSET = 64


def solexa_quality_from_phred(phred_quality: float) -> float:
    """Convert a PHRED quality (range 0 to about 90) to a Solexa quality."""
    if phred_quality is None:
        return None
    elif phred_quality > 0:
        return max(-5.0, 10 * log(10 ** (phred_quality / 10.0) - 1, 10))
    elif phred_quality == 0:
        return -5.0
    else:
        raise ValueError(
            f"PHRED qualities must be positive (or zero), not {phred_quality!r}"
        )


def phred_quality_from_solexa(solexa_quality: float) -> float:
    """Convert a Solexa quality (which can be negative) to a PHRED quality."""
    if solexa_quality is None:
        return None
    if solexa_quality < -5:
        warnings.warn(
            f"Solexa quality less than -5 passed, {solexa_quality!r}", BiopythonWarning
        )
    return 10 * log(10 ** (solexa_quality / 10.0) + 1, 10)


def _get_phred_quality(record: SeqRecord) -> list[float] | list[int]:
    """Extract PHRED qualities from a SeqRecord's letter_annotations (PRIVATE)."""
    try:
        return record.letter_annotations["phred_quality"]
    except KeyError:
        pass
    try:
        return [
            phred_quality_from_solexa(q)
            for q in record.letter_annotations["solexa_quality"]
        ]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None


_phred_to_sanger_quality_str = {
    qp: chr(min(126, qp + SANGER_SCORE_OFFSET)) for qp in range(93 + 1)
}
_solexa_to_sanger_quality_str = {
    qs: chr(min(126, int(round(phred_quality_from_solexa(qs))) + SANGER_SCORE_OFFSET))
    for qs in range(-5, 93 + 1)
}


def _get_sanger_quality_str(record: SeqRecord) -> str:
    """Return a Sanger FASTQ encoded quality string (PRIVATE)."""
    try:
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        pass
    else:
        try:
            return "".join(_phred_to_sanger_quality_str[qp] for qp in qualities)
        except KeyError:
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 93.5:
            warnings.warn(
                "Data loss - max PHRED quality 93 in Sanger FASTQ", BiopythonWarning
            )
        return "".join(
            chr(min(126, int(round(qp)) + SANGER_SCORE_OFFSET)) for qp in qualities
        )
    try:
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    try:
        return "".join(_solexa_to_sanger_quality_str[qs] for qs in qualities)
    except KeyError:
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    if max(qualities) >= 93.5:
        warnings.warn(
            "Data loss - max PHRED quality 93 in Sanger FASTQ", BiopythonWarning
        )
    return "".join(
        chr(min(126, int(round(phred_quality_from_solexa(qs))) + SANGER_SCORE_OFFSET))
        for qs in qualities
    )


assert 62 + SOLEXA_SCORE_OFFSET == 126
_phred_to_illumina_quality_str = {
    qp: chr(qp + SOLEXA_SCORE_OFFSET) for qp in range(62 + 1)
}
_solexa_to_illumina_quality_str = {
    qs: chr(int(round(phred_quality_from_solexa(qs))) + SOLEXA_SCORE_OFFSET)
    for qs in range(-5, 62 + 1)
}


def _get_illumina_quality_str(record: SeqRecord) -> str:
    """Return an Illumina 1.3 to 1.7 FASTQ encoded quality string (PRIVATE)."""
    try:
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        pass
    else:
        try:
            return "".join(_phred_to_illumina_quality_str[qp] for qp in qualities)
        except KeyError:
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 62.5:
            warnings.warn(
                "Data loss - max PHRED quality 62 in Illumina FASTQ", BiopythonWarning
            )
        return "".join(
            chr(min(126, int(round(qp)) + SOLEXA_SCORE_OFFSET)) for qp in qualities
        )
    try:
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    try:
        return "".join(_solexa_to_illumina_quality_str[qs] for qs in qualities)
    except KeyError:
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    if max(qualities) >= 62.5:
        warnings.warn(
            "Data loss - max PHRED quality 62 in Illumina FASTQ", BiopythonWarning
        )
    return "".join(
        chr(min(126, int(round(phred_quality_from_solexa(qs))) + SOLEXA_SCORE_OFFSET))
        for qs in qualities
    )


assert 62 + SOLEXA_SCORE_OFFSET == 126
_solexa_to_solexa_quality_str = {
    qs: chr(min(126, qs + SOLEXA_SCORE_OFFSET)) for qs in range(-5, 62 + 1)
}
_phred_to_solexa_quality_str = {
    qp: chr(min(126, int(round(solexa_quality_from_phred(qp))) + SOLEXA_SCORE_OFFSET))
    for qp in range(62 + 1)
}


def _get_solexa_quality_str(record: SeqRecord) -> str:
    """Return a Solexa FASTQ encoded quality string (PRIVATE)."""
    try:
        qualities = record.letter_annotations["solexa_quality"]
    except KeyError:
        pass
    else:
        try:
            return "".join(_solexa_to_solexa_quality_str[qs] for qs in qualities)
        except KeyError:
            pass
        if None in qualities:
            raise TypeError("A quality value of None was found")
        if max(qualities) >= 62.5:
            warnings.warn(
                "Data loss - max Solexa quality 62 in Solexa FASTQ", BiopythonWarning
            )
        return "".join(
            chr(min(126, int(round(qs)) + SOLEXA_SCORE_OFFSET)) for qs in qualities
        )
    try:
        qualities = record.letter_annotations["phred_quality"]
    except KeyError:
        raise ValueError(
            "No suitable quality scores found in "
            "letter_annotations of SeqRecord (id=%s)." % record.id
        ) from None
    try:
        return "".join(_phred_to_solexa_quality_str[qp] for qp in qualities)
    except KeyError:
        pass
    if None in qualities:
        raise TypeError("A quality value of None was found")
    if max(qualities) >= 62.5:
        warnings.warn(
            "Data loss - max Solexa quality 62 in Solexa FASTQ", BiopythonWarning
        )
    return "".join(
        chr(min(126, int(round(solexa_quality_from_phred(qp))) + SOLEXA_SCORE_OFFSET))
        for qp in qualities
    )
