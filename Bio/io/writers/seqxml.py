# Copyright 2024 by The Biopython Contributors. All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

from ..Interfaces import SequenceWriter


class SeqXmlWriter(SequenceWriter):
    """Writes SeqRecords into seqXML file.

    SeqXML requires the SeqRecord annotations to specify the molecule_type;
    the molecule type is required to contain the term "DNA", "RNA", or
    "protein".
    """

    modes = "b"

    def __init__(
        self, target, source=None, source_version=None, species=None, ncbiTaxId=None
    ):
        """Create Object and start the xml generator.

        Arguments:
         - target - Output stream opened in binary mode, or a path to a file.
         - source - The source program/database of the file, for example
           UniProt.
         - source_version - The version or release number of the source
           program or database from which the data originated.
         - species - The scientific name of the species of origin of all
           entries in the file.
         - ncbiTaxId - The NCBI taxonomy identifier of the species of origin.

        """
        super().__init__(target)
        handle = self.handle
        self.xml_generator = XMLGenerator(handle, "utf-8")
        self.xml_generator.startDocument()
        self.source = source
        self.source_version = source_version
        self.species = species
        self.ncbiTaxId = ncbiTaxId

    def write_records(self, records):
        """Write records to the output file, and return the number of records.

        records - A list or iterator returning SeqRecord objects
        """
        self.write_header()
        count = super().write_records(records)
        self.write_footer()
        return count

    def write_header(self):
        """Write root node with document metadata."""
        attrs = {
            "xmlns:xsi": "http://www.w3.org/2001/XMLSchema-instance",
            "xsi:noNamespaceSchemaLocation": "http://www.seqxml.org/0.4/seqxml.xsd",
            "seqXMLversion": "0.4",
        }

        if self.source is not None:
            attrs["source"] = self.source
        if self.source_version is not None:
            attrs["sourceVersion"] = self.source_version
        if self.species is not None:
            if not isinstance(self.species, str):
                raise TypeError("species should be of type string")
            attrs["speciesName"] = self.species
        if self.ncbiTaxId is not None:
            if not isinstance(self.ncbiTaxId, (str, int)):
                raise TypeError("ncbiTaxID should be of type string or int")
            attrs["ncbiTaxID"] = self.ncbiTaxId

        self.xml_generator.startElement("seqXML", AttributesImpl(attrs))

    def write_record(self, record):
        """Write one record."""
        if not record.id or record.id == "<unknown id>":
            raise ValueError("SeqXML requires identifier")

        if not isinstance(record.id, str):
            raise TypeError("Identifier should be of type string")

        attrb = {"id": record.id}

        if (
            "source" in record.annotations
            and self.source != record.annotations["source"]
        ):
            if not isinstance(record.annotations["source"], str):
                raise TypeError("source should be of type string")
            attrb["source"] = record.annotations["source"]

        self.xml_generator.startElement("entry", AttributesImpl(attrb))
        self._write_species(record)
        self._write_description(record)
        self._write_seq(record)
        self._write_dbxrefs(record)
        self._write_properties(record)
        self.xml_generator.endElement("entry")

    def write_footer(self):
        """Close the root node and finish the XML document."""
        self.xml_generator.endElement("seqXML")
        self.xml_generator.endDocument()

    def _write_species(self, record):
        """Write the species if given (PRIVATE)."""
        local_ncbi_taxid = None
        if "ncbi_taxid" in record.annotations:
            local_ncbi_taxid = record.annotations["ncbi_taxid"]
            if isinstance(local_ncbi_taxid, list):
                # SwissProt parser uses a list (which could cope with chimeras)
                if len(local_ncbi_taxid) == 1:
                    local_ncbi_taxid = local_ncbi_taxid[0]
                elif len(local_ncbi_taxid) == 0:
                    local_ncbi_taxid = None
                else:
                    raise ValueError(
                        "Multiple entries for record.annotations['ncbi_taxid'], %r"
                        % local_ncbi_taxid
                    )
        if "organism" in record.annotations and local_ncbi_taxid:
            local_org = record.annotations["organism"]

            if not isinstance(local_org, str):
                raise TypeError("organism should be of type string")

            if not isinstance(local_ncbi_taxid, (str, int)):
                raise TypeError("ncbiTaxID should be of type string or int")

            # The local species definition is only written if it differs from the global species definition
            if local_org != self.species or local_ncbi_taxid != self.ncbiTaxId:
                attr = {"name": local_org, "ncbiTaxID": str(local_ncbi_taxid)}
                self.xml_generator.startElement("species", AttributesImpl(attr))
                self.xml_generator.endElement("species")

    def _write_description(self, record):
        """Write the description if given (PRIVATE)."""
        if record.description:
            if not isinstance(record.description, str):
                raise TypeError("Description should be of type string")

            description = record.description
            if description == "<unknown description>":
                description = ""

            if len(record.description) > 0:
                self.xml_generator.startElement("description", AttributesImpl({}))
                self.xml_generator.characters(description)
                self.xml_generator.endElement("description")

    def _write_seq(self, record):
        """Write the sequence (PRIVATE).

        Note that SeqXML requires the molecule type to contain the term
        "DNA", "RNA", or "protein".
        """
        seq = bytes(record.seq)

        if not len(seq) > 0:
            raise ValueError("The sequence length should be greater than 0")

        molecule_type = record.annotations.get("molecule_type")
        if molecule_type is None:
            raise ValueError("molecule_type is not defined")
        elif "DNA" in molecule_type:
            seqElem = "DNAseq"
        elif "RNA" in molecule_type:
            seqElem = "RNAseq"
        elif "protein" in molecule_type:
            seqElem = "AAseq"
        else:
            raise ValueError(f"unknown molecule_type '{molecule_type}'")

        self.xml_generator.startElement(seqElem, AttributesImpl({}))
        self.xml_generator.characters(seq)
        self.xml_generator.endElement(seqElem)

    def _write_dbxrefs(self, record):
        """Write all database cross references (PRIVATE)."""
        if record.dbxrefs is not None:
            for dbxref in record.dbxrefs:
                if not isinstance(dbxref, str):
                    raise TypeError("dbxrefs should be of type list of string")
                if dbxref.find(":") < 1:
                    raise ValueError(
                        "dbxrefs should be in the form ['source:id', 'source:id' ]"
                    )

                dbsource, dbid = dbxref.split(":", 1)

                attr = {"source": dbsource, "id": dbid}
                self.xml_generator.startElement("DBRef", AttributesImpl(attr))
                self.xml_generator.endElement("DBRef")

    def _write_properties(self, record):
        """Write all annotations that are key value pairs with values of a primitive type or list of primitive types (PRIVATE)."""
        for key, value in record.annotations.items():
            if key not in ("organism", "ncbi_taxid", "source"):
                if value is None:
                    attr = {"name": key}
                    self.xml_generator.startElement("property", AttributesImpl(attr))
                    self.xml_generator.endElement("property")

                elif isinstance(value, list):
                    for v in value:
                        if v is None:
                            attr = {"name": key}
                        else:
                            attr = {"name": key, "value": str(v)}
                        self.xml_generator.startElement(
                            "property", AttributesImpl(attr)
                        )
                        self.xml_generator.endElement("property")

                elif isinstance(value, (int, float, str)):
                    attr = {"name": key, "value": str(value)}
                    self.xml_generator.startElement("property", AttributesImpl(attr))
                    self.xml_generator.endElement("property")
