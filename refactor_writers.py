import os
import re

io_dir = "Bio/io"
writers_dir = "Bio/io/writers"

# List of files to process
files_to_process = [
    "phd.py",
    "pir.py",
    "seqxml.py",
    "sff.py",
    "tab.py",
    "xdna.py",
    "insdc.py",
]

for filename in files_to_process:
    # Read the content of the parser file
    with open(os.path.join(io_dir, filename), "r") as f:
        content = f.read()

    # Find all writer classes
    writer_classes = re.findall(r"class\s+\w+Writer\(SequenceWriter\):(?s:.)*?if __name__ == \"__main__\":", content)
    if not writer_classes:
        writer_classes = re.findall(r"class\s+\w+Writer\(SequenceWriter\):(?s:.)*$", content)


    if writer_classes:
        # Create the new writer file
        writer_filename = os.path.join(writers_dir, filename)
        with open(writer_filename, "w") as f:
            f.write("# Copyright 2024 by The Biopython Contributors. All rights reserved.\n")
            f.write("#\n")
            f.write("# This file is part of the Biopython distribution and governed by your\n")
            f.write("# choice of the \"Biopython License Agreement\" or the \"BSD 3-Clause License\".\n")
            f.write("# Please see the LICENSE file that should have been included as part of this\n")
            f.write("# package.\n\n")
            # Add necessary imports
            f.write("from ..Interfaces import SequenceWriter\n")
            if filename == "insdc.py":
                f.write("from Bio import SeqFeature\n")
                f.write("from Bio import BiopythonWarning\n")
                f.write("from datetime import datetime, date as datetime_date\n")
                f.write("from string import ascii_letters, digits\n")
            f.write("\n\n")
            f.write("\n".join(writer_classes))

        # Remove the writer classes from the original file
        new_content = content
        for writer_class in writer_classes:
            new_content = new_content.replace(writer_class, "")

        with open(os.path.join(io_dir, filename), "w") as f:
            f.write(new_content)

# Update Bio/io/__init__.py
with open(os.path.join(io_dir, "__init__.py"), "r") as f:
    init_content = f.read()

new_imports = ""
new_writers = ""

for filename in files_to_process:
    module_name = os.path.splitext(filename)[0]
    new_imports += f"from .writers import {module_name} as {module_name}_writer\n"

    if module_name == "insdc":
        new_writers += f'    "gb": {module_name}_writer.GenBankWriter,\n'
        new_writers += f'    "genbank": {module_name}_writer.GenBankWriter,\n'
        new_writers += f'    "embl": {module_name}_writer.EmblWriter,\n'
        new_writers += f'    "imgt": {module_name}_writer.ImgtWriter,\n'
    else:
        writer_class_name = module_name.capitalize() + "Writer"
        if module_name == "seqxml":
            writer_class_name = "SeqXmlWriter"
        new_writers += f'    "{module_name}": {module_name}_writer.{writer_class_name},\n'


init_content = re.sub(r"(# Registry of writers\n_FormatToWriter = {\n)", f"\\1{new_writers}", init_content)

# Add imports
init_content = re.sub(r"(from \.writers import fastq as fastq_writer\n)", f"\\1{new_imports}", init_content)


with open(os.path.join(io_dir, "__init__.py"), "w") as f:
    f.write(init_content)

print("Refactoring complete.")
