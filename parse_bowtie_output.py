__author__ = 'Joy-El R.B. Talbot'
"""Functions to parse Bowtie output for collapse and summary.
Copyright (C) 2014 Joy-El R.B. Talbot

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (LICENSE).
    If not, see <http://www.gnu.org/licenses/>

For details on Bowtie please see their website:
    http://bowtie-bio.sourceforge.net/index.shtml"""

from commonIO import read_chunk


def add_multimapping_tally(open_bowtie_file, chunk_size=2048):
    """Count the number of mappings of a read tag.

    ASSUMES: file is sorted by read tag name
    Either run on unmodified bowtie output or re-sort
    on the first column of the SAM formatted output:
        To preserve the header lines:
            samtools view -SH bowtiefile.sam > sortedbowtiefile.sam
        To add sorted alignment lines:
            samtools view -S bowtiefile.sam | sort -k1,1 >> sortedbowtiefile.sam
    Note that sorting may take a lot of resource to do, that's why it is best
    to add multimapping tally BEFORE any other operations are done to the
    Bowtie output."""
    master_read = ""
    mapping_count = 0
    saved_reads = []
    for alignment in read_chunk(open_bowtie_file, chunk_size):
        columns = alignment.split("\t")
        if columns[0] == master_read:
            mapping_count += 1
            saved_reads.append(alignment)
        else:
            if master_read != "":
                for read in saved_reads:
                    yield "{}\tNH:i:{}".format(read, mapping_count)
            master_read = columns[0]
            saved_reads = [alignment]
            if (int(columns[1]) & 0x4) == 0x4:
                # unmapped read according to SAM 0x4 flag in second column
                mapping_count = 0
            else:
                # mapped read
                mapping_count = 1


