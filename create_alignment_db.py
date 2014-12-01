#!/usr/bin/python
"""create_alignment_db.py creates a set of SQLite3 databases from
an alignment file (currently Bowtie output in SAM format)
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
    If not, see <http://www.gnu.org/licenses/>"""

__author__ = 'Joy-El R.B. Talbot'

from commonIO import read_chunk
from commonIO import CustomParser
import datetime
import sys
import re

CHUNK = 4096  # bytes of input read per IO call with read_chunk

def get_commandline_args():
    """Command-line interface for create_alignment_db.py"""
    parser = CustomParser(
        description='''create_alignment_db.py creates a set of SQLite3 databases from
an alignment file (currently Bowtie output in SAM format)

Copyright (C) 2014 Joy-El R.B. Talbot under the GNU General Public License version 3''')

    parser.add_argument("-i", "--input",
                        help="bowtie file (SAM format), omit to read from commandline",
                        metavar="SAM")
    parser.add_argument("-l", "--library_name",
                        help="Name of library",
                        metavar="NAME",
                        type=str,
                        required=True)
    parser.add_argument("-d", "--database_prefix",
                        help="Prefix for database names; default: current datetime as YMD-HMS",
                        metavar="NAME",
                        type=str,
                        default=datetime.datetime.now().strftime('%y%m%d-%H%M%S'))
    arguments = parser.parse_args()

    if arguments.input is None:
        use_stdin = True
        sys.stderr.write("Reading input from STDIN...\n")
        sys.stderr.flush()
    else:
        use_stdin = False

    return (arguments, use_stdin)


def write_data(string, filename):
    """Write string to file."""
    with open(filename, "a") as output:
        output.write(string + "\n")
    return True


def get_mismatches(flags):
    """Return number of mismatches from MD flag. See SAM format for details."""
    for flag in flags:
        try:
            MD_flag = re.search("MD:Z:([ATCGNatcgn0-9]+)", flag).groups()[0]
        except AttributeError:
            next
        else:
            return len(re.split("[ATCGNatcgn]+", MD_flag))-1


def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence"""
    complement = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N",
                  "a":"T", "t":"A", "g":"C", "c":"G", "n":"N"}  # to account for case
    # creates the complement as a list
    complement_sequence = map(lambda x: complement[x], sequence)
    # reverse and return to a string
    return "".join(complement_sequence[::-1])


def parse_alignment(sam_line):
    """Parse the information in a SAM formatted alignment.

    Returns a tuple of:
        name as string,
        mapped as boolean,
        mismatch_count as int,
        tag as string (sequence),
        position as tab-delimited string of:
            chromosome
            start position (0-based)
            end position (1-based) so the flat file is BED-like
            strand (+, -, or .)"""
    parts = sam_line.strip().split("\t")
    name = parts[0]
    tag_sequence = parts[9]
    if (int(parts[1]) & 0x4) == 0x4:
        # read did not map, can not get info fro mismatch_count or position_string
        mapped = False
        mismatch_count = None
        position_string = None
        strand = None
    else:
        mapped = True
        strand = "+"  # default value
        if (int(parts[1]) & 0x10) == 0x10:
            # aligned in the antisense direction
            strand = "-"
            tag_sequence = reverse_complement(tag_sequence)
        mismatch_count = get_mismatches(parts[11:])
        start = int(parts[3]) - 1  # to make it 0-based
        end = start + len(tag_sequence)  # 1-based
        position_string = "{}\t{}\t{}".format(parts[2], start, end)
    return (name, mapped, mismatch_count, tag_sequence, position_string, strand)


def create_alignment_db(sam_openfile, library_name, database_prefix):
    """Create alignment SQLite3 databases representing alignment data from Bowtie SAM file.

    Database files:
        {prefix}_tagloci.db: partition on (p) Chromosome,
                             index and partition on (ip) Start(0-based),
                             End (0-based),
                             Strand (+,-,or .),
                             Tag Sequence,
                             Number of Mismatches
        {prefix}_{library}.db: (ip) Tag Sequence,
                                Abundance
        {prefix}_tags.db: (ip) Tag Sequence,
                          Total Mappings,
                          Mappings with 0 mismatches (perfect),
                          Mappings with 1 mismatch,
                          Mappings with 2 mismatches"""
    #TODO create_chromosome_db(sam_openfile, database_prefix)
    tagloci = "{}_tagloci".format(database_prefix)
    library = "{}_{}".format(database_prefix, library_name)
    tags = "{}_tags".format(database_prefix)

    # scan through header lines
    header = sam_openfile.readline()
    while header[0] == "@":
        header = sam_openfile.readline()

    mismatch_tally = [0, 0, 0]
    # parse initial alignment
    (read_name, maps, mismatches, tag, position, strand) = parse_alignment(header.strip())
    if maps:
        write_data("{}\t{}\t{}\t{}".format(position, tag, mismatches, strand), "{}.data".format(tagloci))
        mismatch_tally[mismatches] += 1
    last_read_name = read_name
    last_tag = tag

    for alignment in read_chunk(sam_openfile, CHUNK):
        (read_name, maps, mismatches, tag, position, strand) = parse_alignment(alignment)
        if maps:  # don't process unmapped reads
            write_data("{}\t{}\t{}\t{}".format(position, tag, mismatches, strand), "{}.data".format(tagloci))
            if read_name == last_read_name:
                mismatch_tally[mismatches] += 1
            else:
                write_data("{}\t1".format(last_tag), "{}.data".format(library))
                # prepare output for tags database
                mismatch_string = "\t".join(str(m) for m in mismatch_tally)
                write_data("{}\t{}\t{}".format(last_tag, sum(mismatch_tally), mismatch_string),
                           "{}.data".format(tags))
                # reset for next round
                last_read_name = read_name
                last_tag = tag
                mismatch_tally = [0, 0, 0]
                mismatch_tally[mismatches] += 1
    # write out last tag
    write_data("{}\t1".format(last_tag), "{}.data".format(library))
    # prepare output for tags database
    mismatch_string = "\t".join(str(m) for m in mismatch_tally)
    write_data("{}\t{}\t{}".format(last_tag, sum(mismatch_tally), mismatch_string),
               "{}.data".format(tags))


if __name__ == "__main__":
    (args, input_from_stdin) = get_commandline_args()
    try:
        if input_from_stdin:
            in_file = sys.stdin
            create_alignment_db(in_file, args.library_name, args.database_prefix)
        else:
            with open(args.input) as in_file:
                create_alignment_db(in_file, args.library_name, args.database_prefix)
    except IOError as error:
        sys.stderr.write("Could not open file: {}\n".format(error.filename))
        sys.stderr.flush()
        raise IOError(error)
