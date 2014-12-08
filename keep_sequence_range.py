#!/usr/bin/python
"""keep_sequence_range.py return only those fastq reads within a specified
size range. NOT SAFE FOR PAIRED-END DATA!!!
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

from commonIO import read_fastq_chunk
from commonIO import CustomParser
import sys

CHUNK = 4096  # bytes of input read per IO call with read_chunk


def get_commandline_args():
    """Command-line interface for create_alignment_db.py"""
    parser = CustomParser(
        description='''keep_sequence_range.py return only those fastq reads within a specified
size range. NOT SAFE FOR PAIRED-END DATA!!!

Copyright (C) 2014 Joy-El R.B. Talbot under the GNU General Public License version 3''')

    parser.add_argument("-i", "--input",
                        help="FASTQ format reads, omit to read from commandline",
                        metavar="FASTQ")
    parser.add_argument("-o", "--output",
                        help="name of output file, omit to write to commandline",
                        metavar="OUT")
    parser.add_argument("--minimum",
                        help="minimum length to keep, default = 16",
                        metavar="NT",
                        type=int,
			default=16)
    parser.add_argument("--maximum",
                        help="maximum length to keep, default = 35",
                        metavar="NT",
                        type=int,
			default=35)
    arguments = parser.parse_args()

    if arguments.input is None:
        use_stdin = True
        sys.stderr.write("Reading input from STDIN...\n")
        sys.stderr.flush()
    else:
        use_stdin = False

    if arguments.output is None:
        use_stdout = True
        sys.stderr.write("Writing output to STDOUT...\n")
        sys.stderr.flush()
    else:
        use_stdout = False

    return (arguments, use_stdin, use_stdout)


def keep_sequence_range(open_fastq_file, minimum_size, maximum_size, chunk_size):
    """Return chunks of fastq data to keep."""
    output_chunks = ""
    max_reads = chunk_size / 10
    kept_reads = 0
    for read in read_fastq_chunk(open_fastq_file, chunk_size):
        if minimum_size <= len(read[1]) <= maximum_size:
            output_chunks += "\n".join(read) + "\n"
            kept_reads += 1
            if kept_reads > max_reads:
                yield output_chunks
                output_chunks = ""
                kept_reads = 0
    yield output_chunks


if __name__=="__main__":
    (args, input_from_stdin, output_to_stdout) = get_commandline_args()
    try:
        if input_from_stdin:
            in_file = sys.stdin
            for group in keep_sequence_range(in_file, args.minimum, args.maximum, CHUNK):
                if output_to_stdout:
                    print group.strip()
                else:
                    with open(args.output, 'a') as output:
                        output.write(group)
        else:
            with open(args.input) as in_file:
                for group in keep_sequence_range(in_file, args.minimum, args.maximum, CHUNK):
                    if output_to_stdout:
                        print group.strip()
                    else:
                        with open(args.output, 'a') as output:
                            output.write(group)
    except IOError as error:
        sys.stderr.write("Could not open file: {}\n".format(error.filename))
        sys.stderr.flush()
        raise IOError(error)
