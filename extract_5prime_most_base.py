#!/usr/bin/python
"""extract_5prime_most_base.py extracts the chromosomal coordinates of the 5' most base for each read alignment.
Can currently handle SAM and BED formats; will return a BED format file.
Copyright (C) 2015 Joy-El R.B. Talbot

    This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
    later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program (LICENSE).
    If not, see <http://www.gnu.org/licenses/>"""

__author__ = 'Joy-El R.B. Talbot'

import argparse
import re
import sys
from commonIO import CustomParser
from commonIO import read_chunk

CHUNK_SIZE = 4096  # bytes of input read per IO call with read_chunk


class ReadError(Exception):
    """Create Read object error."""
    def __init__(self, message, name='general'):
        self.message = message
        self.name = name

    def __str__(self):
        return 'ReadError - {}: {}'.format(self.name, self.message)


class Read(object):
    """Describes a read based on its chromosomal coordinates"""

    __slots__ = ['chromosome', 'strand', 'name', 'positions']

    input_format = None

    def __init__(self, read_string):
        """Create read object from an alignment.

        self.positions is an array of chromosomal positions from the 5' to 3' ends of the read in a 0-based format"""
        assigned = False
        attempt = 0
        while not assigned and attempt < 3:
            try:
                (self.chromosome,
                 self.strand,
                 self.name,
                 self.positions) = Read.parse_read_string(read_string, Read.input_format)
                assigned = True
            except ReadError as _error:
                if _error.name == 'unmapped':
                    raise ReadError(_error.message, _error.name)
                else:
                    Read.update_input_format(read_string)
                    attempt += 1
            except AttributeError:
                Read.update_input_format(read_string)
                attempt += 1

    def __str__(self):
        return '{}\t{}\t{}\t{}\t0\t{}'.format(self.chromosome,
                                              self.positions[0],
                                              self.positions[-1] + 1,
                                              self.name,
                                              self.strand)

    def print_first_base(self):
        return '{}\t{}\t{}\t{}\t0\t{}'.format(self.chromosome,
                                              self.positions[0],
                                              self.positions[0] + 1,
                                              self.name,
                                              self.strand)

    @classmethod
    def parse_read_string(cls, read_string, read_format):
        """Parse a read string into chromosome, strand, name and positions according to input_format."""
        if read_format == 'BED':
            return cls.parse_BED_read_string(read_string)
        elif read_format == 'SAM':
            return cls.parse_SAM_read_string(read_string)
        else:
            raise ReadError('Can not currently parse read strings of format: {}'.format(read_format))

    @classmethod
    def parse_BED_read_string(cls, read_string):
        """Parses a BED formatted read."""
        BED_data = re.search("^(\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t\S+\t([+-.])", read_string).groups()
        (chromosome, start, end, name, strand) = BED_data
        start = int(start)  # already 0-based
        end = int(end) - 1  # convert from 1-based to 0-based
        if strand == "-":
            positions = range(end, start - 1, -1)  # start - 1 to get all inclusive
        else:  # assume '+' strand if strand is not given ('.')
            positions = range(start, end + 1)  # end + 1 to get all inclusive start to end
        return chromosome, strand, name, positions

    @classmethod
    def parse_SAM_read_string(cls, read_string):
        """Parses a SAM formatted read."""
        # important bitwise tags:
        unmapped = 0x4
        antisense = 0x10
        SAM_data = re.search("^(\S+)\t([0-9]+)\t(\S+)\t([0-9]+)\t[0-9]+\t(\S+)\t\S+\t[0-9]+\t[0-9]+\t(\S+)\t", read_string).groups()
        (name, bitstring, chromosome, start, cigar, sequence) = SAM_data
        bitstring = int(bitstring)
        if (bitstring & unmapped) == unmapped:
            raise ReadError('This read does not map: {}'.format(name), 'unmapped')
        else:
            start = int(start) - 1  # convert from 1-based to 0-based
            length = len(sequence)
            positions = cls.parse_cigar_string(cigar, start, length)
            if (bitstring & antisense) == antisense:
                strand = "-"
                positions.reverse()  # make first position the read start (5' most base)
            else:
                strand = "+"
        return chromosome, strand, name, positions

    @classmethod
    def parse_cigar_string(cls, cigar, read_start, sequence_length):
        """Parses a CIGAR string to return an array of positions covered by the read."""
        cigar_codes = {'M': {'count': True,  'advance': True},   # alignment match (can be either sequence match or mismatch)
                       'I': {'count': False, 'advance': False},  # insertion to the reference
                       'D': {'count': False, 'advance': True},   # deletion from the reference
                       'N': {'count': False, 'advance': True},   # skipped region from the reference
                       'S': {'count': False, 'advance': False},  # soft clipping (clipped sequences present in SEQ)
                       'H': {'count': False, 'advance': False},  # hard clipping (clipped sequences NOT present in SEQ)
                       'P': {'count': False, 'advance': True},   # padding (silent deletion from padded reference)
                       '=': {'count': True,  'advance': True},   # sequence match
                       'X': {'count': True,  'advance': True}}   # sequence mismatch

        # sometimes the CIGAR value is a sole "*", in which case assume a perfect match
        if cigar == "*":
            positions = range(read_start, read_start + sequence_length)  # will give start to end inclusive
        else:
            positions = []
            current_position = read_start
            # separate CIGAR string into nucleotide counts and CIGAR codes
            cigar_entries = re.findall('(\d+)([{}])'.format(''.join(cigar_codes.keys())), cigar)
            for (nucleotide_length, code) in cigar_entries:
                nucleotide_length = int(nucleotide_length)
                if cigar_codes[code]['count']:
                    # add nucleotide_length of positions to the positions array
                    current_end = current_position + nucleotide_length
                    positions += range(current_position, current_end)
                    current_position = current_end
                elif cigar_codes[code]['advance']:
                    # advance the current_position but do not add to the positions array
                    current_position += nucleotide_length
                else:
                    # neither advance the current_position nor add to the positions array
                    pass
        return positions

    @classmethod
    def update_input_format(cls, read_string):
        """Determine the format of a read via a match to a regular expression.  Currently BED vs. SAM"""
        string = r'\S+'
        integer = r'[0-9]+'
        strand = r'[+-.]'
        possible_formats = {'BED': '\t'.join([string, integer, integer, string, string, strand]),
                            'SAM': '\t'.join([string, integer, string, integer, integer, string, string, integer, integer, string, string])}
        unknown_format = True
        for possible_format in possible_formats:
            if re.match(possible_formats[possible_format], read_string) is not None:
                unknown_format = False
                cls.input_format = possible_format
        if unknown_format:
            raise ReadError('''Could not determine the read format of string: {}
            Currently stored formats are:
                {}'''.format(read_string, possible_formats))


def get_arguments():
    """Command-line interface for extract_5prime-most-base.py"""
    parser = CustomParser(
        description='''extract_5prime_most_base.py extracts the chromosomal coordinates of
the 5'-most base for each read alignment. Can currently handle SAM and BED
formats; will return a BED format file.

Copyright (C) 2015 Joy-El R.B. Talbot under
the GNU General Public License version 3

''')

    parser.add_argument('-i', '--input',
                        help='SAM or BED file of read alignments, omit to read from commandline',
                        metavar='ALIGNMENTS')
    parser.add_argument('-o', '--output',
                        help="Name for output BED file of 5'-most bases, omit to write to commandline",
                        metavar="5'-MOST BASES")
    parser.add_argument('--use_stdin',
                        help=argparse.SUPPRESS,
                        default=True)
    parser.add_argument('--use_stdout',
                        help=argparse.SUPPRESS,
                        default=True)

    arguments = parser.parse_args()

    if arguments.input is not None:
        arguments.use_stdin = False
    else:
        sys.stderr.write('Reading input from STDIN...\n')
        sys.stderr.flush()

    if arguments.output is not None:
        arguments.use_stdout = False
    else:
        sys.stderr.write('Writing output to STDOUT...\n')
        sys.stderr.flush()

    return arguments


def confirm_new_file(filename):
    """Confirm that an output filename does not already exist."""
    try:
        with open(filename) as existing_file:
            pass
    except IOError:
        sys.stderr.write('Confirmed that output file {} does not already exist.\n'.format(filename))
        sys.stderr.flush()
    else:
        sys.stderr.write('''Output file {} already exists! \n
Please choose a different filename or delete the existing file.\n'''.format(filename))
        sys.stderr.flush()
        sys.exit(1)


def extract_5prime_most_base(alignments_source, output_to_stdout, output_filename):
    """Extract the 5'-most base from each alignment."""
    for alignment in read_chunk(alignments_source, CHUNK_SIZE):
        if alignment[0] != "@" and alignment[0] != "#":  # skip any header/comment lines
            try:
                read = Read(alignment)
                if output_to_stdout:
                    print read.print_first_base()
                else:
                    with open(output_filename, 'a') as output:
                        output.write('{}\n'.format(read.print_first_base()))
            except ReadError as _error:
                if _error.name != 'unmapped':  # silently skip unmapped reads only
                    raise ReadError(_error.message, _error.name)


if __name__ == '__main__':
    args = get_arguments()
    if not args.use_stdout:
        confirm_new_file(args.output)
    if args.use_stdin:
        input_ = sys.stdin
        extract_5prime_most_base(input_, args.use_stdout, args.output)
    else:
        try:
            with open(args.input) as input_:
                extract_5prime_most_base(input_, args.use_stdout, args.output)
        except IOError as error:
            sys.stderr.write('Could not open alignment file: {}\n'.format(error.filename))
            sys.stderr.flush()
            raise IOError(error)