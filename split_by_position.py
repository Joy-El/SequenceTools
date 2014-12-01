#!/usr/bin/python
"""split_by_position.py splits a BED-like file by chromosome and start position.
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
import sys

CHUNK = 4096  # bytes of input read per IO call with read_chunk
SPLIT = 10000000  # bases per unit


def split_by_position(bed_like_file, base_chunk):
    """Split a file into several subfiles by chromosome and start position.
    """
    for line in read_chunk(bed_like_file, CHUNK):
        parts = line.split("\t")
        outfile_name = "{}_{}".format(parts[0], (int(parts[1]) / base_chunk))
        with open(outfile_name, 'a') as outfile:
            outfile.write(line + "\n")


if __name__ == "__main__":
    infile = sys.stdin
    split_by_position(infile, SPLIT)