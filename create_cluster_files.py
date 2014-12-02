#!/usr/bin/python
"""create_cluster_files.py creates flat data files for clusters and
cluster-tag pairs from merged loci.data files. These files can then
be imported into SQLite.
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
        description='''create_cluster_files.py creates flat data files for clusters and
cluster-tag pairs from merged loci.data files.

Copyright (C) 2014 Joy-El R.B. Talbot under the GNU General Public License version 3''')

    parser.add_argument("-i", "--input",
                        help="merged loci file",
                        metavar="LOCI")
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


def get_unique_tags(all_tags):
    """Return a list of unique tags."""
    return list(set(all_tags))


def create_cluster_files(loci_openfile, library_name, database_prefix):
    """Create cluster and cluster-tag files from merged loci.

    Database files:
        {prefix}_clusters.db: partition on (p) Chromosome,
                             index and partition on (ip) Start(0-based),
                             End (0-based),
                             cluster_name (unique),
                             Number of tags,
                             Strand (+,-,or .)
        {prefix}_clustertags.db: cluster_name,
                                tag_sequence
                                (together the two will be unique)"""
    clusters = "{}_clusters".format(database_prefix)
    clustertags = "{}_clustertags".format(database_prefix)

    cluster_index = 1

    for cluster in read_chunk(loci_openfile, CHUNK):
        parts = cluster.split("\t")
        tags = parts[3].split(";")
        unique_tags = get_unique_tags(tags)
        if len(parts) == 6:
            strand = parts[-1]
        else:
            strand = "."
        with open("{}.data".format(clusters), 'a') as output:
            output.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(parts[0],
                                                           parts[1],
                                                           parts[2],
                                                           "{}_{}".format(library_name, cluster_index),
                                                           len(unique_tags),
                                                           strand))
        with open("{}.data".format(clustertags), 'a') as output:
            for tag in unique_tags:
                output.write("{}\t{}\n".format("{}_{}".format(library_name, cluster_index),
                                               tag))
        cluster_index += 1


if __name__ == "__main__":
    (args, input_from_stdin) = get_commandline_args()
    try:
        if input_from_stdin:
            in_file = sys.stdin
            create_cluster_files(in_file, args.library_name, args.database_prefix)
        else:
            with open(args.input) as in_file:
                create_cluster_files(in_file, args.library_name, args.database_prefix)
    except IOError as error:
        sys.stderr.write("Could not open file: {}\n".format(error.filename))
        sys.stderr.flush()
        raise IOError(error)
