__author__ = 'joy-el'
"""commonIO.py contains extensions of IO functions."""

import sys
import argparse


def read_chunk(open_file_object, chunk_size=1048):
    """Read in file by chunk_size chunks returning one line at a time."""
    # get first chunk
    chunk = open_file_object.read(chunk_size)
    # continue looping until a chunk is just EOF (empty line)
    while chunk:
        chunk_list = chunk.split("\n")
        # yield all but last, potentially incomplete line
        for c in chunk_list[:-1]:
            yield c
        # add incomplete line to beginning of next chunk read
        chunk = chunk_list[-1] + open_file_object.read(chunk_size)


def read_fastq_chunk(fastq_file, chunk_size=1048):
    """Return a tuple representing a fastq read

    DOES NOT WORK WITH MULTI-LINE SEQUENCE FASTQ FILES!!"""
    ##TODO handle comment lines
    ##TODO handle multi-line DNA sequence
    read = []
    for line in read_chunk(fastq_file, chunk_size):
        read.append(line)
        if len(read) > 4:
            yield read[0:4]
            read = read[4:]
    yield read[0:4]


class CustomParser(argparse.ArgumentParser):
    """Custom command line argument parser inheriting from argparse.

    Allows the help menu to print upon error."""
    def __init__(self, description):
        argparse.ArgumentParser.__init__(self, description)

    def error(self, message):
        """Print help message when argparse error occurs.

        Code borrowed from unutbu's answer at:
        http://stackoverflow.com/questions/4042452/display-help-message-with-python-argparse-when-script-is-called-without-any-argu"""
        sys.stderr.write('error: {}\n'.format(message))
        self.print_help()
        sys.exit(2)
