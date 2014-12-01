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
