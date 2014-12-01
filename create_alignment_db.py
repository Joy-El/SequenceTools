__author__ = 'Joy-El R.B. Talbot'
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

from commonIO import read_chunk
from commonIO import CustomParser
from parse_bowtie_output import add_multimapping_tally


