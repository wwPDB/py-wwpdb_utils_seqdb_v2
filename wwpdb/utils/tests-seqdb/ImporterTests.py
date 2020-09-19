##
# File: ImportTests.py
# Date:  06-Oct-2018  E. Peisach
#
# Updates:
##
"""Test cases for seqdb_v2"""

__docformat__ = "restructuredtext en"
__author__ = "Ezra Peisach"
__email__ = "peisach@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.01"

import unittest

import wwpdb.utils.seqdb_v2.BlastProcess  # noqa: F401
import wwpdb.utils.seqdb_v2.FastaUtil  # noqa: F401
import wwpdb.utils.seqdb_v2.FetchNcbiXml  # noqa: F401
import wwpdb.utils.seqdb_v2.FetchUniProtEntry  # noqa: F401
import wwpdb.utils.seqdb_v2.FetchUnpXml  # noqa: F401
import wwpdb.utils.seqdb_v2.NcbiBlastService  # noqa: F401
import wwpdb.utils.seqdb_v2.ReadNcbiBlastXml  # noqa: F401
import wwpdb.utils.seqdb_v2.ReadNcbiSummary  # noqa: F401
import wwpdb.utils.seqdb_v2.ReadNcbiXml  # noqa: F401
import wwpdb.utils.seqdb_v2.ReadUnpBlastXml  # noqa: F401
import wwpdb.utils.seqdb_v2.ReadUnpXml  # noqa: F401
import wwpdb.utils.seqdb_v2.UnpBlastService  # noqa: F401
import wwpdb.utils.seqdb_v2.mmCIFUtil  # noqa: F401


class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        """Tests simple instantiation"""
        pass
