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


import wwpdb.utils.seqdb_v2.BlastProcess
import wwpdb.utils.seqdb_v2.FastaUtil
import wwpdb.utils.seqdb_v2.FetchNcbiXml
import wwpdb.utils.seqdb_v2.FetchUniProtEntry
import wwpdb.utils.seqdb_v2.FetchUnpXml
import wwpdb.utils.seqdb_v2.NcbiBlastService
import wwpdb.utils.seqdb_v2.ReadNcbiBlastXml
import wwpdb.utils.seqdb_v2.ReadNcbiSummary
import wwpdb.utils.seqdb_v2.ReadNcbiXml
import wwpdb.utils.seqdb_v2.ReadUnpBlastXml
import wwpdb.utils.seqdb_v2.ReadUnpXml
import wwpdb.utils.seqdb_v2.UnpBlastService
import wwpdb.utils.seqdb_v2.mmCIFUtil

class ImportTests(unittest.TestCase):
    def setUp(self):
        pass

    def testInstantiate(self):
        """Tests simple instantiation"""
        pass
