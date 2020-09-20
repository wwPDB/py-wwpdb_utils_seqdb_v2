##
# File:    FetchUniProtEntryTests.py
# Author:  j. westbrook
# Date:    29-July-2014
# Version: 0.001
#
# Update:
##
"""
Test cases for wrapper for individual and batch fetch of UniProt sequence entries (including variants)

This version extracts variant sequence from UniProt Variant FASTA file to avoid any inconsitency with
Blast search results --

"""

import sys
import unittest
import os
import os.path
import string
import traceback
import platform

from wwpdb.utils.config.ConfigInfo import getSiteId
from wwpdb.utils.seqdb_v2.FetchUniProtEntry import FetchUniProtEntry

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, 'test-output', platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)


@unittest.skip("Until have downloaded uniprot variant file")
class FetchUniProtEntryTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__siteId = getSiteId(defaultSiteId='WWPDB_DEPLOY_TEST')
        self.__lfh.write("\nTesting with site environment for:  %s\n" % self.__siteId)
        # self.__unpIdListV=['P42284','P42284-1','P42284-2','P42284-3','P29994-1','P29994-2','P29994-3','P29994-4','P29994-5','P29994-6','P29994-7','Q9NYG8-2']
        # self.__unpIdListV = ['P01116-2', 'P01116-2', 'Q9H8M2-1']
        self.__unpIdListV = ['Q9D6K9', 'Q9D6K9-1', 'Q9D6K9-2']
        #

    def tearDown(self):
        pass

    def testFetchVariantIds(self):
        """
        """
        self.__lfh.write("\nStarting FetchUniProtEntryTests testFetchVariantIds\n")
        try:
            fobj = FetchUniProtEntry(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for unpid in self.__unpIdListV:
                ok = fobj.fetchList([unpid])
                if ok:
                    fobj.writeUnpXml(os.path.join(TESTOUTPUT, unpid + '.xml'))
                    rdict = fobj.getResult()
                    for (eId, eDict) in rdict.items():
                        for k, v in eDict.items():
                            self.__lfh.write(" eId = %s   key %r : value: %r\n" % (eId, k, v))
                    for (eId, eDict) in rdict.items():
                        if 'db_isoform' in eDict and eId == unpid:
                            self.__lfh.write("------ sequence database code  %s has key db_isoform:  %r\n" % (eId, eDict['db_isoform']))
                            self.__lfh.write("------ sequence database code  %s sequence length %d\n" % (eId, len(self.__cleanString(eDict['sequence']))))
                            self.__lfh.write("%s\n" % self.__cleanString(eDict['sequence']))
                        elif eId == unpid:
                            self.__lfh.write("------ No matching isoform for %s\n" % unpid)
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % unpid)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __cleanString(self, strIn):
        sL = []
        for ss in strIn:
            if ss in string.whitespace:
                continue
            sL.append(ss)
        return ''.join(sL)


def suiteFetchVariantTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FetchUniProtEntryTests("testFetchVariantIds"))
    #
    return suiteSelect


if __name__ == '__main__':
    mySuite = suiteFetchVariantTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
