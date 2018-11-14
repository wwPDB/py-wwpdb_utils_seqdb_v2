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

from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.utils.seqdb_v2.FetchUniProtEntry import FetchUniProtEntry


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
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            fobj = FetchUniProtEntry(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for id in self.__unpIdListV:
                ok = fobj.fetchList([id])
                if ok:
                    fobj.writeUnpXml(id + '.xml')
                    dict = fobj.getResult()
                    for (eId, eDict) in dict.items():
                        for k, v in eDict.items():
                            self.__lfh.write(" eId = %s   key %r : value: %r\n" % (eId, k, v))
                    for (eId, eDict) in dict.items():
                        if 'db_isoform' in eDict and eId == id:
                            self.__lfh.write("------ sequence database code  %s has key db_isoform:  %r\n" % (eId, eDict['db_isoform']))
                            self.__lfh.write("------ sequence database code  %s sequence length %d\n" % (eId, len(self.__cleanString(eDict['sequence']))))
                            self.__lfh.write("%s\n" % self.__cleanString(eDict['sequence']))
                        elif eId == id:
                            self.__lfh.write("------ No matching isoform for %s\n" % id)
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)
        except:
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
    #
    if True:
        mySuite = suiteFetchVariantTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
