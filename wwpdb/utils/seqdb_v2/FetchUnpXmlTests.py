##
# File:    FetchUnpXmlTests.py
# Author:  j. westbrook
# Date:    2-Jan-2013
# Version: 0.001
#
# Update:
#     3-Jan-2013 jdw add variant tests
##
"""
Test cases for individual and batch fetch of UniProt sequence entries.

"""

import sys
import unittest
import os
import os.path
import string
import traceback

from wwpdb.api.facade.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml


class FetchUnpXmlTests(unittest.TestCase):

    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        # Pick up site information from the environment or failover to the development site id.
        self.__siteId = getSiteId(defaultSiteId='WWPDB_DEPLOY_TEST')
        self.__lfh.write("\nTesting with site environment for:  %s\n" % self.__siteId)
        #
        cI = ConfigInfo(self.__siteId)
        self.__testFilePath = './data'
        self.__unpIdList1 = ['P20937', 'P21877', 'P22868', 'P23832', 'P25665', 'P26562', 'P27614']
        self.__unpIdList2 = ['P29490', 'P29496', 'P29498', 'P29499', 'P29503', 'P29506', 'P29508', 'P29509', 'P29525', 'P29533', 'P29534', 'P29547',
                             'P29549', 'P29555', 'P29557', 'P29558', 'P29559', 'P29563', 'P29588', 'P29589', 'P29590', 'P29597', 'P29599', 'P29600',
                             'P29602', 'P29603', 'P29617', 'P29678', 'P29691', 'P29715', 'P29717', 'P29723', 'P29724', 'P29736', 'P29741', 'P29745',
                             'P29748', 'P29749', 'P29752', 'P29758', 'P29768', 'P29803', 'P29808', 'P29813', 'P29827', 'P29830', 'P29837', 'P29838',
                             'P29846', 'P29848', 'P29882', 'P29894', 'P29898', 'P29899', 'P29929', 'P29946', 'P29957', 'P29960', 'P29965', 'P29966',
                             'P29972', 'P29978', 'P29986', 'P29987', 'P29988', 'P29989', 'P29990', 'P29991', 'P29994']

        self.__unpIdListV = ['P42284', 'P42284-1', 'P42284-2', 'P42284-3', 'P29994-1', 'P29994-2', 'P29994-3', 'P29994-4', 'P29994-5', 'P29994-6', 'P29994-7']

    def tearDown(self):
        pass

    def testFetchIds(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            for id in self.__unpIdList1:
                ok = fobj.fetchList([id])
                if ok:
                    fobj.writeUnpXml(id + '.xml')
                    dict = fobj.getResult()
                    for (k, v) in dict.items():
                        self.__lfh.write("%s=%s" % (k, v))
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBatchFetch(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            ok = fobj.fetchList(self.__unpIdList2)
            if ok:
                fobj.writeUnpXml('batch-fetch.xml')
                dict = fobj.getResult()
                for (k, v) in dict.items():
                    self.__lfh.write("%s=%s" % (k, v))
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

    def testFetchVariantIds(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            for id in self.__unpIdListV:
                ok = fobj.fetchList([id])
                if ok:
                    fobj.writeUnpXml(id + '.xml')
                    dict = fobj.getResult()
                    for (eId, eDict) in dict.items():
                        if 'db_isoform' in eDict and eId == id:
                            self.__lfh.write("------ sequence database code  %s has key db_isoform:  %r\n" % (eId, eDict['db_isoform']))
                            self.__lfh.write("------ sequence database code  %s sequence length %d\n" % (eId, len(self.__cleanString(eDict['sequence']))))
                            #self.__lfh.write("------ sequence database code  %s keys %r\n" % (eId,eDict.keys()))
                            self.__lfh.write("%s\n" % self.__cleanString(eDict['sequence']))
                        elif eId == id:
                            self.__lfh.write("------ No matching isoform for %s\n" % id)
                        # for k,v in eDict.items():
                        #    self.__lfh.write("%-25s = %s\n" % (k, v))
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)
        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBatchFetchVariants(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            ok = fobj.fetchList(self.__unpIdListV)
            if ok:
                fobj.writeUnpXml('variant-batch-fetch.xml')
                dict = fobj.getResult()
                for (eId, eDict) in dict.items():
                    self.__lfh.write("\n\n------ Entry id %s\n" % eId)
                    for k, v in eDict.items():
                        self.__lfh.write("%-25s = %s\n" % (k, v))
            else:
                self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)

        except:
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteFetchTests():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(FetchUnpXmlTests("testFetchIds"))
    suiteSelect.addTest(FetchUnpXmlTests("testBatchFetch"))
    #
    return suiteSelect


def suiteFetchVariantTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FetchUnpXmlTests("testFetchVariantIds"))
    # suiteSelect.addTest(FetchUnpXmlTests("testBatchFetchVariants"))
    #
    return suiteSelect


if __name__ == '__main__':
    # Run all tests --
    # unittest.main()
    #
    if False:
        mySuite = suiteFetchTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
        #
    if True:
        mySuite = suiteFetchVariantTests()
        unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
    #
