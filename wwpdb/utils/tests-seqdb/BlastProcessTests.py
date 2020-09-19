##
# File:    BlastProcessTests.py
# Author:  j. westbrook
# Date:    2-Jan-2013
# Version: 0.001
#
# Update:

##
"""
Test cases from funning reference sequence database blast searches and processing the search results.

"""

import sys
import unittest
import os
import os.path
import traceback
import platform

from wwpdb.utils.seqdb_v2.BlastProcess import BlastProcess

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, 'test-output', platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)
mockTopPath = os.path.join(TOPDIR, 'wwpdb', 'mock-data')


class BlastProcessTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__testModelPath = os.path.join(mockTopPath, 'MODELS')
        self.__testFileCif = '4ec0.cif'
        self.__testFileFragmentsCif = '3l2j.cif'
        self.__testTaxPath = os.path.join(mockTopPath, 'TAXONOMY')
        self.__taxonomyDataFile = 'nodes.dmp.gz'

    def tearDown(self):
        pass

    def testGetPolymerEntityDetails(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)

                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)

                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails['seq'])
                    self.__lfh.write("length          = %d\n" % len(polyDetails['seq']))
                    self.__lfh.write("type            = %s\n" % polyDetails['type'])
                    self.__lfh.write("fragments       = %r\n" % polyDetails['fragments'])
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetPolymerEntityDetailsFragments(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)
                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails['seq'])
                    self.__lfh.write("length          = %d\n" % len(polyDetails['seq']))
                    self.__lfh.write("type            = %s\n" % polyDetails['type'])
                    self.__lfh.write("fragments       = %r\n" % polyDetails['fragments'])
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testPolymerSearch1(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath,
                                                self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath,
                                  taxonomyFilePath=taxonomyFilePath,
                                  verbose=self.__verbose, log=self.__lfh)
                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    resultlist = bp.Run(entityId=entityId)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    for ii, rL in enumerate(resultlist):
                        for r in rL:
                            for k in sorted(r.keys()):
                                self.__lfh.write("%d: %s:  %s\n" % (ii, k, r[k]))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testPolymerSearchAndStore(self):
        """
        """
        self.__lfh.write("\nStarting %s %s\n" % (self.__class__.__name__, sys._getframe().f_code.co_name))
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                entryId, fExt = os.path.splitext(fn)
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath,
                                                self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath,
                                  taxonomyFilePath=taxonomyFilePath,
                                  verbose=self.__verbose, log=self.__lfh)
                bp.saveBlastResults(blastPath=TESTOUTPUT, blastFileNamePrefix=entryId)
                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    ofn = os.path.join(TESTOUTPUT,
                                       entryId + "_seqdb-match_P" + str(entityId) + ".cif")
                    ok = bp.RunAndSave(entityId=entityId, fName=ofn)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s ofn = %s status = %r -------  \n" % (fn, entityId, ofn, ok))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteSearchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(BlastProcessTests("testGetPolymerEntityDetails"))
    # suiteSelect.addTest(BlastProcessTests("testPolymerSearch1"))e
    suiteSelect.addTest(BlastProcessTests("testPolymerSearchAndStore"))
    #
    return suiteSelect


if __name__ == '__main__':
    # Run all tests --
    # unittest.main()
    #
    mySuite = suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
