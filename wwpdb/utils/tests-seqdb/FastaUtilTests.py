##
# File:    FastaUtilTests.py
# Author:  j. westbrook
# Date:    29-Jul-2014
# Version: 0.001
#
# Update:
##
"""
Test cases for reading UniProt Fasta files -

"""

import sys
import unittest
import os
import os.path
import string
import traceback

from wwpdb.utils.config.ConfigInfo import ConfigInfo, getSiteId
from wwpdb.utils.seqdb_v2.FastaUtil import FastaUtil


@unittest.skip("Cannot test without sequence fasta files")
class FastaUtilTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        # Pick up site information from the environment or failover to the development site id.
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
        self.__lfh.write("\nTesting with site environment for:  %s\n" % self.__siteId)
        #
        cI = ConfigInfo(self.__siteId)
        self.__fastaPath = cI.get("SITE_REFDATA_SEQUENCE_DB_PATH")
        self.__variantFastaFilePath = os.path.join(self.__fastaPath, "uniprot_sprot_varsplic.fasta")
        self.__unpIdListV = ["P42284", "P42284-1", "P42284-3", "P29994-1", "P29994-2", "P29994-3", "P29994-4", "P29994-5", "P29994-6", "P29994-7"]

    def tearDown(self):
        pass

    def __cleanString(self, strIn):
        sL = []
        for ss in strIn:
            if ss in string.whitespace:
                continue
            sL.append(ss)
        return "".join(sL)

    def testReadFasta(self):
        """"""
        self.__lfh.write("\nStarting FastaUtilTests testReadFasts\n")
        try:
            fau = FastaUtil(verbose=self.__verbose, log=self.__lfh)
            _sL, sD = fau.loadFastaUniProt(fastaFilePath=self.__variantFastaFilePath)
            for uid in self.__unpIdListV:
                if uid in sD:
                    self.__lfh.write(
                        "id %s accession %s isoform %s description %s gene %s org %s length sequence %d\n"
                        % (
                            uid,
                            sD[uid]["db_accession"],
                            sD[uid]["db_isoform"],
                            sD[uid]["description"],
                            sD[uid]["gene_name"],
                            sD[uid]["org"],
                            len(self.__cleanString(sD[uid]["sequence"])),
                        )
                    )
                    self.__lfh.write("%s\n" % self.__cleanString(sD[uid]["sequence"]))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()


def suiteReadTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FastaUtilTests("testReadFasta"))
    #
    return suiteSelect


if __name__ == "__main__":

    mySuite = suiteReadTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
