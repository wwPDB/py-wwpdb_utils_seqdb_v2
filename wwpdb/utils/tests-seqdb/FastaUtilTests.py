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

from wwpdb.utils.seqdb_v2.FastaUtil import FastaUtil


class FastaUtilTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        #
        here = os.path.abspath(os.path.dirname(__file__))
        fastaPath = os.path.join(here, "refdata")
        self.__variantFastaFilePath = os.path.join(fastaPath, "uniprot_sprot_varsplic.fasta")
        self.__unpIdListV = ["P42284-1", "P42284-3", "P29994-2", "P29994-3", "P29994-4", "P29994-5", "P29994-6", "P29994-7"]
        self.__unpIdNotV = ["P42284"]

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
                else:
                    self.__lfh.write("id %s not found\n" % uid)
                    self.fail()

            # Test for entries that should not be there
            for uid in self.__unpIdNotV:
                if uid in sD:
                    self.__lfh.write("id %s found in variant file when it should not\n" % uid)
                    self.fail()

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
