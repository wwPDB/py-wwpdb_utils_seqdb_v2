##
# File:    FetchUniprotEntry.py
# Author:  John Westbrook
# Update:  29-July-2014
# Version: 001  Initial version
#
# Update:
#  17-Sep-2014  jdw add attribute to reflect isoform sequence replaced from FASTA resource file
#
##


import sys
import os
import copy
import traceback
from wwpdb.utils.config.ConfigInfo import ConfigInfo
from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml
from wwpdb.utils.seqdb_v2.FastaUtil import FastaUtil


class FetchUniProtEntry:
    """
    Manage execution of fetch queries for UniProt entries.

    This class wraps fetch data for UniProt entries including capturing
    correct variant sequences from UniProt distributed Fasta file.
    This avoids any inconsitency with Blast search results --

    """

    def __init__(self, siteId=None, maxLength=100, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__siteId = siteId
        self.__lfh = log
        self.__debug = False
        self.__maxLength = maxLength
        #
        self.__fXml = FetchUnpXml(maxLength=self.__maxLength, verbose=self.__verbose, log=self.__lfh)
        self.__fFa = FastaUtil(verbose=self.__verbose, log=self.__lfh)
        self.__vList = []
        self.__vD = {}
        self.__faLoaded = False
        #
        cI = ConfigInfo(self.__siteId)
        self.__fastaPath = cI.get("SITE_REFDATA_SEQUENCE_DB_PATH")
        self.__variantFastaFilePath = os.path.join(self.__fastaPath, "uniprot_sprot_varsplic.fasta")
        #

        self.__uidList = []

    def __loadVariants(self):
        try:
            self.__vList, self.__vD = self.__fFa.loadFastaUniProt(fastaFilePath=self.__variantFastaFilePath)
            self.__faLoaded = True
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+FetchUniprotEntry.__loadVariants() failed for file %s\n" % self.__variantFastaFilePath)
                traceback.print_exc(file=self.__lfh)

    def fetchList(self, idList):
        self.__uidList = copy.deepcopy(idList)
        return self.__fXml.fetchList(idList)

    def writeUnpXml(self, filename):
        return self.__fXml.writeUnpXml(filename)

    def getResult(self):
        """Get the parsed UniProt entry data in a dictionary structure.

        Return dictionary has a key of input accession code which references
        a dictionary containing the entry data.

        """
        try:
            rD = self.__fXml.getResult()
            for uid, uD in rD.items():
                if "db_isoform" in uD:
                    if self.__verbose:
                        self.__lfh.write("+FetchUniProtEntry.getResult() searching isoform %s\n" % uD["db_isoform"])
                    if not self.__faLoaded:
                        self.__loadVariants()
                    if uid in self.__vD:
                        uD["sequence"] = self.__vD[uid]["sequence"]
                        uD["db_isoform_description"] = self.__vD[uid]["description"]
                        uD["isoform_sequence_updated_from_fasta"] = "Y"
                        if self.__verbose:
                            self.__lfh.write("+FetchUniProtEntry.getResult() Updating isoform sequences from FASTA data for %s\n" % uid)
                            self.__lfh.write("+FetchUniProtEntry.getResult() Updating isoform description from FASTA  %s\n" % uD["db_isoform_description"])
            return rD
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)

    #


def main():  # pragma: no cover
    _funp = FetchUniProtEntry()  # noqa: F841
