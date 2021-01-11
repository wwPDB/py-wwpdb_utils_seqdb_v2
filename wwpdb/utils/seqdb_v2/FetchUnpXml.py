##
# File:    FetchUnpXml.py
# Author:  Zukang Feng
# Update:  24-May-2010
# Version: 001  Initial version
#
# Update:
#
#  02-Jan-2013 jdw revise to support lists of accession codes as input.
#  03-Jan-2013 jdw revise variant handling.
#  13-Mar-2013 jdw Fix initialization issue
#
#  29-Jul-2014 jdw WARNING -- Handling of variant seqeuences is incomplete --
#  02-Sep-2020 zf  added a work-around for isoforms sequence.
#  29-Dec-2020 zf  added getMultipleResultDict()
##

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import sys
import getopt
import requests
import urllib3
import traceback

try:
    # Python 3
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest

from wwpdb.utils.seqdb_v2.ReadUnpXml import ReadUnpXmlString
from wwpdb.utils.config.ConfigInfo import ConfigInfo
import logging

logger = logging.getLogger()

# To disable warning about not checking ssl certificates. Still needed?
urllib3.disable_warnings()


class FetchUnpXml:
    """
    Manage execution of fetch queries for UniProt entries.

    XML entry data is parsed into a dictionary structure.

    """

    def __init__(self, maxLength=100, verbose=False, log=sys.stderr):

        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        self._baseUrl = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch"
        self._baseUrlUnp = "https://www.ebi.ac.uk/proteins/api/proteins"
        #
        self.__maxLength = maxLength
        #
        self.__dataList = []
        self.__result = {}
        self.__multipleResult = {}
        #
        self.__idList = []
        self.__searchIdList = []
        self.__variantD = {}
        #
        self.__isoformIdList = []
        #
        cI = ConfigInfo()
        self.__forcefallback = cI.get("FETCH_UNP_FORCE_FALLBACK", None)
        #

    def fetchList(self, idList):
        """Execute a fetch query for the input id list.
        The input list is filtered for variants (e.g. ids with appended '-#').

        Divide the input list into manageable chunks, fetch each chunk,
        and concatenate the result.

        Return True for success or False otherwise.

        """
        try:
            if not idList:
                return False
            #
            self.__result = {}
            self.__multipleResult = {}
            self.__dataList = []
            self.__idList = []
            self.__isoformIdList = []

            for accId in idList:
                splitList = accId.split("-")
                if len(splitList) == 2:
                    if splitList[0] not in self.__isoformIdList:
                        self.__isoformIdList.append(splitList[0])
                    #
                else:
                    if accId not in self.__idList:
                        self.__idList.append(accId)
                    #
                #
            #

            num = self.__processIdList()
            if self.__verbose:
                self.__lfh.write("+FetchUnpXml.fetchList() input id list len %d\n" % len(idList))
            if self.__debug:
                self.__lfh.write("+FetchUnpXml.fetchList() input id list %s\n" % idList)
                self.__lfh.write("+FetchUnpXml.fetchList() search   list %s\n" % self.__searchIdList)
                self.__lfh.write("+FetchUnpXml.fetchList() variants      %s\n" % self.__variantD.items())
            #
            if (num == 0) and (not self.__isoformIdList):
                return False
            #
            if num:
                subLists = self.__makeSubLists(self.__maxLength, self.__searchIdList)

                for subList in subLists:
                    idString = ",".join(subList)
                    if self.__debug:
                        self.__lfh.write("+FetchUnpXml.fetchList() subList %s string %s\n" % (subList, idString))
                    #
                    # If primary site has issues, automatic fallback
                    if self.__forcefallback:
                        xmlText = self.__RequestUnpXml(idString, fallback=True)
                    else:
                        try:
                            xmlText = self.__RequestUnpXml(idString)
                        except requests.exceptions.HTTPError:
                            xmlText = self.__RequestUnpXml(idString, fallback=True)
                        #
                    #
                    # filter possible simple text error messages from the failed queries.
                    if (xmlText is not None) and not xmlText.startswith("ERROR"):
                        self.__dataList.append(xmlText)
                    #
                #
            #
            for accId in self.__isoformIdList:
                xmlText = self.__RequestIsoformsUnpXml(accId)
                if xmlText:
                    self.__dataList.append(xmlText)
                #
            #
            if len(self.__dataList) > 0:
                ok = self.__ParseUnpXmlData()
            else:
                ok = False
            #
            return ok
        except Exception as e:
            if self.__verbose:
                self.__lfh.write("+FetchUnpXml.fetchList() exception %s\n" % str(e))
                traceback.print_exc(file=self.__lfh)
            #
        #
        return False

    def writeUnpXml(self, filename):
        wfile = open(filename, "w")
        for data in self.__dataList:
            wfile.write(data)
        #
        wfile.close()

    def getResult(self):
        """Get the parsed UniProt entry data in a dictionary structure.

        Return dictionary has a key of input accession code which references
        a dictionary containing the entry data.

        """
        return self.__result

    def getMultipleResultDict(self):
        return self.__multipleResult

    def __processIdList(self):
        """
        Filter the input id list for variants and create the list of unique
        searchable id codes.   Create a dictionary of variant identifiers and
        their corresponding searchable id codes.

        """
        self.__searchIdList = []
        self.__variantD = {}
        tList = []
        #
        for vid in self.__idList:
            # check for variant id
            idx = vid.find("-")
            if idx == -1:
                sId = vid
            else:
                sId = vid[0:idx]
                self.__variantD[vid] = sId
            #
            tList.append(sId)
        #
        # unique list of searchable accessions
        #
        self.__searchIdList = list(set(tList))
        #
        return len(self.__searchIdList)

    def __makeSubLists(self, n, iterable):
        args = [iter(iterable)] * n
        return ([e for e in t if e is not None] for t in zip_longest(*args))

    def __makeSubListsWithPadding(self, n, iterable, padvalue=None):
        "__sublist(3, 'abcdefg', 'x') --> ('a','b','c'), ('d','e','f'), ('g','x','x')"
        return zip_longest(*[iter(iterable)] * n, fillvalue=padvalue)

    def __RequestUnpXml(self, idString, fallback=False):
        """Execute fetch Request for the input comma separated accession list.

        Return xml text for the corresentry  UniProt entries
        """
        if not fallback:
            params = {}
            params["db"] = "uniprotkb"
            params["id"] = idString

            params["format"] = "uniprotxml"
            params["style"] = "raw"
            # params['Retrieve'] = 'Retrieve'
            logger.debug("Request %s with data %s", self._baseUrl, params)
            reqH = requests.post(self._baseUrl, data=params, verify=False)
            reqH.raise_for_status()
        else:
            params = {}
            params["size"] = "-1"
            params["accession"] = idString

            # Need to do this as UNP service will not take POST - so force GET
            logger.debug("Request with data %s", params)
            reqH = requests.get(self._baseUrlUnp, params=params, headers={"Accept": "application/xml"}, verify=False)
            reqH.raise_for_status()
        #
        # data = reqH.read()
        data = reqH.text
        return data

    def __RequestIsoformsUnpXml(self, accId):
        """Using URL https://www.ebi.ac.uk/proteins/api/proteins/{accession}/isoforms"""
        isoformUrl = self._baseUrlUnp + "/" + accId + "/isoforms"
        params = {}
        reqH = requests.get(isoformUrl, params=params, headers={"Accept": "application/xml"}, verify=False)
        reqH.raise_for_status()
        data = reqH.text
        if data.find("errorMessages") >= 0:
            return ""
        #
        return data

    def __ParseUnpXmlData(self):
        """
        Parse the accumulated xml text for each chunk and store the parsed data in
        the internal result dictionary.

        variant id codes are regisistered with the parser so that the returned dictionary
        contains keys corresponding to the original input id list.

        """
        try:
            self.__result = {}
            self.__multipleResult = {}
            for data in self.__dataList:
                readxml = ReadUnpXmlString(data, verbose=self.__verbose, log=self.__lfh)
                for vId, aId in self.__variantD.items():
                    readxml.addVariant(aId, vId)
                #
                self.__result.update(readxml.getResult())
                self.__multipleResult.update(readxml.getMultipleResultDict())
            #
            return True
        except Exception as e:
            if self.__verbose:
                self.__lfh.write("+FetchUnpXml.__ParseUnpXmlData() exception %s\n" % str(e))
                traceback.print_exc(file=self.__lfh)
            #
            return False
        #


def main(argv):  # pragma: no cover
    opts, _args = getopt.getopt(argv, "i:", ["id="])
    for opt, arg in opts:
        if opt in ("-i", "--id"):
            uid = arg
            fobj = FetchUnpXml(verbose=True)
            fobj.fetchList([uid])
            fobj.writeUnpXml(uid + ".xml")
            rdict = fobj.getResult()
            for (k, v) in rdict.items():
                sys.stdout.write("%-30s = %s\n" % (k, v))
            #
        #
    #


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        sys.stderr.write(exc)
        sys.exit(1)
    #
