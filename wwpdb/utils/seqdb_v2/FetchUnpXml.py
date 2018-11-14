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
##

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import sys
import getopt
try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError
import ssl
import traceback
from xml.dom import minidom
try:
    # Python 3
    from itertools import zip_longest
except ImportError:
    # Python 2
    from itertools import izip_longest as zip_longest

from wwpdb.utils.seqdb_v2.ReadUnpXml import ReadUnpXmlString


class FetchUnpXml:
    """
    Manage execution of fetch queries for UniProt entries.

    XML entry data is parsed into a dictionary structure.

    """
    def __init__(self, maxLength=100, verbose=False, log=sys.stderr):

        self.__verbose = verbose
        self.__debug = False
        self.__lfh = log
        self._baseUrl = 'https://www.ebi.ac.uk/Tools/dbfetch/dbfetch'
        self._baseUrlUnp = 'https://www.ebi.ac.uk/proteins/api/proteins'
        #

        self.__maxLength = maxLength
        #
        self.__dataList = []
        self.__result = {}
        #
        self.__idList = []
        self.__searchIdList = []
        self.__variantD = {}
        #

    def fetchList(self, idList):
        """     Execute a fetch query for the input id list.
                The input list is filtered for variants (e.g. ids with appended '-#').

                Divide the input list into manageable chunks, fetch each chunk,
                and concatenate the result.

                Return True for success or False otherwise.

        """
        try:
            self.__result = {}
            self.__dataList = []

            self.__idList = idList
            #
            num = self.__processIdList()
            if self.__verbose:
                self.__lfh.write("+FetchUnpXml.fetchList() input id list len %d\n" % len(idList))
            if self.__debug:
                self.__lfh.write("+FetchUnpXml.fetchList() input id list %s\n" % idList)
                self.__lfh.write("+FetchUnpXml.fetchList() search   list %s\n" % self.__searchIdList)
                self.__lfh.write("+FetchUnpXml.fetchList() variants      %s\n" % self.__variantD.items())

            if num == 0:
                return False
            #
            self.__subLists = self.__makeSubLists(self.__maxLength, self.__searchIdList)

            for subList in self.__subLists:
                idString = ','.join(subList)
                if (self.__debug):
                    self.__lfh.write("+FetchUnpXml.fetchList() subList %s string %s\n" % (subList, idString))
                xmlText = self.__RequestUnpXml(idString)
                # filter possible simple text error messages from the failed queries.
                if ((xmlText is not None) and not xmlText.startswith("ERROR")):
                    self.__dataList.append(xmlText)
            if len(self.__dataList) > 0:
                ok = self.__ParseUnpXmlData()
            else:
                ok = False

            return ok
        except:
            if (self.__verbose):
                traceback.print_exc(file=self.__lfh)
        return False

    def writeUnpXml(self, filename):
        file = open(filename, 'w')
        for data in self.__dataList:
            file.write(data)
        file.close()

    def getResult(self):
        """ Get the parsed UniProt entry data in a dictionary structure.

            Return dictionary has a key of input accession code which references
            a dictionary containing the entry data.

        """
        return self.__result
    #

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
        for id in self.__idList:
            # check for variant id
            idx = id.find('-')
            if idx == -1:
                sId = id
            else:
                sId = id[0:idx]
                self.__variantD[id] = sId
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

    def __RequestUnpXml(self, idString):
        """Execute fetch Request for the input comma separated accession list.

           Return xml text for the corresentry  UniProt entries
           """
        #
        gcontext = ssl._create_unverified_context()

        if False:
            params = {}
            params['db'] = 'uniprotkb'
            params['id'] = idString

            params['format'] = 'uniprotxml'
            params['style'] = 'raw'
            #params['Retrieve'] = 'Retrieve'
            requestData = urlencode(params)
            reqH = urlopen(self._baseUrl, requestData, context=gcontext)
        else:
            params = {}
            params['size'] = '-1'
            params['accession'] = idString

            requestData = urlencode(params)
            # Need to do this as UNP service will not take POST - so forect GET
            request = Request('%s?%s' % (self._baseUrlUnp, requestData))
            request.add_header("Accept", "application/xml")
            reqH = urlopen(request, context=gcontext)

        data = reqH.read()
        reqH.close()
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
            for data in self.__dataList:
                readxml = ReadUnpXmlString(data, verbose=self.__verbose, log=self.__lfh)
                for vId, aId in self.__variantD.items():
                    readxml.addVariant(aId, vId)
                self.__result.update(readxml.getResult())
            return True
        except:
            if (self.__verbose):
                traceback.print_exc(file=self.__lfh)
            return False


def main(argv):
    opts, args = getopt.getopt(argv, "i:", ["id="])
    for opt, arg in opts:
        if opt in ("-i", "--id"):
            id = arg
            fobj = FetchUnpXml()
            fobj.fetchList([id])
            fobj.writeUnpXml(id + '.xml')
            dict = fobj.getResult()
            for (k, v) in dict.items():
                sys.stdout.write("%-30s = %s\n" % (k, v))
if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        sys.stderr.write(exc)
        sys.exit(1)
