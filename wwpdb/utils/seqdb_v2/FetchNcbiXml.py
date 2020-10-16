"""
File:    FetchNcbiXml.py
Author:  Zukang Feng
Update:  24-May-2010
Version: 001  Initial version

Updates:

25-Feb-2013  jdw create a class version to return full XML entries -
20-Mar-2013  jdw Fix error in import path

"""

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import sys
import requests
import getopt
import urllib3

from wwpdb.utils.seqdb_v2.ReadNcbiXml import ReadNcbiXmlString
from wwpdb.utils.seqdb_v2.ReadNcbiSummary import ReadNcbiSummaryString

# To disable warning about not checking ssl certificates. Still needed?
urllib3.disable_warnings()


class FetchNcbiXml:
    """Using esummary.fcgi utility to get entry summary xml file from NCBI site and
    using ReadNcbiSummary class to parse the result
    """

    def __init__(self, qid, database, apikey=None):
        self._id = qid
        self._database = database
        self._apikey = apikey
        self._baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        self._data = self._RequestNcbiXml()

    def _RequestNcbiXml(self):
        """Request summary xml from NCBI site"""
        params = {}
        params["db"] = self._database
        params["id"] = self._id
        params["retmode"] = "xml"
        if self._apikey:
            params["api_key"] = self._apikey
        #
        try:
            reqH = requests.get(self._baseUrl, params=params, verify=False)
            reqH.raise_for_status()
            data = reqH.text
            return data
        except:  # noqa: E722 pylint: disable=bare-except
            return ""
        #

    def WriteNcbiXml(self, filename):
        try:
            file = open(filename, "w")
            file.write(self._data)
            file.close()
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False

    def ParseNcbiXmlData(self):
        """Using ReadNcbiSummary class to parse the result"""
        readxml = ReadNcbiSummaryString(self._data)
        return readxml.GetResult()


class FetchFullNcbiXml:
    """Using efetch.fcgi utility to get and parse the full entry xml file from NCBI site."""

    def __init__(self, qid, database, apikey=None):
        self._id = qid
        self._database = database
        self._apikey = apikey
        self._baseUrl = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self._data = self._RequestNcbiXml()

    def _RequestNcbiXml(self):
        """Request summary xml from NCBI site"""
        params = {}
        params["db"] = self._database
        params["id"] = self._id
        params["retmode"] = "xml"
        if self._apikey:
            params["api_key"] = self._apikey
        #
        try:
            reqH = requests.get(self._baseUrl, params=params, verify=False)
            reqH.raise_for_status()
            data = reqH.text
            return data
        except:  # noqa: E722 pylint: disable=bare-except
            return ""
        #

    def WriteNcbiXml(self, filename):
        file = open(filename, "w")
        file.write(self._data)
        file.close()

    def ParseNcbiXmlData(self):
        """Parse xml result"""
        readxml = ReadNcbiXmlString(self._data)
        return readxml.GetResult()


def main(argv):  # pragma: no cover
    opts, _args = getopt.getopt(argv, "i:d:a:", ["id=", "db=", "apikey="])

    tid = None
    db = None
    apikey = None
    for opt, arg in opts:
        if opt in ("-i", "--id"):
            tid = arg
        elif opt in ("-d", "--db"):
            db = arg
        elif opt in ("-a", "--apikey"):
            apikey = arg

    if tid and db:
        fetchobj = FetchNcbiXml(tid, db, apikey)
        # fetchobj = FetchFullNcbiXml(tid, db, apikey)
        fetchobj.WriteNcbiXml(tid + ".xml")
        pdict = fetchobj.ParseNcbiXmlData()
        for (k, v) in pdict.items():
            print("%s=%s" % (k, v))


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)
