"""
File:    FetchNcbiXml.py
Author:  Zukang Feng
Update:  24-May-2010
Version: 001  Initial version

Updates:

25-Feb-2013  jdw create a class version to return full XML entries -
20-Mar-2013  jdw Fix error in import path

"""

__author__  = "Zukang Feng"
__email__   = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

from xml.dom import minidom
import urllib, urllib2, ssl
from wwpdb.utils.seqdb_v2.ReadNcbiXml     import ReadNcbiXmlString
from wwpdb.utils.seqdb_v2.ReadNcbiSummary import ReadNcbiSummaryString
import sys
import getopt



class FetchNcbiXml:
    """Using esummary.fcgi utility to get entry summary xml file from NCBI site and
       using ReadNcbiSummary class to parse the result
    """
    def __init__(self, id, database):
        self._id = id;
        self._database = database
#       self._baseUrl = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        self._baseUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        self._data = self._RequestNcbiXml()

    def _RequestNcbiXml(self):
        """Request summary xml from NCBI site"""
        params = {}
        params['db'] = self._database
        params['id'] = self._id;
        params['retmode'] = 'xml'
        gcontext = ssl._create_unverified_context()
        requestData = urllib.urlencode(params)
        reqH = urllib2.urlopen(self._baseUrl, requestData, context=gcontext)
        data = reqH.read()
        reqH.close()
        return data

    def WriteNcbiXml(self, filename):
        try:
            file = open(filename, 'w')
            file.write(self._data)
            file.close()
            return True
        except:
            return False

    def ParseNcbiXmlData(self):
        """Using ReadNcbiSummary class to parse the result"""
        readxml = ReadNcbiSummaryString(self._data)
        return readxml.GetResult()


class FetchFullNcbiXml:
    """Using efetch.fcgi utility to get and parse the full entry xml file from NCBI site.
    """
    def __init__(self, id, database):
        self._id = id;
        self._database = database
        self._baseUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
        self._data = self._RequestNcbiXml()

    def _RequestNcbiXml(self):
        """Request summary xml from NCBI site"""
        params = {}
        params['db'] = self._database
        params['id'] = self._id;
        params['retmode'] = 'xml'
        requestData = urllib.urlencode(params)
        gcontext = ssl._create_unverified_context()
        reqH = urllib2.urlopen(self._baseUrl, requestData, context=gcontext)
        data = reqH.read()
        reqH.close()
        return data

    def WriteNcbiXml(self, filename):
        file = open(filename, 'w')
        file.write(self._data)
        file.close()

    def ParseNcbiXmlData(self):
        """Parse xml result
        """
        readxml = ReadNcbiXmlString(self._data)
        return readxml.GetResult()

        
def main(argv):
    opts, args = getopt.getopt(argv, "i:d:", ["id=", "db="])

    id = None
    db = None
    for opt, arg in opts:
        if opt in ("-i", "--id"):
            id = arg
        elif opt in ("-d", "--db"):
            db = arg

    if id and db:
        fetchobj = FetchNcbiXml(id, db)
        fetchobj.WriteNcbiXml(id + '.xml')
        dict = fetchobj.ParseNcbiXmlData()
        for (k, v) in dict.items():
            print "%s=%s" % (k, v)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception, exc:
        print exc
        sys.exit(1)
