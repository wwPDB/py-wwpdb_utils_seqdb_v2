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
import platform
from xml.dom import minidom

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch
import requests
import requests_mock

try:
    from urllib.parse import parse_qs
except ImportError:
    from urlparse import parse_qs

from wwpdb.utils.config.ConfigInfo import getSiteId
from wwpdb.utils.config.ConfigInfoApp import ConfigInfoAppCommon
from wwpdb.utils.seqdb_v2.FetchUniProtEntry import FetchUniProtEntry

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, "test-output", platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)


class MyConfigInfo(ConfigInfoAppCommon):
    """A class to bypass setting of refdata"""

    def __init__(self, siteId=None, verbose=True, log=sys.stderr):
        super(MyConfigInfo, self).__init__(siteId=siteId, verbose=verbose, log=log)

    def get_site_refdata_sequence_db_path(self):
        return os.path.join(HERE, "refdata")


# Mock request handlers
# dbfetch parses the request and returns data
def dbfetchTextCallBack(request, context):
    """Handles dbfetch requests"""
    dat = parse_qs(request.body)
    accs = dat["id"][0]

    # Create a combined file of contents
    node = None
    root = None

    for acc in accs.split(","):
        fpath = os.path.join(HERE, "refdata", "dbfetch", acc + ".xml")
        if not os.path.exists(fpath):
            context.status_code = 404
            print("XXX COULD NOT FIND", fpath)
            return ""

        doc = minidom.parse(fpath)
        if not node:
            node = doc
            root = node.documentElement
        else:
            ent = doc.getElementsByTagName("entry")
            for e in ent:
                root.appendChild(e)

    if node:
        ret = root.toxml()
    else:
        ret = ""
    return ret


def isoformMatcher(request):
    """Isoform search matcher of the form
    https://www.ebi.ac.uk/proteins/api/proteins/P29994/isoforms
    """
    if "https://www.ebi.ac.uk/proteins/api" in request.url:
        resp = requests.Response()

        # Parse the request to get the id
        preq = request.path_url.split("/")
        if len(preq) < 6:
            return None
        iso = preq[5]
        acc = preq[4]
        if iso != "isoforms":
            return None

        fpath = os.path.join(HERE, "refdata", "ebi_proteins_isoforms", acc + ".xml")
        if not os.path.exists(fpath):
            resp.status_code = 404
            return resp

        resp.status_code = 200
        with open(fpath, "rb") as fin:
            resp._content = fin.read()  # pylint: disable=protected-access
        return resp
    # Error - not found
    return None


def isoformMatcher_dbfetch(request):
    """Isoform search matcher of the form
    https://www.ebi.ac.uk/proteins/api/proteins/P29994/isoforms
    https://www.ebi.ac.uk/Tools/dbfetch/dbfetch
    """
    if "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch" in request.url:
        print("CKP1", request.body)
        resp = requests.Response()

        # Parse the request to get the id
        dat = parse_qs(request.body)
        acc = dat["id"][0]

        fpath = os.path.join(HERE, "refdata", "dbfetch", acc + ".xml")
        if not os.path.exists(fpath):
            resp.status_code = 404
            return resp

        resp.status_code = 200
        with open(fpath, "rb") as fin:
            resp._content = fin.read()  # pylint: disable=protected-access
        return resp
    # Error - not found
    return None


class FetchUniProtEntryTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__siteId = getSiteId(defaultSiteId="WWPDB_DEPLOY_TEST")
        self.__lfh.write("\nTesting with site environment for:  %s\n" % self.__siteId)
        # self.__unpIdListV=['P42284','P42284-1','P42284-2','P42284-3','P29994-1','P29994-2','P29994-3','P29994-4','P29994-5','P29994-6','P29994-7','Q9NYG8-2']
        # self.__unpIdListV = ['P01116-2', 'P01116-2', 'Q9H8M2-1']
        self.__unpIdListV = ["Q9D6K9", "Q9D6K9-1", "Q9D6K9-2"]
        if "MOCKREQUESTS" in os.environ:
            self.__mock = requests_mock.Mocker()
            self.__mock.post("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch", text=dbfetchTextCallBack)
            self.__mock.add_matcher(isoformMatcher_dbfetch)
            self.__mock.start()
        else:
            self.__mock = None

    def tearDown(self):
        if self.__mock is not None:
            self.__mock.stop()

    @patch("wwpdb.utils.seqdb_v2.FetchUniProtEntry.ConfigInfoAppCommon", side_effect=MyConfigInfo)
    def testFetchVariantIds(self, mock1):  # pylint: disable=unused-argument
        """"""
        self.__lfh.write("\nStarting FetchUniProtEntryTests testFetchVariantIds\n")
        fobj = FetchUniProtEntry(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
        try:
            fobj = FetchUniProtEntry(siteId=self.__siteId, verbose=self.__verbose, log=self.__lfh)
            for unpid in self.__unpIdListV:
                ok = fobj.fetchList([unpid])
                if ok:
                    fobj.writeUnpXml(os.path.join(TESTOUTPUT, unpid + ".xml"))
                    rdict = fobj.getResult()
                    for (eId, eDict) in rdict.items():
                        for k, v in eDict.items():
                            self.__lfh.write(" eId = %s   key %r : value: %r\n" % (eId, k, v))
                    for (eId, eDict) in rdict.items():
                        if "db_isoform" in eDict and eId == unpid:
                            self.__lfh.write("------ sequence database code  %s has key db_isoform:  %r\n" % (eId, eDict["db_isoform"]))
                            self.__lfh.write("------ sequence database code  %s sequence length %d\n" % (eId, len(self.__cleanString(eDict["sequence"]))))
                            self.__lfh.write("%s\n" % self.__cleanString(eDict["sequence"]))
                        elif eId == unpid:
                            self.__lfh.write("------ No matching isoform for %s\n" % unpid)
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % unpid)
                    self.fail()
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __cleanString(self, strIn):
        sL = []
        for ss in strIn:
            if ss in string.whitespace:
                continue
            sL.append(ss)
        return "".join(sL)


def suiteFetchVariantTests():  # pragma: no cover
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FetchUniProtEntryTests("testFetchVariantIds"))
    #
    return suiteSelect


if __name__ == "__main__":  # pragma: no cover
    mySuite = suiteFetchVariantTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
