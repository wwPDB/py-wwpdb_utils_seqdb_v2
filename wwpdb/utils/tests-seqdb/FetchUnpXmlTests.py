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
import platform
import traceback
import requests
import requests_mock

try:
    from urllib.parse import parse_qs
except ImportError:
    from urlparse import parse_qs
from xml.dom import minidom

from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, "test-output", platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)


# Mock request handler
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
            # print("XXX COULD NOT FIND", fpath)
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


class FetchUnpXmlTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        # Pick up site information from the environment or failover to the development site id.
        #
        self.__unpIdList1 = ["P20937", "P21877", "P22868", "P23832", "P25665", "P26562", "P27614"]
        self.__unpIdList2 = [
            "P29490",
            "P29496",
            "P29498",
            "P29499",
            "P29503",
            "P29506",
            "P29508",
            "P29509",
            "P29525",
            "P29533",
            "P29534",
            "P29547",
            "P29549",
            "P29555",
            "P29557",
            "P29558",
            "P29559",
            "P29563",
            "P29588",
            "P29589",
            "P29590",
            "P29597",
            "P29599",
            "P29600",
            "P29602",
            "P29603",
            "P29617",
            "P29678",
            "P29691",
            "P29715",
            "P29717",
            "P29723",
            "P29724",
            "P29736",
            "P29741",
            "P29745",
            "P29748",
            "P29749",
            "P29752",
            "P29758",
            "P29768",
            "P29803",
            "P29808",
            "P29813",
            "P29827",
            "P29830",
            "P29837",
            "P29838",
            "P29846",
            "P29848",
            "P29882",
            "P29894",
            "P29898",
            "P29899",
            "P29929",
            "P29946",
            "P29957",
            "P29960",
            "P29965",
            "P29966",
            "P29972",
            "P29978",
            "P29986",
            "P29987",
            "P29988",
            "P29989",
            "P29990",
            "P29991",
            "P29994",
        ]

        self.__unpIdListV = ["P42284", "P42284-1", "P42284-2", "P42284-3", "P29994-1", "P29994-2", "P29994-3", "P29994-4", "P29994-5", "P29994-6", "P29994-7"]
        self.__mock = None
        if "MOCKREQUESTS" in os.environ:
            self.__mock = requests_mock.Mocker()
            self.__mock.post("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch", text=dbfetchTextCallBack)
            self.__mock.add_matcher(isoformMatcher)
            self.__mock.start()

    def tearDown(self):
        if self.__mock is not None:
            self.__mock.stop()

    def testFetchIds(self):
        """"""
        self.__lfh.write("\nStarting FetchUnpXmlTests testFetchIds\n")
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            for unpid in self.__unpIdList1:
                ok = fobj.fetchList([unpid])
                if ok:
                    fobj.writeUnpXml(os.path.join(TESTOUTPUT, unpid + ".xml"))
                    rdict = fobj.getResult()
                    for (k, v) in rdict.items():
                        self.__lfh.write("%s=%s" % (k, v))
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % unpid)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBatchFetch(self):
        """"""
        self.__lfh.write("\nStarting FetchUnpXmlTests testBatchFetch\n")
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            ok = fobj.fetchList(self.__unpIdList2)
            if ok:
                fobj.writeUnpXml(os.path.join(TESTOUTPUT, "batch-fetch.xml"))
                rdict = fobj.getResult()
                for (k, v) in rdict.items():
                    self.__lfh.write("%s=%s" % (k, v))
            else:
                self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)
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

    def testFetchVariantIds(self):
        """"""
        self.__lfh.write("\nStarting FetchUnpXmlTests testFetchVariantIds\n")
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            for unpid in self.__unpIdListV:
                ok = fobj.fetchList([unpid])
                if ok:
                    fobj.writeUnpXml(os.path.join(TESTOUTPUT, unpid + ".xml"))
                    rdict = fobj.getResult()
                    for (eId, eDict) in rdict.items():
                        if eId == unpid and ("db_isoform" in eDict or len(eDict.get("sequence", "")) > 0):
                            if "db_isoform" in eDict:
                                self.__lfh.write("------ sequence database code  %s has key db_isoform:  %r\n" % (eId, eDict["db_isoform"]))
                            self.__lfh.write("------ sequence database code  %s sequence length %d\n" % (eId, len(self.__cleanString(eDict["sequence"]))))
                            # self.__lfh.write("------ sequence database code  %s keys %r\n" % (eId,eDict.keys()))
                            self.__lfh.write("%s\n" % self.__cleanString(eDict["sequence"]))
                        elif eId == unpid:
                            self.__lfh.write("------ No matching isoform for %s\n" % unpid)
                        # for k,v in eDict.items():
                        #    self.__lfh.write("%-25s = %s\n" % (k, v))
                else:
                    self.__lfh.write("+WARNING - Fetch failed for id %s\n" % unpid)
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testBatchFetchVariants(self):
        """"""
        self.__lfh.write("\nStarting FetchUnpXmlTests testBatchFetchVariants\n")
        fail = False
        try:
            fobj = FetchUnpXml(verbose=self.__verbose, log=self.__lfh)
            ok = fobj.fetchList(self.__unpIdListV)
            if ok:
                fobj.writeUnpXml(os.path.join(TESTOUTPUT, "variant-batch-fetch.xml"))
                rdict = fobj.getResult()
                for (eId, eDict) in rdict.items():
                    self.__lfh.write("\n\n------ Entry id %s\n" % eId)
                    for k, v in eDict.items():
                        self.__lfh.write("%-25s = %s\n" % (k, v))
            else:
                self.__lfh.write("+WARNING - Fetch failed for id %s\n" % id)
                fail = True

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

        if fail:
            self.fail()


def suiteFetchTests():
    suiteSelect = unittest.TestSuite()
    # suiteSelect.addTest(FetchUnpXmlTests("testFetchIds"))
    suiteSelect.addTest(FetchUnpXmlTests("testBatchFetch"))
    #
    return suiteSelect


def suiteFetchVariantTests():  # pragma: no cover
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FetchUnpXmlTests("testFetchVariantIds"))
    # suiteSelect.addTest(FetchUnpXmlTests("testBatchFetchVariants"))
    #
    return suiteSelect


if __name__ == "__main__":  # pragma: no cover
    mySuite = suiteFetchVariantTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
    #
