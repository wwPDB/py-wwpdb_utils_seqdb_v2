##
# File:    FetchUnpXmlTests.py
# Author:  E. Peisach
# Date:    15-Oct-2020
# Version: 0.001
#
##
"""
Test cases for fetching using NCBI XML

"""

import sys
import unittest
import os
import os.path
import platform
import requests_mock
import time

from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchNcbiXml, FetchFullNcbiXml

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, "test-output", platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)


# dbfetch parses the request and returns data
def dbfetchSummaryCallBack(request, context):
    """Handles dbfetch requests"""
    dat = request.qs
    accs = dat["id"][0]

    fpath = os.path.join(HERE, "refdata", "ncbiquery", accs + "-summary.xml")
    if not os.path.exists(fpath):
        context.status_code = 404
        # print("XXX COULD NOT FIND", fpath)
        return ""

    context.status_code = 200
    with open(fpath, "r") as fin:
        resp = fin.read()  # pylint: disable=protected-access
    return resp


def dbfetchFullCallBack(request, context):
    """Handles dbfetch requests"""
    dat = request.qs
    accs = dat["id"][0]

    fpath = os.path.join(HERE, "refdata", "ncbiquery", accs + "-full-summary.xml")
    if not os.path.exists(fpath):
        context.status_code = 404
        # print("XXX COULD NOT FIND", fpath)
        return ""

    context.status_code = 200
    with open(fpath, "r") as fin:
        resp = fin.read()  # pylint: disable=protected-access
    return resp


class FetchNcbiXmlTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__nucList = ["21614549"]
        self.__mock = None
        if "MOCKREQUESTS" in os.environ:
            self.__mock = requests_mock.Mocker()
            self.__mock.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi", text=dbfetchSummaryCallBack)
            self.__mock.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", text=dbfetchFullCallBack)
            self.__mock.start()

    def tearDown(self):
        if self.__mock is not None:
            self.__mock.stop()

    def testFetchSummaries(self):
        """Fetch NCBI Summaries"""
        self.__lfh.write("\nStarting FetchNcbiXmlTests testFetchIds\n")
        for gid in self.__nucList:
            fobj = FetchNcbiXml(database="nuccore", qid=gid)
            pd = fobj.ParseNcbiXmlData()
            # print("PD is %s" % pd)
            self.assertIn("name", pd)
            self.assertIn("taxonomy_id", pd)
            ok = fobj.WriteNcbiXml(os.path.join(TESTOUTPUT, "%s-summary.xml" % gid))
            self.assertTrue(ok)
            # Limit no more than 3 requests/sec
            if self.__mock is None:
                time.sleep(0.5)

    def testFetchFullNcbiXml(self):
        """Fetch Full NCBI XML Summaries"""
        self.__lfh.write("\nStarting FetchNcbiXmlTests testFetchFullNcbiXml\n")
        for qid in self.__nucList:
            self.__lfh.write("About to do %s\n" % qid)
            fobj = FetchFullNcbiXml(database="nuccore", qid=qid)
            fobj.WriteNcbiXml(os.path.join(TESTOUTPUT, "%s-full-summary.xml" % qid))
            # No status
            pd = fobj.ParseNcbiXmlData()
            self.assertIn("db_code", pd)
            self.assertIn("taxonomy_id", pd)
            if self.__mock is None:
                time.sleep(0.5)


def suiteFetchTests():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(FetchNcbiXmlTests("testFetchSummaries"))
    suiteSelect.addTest(FetchNcbiXmlTests("testFetchFullNcbiXml"))
    #
    return suiteSelect


if __name__ == "__main__":
    #
    mySuite = suiteFetchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
