##
# File:    BlastProcessTests.py
# Author:  j. westbrook
# Date:    2-Jan-2013
# Version: 0.001
#
# Update:

##
"""
Test cases from running reference sequence database blast searches and processing the search results.

"""

import sys
import unittest
import os
import os.path
import traceback
import platform

try:
    from unittest.mock import patch
except ImportError:
    from mock import patch
import requests_mock


try:
    from urllib.parse import parse_qs
except ImportError:
    from urlparse import parse_qs
from xml.dom import minidom

from wwpdb.utils.seqdb_v2.BlastProcess import BlastProcess
from wwpdb.utils.seqdb_v2.UnpBlastService import UnpBlastService

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TESTOUTPUT = os.path.join(HERE, "test-output", platform.python_version())
if not os.path.exists(TESTOUTPUT):
    os.makedirs(TESTOUTPUT)
mockTopPath = os.path.join(TOPDIR, "wwpdb", "mock-data")


class MyUnpBlastService(UnpBlastService):
    """A class to allow setting the check interval to speed up tests"""

    def __init__(self, sequence):
        super(MyUnpBlastService, self).__init__(sequence)
        if "MOCKREQUESTS" in os.environ:
            self._checkInterval = 0.001
            self._initWait = 0.001


def PostCallBack(request, context):
    """Simulates blast post request"""
    dat = parse_qs(request.body)

    context.status_code = 404
    if dat["program"] != ["blastp"] or dat["database"] != ["uniprotkb"] or dat["stype"] != ["protein"]:
        return ""

    smap = {
        "HMPNYKLTYFNMRGRAEIIRYIFAYLDIQYEDHRIEQADWPEIKSTLPFGKIPILEVDGLTLHQSLAIARYLTKNTDLAGNTEMEQCHVDAIVDTLDDFMSCFPWAEKKQDVKEQMFNELLTYNAPHLMQDLDTYLGGREWLIGNSVTWADFYWEICSTTLLVFKPDLLDNHPRLVTLRKKVQAIPAVANWIKRRPQTKL": "requestblast-1",  # noqa: E501
        "EEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQT": "requestblast-2",  # noqa: E501
        "DDVMTKEEQIFLLHRAQAQCEKRLKEVLQRPASIMESDKGWTSASTSGKPRKDKASGKLYPESEEDKEAPTGSRYRGRPCLPEWDHILCWPLGAPGEVVAVPCPDYIYDFNHKGHAYRRCDRNGSWELVPGHNRTWANYSECVKFLTNETREREVFDRL": "requestblast-3",  # noqa: E501
    }

    seq = dat["sequence"][0]
    if seq in smap:
        context.status_code = 200
        return smap[seq]

    return ""


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


def ResultCallBack1(_request, context):
    return commonResult("requestblast-1", context)


def ResultCallBack2(_request, context):
    return commonResult("requestblast-2", context)


def ResultCallBack3(_request, context):
    return commonResult("requestblast-3", context)


def commonResult(req, context):
    fpath = os.path.join(HERE, "refdata", "ncbiblast", req + ".xml")
    if not os.path.exists(fpath):
        context.status_code = 404
        return ""
    else:
        with open(fpath, "r") as fin:
            ret = fin.read()
        return ret


class BlastProcessTests(unittest.TestCase):
    def setUp(self):
        self.__verbose = True
        self.__lfh = sys.stderr
        self.__testModelPath = os.path.join(mockTopPath, "MODELS")
        self.__testFileCif = "4ec0.cif"
        self.__testFileFragmentsCif = "3l2j.cif"
        self.__testTaxPath = os.path.join(mockTopPath, "TAXONOMY")
        self.__taxonomyDataFile = "nodes.dmp.gz"

    def setup_mock(self, m):
        """Sets up the mock m"""

        url = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
        m.post(url + "/run/", text=PostCallBack)
        m.get(url + "/status/requestblast-1", [{"text": "PENDING"}, {"text": "RUNNING"}, {"text": "FINISHED"}])
        m.get(url + "/result/requestblast-1/xml", text=ResultCallBack1)
        m.get(url + "/status/requestblast-2", [{"text": "FINISHED"}])
        m.get(url + "/result/requestblast-2/xml", text=ResultCallBack2)
        m.get(url + "/status/requestblast-3", [{"text": "FINISHED"}])
        m.get(url + "/result/requestblast-3/xml", text=ResultCallBack3)
        m.post("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch", text=dbfetchTextCallBack)

    def testGetPolymerEntityDetails(self):
        """"""
        self.__lfh.write("\nStarting BlastProcessTests testGetPolymerEntityDetails\n")
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)

                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)

                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails["seq"])
                    self.__lfh.write("length          = %d\n" % len(polyDetails["seq"]))
                    self.__lfh.write("type            = %s\n" % polyDetails["type"])
                    self.__lfh.write("fragments       = %r\n" % polyDetails["fragments"])
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def testGetPolymerEntityDetailsFragments(self):
        """"""
        self.__lfh.write("\nStarting BlastProcessTests testGetPolymerEntityDetailsFragments\n")
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)
                ped = bp.getPolymerEntityDetails()
                for (entityId, polyDetails) in ped.items():
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    self.__lfh.write("one-letter-code = %s\n" % polyDetails["seq"])
                    self.__lfh.write("length          = %d\n" % len(polyDetails["seq"]))
                    self.__lfh.write("type            = %s\n" % polyDetails["type"])
                    self.__lfh.write("fragments       = %r\n" % polyDetails["fragments"])
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    def __polymersearch1(self):
        """The real workings of PolymerSearch1"""
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)
                ped = bp.getPolymerEntityDetails()
                for (entityId, _polyDetails) in ped.items():
                    resultlist = bp.Run(entityId=entityId)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s -------  \n" % (fn, entityId))
                    for ii, rL in enumerate(resultlist):
                        for r in rL:
                            for k in sorted(r.keys()):
                                self.__lfh.write("%d: %s:  %s\n" % (ii, k, r[k]))
        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    @patch("wwpdb.utils.seqdb_v2.BlastProcess.UnpBlastService", side_effect=MyUnpBlastService)
    def testPolymerSearch1(self, mock1):  # pylint: disable=unused-argument
        """"""
        sys.stderr.write("ABOUT TO SLEEP\n")
        self.__lfh.write("\nStarting BlastProcessTests testPolymerSearch1\n")

        if "MOCKREQUESTS" in os.environ:
            with requests_mock.Mocker(real_http=False) as m:
                self.setup_mock(m)
                self.__polymersearch1()
        else:
            self.__polymersearch1()

    def __polymerSearchAndStore(self):
        """"""
        try:
            for fn in [self.__testFileCif, self.__testFileFragmentsCif]:
                entryId, _fExt = os.path.splitext(fn)
                cifFilePath = os.path.join(self.__testModelPath, fn)
                taxonomyFilePath = os.path.join(self.__testTaxPath, self.__taxonomyDataFile)
                bp = BlastProcess(cifFilePath=cifFilePath, taxonomyFilePath=taxonomyFilePath, verbose=self.__verbose, log=self.__lfh)
                bp.saveBlastResults(blastPath=TESTOUTPUT, blastFileNamePrefix=entryId)
                ped = bp.getPolymerEntityDetails()
                for (entityId, _polyDetails) in ped.items():
                    ofn = os.path.join(TESTOUTPUT, entryId + "_seqdb-match_P" + str(entityId) + ".cif")
                    ok = bp.RunAndSave(entityId=entityId, fName=ofn)
                    self.__lfh.write("\n----Entry %s  Entity Id = %s ofn = %s status = %r -------  \n" % (fn, entityId, ofn, ok))

        except:  # noqa: E722 pylint: disable=bare-except
            traceback.print_exc(file=self.__lfh)
            self.fail()

    @patch("wwpdb.utils.seqdb_v2.BlastProcess.UnpBlastService", side_effect=MyUnpBlastService)
    def testPolymerSearchAndStore(self, mock1):  # pylint: disable=unused-argument
        """"""
        self.__lfh.write("\nStarting BlastProcessTests testPolymerSearchAndStore\n")

        if "MOCKREQUESTS" in os.environ:
            with requests_mock.Mocker(real_http=False) as m:
                self.setup_mock(m)
                self.__polymerSearchAndStore()
        else:
            self.__polymerSearchAndStore()


def suiteSearchTests():  # pragma: no cover
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(BlastProcessTests("testGetPolymerEntityDetails"))
    suiteSelect.addTest(BlastProcessTests("testGetPolymerEntityDetailsFragments"))
    suiteSelect.addTest(BlastProcessTests("testPolymerSearch1"))
    suiteSelect.addTest(BlastProcessTests("testPolymerSearchAndStore"))
    #
    return suiteSelect


if __name__ == "__main__":  # pragma: no cover
    # Run all tests --
    # unittest.main()
    #
    mySuite = suiteSearchTests()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
    #
