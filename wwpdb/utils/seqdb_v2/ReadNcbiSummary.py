##
# File:    ReadNcbiSummary.py
# Author:  Zukang Feng
# Update:  24-May-2010
# Version: 001  Initial version
#
# Updates -
#  20-Apr-2013 jdw  Capture title attribute in summary as key name
#
##

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

from xml.dom import minidom
import sys
import getopt


class ReadNcbiSummary:
    """Parse eSummaryResult xml file from NCBI and put source organism and
    taxonomy information into dictionary:
    dict['source_scientific']
    dict['source_common']
    dict['taxonomy_id']
    dict['name']
    """

    def __init__(self, doc):
        self._dict = {}
        if doc:
            self._dict = self._parse(doc)

    def GetResult(self):
        return self._dict

    def _parse(self, doc):
        """Parse the eSummaryResult xml file:

            <eSummaryResult>
            <DocSum>
                    <Id>300852</Id>
                    <Item Name="Title" Type="String">Thermus thermophilus HB8 xxxxxx  xxxxxx </Item>
                    <Item Name="ScientificName" Type="String">Thermus thermophilus HB8</Item>
                    <Item Name="CommonName" Type="String"></Item>
                    <Item Name="TaxId" Type="Integer">300852</Item>
            /DocSum>
            </eSummaryResult>

        and put values into dictionary:

        dict['source_scientific']
        dict['source_common']
        dict['taxonomy_id']
        dict['name']
        """
        rdict = {}

        entryList = doc.getElementsByTagName("Item")
        if entryList:
            for node in entryList:
                if node.nodeType != node.ELEMENT_NODE:
                    continue

                name = node.attributes["Name"]
                if not name:
                    continue
                if name.value == "TaxId":
                    if node.firstChild and node.firstChild.data != "0":
                        rdict["taxonomy_id"] = node.firstChild.data
                elif name.value == "ScientificName":
                    if node.firstChild:
                        rdict["source_scientific"] = node.firstChild.data
                elif name.value == "CommonName":
                    if node.firstChild:
                        rdict["source_common"] = node.firstChild.data
                elif name.value == "Title":
                    if node.firstChild:
                        rdict["name"] = node.firstChild.data

        return rdict


class ReadNcbiSummaryFile(ReadNcbiSummary):
    """"""

    def __init__(self, fileName):
        self._fileName = fileName
        self._doc = ""
        try:
            self._doc = minidom.parse(self._fileName)
        except Exception as exc:  # noqa: F841 pylint: disable=unused-variable
            pass
        ReadNcbiSummary.__init__(self, self._doc)


class ReadNcbiSummaryString(ReadNcbiSummary):
    """"""

    def __init__(self, data):
        self._data = data
        self._doc = ""
        try:
            self._doc = minidom.parseString(self._data)
        except Exception as exc:  # noqa: F841 pylint: disable=unused-variable
            pass
        ReadNcbiSummary.__init__(self, self._doc)


def main(argv):  # pragma: no cover
    opts, _args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadNcbiSummaryFile(arg)
            rdict = obj.GetResult()
            for (k, v) in rdict.items():
                print("%s=%s" % (k, v))


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)
