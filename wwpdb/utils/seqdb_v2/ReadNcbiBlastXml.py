"""
File:    ReadNcbiBlastXml.py
Author:  Zukang Feng
Update:  24-May-2010
Version: 001  Initial version

"""

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

from xml.dom import minidom
import sys
import getopt
from wwpdb.utils.seqdb_v2.FetchNcbiXml import FetchNcbiXml


class ReadNcbiBlastXml:
    """Read Uniprot blast result xml file and return the list of hits which contains:

    dict['db_name']
    dict['db_code']
    dict['db_accession']
    dict['identity']
    dict['positive']
    dict['gaps']
    dict['midline']
    dict['query']
    dict['queryFrom']
    dict['queryTo']
    dict['subject']
    dict['hitFrom']
    dict['hitTo']
    dict['alignLen']
    dict['match_length']
    dict['source_scientific']
    dict['source_common']
    dict['taxonomy_id']
    """

    def __init__(self, doc):
        self._result = self._parse(doc)

    def GetResult(self):
        return self._result

    def _parse(self, doc):
        resultlist = []
        hlist = doc.getElementsByTagName("Hit")
        if not hlist:
            return resultlist

        length = None
        for node in doc.getElementsByTagName("BlastOutput_query-len"):
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == "BlastOutput_query-len":
                length = node.firstChild.data
                break

        for node in hlist:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            alignlist = self._ProcessHitTag(node.childNodes, length)
            if alignlist:
                for align in alignlist:
                    if length:
                        align["query_length"] = length
                    resultlist.append(align)

        return resultlist

    def _ProcessHitTag(self, nodelist, length):
        resultlist = []
        alignlist = []
        description = ""
        length = ""
        gi = ""
        code = ""
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == "Hit_id":
                gi = self._parseID(node.firstChild.data)
                if not gi:
                    return resultlist
            elif node.tagName == "Hit_accession":
                code = node.firstChild.data
            elif node.tagName == "Hit_def":
                description = node.firstChild.data
            elif node.tagName == "Hit_len":
                length = node.firstChild.data
            elif node.tagName == "Hit_hsps":
                hlist = self._ProcessHit_hspsTag(node.childNodes, length)
                if hlist:
                    for li in hlist:
                        alignlist.append(li)

        if not gi or not alignlist:
            return resultlist

        # Get genbank sequence information
        fetchobj = FetchNcbiXml(gi, "Nucleotide")
        pdict = fetchobj.ParseNcbiXmlData()
        if pdict and "taxonomy_id" in pdict:
            fetchobj = FetchNcbiXml(pdict["taxonomy_id"], "taxonomy")
            pdict = fetchobj.ParseNcbiXmlData()

        for align in alignlist:
            align["db_name"] = "GB"
            align["db_code"] = code
            align["db_accession"] = gi
            align["db_description"] = description
            align["db_length"] = length

            # add genbank sequence information
            if pdict:
                for (k, v) in pdict.items():
                    if k not in align:
                        align[k] = v

            resultlist.append(align)

        return resultlist

    def _parseID(self, data):
        # Parses NCBI fasta descriptor. Ignore PDB references
        gi = ""
        dlist = data.split("|")
        if len(dlist) >= 2 and dlist[2] != "pdb":
            gi = dlist[1]
        return gi

    def _ProcessHit_hspsTag(self, nodelist, length):
        resultlist = []
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            if node.tagName != "Hsp":
                continue

            rdict = self._GetMatchAlignment(node.childNodes, length)
            if rdict:
                resultlist.append(rdict)
        return resultlist

    def _GetMatchAlignment(self, nodelist, length):
        rdict = {}
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == "Hsp_identity":
                rdict["identity"] = node.firstChild.data
            elif node.tagName == "Hsp_positive":
                rdict["positive"] = node.firstChild.data
            elif node.tagName == "Hsp_gaps":
                rdict["gaps"] = node.firstChild.data
            elif node.tagName == "Hsp_midline":
                rdict["midline"] = node.firstChild.data.replace("T", "U")
            elif node.tagName == "Hsp_qseq":
                rdict["query"] = node.firstChild.data.replace("T", "U")
            elif node.tagName == "Hsp_query-from":
                rdict["queryFrom"] = node.firstChild.data
            elif node.tagName == "Hsp_query-to":
                rdict["queryTo"] = node.firstChild.data
            elif node.tagName == "Hsp_hseq":
                rdict["subject"] = node.firstChild.data.replace("T", "U")
            elif node.tagName == "Hsp_hit-from":
                rdict["hitFrom"] = node.firstChild.data
            elif node.tagName == "Hsp_hit-to":
                rdict["hitTo"] = node.firstChild.data
            elif node.tagName == "Hsp_align-len":
                rdict["alignLen"] = node.firstChild.data
                rdict["match_length"] = node.firstChild.data

        if "identity" in rdict:
            if length:
                identity = int(rdict["identity"]) * 100 / int(length)
                if identity < 70:
                    rdict.clear()
            elif "alignLen" in rdict:
                identity = int(rdict["identity"]) * 100 / int(rdict["alignLen"])
                if identity < 70:
                    rdict.clear()

        return rdict


class ReadNcbiBlastXmlFile(ReadNcbiBlastXml):
    """"""

    def __init__(self, fileName):
        self._fileName = fileName
        self._doc = minidom.parse(self._fileName)
        ReadNcbiBlastXml.__init__(self, self._doc)


class ReadNcbiBlastXmlString(ReadNcbiBlastXml):
    """"""

    def __init__(self, data):
        self._data = data
        self._doc = minidom.parseString(self._data)
        ReadNcbiBlastXml.__init__(self, self._doc)


def main(argv):  # pragma: no cover
    opts, _args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadNcbiBlastXmlFile(arg)
            result = obj.GetResult()
            for align in result:
                for (k, v) in align.items():
                    print("%s=%s" % (k, v))
                print("\n")


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)
