##
# File:    ReadUnpXml.py
# Author:  Zukang Feng
# Update:  24-May-2010
# Version: 001  Initial version
#
# Update:
#  31-Dec-2012  jdw change parser to manage files containing multiple entries
#  02-Jan-2013  jdw provide for variant registration
#  03-Jan-2013  jdw overhaul variant processing.
#  22-May-2013  jdw add support for parsing submittedName's as protein names.
#  07-Mar-2014  jdw fix concatenation problem with EC numbers --
#   3-Jun-2014  jdw pull the dataset attribute value from the entry element as db_name
#  17-Sep-2014  jdw return variant details even if these are incomplete.
#  31-Jul-2015  jdw adjust diagnostic output
#  26-Aug-2015  jdw trap what looks like a data corruption or missing data issue for isoforms
#  12-Dec-2015  jdw remove requirement for referenced isoforms - take any found c
#  12-Dec-2015  jdw update constructor for ReadUnpXmlString
#  13-Feb-2019  ep  strip newlines from sequences returned. dbfetch service started adding them.
#  29-Dec-2020  zf  add reading protein names under <protein><component>...</component></protein> tags
#
##

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

from xml.dom import minidom
import sys
import getopt
import copy
import traceback


class ReadUnpXml(object):

    """Read Uniprot entry xml file and put the following information into
    dictionary:
        dict['db_code']           - code
        dict['db_accession']      - first accession code
        dict['sequence']          - sequence
        dict['ec']                - EC number(s)
        dict['keyword']           - keywords
        dict['name']              - protein name
        dict['synonyms']          - protein synonyms
        dict['gene']              - gene
        dict['source_scientific'] - source scientific name
        dict['source_common']     - source common name
        dict['taxonomy_id']       - source taxonomy ID
        dict['comments']          - Uniprot comments

     If there is a registered variant, <isoform> tags are parsed:

        <isoform>
          <id>P42284-2</id>
          <name>V</name>
          <name>Ohsako-G</name>
          <sequence type="displayed"/>
        </isoform>
        <isoform>
          <id>P42284-3</id>
          <name>H</name>
          <name>Ohsako-M</name>
          <sequence type="described" ref="VSP_015404 VSP_015406"/>
        </isoform>

    and <feature type="splice variant"> tags:

        <feature type="splice variant" id="VSP_015404" description="(in isoform H)">
          <original>DVSTNQTVVLPHYSIYHYYSNIYYLLSHTTIYEADRTVSVSCPGKLNCLPQRNDLQETKSVTVL</original>
          <variation>DEAGQNEGGESRIRVRNWLMLADKSIIGKSSDEPSVLHIVLLLSTHRHIISFLLIIQSFIDKIY</variation>
          <location>
            <begin position="455"/>
            <end position="518"/>
          </location>
        </feature>
        <feature type="splice variant" id="VSP_015406" description="(in isoform H)">
          <location>
            <begin position="519"/>
            <end position="549"/>
          </location>
        </feature>

    to find the isoform sequence. If no match found, the default sequence from <sequence> tag
    will be used.
    """

    def __init__(self, doc, verbose=True, log=sys.stderr):
        self.__lfh = log
        self.__verbose = verbose
        self.__debug = False
        self.__variantD = {}
        self.__accessionD = {}
        self.__doc = doc
        self.__multipleResultDict = {}
        #
        self.__dictList = []
        self.__commentsDict = {}
        self.__rangesDict = {}
        self.__variantList = []

    def addVariant(self, accessionId, varId):
        """Register a variant id with the input accession code."""
        try:
            self.__variantD[varId] = accessionId
            return True
        except:  # noqa: E722 pylint: disable=bare-except
            return False
        #

    def __updateAccessionDict(self):
        """Update the list of registered variants for each accession code."""
        self.__accessionD = {}
        for vId, aId in self.__variantD.items():
            if aId not in self.__accessionD:
                self.__accessionD[aId] = []
            #
            self.__accessionD[aId].append(vId)
        #

    def __getVariantList(self, accessionId):
        try:
            return self.__accessionD[accessionId]
        except:  # noqa: E722 pylint: disable=bare-except
            return []
        #

    def getResult(self):
        self.__updateAccessionDict()
        return self._parse(self.__doc)

    def getMultipleResultDict(self):
        return self.__multipleResultDict

    def _parse(self, doc):
        entryDict = {}
        entryList = doc.getElementsByTagName("entry")
        for entry in entryList:
            # JDW 3-June-2014 Return the data set as db_name -- using PDB conventions --
            if entry.nodeType != entry.ELEMENT_NODE:
                continue
            #
            self.__dictList = []
            self.__commentsDict = {}
            self.__rangesDict = {}
            self.__variantList = []
            rdict = {}
            if entry.attributes["dataset"].value == "Swiss-Prot":
                rdict["db_name"] = "SP"
            elif entry.attributes["dataset"].value == "TrEMBL":
                rdict["db_name"] = "TR"
            else:
                rdict["db_name"] = str(entry.attributes["dataset"].value)
            #
            for node in entry.childNodes:
                if node.nodeType != node.ELEMENT_NODE:
                    continue
                #
                if node.tagName == "name":
                    # Get entry code
                    rdict["db_code"] = node.firstChild.data
                elif node.tagName == "accession":
                    # Get entry first accession
                    if "db_accession" in rdict:
                        pass
                    else:
                        rdict["db_accession"] = node.firstChild.data
                    #
                elif node.tagName == "sequence":
                    # Get sequence
                    # Sequence must have newlines removed
                    rdict["sequence"] = node.firstChild.data.replace("\n", "")
                elif node.tagName == "protein":
                    self._findProteinName(node.childNodes, rdict)
                    self._findMultipleProteinName(node.childNodes, True)
                elif node.tagName == "gene":
                    self._findGeneName(node.childNodes, rdict)
                elif node.tagName == "organism":
                    self._findSourceOrganism(node.childNodes, rdict)
                elif node.tagName == "dbReference":
                    # Get EC number from <dbReference type="EC" key="1" id="3.1.-.-"/>
                    # and concatenate them using comma separator
                    ntype = node.attributes["type"]
                    if ntype and ntype.value == "EC":
                        eid = node.attributes["id"]
                        if eid:
                            if "ec" in rdict:
                                rdict["ec"] = rdict["ec"] + ", " + eid.value
                            else:
                                rdict["ec"] = eid.value
                            #
                        #
                    #
                elif node.tagName == "keyword":
                    # Get keyword from <keyword id="KW-0181">Complete proteome</keyword>
                    # and concatenate them using comma separator
                    if "keyword" in rdict:
                        rdict["keyword"] = rdict["keyword"] + ", " + node.firstChild.data
                    else:
                        rdict["keyword"] = node.firstChild.data
                    #
                elif node.tagName == "comment":
                    self._findComments(node, rdict)
                elif node.tagName == "feature":
                    self._findFeatures(node)
                #
            #
            # This is an improbable situation of entry lacking an accession code.
            #
            if "db_accession" not in rdict:
                continue
            #
            if self.__variantList:
                rdict["variant"] = ",".join(self.__variantList)
            #
            dbAccession = rdict["db_accession"]
            #
            # Add variants if these have been specified --
            #
            vList = self.__getVariantList(dbAccession)
            if len(vList) > 0 and "sequence" in rdict:
                for vId in vList:
                    vDict = copy.deepcopy(rdict)
                    ok, seqUpdated = self._FindIsoFormSeq(doc, vId, vDict)
                    if seqUpdated:
                        vDict["isoform_sequence_updated"] = "Y"
                    else:
                        vDict["isoform_sequence_updated"] = "N"
                    if ok:
                        vDict["db_isoform"] = vId
                        entryDict[vId] = vDict
                    #
                #
            #
            entryDict[rdict["db_accession"]] = rdict
            if len(self.__dictList) > 1:
                foundComment = False
                foundError = False
                for dic in self.__dictList:
                    if ("name" not in dic) or (not dic["name"]):
                        foundError = True
                        continue
                    #
                    if (dic["name"] in self.__commentsDict) and self.__commentsDict[dic["name"]]:
                        foundComment = True
                        dic["comments"] = "\n".join(self.__commentsDict[dic["name"]])
                    #
                    if (dic["name"] in self.__rangesDict) and self.__rangesDict[dic["name"]]:
                        dic.update(self.__rangesDict[dic["name"]])
                    elif ("name" in rdict) and (dic["name"] != rdict["name"]):
                        foundError = True
                    #
                #
                if foundError:
                    continue
                #
                copyItems = ("db_code", "db_accession", "gene", "source_scientific", "source_common", "taxonomy_id", "sequence")
                if not foundComment:
                    copyItems = ("db_code", "db_accession", "gene", "source_scientific", "source_common", "taxonomy_id", "comments", "sequence")
                #
                for dic in self.__dictList:
                    if self.__variantList:
                        dic["variant"] = ",".join(self.__variantList)
                    #
                    for item in copyItems:
                        if (item in rdict) and rdict[item]:
                            dic[item] = rdict[item]
                        #
                    #
                #
                if ("keyword" in rdict) and rdict["keyword"]:
                    self.__dictList[-1]["keyword"] = rdict["keyword"]
                #
                self.__multipleResultDict[rdict["db_accession"]] = self.__dictList
            #
        #
        return entryDict

    def _findProteinName(self, nodeList, rdict):
        """In content:
          <recommendedName>
            <fullName>Platelet-derived growth factor subunit B</fullName>
            <shortName>PDGF subunit B</shortName>
          </recommendedName>
          <alternativeName>
            <fullName>Platelet-derived growth factor B chain</fullName>
          </alternativeName>
          <alternativeName>
            <fullName>Platelet-derived growth factor beta polypeptide</fullName>
          </alternativeName>
          .....
        Get protein name from <recommendedName><fullName>...</fullName></recommendedName>
        and put rest names to synonyms using comma separator
        """
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName in ("recommendedName", "alternativeName", "submittedName"):
                namelist = self._findName(node.childNodes)
                for k, v in namelist.items():
                    if k == "fullName":
                        if "name" not in rdict:
                            rdict["name"] = v
                        elif "synonyms" in rdict:
                            rdict["synonyms"] = rdict["synonyms"] + ", " + v
                        else:
                            rdict["synonyms"] = v
                        #
                    elif k == "shortName":
                        if "synonyms" in rdict:
                            rdict["synonyms"] = rdict["synonyms"] + ", " + v
                        else:
                            rdict["synonyms"] = v
                        #
                    #
                #
            #
        #

    def _findMultipleProteinName(self, nodeList, allFlag):
        rdict = {}
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName in ("recommendedName", "alternativeName", "submittedName"):
                namelist = self._findName(node.childNodes)
                for k, v in namelist.items():
                    if k == "fullName":
                        if "name" not in rdict:
                            rdict["name"] = v
                        elif "synonyms" in rdict:
                            rdict["synonyms"] = rdict["synonyms"] + ", " + v
                        else:
                            rdict["synonyms"] = v
                        #
                    elif k == "shortName":
                        if "synonyms" in rdict:
                            rdict["synonyms"] = rdict["synonyms"] + ", " + v
                        else:
                            rdict["synonyms"] = v
                        #
                    elif k == "ecNumber":
                        if "ec" in rdict:
                            rdict["ec"] = rdict["ec"] + ", " + v
                        else:
                            rdict["ec"] = v
                        #
                    #
                #
            elif node.tagName == "component":
                self._findMultipleProteinName(node.childNodes, False)
            #
        #
        if rdict:
            if allFlag:
                rdict["all"] = "yes"
            #
            self.__dictList.append(rdict)
        #

    def _findName(self, nodeList):
        """Get names from <fullName> & <shortName> tags:

        <fullName>Platelet-derived growth factor subunit B</fullName>
        <shortName>PDGF subunit B</shortName>
        """
        d = {}
        tagName = ("fullName", "shortName", "ecNumber")
        for node in nodeList:
            if (node.nodeType != node.ELEMENT_NODE) or (node.tagName not in tagName):
                continue
            #
            if node.tagName in d:
                d[node.tagName] = d[node.tagName] + ", " + node.firstChild.data
            else:
                d[node.tagName] = node.firstChild.data
            #
        #
        return d

    def _findGeneName(self, nodeList, rdict):
        """Get genes from
           <gene>
             <name type="primary">PDGFB</name>
             <name type="synonym">PDGF2</name>
             <name type="synonym">SIS</name>
           </gene>
        and concatenate them using comma separator
        """
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName == "name":
                if "gene" in rdict:
                    rdict["gene"] = rdict["gene"] + ", " + node.firstChild.data
                else:
                    rdict["gene"] = node.firstChild.data
                #
            #
        #

    def _findSourceOrganism(self, nodeList, rdict):
        """Get organism's scientific name, common name and NCBI Taxonomy ID from
        <name type="scientific">Homo sapiens</name>
        <name type="common">Human</name>
        <dbReference type="NCBI Taxonomy" key="1" id="9606"/>
        """
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName == "name":
                ntype = node.attributes["type"]
                if ntype:
                    if ntype.value == "scientific":
                        rdict["source_scientific"] = node.firstChild.data
                    elif (ntype.value == "common") or (ntype.value == "synonym"):
                        if "source_common" not in rdict:
                            rdict["source_common"] = node.firstChild.data
                        else:
                            rdict["source_common"] = rdict["source_common"] + ", " + node.firstChild.data
                        #
                    #
                #
            elif node.tagName == "dbReference":
                ntype = node.attributes["type"]
                if ntype and ntype.value == "NCBI Taxonomy":
                    tid = node.attributes["id"]
                    if tid:
                        rdict["taxonomy_id"] = tid.value
                    #
                #
            #
        #

    def _findComments(self, node, cDict):
        """From
           <comment type="function">
             <text>Platelet-derived .... </text>
           </comment>
           <comment type="subunit" evidence="EC1">
             <text status="by similarity">Antiparallel disulfide-linked .... </text>
           </comment>
           <comment type="miscellaneous">
             <text>A-A and B-B, as well as A-B, dimers can bind to the PDGF receptor.</text>
           </comment>
           <comment type="similarity">
             <text>Belongs to the PDGF/VEGF growth factor family.</text>
           </comment>
           .....
           <comment type="online information" name="Regranex">
             <link uri="http://www.regranex.com/"/>
             <text>Clinical information on Regranex</text>
           </comment>
        Get "type": "text" content and concatenate them using newline separator
        Comments from <comment type="online information"> will be ignored.
        """

        ntype = node.attributes["type"]
        if ntype and ntype.value != "online information":
            text, molecule = self._findText(node.childNodes)
            if text:
                if "comments" in cDict:
                    cDict["comments"] = cDict["comments"] + "\n" + ntype.value + ": " + text
                else:
                    cDict["comments"] = ntype.value + ": " + text
                #
                if molecule:
                    if molecule in self.__commentsDict:
                        self.__commentsDict[molecule].append(ntype.value + ": " + text)
                    else:
                        self.__commentsDict[molecule] = [ntype.value + ": " + text]
                    #
                #
            #
        #

    def _findText(self, nodeList):
        """Get text value from
        <text status="by similarity">Antiparallel disulfide-linked .... </text>
        """
        text = ""
        molecule = ""
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            #
            if node.tagName == "text":
                text = node.firstChild.data
            elif node.tagName == "molecule":
                molecule = node.firstChild.data
            #
        #
        return text, molecule

    def _findFeatures(self, node):
        """"""
        if node.hasAttribute("type") is False:
            return
        #
        ntype = node.attributes["type"]
        if not ntype:
            return
        #
        if (ntype.value == "chain") or (ntype.value == "peptide"):
            if node.hasAttribute("description") is False:
                return
            #
            ndescription = node.attributes["description"]
            if not ndescription:
                return
            #
            rDict = {}
            for node1 in node.childNodes:
                if (node1.nodeType != node1.ELEMENT_NODE) or (node1.tagName != "location"):
                    continue
                #
                for node2 in node1.childNodes:
                    if node2.nodeType != node2.ELEMENT_NODE:
                        continue
                    #
                    if node2.tagName in ("begin", "end"):
                        try:
                            rDict[node2.tagName] = int(node2.attributes["position"].value)
                        except:  # noqa: E722 pylint: disable=bare-except
                            traceback.print_exc(file=self.__lfh)
                        #
                    #
                #
            #
            if len(rDict) == 2:
                self.__rangesDict[ndescription.value] = rDict
            #
        elif ntype.value == "sequence variant":
            original = ""
            variation = ""
            position = ""
            for node1 in node.childNodes:
                if node1.nodeType != node1.ELEMENT_NODE:
                    continue
                #
                if node1.tagName == "original":
                    original = node1.firstChild.data
                elif node1.tagName == "variation":
                    variation = node1.firstChild.data
                elif node1.tagName == "location":
                    for node2 in node1.childNodes:
                        if node2.nodeType != node2.ELEMENT_NODE:
                            continue
                        #
                        if node2.tagName == "position":
                            try:
                                position = str(node2.attributes["position"].value)
                            except:  # noqa: E722 pylint: disable=bare-except
                                traceback.print_exc(file=self.__lfh)
                            #
                        #
                    #
                #
            #
            if original and variation and position:
                variant = original + position + variation
                if variant not in self.__variantList:
                    self.__variantList.append(variant)
                #
            #
        #

    def _FindIsoFormSeq(self, doc, vId, vDict):
        """Get isoform sequence for vId if it exists  -  """
        if self.__debug:
            self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - starting vId %s dict %r\n" % (vId, vDict.items()))
        try:
            isoformdic = self._FindIsoFormIds(doc)

            if self.__debug:
                self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - vId %s isoformdic  %r\n" % (vId, isoformdic.items()))

            if not isoformdic:
                return False, False

            if vId not in isoformdic:
                return False, False

            # JDW 12-DEC-2015  Remove this test for a reference - Finding isoform data in comments lacking this
            if isoformdic[vId]["type"] == "displayed" or "ref" not in isoformdic[vId]:
                # return True, False
                return True, False

            refdic = self._FindIsoFormRefs(doc)
            if self.__debug:
                self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - vId %s refdic  %r\n" % (vId, refdic.items()))

            if not refdic:
                # return with sequence updated = False
                return True, False

            reflist = isoformdic[vId]["ref"].split(" ")
            # Reverse the ref list order so that sequence manipulation starts from C-terminal
            reflist.reverse()
            for ref in reflist:
                if ref in refdic:
                    vDict["sequence"] = self._ProcessIsoFormSeq(vDict["sequence"], refdic[ref])
            # return with seqquence updated = True
            return True, True
        except Exception as e:
            self.__lfh.write("_findIsoFormSeq() failing -- %s %s\n" % (vId, str(e)))
            traceback.print_exc(file=self.__lfh)
        return False, False

    def _FindIsoFormIds(self, doc):
        """Get isoform information from:
            <isoform>
              <id>P42284-2</id>
              <name>V</name>
              <name>Ohsako-G</name>
              <sequence type="displayed"/>
            </isoform>
            <isoform>
              <id>P42284-3</id>
              <name>H</name>
              <name>Ohsako-M</name>
              <sequence type="described" ref="VSP_015404 VSP_015406"/>
            </isoform>

        and put them into dictionary:

            { 'P42284-2' : { 'type' : 'displayed'},
              'P42284-3' : { 'type' : 'described', 'ref' : 'VSP_015404 VSP_015406' } }
        """
        dic = {}
        entryList = doc.getElementsByTagName("isoform")
        if not entryList:
            return dic

        for node in entryList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            # id = None
            lid = []
            stype = None
            ref = None
            for node1 in node.childNodes:
                if node1.nodeType != node1.ELEMENT_NODE:
                    continue
                if node1.tagName == "id":
                    lid.append(node1.firstChild.data)

                elif node1.tagName == "sequence":
                    stype = node1.attributes["type"].value
                    # JDW Aug-26 The following is behaving badly  --
                    try:
                        if "ref" in node1.attributes:
                            ref = node1.attributes["ref"].value
                    except:  # noqa: E722 pylint: disable=bare-except
                        pass

            if len(lid) < 1 or not stype:
                continue
            d = {}
            d["type"] = stype
            if ref:
                d["ref"] = ref
            dic[lid[0]] = d

        return dic

    def _FindIsoFormRefs(self, doc):
        """Get variant information from

            <feature type="splice variant" id="VSP_015404" description="(in isoform H)">
              <original>DVSTNQTVVLPHYSIYHYYSNIYYLLSHTTIYEADRTVSVSCPGKLNCLPQRNDLQETKSVTVL</original>
              <variation>DEAGQNEGGESRIRVRNWLMLADKSIIGKSSDEPSVLHIVLLLSTHRHIISFLLIIQSFIDKIY</variation>
              <location>
                <begin position="455"/>
                <end position="518"/>
              </location>
            </feature>
            <feature type="splice variant" id="VSP_015406" description="(in isoform H)">
              <location>
                <begin position="519"/>
                <end position="549"/>
              </location>
            </feature>

        and put them into dictionary:

            { 'VSP_015404' : { 'begin' : '455', 'end' : '518', 'variation' : 'DEAGQNEGG....' },
              'VSP_015406' : { 'begin' : '519', 'end' : '549' } }
        """
        dic = {}
        entryList = doc.getElementsByTagName("feature")
        if not entryList:
            return dic

        for node in entryList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.attributes["type"].value != "splice variant":
                continue

            if "id" not in node.attributes:
                continue

            aid = node.attributes["id"].value
            begin = None
            end = None
            variation = None
            for node1 in node.childNodes:
                if node1.nodeType != node1.ELEMENT_NODE:
                    continue
                if node1.tagName == "variation":
                    variation = node1.firstChild.data.replace("\n", "")
                elif node1.tagName == "location":
                    for node2 in node1.childNodes:
                        if node2.nodeType != node2.ELEMENT_NODE:
                            continue
                        if node2.tagName == "begin":
                            begin = node2.attributes["position"].value
                        elif node2.tagName == "end":
                            end = node2.attributes["position"].value

            if not begin or not end:
                continue
            d = {}
            d["begin"] = begin
            d["end"] = end
            if variation:
                d["variation"] = variation
            dic[aid] = d

        return dic

    def _ProcessIsoFormSeq(self, seq, ref):
        """Manipulate sequence using information from dictionary ref:

        { 'begin' : '455', 'end' : '518', 'variation' : 'DEAGQNEGG....' }
        """
        begin = int(ref["begin"]) - 1
        end = int(ref["end"])
        seq1 = seq[0:begin]
        if "variation" in ref:
            seq1 += ref["variation"]
        seq1 += seq[end:]
        return seq1


class ReadUnpXmlFile(ReadUnpXml):

    """"""

    def __init__(self, fileName):
        self._fileName = fileName
        self._doc = minidom.parse(self._fileName)
        ReadUnpXml.__init__(self, self._doc)


class ReadUnpXmlString(ReadUnpXml):

    """"""

    def __init__(self, data, verbose=False, log=sys.stderr):
        self._data = data
        self._doc = minidom.parseString(self._data)
        ReadUnpXml.__init__(self, self._doc, verbose=verbose, log=log)


def main(argv):  # pragma: no cover
    opts, _args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadUnpXmlFile(arg)
            rdict = obj.getResult()
            for (k, v) in rdict.items():
                print("%s=%s" % (k, v))


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)
