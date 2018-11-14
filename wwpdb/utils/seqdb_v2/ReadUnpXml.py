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


class ReadUnpXml:

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

    def addVariant(self, accessionId, varId):
        """ Register a variant id with the input accession code.
        """
        try:
            self.__variantD[varId] = accessionId
            return True
        except:
            return False

    def __updateAccessionDict(self):
        """ Update the list of registered variants for each accession code.
        """
        self.__accessionD = {}
        for vId, aId in self.__variantD.items():
            if aId not in self.__accessionD:
                self.__accessionD[aId] = []
            self.__accessionD[aId].append(vId)

    def __getVariantList(self, accessionId):
        try:
            return self.__accessionD[accessionId]
        except:
            return []

    def getResult(self):
        self.__updateAccessionDict()
        self._entryDict = self._parse(self.__doc)
        return self._entryDict

    def _parse(self, doc):
        entryDict = {}
        entryList = doc.getElementsByTagName('entry')

        entryDict = {}
        for entry in entryList:

            # JDW 3-June-2014 Return the data set as db_name -- using PDB conventions --
            if entry.nodeType != entry.ELEMENT_NODE:
                continue

            dict = {}
            if entry.attributes['dataset'].value == 'Swiss-Prot':
                dict['db_name'] = 'SP'
            elif entry.attributes['dataset'].value == 'TrEMBL':
                dict['db_name'] = 'TR'
            else:
                dict['db_name'] = str(entry.attributes['dataset'].value)

            for node in entry.childNodes:
                if node.nodeType != node.ELEMENT_NODE:
                    continue

                if node.tagName == 'name':
                    'Get entry code'
                    dict['db_code'] = node.firstChild.data

                elif node.tagName == 'accession':
                    'Get entry first accession'
                    if 'db_accession' in dict:
                        pass
                    else:
                        dict['db_accession'] = node.firstChild.data

                elif node.tagName == 'sequence':
                    'Get sequence'
                    dict['sequence'] = node.firstChild.data

                elif node.tagName == 'protein':
                    self._findProteinName(node.childNodes, dict)

                elif node.tagName == 'gene':
                    self._findGeneName(node.childNodes, dict)

                elif node.tagName == 'organism':
                    self._findSourceOrganism(node.childNodes, dict)

                elif node.tagName == 'dbReference':
                    """Get EC number from <dbReference type="EC" key="1" id="3.1.-.-"/>
                        and concatenate them using comma separator
                    """
                    type = node.attributes['type']
                    if type and type.value == 'EC':
                        id = node.attributes['id']
                        if id:
                            if 'ec' in dict:
                                dict['ec'] = dict['ec'] + ', ' + id.value
                            else:
                                dict['ec'] = id.value

                elif node.tagName == 'keyword':
                    """Get keyword from <keyword id="KW-0181">Complete proteome</keyword>
                        and concatenate them using comma separator
                    """
                    if 'keyword' in dict:
                        dict['keyword'] = dict['keyword'] + ', ' + node.firstChild.data
                    else:
                        dict['keyword'] = node.firstChild.data

                elif node.tagName == 'comment':
                    self._findComments(node, dict)

            #
            # This is an improbable situation of entry lacking an accession code.
            #
            if 'db_accession' not in dict:
                continue

            dbAccession = dict['db_accession']
            #
            # Add variants if these have been specified --
            #
            vList = self.__getVariantList(dbAccession)
            if len(vList) > 0 and 'sequence' in dict:
                for vId in vList:
                    vDict = copy.deepcopy(dict)
                    ok, seqUpdated = self._FindIsoFormSeq(doc, vId, vDict)
                    if seqUpdated:
                        vDict['isoform_sequence_updated'] = 'Y'
                    else:
                        vDict['isoform_sequence_updated'] = 'N'
                    if ok:
                        vDict['db_isoform'] = vId
                        entryDict[vId] = vDict

            entryDict[dict['db_accession']] = dict

        return entryDict

    def _findProteinName(self, nodeList, dict):
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

            if node.tagName == 'recommendedName':
                namelist = self._findName(node.childNodes)
                for k, v in namelist.items():
                    if k == 'fullName':
                        dict['name'] = v
                    elif k == 'shortName':
                        if 'synonyms' in dict:
                            dict['synonyms'] = dict['synonyms'] + ', ' + v
                        else:
                            dict['synonyms'] = v
            elif node.tagName == 'alternativeName':
                namelist = self._findName(node.childNodes)
                for v in namelist.values():
                    if 'synonyms' in dict:
                        dict['synonyms'] = dict['synonyms'] + ', ' + v
                    else:
                        dict['synonyms'] = v
            elif node.tagName == 'submittedName':
                namelist = self._findName(node.childNodes)
                for k, v in namelist.items():
                    if k == 'fullName' and 'name' not in dict:
                        dict['name'] = v
                    elif 'synonyms' in dict:
                        dict['synonyms'] = dict['synonyms'] + ', ' + v
                    else:
                        dict['synonyms'] = v

    def _findName(self, nodeList):
        """Get names from <fullName> & <shortName> tags:

                <fullName>Platelet-derived growth factor subunit B</fullName>
                <shortName>PDGF subunit B</shortName>
        """
        d = {}
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'fullName':
                d['fullName'] = node.firstChild.data
            elif node.tagName == 'shortName':
                d['shortName'] = node.firstChild.data
        return d

    def _findGeneName(self, nodeList, dict):
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

            if node.tagName == 'name':
                if 'gene' in dict:
                    dict['gene'] = dict['gene'] + ', ' + node.firstChild.data
                else:
                    dict['gene'] = node.firstChild.data

    def _findSourceOrganism(self, nodeList, dict):
        """Get organism's scientific name, common name and NCBI Taxonomy ID from
               <name type="scientific">Homo sapiens</name>
               <name type="common">Human</name>
               <dbReference type="NCBI Taxonomy" key="1" id="9606"/>
        """
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'name':
                type = node.attributes['type']
                if type:
                    if type.value == 'scientific':
                        dict['source_scientific'] = node.firstChild.data
                    elif type.value == 'common':
                        dict['source_common'] = node.firstChild.data

            elif node.tagName == 'dbReference':
                type = node.attributes['type']
                if type and type.value == 'NCBI Taxonomy':
                    id = node.attributes['id']
                    if id:
                        dict['taxonomy_id'] = id.value

    def _findComments(self, node, dict):
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

        type = node.attributes['type']
        if type and type.value != 'online information':
            text = self._findText(node.childNodes)
            if text is not None:
                if 'comments' in dict:
                    dict['comments'] = dict['comments'] + '\n' + type.value + ': ' + text
                else:
                    dict['comments'] = type.value + ': ' + text

    def _findText(self, nodeList):
        """Get text value from
           <text status="by similarity">Antiparallel disulfide-linked .... </text>
        """
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            if node.tagName == 'text':
                return node.firstChild.data
        return None

    def _FindIsoFormSeq(self, doc, vId, dict):
        """Get isoform sequence for vId if it exists  -  """
        if (self.__debug):
            self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - starting vId %s dict %r\n" % (vId, dict.items()))
        try:
            isoformdic = self._FindIsoFormIds(doc)

            if self.__debug:
                self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - vId %s isoformdic  %r\n" % (vId, isoformdic.items()))

            if not isoformdic:
                return False, False

            if vId not in isoformdic:
                return False, False

            # JDW 12-DEC-2015  Remove this test for a reference - Finding isoform data in comments lacking this
            if isoformdic[vId]['type'] == 'displayed' or 'ref' not in isoformdic[vId]:
                # return True, False
                return True, False

            refdic = self._FindIsoFormRefs(doc)
            if (self.__debug):
                self.__lfh.write("+ReadUnpXML._FindIsoFormSeq - vId %s refdic  %r\n" % (vId, refdic.items()))

            if not refdic:
                # return with sequence updated = False
                return True, False

            reflist = isoformdic[vId]['ref'].split(' ')
            # Reverse the ref list order so that sequence manipulation starts from C-terminal
            reflist.reverse()
            for ref in reflist:
                if ref in refdic:
                    dict['sequence'] = self._ProcessIsoFormSeq(dict['sequence'], refdic[ref])
            # return with seqquence updated = True
            return True, True
        except:
            self.__lfh.write("_findIsoFormSeq() failing -- %s\n" % vId)
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
        entryList = doc.getElementsByTagName('isoform')
        if not entryList:
            return dic

        for node in entryList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            #id = None
            id = []
            type = None
            ref = None
            for node1 in node.childNodes:
                if node1.nodeType != node1.ELEMENT_NODE:
                    continue
                if node1.tagName == 'id':
                    id.append(node1.firstChild.data)

                elif node1.tagName == 'sequence':
                    type = node1.attributes['type'].value
                    # JDW Aug-26 The following is behaving badly  --
                    try:
                        if 'ref' in node1.attributes:
                            ref = node1.attributes['ref'].value
                    except:
                        pass

            if len(id) < 1 or not type:
                continue
            d = {}
            d['type'] = type
            if ref:
                d['ref'] = ref
            dic[id[0]] = d

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
        entryList = doc.getElementsByTagName('feature')
        if not entryList:
            return dic

        for node in entryList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.attributes['type'].value != 'splice variant':
                continue

            if 'id' not in node.attributes:
                continue

            id = node.attributes['id'].value
            begin = None
            end = None
            variation = None
            for node1 in node.childNodes:
                if node1.nodeType != node1.ELEMENT_NODE:
                    continue
                if node1.tagName == 'variation':
                    variation = node1.firstChild.data
                elif node1.tagName == 'location':
                    for node2 in node1.childNodes:
                        if node2.nodeType != node2.ELEMENT_NODE:
                            continue
                        if node2.tagName == 'begin':
                            begin = node2.attributes['position'].value
                        elif node2.tagName == 'end':
                            end = node2.attributes['position'].value

            if not begin or not end:
                continue
            d = {}
            d['begin'] = begin
            d['end'] = end
            if variation:
                d['variation'] = variation
            dic[id] = d

        return dic

    def _ProcessIsoFormSeq(self, seq, ref):
        """Manipulate sequence using information from dictionary ref:

               { 'begin' : '455', 'end' : '518', 'variation' : 'DEAGQNEGG....' }
        """
        begin = int(ref['begin']) - 1
        end = int(ref['end'])
        seq1 = seq[0:begin]
        if 'variation' in ref:
            seq1 += ref['variation']
        seq1 += seq[end:]
        return seq1


class ReadUnpXmlFile(ReadUnpXml):

    """
    """

    def __init__(self, fileName):
        self._fileName = fileName
        self._doc = minidom.parse(self._fileName)
        ReadUnpXml.__init__(self, self._doc)


class ReadUnpXmlString(ReadUnpXml):

    """
    """

    def __init__(self, data, verbose=False, log=sys.stderr):
        self._data = data
        self._doc = minidom.parseString(self._data)
        ReadUnpXml.__init__(self, self._doc, verbose=verbose, log=log)


def main(argv):
    opts, args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadUnpXmlFile(arg)
            dict = obj.GetResult()
            for (k, v) in dict.items():
                print("%s=%s" % (k, v))

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)