"""
File:    ReadNcbiXml.py
Author:  Zukang Feng
Update:  24-May-2010
Version: 001  Initial version

Updates -- 

    20-Apr-2013  jdw  -  provide options to avoid saving full sequence.
"""

__author__  = "Zukang Feng"
__email__   = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

from xml.dom import minidom
import sys
import getopt

class ReadNcbiXml:
    """Parse GeneBank entry xml file and put seqeunce, source organism, and
       taxonomy information into dictionary:
       dict['sequence']
       dict['source_scientific']
       dict['taxonomy_id']
    """
    def __init__(self, doc):
        self._dict = {}
        self.__skipSequence=False
        if doc:
            self._dict = self._parse(doc)

    def skipSequence(self):
        self.__skipSequence=True

    def GetResult(self):
        return self._dict

    def _parse(self, doc):
        dict = {}
        db_code=None
        entryList = doc.getElementsByTagName('GBSeq_accession-version')
        if entryList:
            for node in entryList:
                if node.nodeType != node.ELEMENT_NODE:
                    continue
                if node.tagName == 'GBSeq_accession-version':
                    db_code = node.firstChild.data
                    break

        # Get sequence:
        sequence = None
        if not self.__skipSequence:
            entryList = doc.getElementsByTagName('GBSeq_sequence')
            if entryList:
                for node in entryList:
                    if node.nodeType != node.ELEMENT_NODE:
                        continue
                    if node.tagName == 'GBSeq_sequence':
                        sequence = node.firstChild.data
                        break

            if sequence:
                sequence = self._ProcessRNASequence(sequence)

        entry = None
        entryList = doc.getElementsByTagName('GBFeature')
        if entryList:
            for node in entryList:
                if node.nodeType != node.ELEMENT_NODE:
                    continue

                dict = self._ProcessGBFeatureTag(node.childNodes)
                if dict:
                    break;

        if sequence:
            dict['sequence'] = sequence

        if db_code:
            dict['db_code'] = db_code

        return dict

    def _ProcessRNASequence(self, sequence):
        seq = sequence.upper()
        seq = seq.replace('T', 'U')
        return seq;

    def _ProcessGBFeatureTag(self, nodeList):
        """Get source information from
              <GBFeature>
                <GBFeature_key>source</GBFeature_key>
                ......
                <GBFeature_quals>
                  <GBQualifier>
                    <GBQualifier_name>organism</GBQualifier_name>
                    <GBQualifier_value>Thermus thermophilus HB8</GBQualifier_value>
                  </GBQualifier>
                  ......
                  <GBQualifier>
                    <GBQualifier_name>db_xref</GBQualifier_name>
                    <GBQualifier_value>taxon:300852</GBQualifier_value>
                  </GBQualifier>
                  ......
                </GBFeature_quals>
              </GBFeature>
        """
        dict = {}
        type = ''
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'GBFeature_key':
                type = node.firstChild.data
                if type != 'source':
                    dict.clear()
                    return dict;
            elif node.tagName == 'GBFeature_quals':
                dict = self._ProcessGBQualifierTag(node.childNodes)

        if type != 'source':
            dict.clear()

        return dict;

    def _ProcessGBQualifierTag(self, nodeList):
        """Get source_scientific from
                  <GBQualifier>
                    <GBQualifier_name>organism</GBQualifier_name>
                    <GBQualifier_value>Thermus thermophilus HB8</GBQualifier_value>
                  </GBQualifier>

           and taxonomy_id from
                  <GBQualifier>
                    <GBQualifier_name>db_xref</GBQualifier_name>
                    <GBQualifier_value>taxon:300852</GBQualifier_value>
                  </GBQualifier>
        """
        dict = {}
        for node in nodeList:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            name = None
            value = None
            for childnode in node.childNodes:
                if childnode.nodeType != node.ELEMENT_NODE:
                    continue

                if childnode.tagName == 'GBQualifier_name':
                    name = childnode.firstChild.data 
                elif childnode.tagName == 'GBQualifier_value':
                    value = childnode.firstChild.data

            if not name or not value:
                continue

            if name == 'organism':
                dict['source_scientific'] = value
            elif name == 'db_xref' and value.find('taxon:') >= 0:
                list = value.split(':')
                dict['taxonomy_id'] = list[1]

        return dict

class ReadNcbiXmlFile(ReadNcbiXml):
    """
    """
    def __init__(self, fileName): 
        self._fileName = fileName
        self._doc = ""
        try:
            self._doc = minidom.parse(self._fileName)
        except Exception, exc:
            pass
        ReadNcbiXml.__init__(self, self._doc)

class ReadNcbiXmlString(ReadNcbiXml):
    """
    """
    def __init__(self, data):
        self._data = data
        self._doc = ""
        try:
            self._doc = minidom.parseString(self._data)
        except Exception, exc:
            pass
        ReadNcbiXml.__init__(self, self._doc)

def main(argv):
    opts, args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadNcbiXmlFile(arg)
            dict = obj.GetResult()
            for (k, v) in dict.items():
                print "%s=%s" % (k, v)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception, exc:
        print exc
        sys.exit(1)
