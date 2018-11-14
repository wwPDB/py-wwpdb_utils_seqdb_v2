##
# File:    ReadUnpBlastXml.py
# Author:  Zukang Feng
# Update:  24-May-2010
# Version: 001  Initial version
#
# Updates:
#
# 03-Jan-2013 jdw adopt new FetchUnpXml class and fetch protocol.
# 05-Apr-2013 jdw add description and length details from the hit attributes
# 07-Apr-2013 jdw return isoform variant as a separate item.
# 19-Apr-2013 jdw restore db_code from the full uniprot entry. 
##

__author__  = "Zukang Feng"
__email__   = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import sys, math, getopt
from xml.dom import minidom
from wwpdb.utils.seqdb_v2.FetchUnpXml import FetchUnpXml

class ReadUnpBlastXml:
    """Read Uniprot blast result xml file and return the list of hits which contains:

           dict['db_name']
           dict['db_code']
           dict['db_accession']
           dict['db_length']
           dict['db_description']
           dict['db_isoform']
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
           dict['sequence']
           dict['ec']
           dict['keyword']
           dict['name']
           dict['synonyms']
           dict['gene']
           dict['source_scientific']
           dict['source_common']
           dict['taxonomy_id']
           dict['comments']
    """

    def __init__(self, doc, verbose=True, log=sys.stderr):
        self.__verbose=verbose
        self.__lfh=log
        self._result = self._parse(doc)

    def GetResult(self):
        return self._result

    def _parse(self, doc):
        length = None
        for node in doc.getElementsByTagName('sequence'):
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName == 'sequence' and node.attributes['length']:
                length = node.attributes['length'].value
                break

        mapping = {}
        mapping['database']    = 'db_name'
        mapping['id']          = 'db_code'
        mapping['ac']          = 'db_accession'
        mapping['description'] = 'db_description'
        mapping['length']      = 'db_length'
        #
        resultlist = []
        entry = None
        entryList = doc.getElementsByTagName('hits')
        if len(entryList) > 0:
            entry = entryList[0]
        else:
            return resultlist

        total = entry.attributes['total']
        if total:
            if total.value == '0':
                return resultlist

        fetchList = []
        for node in entry.childNodes:
            if node.nodeType != node.ELEMENT_NODE:
                continue

            if node.tagName != 'hit':
                continue


            alignlist = self._ProcessAlignmentsTag(node.childNodes, length)
            if alignlist:
                for align in alignlist:
                    for (k, v) in mapping.items():
                        align[v] = node.attributes[k].value
                    if length:
                        align['query_length'] = length
                    #
                    isoForm=''
                    if (align.has_key('db_code') and (align['db_code'].find('-') != -1)):
                        tL=align['db_code'].split('-')
                        if len(tL)>1 and len(tL[1])>0:
                            isoForm=str(tL[1]).strip()
                            dbCode=str(tL[0]).strip()
                            align['db_code']=dbCode
                    if (align.has_key('db_accession') and (align['db_accession'].find('-') != -1)):
                        tL=align['db_accession'].split('-')
                        if len(tL)>1 and len(tL[1])>0:
                            isoForm=str(tL[1]).strip()
                            dbAcc=str(tL[0]).strip()
                            align['db_accession']=dbAcc

                    align['db_isoform']=isoForm
                    # Get uniprot sequence information

                    if align.has_key('db_accession'):
                        fetchList.append(align['db_accession'])
                        

                    resultlist.append(align)
        #
        # incorporate additional details about each matching sequence entry.
        #  --  DO NOT OVERWRITE EXISTING ANNOTATION other than db_code -- 
        #
        if len(fetchList) > 0:
            fobj = FetchUnpXml(verbose=self.__verbose,log=self.__lfh)
            fobj.fetchList(fetchList)
            eD = fobj.getResult()
            for align in resultlist:
                if align.has_key('db_accession'):
                    acId=align['db_accession']
                    if eD.has_key(acId):
                        dict=eD[acId]
                        for (k, v) in dict.items():
                            if (k in ['db_code']):
                                align[k] = v                            
                            elif not align.has_key(k):
                                align[k] = v


        return resultlist

    def _ProcessAlignmentsTag(self, nodelist, length):
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            if node.tagName != 'alignments':
                continue
            total = node.attributes['total']
            if total:
                if total.value == '0':
                    continue
    
            return self._ProcessAlignmentTag(node.childNodes, length)

    def _ProcessAlignmentTag(self, nodelist, length):
        resultlist = []
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
            if node.tagName != 'alignment':
                continue

            dict = self._GetMatchAlignment(node.childNodes, length)
            if dict:
                resultlist.append(dict)
        return resultlist

    def _GetMatchAlignment(self, nodelist, query_length):
        dict = {}
        for node in nodelist:
            if node.nodeType != node.ELEMENT_NODE:
                continue
    
            if node.tagName == 'identity':
                dict['identity'] = node.firstChild.data
            elif node.tagName == 'positives':
                dict['positive'] = node.firstChild.data
            elif node.tagName == 'gaps':
                dict['gaps'] = node.firstChild.data 
            elif node.tagName == 'pattern':
                dict['midline'] = node.firstChild.data
            elif node.tagName == 'querySeq':
                dict['queryFrom'] = node.attributes['start'].value
                dict['queryTo'] = node.attributes['end'].value
                dict['query'] = node.firstChild.data
            elif node.tagName == 'matchSeq':
                dict['hitFrom'] = node.attributes['start'].value
                dict['hitTo'] = node.attributes['end'].value
                dict['subject'] = node.firstChild.data

        if float(dict['identity']) < 70:
            dict.clear()
            return dict

        if dict:
            if dict.has_key('query'):
                dict['alignLen'] = str(len(dict['query']))
            if dict.has_key('queryFrom') and dict.has_key('queryTo'):
                length = int(dict['queryTo']) - int(dict['queryFrom']) + 1
                dict['match_length'] = str(length)
                if dict.has_key('identity'):
                    identity = int(math.ceil(float(length) * float(dict['identity']) / 100.0))
                    if query_length:
                        percent = identity * 100 / int(query_length)
                        if percent < 70:
                            dict.clear()
                            return dict
                    dict['identity'] = str(identity)
                if dict.has_key('positive'):
                    positive = int(math.ceil(float(length) * float(dict['positive']) / 100.0))
                    dict['positive'] = str(positive)
                if dict.has_key('gaps'):
                    gaps = int(math.ceil(float(length) * float(dict['gaps']) / 100.0))
                    dict['gaps'] = str(gaps)
            if not dict.has_key('gaps'):
                dict['gaps'] = '0'
        return dict

class ReadUnpBlastXmlFile(ReadUnpBlastXml):
    """
    """
    def __init__(self, fileName,verbose=True,log=sys.stderr):
        self._fileName = fileName
        self._doc = minidom.parse(self._fileName)
        ReadUnpBlastXml.__init__(self, self._doc,verbose=verbose,log=log)

class ReadUnpBlastXmlString(ReadUnpBlastXml):
    """
    """
    def __init__(self, data,verbose=True,log=sys.stderr):
        self._data = data
        self._doc = minidom.parseString(self._data)
        ReadUnpBlastXml.__init__(self, self._doc, verbose=verbose,log=log)


def main(argv):
    opts, args = getopt.getopt(argv, "x:", ["xml="])
    for opt, arg in opts:
        if opt in ("-x", "--xml"):
            obj = ReadUnpBlastXmlFile(arg)
            result = obj.GetResult()
            for align in result:
                for (k, v) in align.items():
                    sys.stdout.write('%s=%s\n' % (k, v))

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        sys.stderr.write( exc )
        sys.exit(1)
