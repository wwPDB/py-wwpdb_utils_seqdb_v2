##
# File: BlastProcess.py
# Date: 25-April-2010
# Author:  Zukang Feng
#
# Updates:
# 27-Apr-2010 jdw add method to run and save by entity to a named output file.
# 31-Dec-2012 jdw migrate taxonomy extensions to V2 and consolidate batch fetch of UniProt entries.
#  2-Jan-2013 jdw  add additional fragment tracking items -
##

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import os
import sys
import time
import getopt
import re
import gzip
import traceback
from wwpdb.utils.seqdb_v2.UnpBlastService import UnpBlastService
from wwpdb.utils.seqdb_v2.ReadUnpBlastXml import ReadUnpBlastXmlString
from wwpdb.utils.seqdb_v2.NcbiBlastService import NcbiBlastService
from wwpdb.utils.seqdb_v2.ReadNcbiBlastXml import ReadNcbiBlastXmlString
from wwpdb.utils.seqdb_v2.mmCIFUtil import mmCIFUtil
from mmcif.api.PdbxContainers import DataContainer
from mmcif.api.DataCategory import DataCategory
from mmcif.io.PdbxWriter import PdbxWriter

#


class RunBlastPerSeq(object):
    """Run blast search and save result into cif file"""

    def __init__(self, entityId=None, entityInfo=None, taxonomyData=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__entity = entityId
        self.__sequence = ""
        self.__type = ""
        self.__taxonomyData = taxonomyData
        self.__fragments = []
        self.__result = []
        #
        self.__saveBlastResults = False
        self.__blastFileNamePrefix = None
        #
        self.__blastPath = "."
        #
        if entityInfo:
            self.__sequence = entityInfo["seq"]
            self.__type = entityInfo["type"]
            self.__fragmentList = entityInfo["fragments"]

    def saveBlastResults(self, blastPath=".", blastFileNamePrefix="test"):
        self.__saveBlastResults = True
        self.__blastPath = blastPath
        self.__blastFileNamePrefix = blastFileNamePrefix

    def GetResult(self):
        return self.__result

    def Run(self):
        """Manage the reference sequence database search for each fragment of the input
        sequence.

        Return sorted results using taxonomy data if this is avalable.
        """
        if (not self.__sequence) or (not self.__fragmentList):
            return

        for (fragId, fragment) in enumerate(self.__fragmentList):
            indFrag = fragId + 1
            start = int(fragment["beg"]) - 1
            end = int(fragment["end"])
            sequence = self.__sequence[start:end]

            if self.__saveBlastResults:
                bPath = os.path.join(self.__blastPath, self.__blastFileNamePrefix + "_seqdb-blast_E" + self.__entity + "_F " + str(indFrag) + ".xml")
            else:
                bPath = None

            result = self._Run(sequence, self.__type, blastFilePath=bPath)
            if not result:
                continue
            #
            for match in result:
                match["beg_seq_num"] = fragment["beg"]
                match["end_seq_num"] = fragment["end"]
                match["fragment_id"] = indFrag

            #
            # re-sorting result
            #
            taxids = {}
            if "taxid" in fragment:
                taxids = self._GetTaxonomyTree(fragment["taxid"])
            #
            auth_accession_code = None
            if "accession" in fragment:
                auth_accession_code = fragment["accession"]
            self._SortResult(result, taxids, auth_accession_code)

            self.__result.append(result)

    def _Run(self, sequence, rtype, blastFilePath=None):
        """Internal method to execute the sequence search according to input polymer type."""
        if self.__verbose:
            self.__lfh.write("+INFO (RunBlastPerSeq._Run) Launch search for sequence      = %s\n" % sequence)
            self.__lfh.write("+INFO (RunBlastPerSeq._Run) Launch search for sequence type = %s\n" % rtype)
        #
        timeBegin = time.time()
        blast_match_result = []

        # run uniprot blast service for protein sequence
        if rtype == "polypeptide":
            service = UnpBlastService(sequence)
            service.RunService()
            # fetch the raw XML result from the Blast search
            xmlResult = service.GetResult()
            if blastFilePath is not None:
                service.WriteResultFile(blastFilePath)
            #
            timeBlast = time.time()
            if self.__verbose:
                self.__lfh.write("+INFO (RunBlastPerSeq._Run) Blast search completed  in %d seconds\n" % (timeBlast - timeBegin))
            #
            if xmlResult:
                blastresult = ReadUnpBlastXmlString(xmlResult, verbose=self.__verbose, log=self.__lfh)
                blast_match_result = blastresult.GetResult()

        # run ncbi blast service for RNA sequence
        elif rtype == "polyribonucleotide":
            service = NcbiBlastService(sequence)
            service.RunService()
            # fetch the raw XML result from the Blast search
            xmlResult = service.GetResult()
            if blastFilePath is not None:
                service.WriteResultFile(blastFilePath)
            timeBlast = time.time()
            if self.__verbose:
                self.__lfh.write("+INFO (RunBlastPerSeq._Run) Blast search completed  in %d seconds\n" % (timeBlast - timeBegin))
            if xmlResult:
                blastresult = ReadNcbiBlastXmlString(xmlResult)
                blast_match_result = blastresult.GetResult()

        else:
            self.__lfh.write("+INFO (RunBlastPerSeq._Run) Search failed for unknown type    = %s\n" % rtype)

        timeProc = time.time()
        if self.__verbose:
            if self.__verbose:
                self.__lfh.write("+INFO (RunBlastPerSeq._Run) Search processing completed in %d seconds\n" % (timeProc - timeBegin))

            self.__lfh.write("+INFO (RunBlastPerSeq._Run) Result length       = %d\n" % len(blast_match_result))

        return blast_match_result

    def _GetTaxonomyTree(self, taxid):
        rdict = {}
        rdict["id"] = taxid
        if self.__taxonomyData:
            if taxid in self.__taxonomyData:
                parent_id = self.__taxonomyData[taxid]
                rdict["p_id"] = parent_id
                if parent_id in self.__taxonomyData:
                    gparent_id = self.__taxonomyData[parent_id]
                    rdict["gp_id"] = gparent_id
                #
            #
        #
        return rdict

    def _SortResult(self, result, taxids, auth_accession_code):
        """Add a sorting index to the dictionary of sequence correspondence results obtained from
        the Blast search.   The sort index is based on heurestic which includes sequence matching
        metrics and taxonomy data.

        """
        if not result:
            return

        sorting_list = []
        for i in range(0, len(result)):
            identity = (int(result[i]["identity"]) - int(result[i]["gaps"])) * 100 / int(result[i]["query_length"])
            if result[i]["db_name"] == "SP":
                identity += 2
            if auth_accession_code and (result[i]["db_code"] == auth_accession_code or result[i]["db_accession"] == auth_accession_code):
                identity += 1
            #
            target_taxids = {}
            if "taxonomy_id" in result[i]:
                target_taxids = self._GetTaxonomyTree(result[i]["taxonomy_id"])
            #
            taxid_match = 0
            if taxids and target_taxids:
                if taxids["id"] == target_taxids["id"]:
                    taxid_match = 3
                elif "p_id" in taxids and taxids["p_id"] == target_taxids["id"]:
                    taxid_match = 2
                elif "p_id" in target_taxids and target_taxids["p_id"] == taxids["id"]:
                    taxid_match = 2
                elif "gp_id" in taxids and taxids["gp_id"] == target_taxids["id"]:
                    taxid_match = 1
                elif "gp_id" in target_taxids and target_taxids["gp_id"] == taxids["id"]:
                    taxid_match = 1
            #
            identity = identity * 4 + taxid_match
            ilist = [identity, i]
            sorting_list.append(ilist)
            #
        #
        sorting_list.sort(key=lambda val: val[0], reverse=True)
        #
        for i in range(0, len(sorting_list)):
            result[sorting_list[i][1]]["sort_order"] = str(i + 1)
        #

    def WriteCifFile(self, EntryID, fileName):
        """Export the reference sequence matching data to a CIF data file."""
        if not self.__result:
            return

        curContainer = DataContainer("match_entity")
        t = DataCategory("info")
        t.appendAttribute("struct_id")
        t.appendAttribute("entity_id")
        t.appendAttribute("sequence")
        t.appendAttribute("fragment_count")
        #
        t.setValue(EntryID, "struct_id", 0)
        t.setValue(self.__entity, "entity_id", 0)
        t.setValue(self._FormatSequence(self.__sequence), "sequence", 0)
        t.setValue(len(self.__fragmentList), "fragment_count", 0)

        curContainer.append(t)

        _table_items = []
        _table_items.append("id")
        _table_items.append("beg_seq_num")
        _table_items.append("end_seq_num")
        _table_items.append("sort_order")
        _table_items.append("db_name")
        _table_items.append("db_code")
        _table_items.append("db_accession")
        _table_items.append("match_length")
        _table_items.append("queryFrom")
        _table_items.append("queryTo")
        _table_items.append("hitFrom")
        _table_items.append("hitTo")
        _table_items.append("identity")
        _table_items.append("positive")
        _table_items.append("gaps")
        _table_items.append("alignLen")
        _table_items.append("query")
        _table_items.append("subject")
        _table_items.append("midline")
        _table_items.append("query_length")
        _table_items.append("name")
        _table_items.append("source_scientific")
        _table_items.append("source_common")
        _table_items.append("taxonomy_id")
        _table_items.append("gene")
        _table_items.append("synonyms")
        _table_items.append("comments")
        _table_items.append("keyword")
        _table_items.append("ec")
        _table_items.append("fragment_id")

        t = DataCategory("match_entity")
        for item in _table_items:
            t.appendAttribute(item)

        t1 = DataCategory("org_sequence")
        t1.appendAttribute("id")
        t1.appendAttribute("sequence")

        row = 0
        for fragment_result in self.__result:
            for match in fragment_result:
                t.setValue(str(row + 1), "id", row)
                for item in _table_items:
                    if item in match:
                        t.setValue(str(match[item]), item, row)

                t1.setValue(str(row + 1), "id", row)
                if "sequence" in match:
                    seq = re.sub("[\t \n]", "", match["sequence"])
                    t1.setValue(self._FormatSequence(str(seq)), "sequence", row)

                row = row + 1
            #
        #

        curContainer.append(t)
        curContainer.append(t1)

        myDataList = []
        myDataList.append(curContainer)
        ofh = open(fileName, "w")
        pdbxW = PdbxWriter(ofh)
        pdbxW.write(myDataList)
        ofh.close()

    def _FormatSequence(self, sequence):
        num_per_line = 60
        lines = int(len(sequence) / num_per_line)
        x = len(sequence) % num_per_line
        m = lines
        if x:
            m = lines + 1

        seq = ""
        for i in range(m):
            n = num_per_line
            if i == lines:
                n = x
            seq += sequence[i * num_per_line : i * num_per_line + n]
            if i != (m - 1):
                seq += "\n"

        return seq


class BlastProcess(object):
    """Search reference sequence database for sequence correspondences for each polymer entity.

    Obtain sequence search target from the the one-letter-code sequences stored in category
    entity_poly.

    Launch the Blast search for each entity and store processed results in CIF format data
    file for each entity.

    """

    def __init__(self, cifFilePath=None, taxonomyFilePath=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        #
        self.__saveBlastResults = False
        self.__blastFileNamePrefix = None
        #
        self.__filePath = cifFilePath
        self.__taxonomyFilePath = taxonomyFilePath
        self.__taxonomyData = {}
        #
        self.__blastPath = "."
        #
        # Read input file and create dictionary of polymer sequence details
        self.__cifObj = mmCIFUtil(filePath=self.__filePath)
        self.__BlockID = self.__cifObj.GetBlockID()
        self.__entitySeq = self._GetEntitySeq()
        #
        self.__taxFlag = True
        if self.__taxFlag:
            self.__taxonomyData = self._readTaxonomyData()
            if self.__verbose:
                self.__lfh.write("+INFO BlastProcess() taxonomy file path   %s\n" % self.__taxonomyFilePath)
                self.__lfh.write("+INFO BlastProcess() taxonomy data length %d\n" % len(self.__taxonomyData))

        #
        if self.__verbose:
            self.__lfh.write("+INFO BlastProcess() cifFilePath      %s\n" % self.__filePath)
            self.__lfh.write("+INFO BlastProcess() taxonomyFilePath %s\n" % self.__taxonomyFilePath)
            self.__lfh.write("+INFO BlastProcess() entity sequence dictionary  %r\n" % self.__entitySeq)

    def saveBlastResults(self, blastPath=".", blastFileNamePrefix="test"):
        self.__saveBlastResults = True
        self.__blastPath = blastPath
        self.__blastFileNamePrefix = blastFileNamePrefix

    def getPolymerEntityDetails(self):
        """Return a referenc to the internal dictionary containing polymer entity details."""
        return self.__entitySeq

    def _GetEntitySeq(self):
        """Read the entity_poly category and

        Return a dictionary containing for each polymer entity:

        ['seq']       = one-letter-code sequence
        ['type']      = sequence type --
        ['fragments'] = list of dictionaries for each sequece fragment  containing keys - beg (one based), end, & tax_id

        """
        rdict = {}

        dList = self.__cifObj.GetValue("entity_poly")
        if not dList:
            return rdict

        for d in dList:
            if "entity_id" not in d:
                continue

            sequence = ""
            if "pdbx_seq_one_letter_code" in d:
                sequence = d["pdbx_seq_one_letter_code"]
            elif "ndb_seq_one_letter_code" in d:
                sequence = d["ndb_seq_one_letter_code"]

            sequence = self._CleanSequence(sequence)
            if len(sequence) < 5:
                continue

            rtype = None
            if "type" in d:
                rtype = d["type"]
                rtype = self._FindCorrectType(rtype)

            if not rtype:
                rtype = self._FindTypeFromSeq(sequence)

            if rtype not in ["polyribonucleotide", "polypeptide"]:
                continue

            if rtype == "polyribonucleotide" and len(sequence) < 50:
                continue

            rdict[d["entity_id"]] = {}
            rdict[d["entity_id"]]["seq"] = sequence
            rdict[d["entity_id"]]["type"] = rtype
            fragmentList = self._GetEntityProperty(d["entity_id"], len(sequence))
            if fragmentList:
                rdict[d["entity_id"]]["fragments"] = fragmentList

        return rdict

    def _GetEntityProperty(self, entityid, length):
        pList = self._GetTaxID(entityid, length, "entity_src_gen", "pdbx_gene_src_ncbi_taxonomy_id")
        if not pList:
            pList = self._GetTaxID(entityid, length, "entity_src_nat", "pdbx_ncbi_taxonomy_id")
        if not pList:
            pList = self._GetTaxID(entityid, length, "pdbx_entity_src_syn", "ncbi_taxonomy_id")

        if not pList:
            rdict = {}
            rdict["beg"] = "1"
            rdict["end"] = str(length)
            pList.append(rdict)

        return pList

    def _GetTaxID(self, entityid, length, category, item):
        pList = []

        dList = self.__cifObj.GetValue(category)
        if not dList:
            return pList

        for d in dList:
            if "entity_id" not in d:
                continue

            if d["entity_id"] != entityid:
                continue

            rdict = {}
            if item in d:
                rdict["taxid"] = d[item]

            if "pdbx_beg_seq_num" in d and "pdbx_end_seq_num" in d:
                rdict["beg"] = d["pdbx_beg_seq_num"]
                rdict["end"] = d["pdbx_end_seq_num"]
            else:
                rdict["beg"] = "1"
                rdict["end"] = str(length)

            pList.append(rdict)
        #

        if pList:
            if not self._checkFragments(length, pList):
                pList = pList[0:1]
                pList[0]["beg"] = "1"
                pList[0]["end"] = str(length)

            for fragment in pList:
                if "taxid" in fragment:
                    self.__taxFlag = True
                    break
                #
            #
        #

        return pList

    def _checkFragments(self, length, pList):
        if not pList:
            return True

        # check start number in first fragment >= 1
        start = int(pList[0]["beg"])
        if start < 1:
            return False

        # check end number in last fragment <= length of whole sequence
        n = len(pList)
        iend = int(pList[n - 1]["end"])
        if iend > length:
            return False

        # check start number < end number in same fragment
        for fragment in pList:
            start = int(fragment["beg"])
            end = int(fragment["end"])
            if start >= end:
                return False

        if n < 2:
            return True

        # check start number in current fragment > end number in previous fragment
        for i in range(1, n):
            end = int(pList[i - 1]["end"])
            start = int(pList[i]["beg"])
            if end >= start:
                return False

        return True

    def _readTaxonomyData(self):
        """Read the NCBI taxonomy data file 'nodes.dmp' and
        return a dictionary of parent taxonomy id's.

        d[tax_id]=parent_tax_id

        """
        if not self.__taxonomyFilePath:
            return

        d = {}
        try:
            if self.__taxonomyFilePath[-3:] == ".gz":
                f = gzip.open(self.__taxonomyFilePath, "r")
            else:
                f = open(self.__taxonomyFilePath, "r")
            data = f.read()
            f.close()
            #
            # Python3 coming in as bytes from compressed file
            tlist = data.decode("utf-8").split("\t|\n")
            for line in tlist:
                if not line:
                    continue
                #
                list1 = line.split("\t|\t")
                d[list1[0]] = list1[1]
        except Exception as e:
            traceback.print_exc(file=self.__lfh)
            self.__lfh.write("+ERROR BlastProcess() failed to read taxonomy data file %s %s\n" % (str(e), self.__taxonomyFilePath))
        return d

    def _CleanSequence(self, seq):
        seq = seq.upper()
        seq = re.sub("[\t \n]", "", seq)
        seq = re.sub("\\(MSE\\)", "M", seq)
        seq = re.sub("\\([A-Z]{2,3}\\)", "X", seq)
        return seq

    def _FindCorrectType(self, ptype):
        ptype = ptype.lower()
        ptype = ptype.strip()
        if not ptype:
            return ptype
        if ptype.find("polydeoxyribonucleotide") >= 0 or ptype.find("dna") >= 0:
            return "polydeoxyribonucleotide"
        elif ptype == "polyribonucleotide" or ptype == "rna":
            return "polyribonucleotide"
        elif ptype.find("polypeptide") >= 0 or ptype == "protein":
            return "polypeptide"
        elif ptype == "?" or ptype == ".":
            return ""
        else:
            return ptype

    def _FindTypeFromSeq(self, seq):
        if re.search("[DEFHIKLMNPQRSVWY]", seq):
            return "polypeptide"
        elif re.search("[ACGU]", seq):
            return "polyribonucleotide"
        else:
            return ""

    def Run(self, entityId=None):
        """Execute the reference sequence search and return processed correspondence results -

        If entityId is specified then results are returned in a dictionary for each polymer
        entity.

        Otherwise, the all polymer entities are processed and the results are written
        to CIF data files for each entity.

        The output file name convention is <Entry_id>.<entity_id>.info.cif

        """
        hitlist = {}
        if not self.__entitySeq:
            return hitlist

        if entityId:
            if entityId in self.__entitySeq:
                perseq = RunBlastPerSeq(entityId=entityId, entityInfo=self.__entitySeq[entityId], taxonomyData=self.__taxonomyData, verbose=self.__verbose, log=self.__lfh)
                if self.__saveBlastResults:
                    perseq.saveBlastResults(blastPath=self.__blastPath, blastFileNamePrefix=self.__blastFileNamePrefix)
                perseq.Run()
                hitlist = perseq.GetResult()
            return hitlist

        for (k, v) in self.__entitySeq.items():
            perseq = RunBlastPerSeq(entityId=k, entityInfo=v, taxonomyData=self.__taxonomyData, verbose=self.__verbose, log=self.__lfh)
            if self.__saveBlastResults:
                perseq.saveBlastResults(blastPath=self.__blastPath, blastFileNamePrefix=self.__blastFileNamePrefix)

            perseq.Run()
            perseq.WriteCifFile(self.__BlockID, self.__BlockID + "." + k + ".info.cif")
        return hitlist

    def RunAndSave(self, entityId=None, fName=None):
        """Execute the reference sequence search and return processed correspondence results in
        the named file.

        This is the  entry point for wwDPB API plugin.

        """
        if entityId is None or fName is None:
            return False

        eId = str(entityId)
        if eId not in self.__entitySeq:
            return False

        perseq = RunBlastPerSeq(entityId=eId, entityInfo=self.__entitySeq[eId], taxonomyData=self.__taxonomyData, verbose=self.__verbose, log=self.__lfh)
        if self.__saveBlastResults:
            perseq.saveBlastResults(blastPath=self.__blastPath, blastFileNamePrefix=self.__blastFileNamePrefix)

        perseq.Run()
        perseq.WriteCifFile(self.__BlockID, str(fName))
        return True


def main(argv):  # pragma: no cover

    opts, _args = getopt.getopt(argv, "i:e:t:", ["ciffile=", "entity=", "taxonomy="])

    ciffile = None
    entity = None
    taxfile = None
    for opt, arg in opts:
        if opt in ("-i", "--ciffile"):
            ciffile = arg
        elif opt in ("-e", "--entity"):
            entity = arg
        elif opt in ("-t", "--taxonomy"):
            taxfile = arg

    if ciffile:
        try:
            process = BlastProcess(cifFilePath=ciffile, taxonomyFilePath=taxfile)
            resultlist = process.Run(entityId=entity)
            for ii, rst in enumerate(resultlist):
                sys.stdout.write("%d: %r\n" % (ii, rst))

        except Exception as exc:
            sys.stderr.write(exc)
            traceback.print_exc(file=sys.stderr)


if __name__ == "__main__":  # pragma: no cover
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        sys.stderr.write(exc)
        sys.exit(1)
