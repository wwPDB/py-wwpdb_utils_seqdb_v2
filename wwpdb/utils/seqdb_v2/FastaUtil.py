##
# File: FastaUtil.py
# Date:  30-Jan-2014  Jdw
#
import sys
from gzip import GzipFile


class FastaUtil(object):
    """Simple FASTA reader and writer with methods adjusted  UniProt Fasta files -"""

    def __init__(self, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__lfh = log
        self.__debug = False

    def loadFastaUniProt(self, fastaFilePath):
        """Read the the input FASTA file and return a dictionary sequence using the
        first white-space delimited token on each comment line as dicitonary key.

        A list of sequence id's is also returned.
        """
        seqDict = {}
        seqIdList = []

        if fastaFilePath[-3:] == ".gz":
            ifh = GzipFile(fastaFilePath)
        else:
            ifh = open(fastaFilePath, "r")

        for cmtLine, sequence in self.__read_record_fasta(ifh):
            try:
                org = ""
                geneName = ""
                dbName = ""
                dbAccession = ""
                dbIsoform = ""
                seqId = ""
                description = ""
                #
                ff = cmtLine[1:].split("|")
                dbName = ff[0].upper()
                seqId = ff[1]
                tS = ff[2]
                #
                tt = seqId.split("-")
                dbAccession = seqId
                dbIsoform = ""
                if len(tt) > 1:
                    dbAccession = tt[0]
                    dbIsoform = tt[1]

                #
                idx = tS.find("OS=")
                if idx > 0:
                    description = tS[: idx - 1]
                    ff = tS[idx:].split("=")
                    if len(ff) == 3:
                        if ff[0] == "OS":
                            org = ff[1][:-2]
                            geneName = ff[2]
                    elif len(ff) > 1:
                        if ff[0] == "OS":
                            org = ff[1]
            except Exception as e:
                if self.__verbose:
                    self.__lfh.write("+FastaUtil.loadFastaUniProt() Read failing for file %s at %s err %s\n" % (fastaFilePath, cmtLine, str(e)))
                ifh.close()
                return seqIdList, seqDict
                #
            if seqId not in seqDict:
                seqIdList.append(seqId)
                seqDict[seqId] = {
                    "sequence": sequence.upper(),
                    "description": description,
                    "org": org,
                    "gene_name": geneName,
                    "db_accession": dbAccession,
                    "db_name": dbName,
                    "db_isoform": dbIsoform,
                }
            else:
                if self.__debug:
                    self.__lfh.write("+FastaUtil.loadFastaUniProt() Duplicate sequence identifier %s\n%s\n" % (seqId, sequence))

        if self.__verbose:
            self.__lfh.write("+FastaUtil.loadFastaUniProt() Read %d sequences in %s\n" % (len(seqIdList), fastaFilePath))

        ifh.close()
        return seqIdList, seqDict

    def __read_record_fasta(self, ifh):
        """Returen the next FASTA record in the input file handle."""
        comment, sequence = None, []
        for line in ifh:
            line = line.rstrip()
            if line.startswith(">"):
                if comment:
                    yield (comment, "".join(sequence))
                comment, sequence = line, []
            else:
                sequence.append(line)
        if comment:
            yield (comment, "".join(sequence))

    def __writeFasta(self, ofh, seqDict, seqIdList, maxLineLength=70):
        """"""
        lCount = 0
        for seqId in seqIdList:
            sA = []
            sequence = seqDict[seqId]
            sA.append(">%s\n" % seqId)
            lCount += 1
            for i in range(0, len(sequence), maxLineLength):
                sA.append(sequence[i : i + maxLineLength] + "\n")
                result = "".join(sA)
                lCount += 1
            ofh.write(result)
        self.__lfh.write("+FastaUtil.__writeFasta() Formatted %d sequences in %d text lines\n" % (len(seqIdList), lCount))

    def __formatFasta(self, seqDict, seqIdList, maxLineLength=70):
        """"""
        sA = []
        lCount = 0
        for seqId in seqIdList:
            sequence = seqDict[seqId]
            sA.append(">%s\n" % seqId)
            lCount += 1
            for i in range(0, len(sequence), maxLineLength):
                sA.append(sequence[i : i + maxLineLength] + "\n")
                result = "".join(sA)
                lCount += 1
        self.__lfh.write("+FastaUtil.__formatFasta() Formatted %d sequences in %d text lines\n" % (len(seqIdList), lCount))
        return result


if __name__ == "__main__":  # pragma: no cover
    fa = FastaUtil(verbose=True, log=sys.stderr)
    sList, sDict = fa.loadFastaUniProt("uniprot_sprot_varsplic.fasta")
