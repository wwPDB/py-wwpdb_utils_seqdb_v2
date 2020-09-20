"""
File:    mmCIFUtil.py
Author:  Zukang Feng
Update:  21-August-2012
Version: 001  Initial version

"""

__author__ = "Zukang Feng"
__email__ = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import sys
from mmcif_utils.persist.PdbxPyIoAdapter import PdbxPyIoAdapter as PdbxIoAdapter


class mmCIFUtil(object):
    """Using pdbx mmCIF utility to parse mmCIF file"""

    def __init__(self, verbose=False, log=sys.stderr, filePath=None):
        self.__verbose = verbose
        self.__lfh = log
        self.__filePath = filePath
        self.__container = None
        self.__blockID = None
        self.__read()
        #

    def __read(self):
        myReader = PdbxIoAdapter(self.__verbose, self.__lfh)
        _ok = myReader.read(pdbxFilePath=self.__filePath)  # noqa: F841
        containerNameList = myReader.getContainerNameList()
        if not containerNameList:
            return
        #
        self.__blockID = containerNameList[0]
        self.__container = myReader.getContainer(containerNameList[0])
        #

    def GetBlockID(self):
        return self.__blockID
        #

    def GetValue(self, catName):
        """Get category values based on category name 'catName'. The results are stored
        in a list of dictionaries with item name as key
        """
        dList = []
        if not self.__container:
            return dList
        #
        catObj = self.__container.getObj(catName)
        if not catObj:
            return dList
        #
        # Get column name index
        #
        itNameList = catObj.getItemNameList()
        #
        rowList = catObj.getRowList()
        for row in rowList:
            tD = {}
            for idxIt, itName in enumerate(itNameList):
                if row[idxIt] != "?" and row[idxIt] != ".":
                    tlist = itName.split(".")
                    tD[tlist[1]] = row[idxIt]
            #
            if tD:
                dList.append(tD)
        #
        return dList
        #
