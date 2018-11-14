##
# File:    UnpBlastService.py
# Author:  Zukang Feng
# Update:  24-May-2010
# Version: 001  Initial version
#
# Updates:
#
# Apr 13, 2103  jdw Adjust blast service parameters - 
#                   Add logging parameters
##

__author__  = "Zukang Feng"
__email__   = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

try:
    from urllib.parse import urlparse, urlencode
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError
except ImportError:
    from urlparse import urlparse
    from urllib import urlencode
    from urllib2 import urlopen, Request, HTTPError

import ssl
import os, sys, time
import getopt

class UnpBlastService:
    """ Utility class for running blastp service using uniprotkb database from Uniprot site
    """
    def __init__(self, sequence, verbose=False, log=sys.stderr):
        self._sequence = sequence;
        self._baseUrl = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast'
        self._result = None
        self._checkInterval = 10
        self.__verbose=verbose
        self.__lfh=log

    def RunService(self):
        """Main function for running blastp service"""
        params = {}
        params['email'] = 'zfeng@rcsb.rutgers.edu'
        params['sequence'] = self._sequence
        #params['gapalign'] = False
        params['program'] = 'blastp'
        params['database'] = 'uniprotkb'
        params['stype'] = 'protein'
        # Adjusted to match the defaults on the UniProt web -
        params['gapalign'] = True
        params['ext'] = 10
        # Submit the job
        jobid = self._serviceSubmit(params)
        if self.__verbose:
            self.__lfh.write("+UnpBlastService.RunService() Blast service job %s started at %s\n" % 
                             (jobid,time.strftime("%Y %m %d %H:%M:%S", time.localtime())))
        time.sleep(5)
        self._result = self._getResultfromServer(jobid)

    def GetResult(self):
        return self._result

    def WriteResultFile(self, filename):
        fh = open(filename, 'w')
        fh.write(self._result)
        fh.close()

    def _serviceSubmit(self, params):
        """Submit job"""
        try:
            requestUrl = self._baseUrl + '/run/'
            # Get the data for the options
            requestData = urlencode(params)
            gcontext = ssl._create_unverified_context()
            reqH = urlopen(requestUrl, requestData, context=gcontext)
            jobId = reqH.read()
            reqH.close()    
            return jobId
        except Exception as exc:
            return ''

    def _getResultfromServer(self, jobId):
        """Get result from server"""
        self._clientPoll(jobId)
        requestUrl = self._baseUrl + '/result/' + jobId + '/xml'
        return self._restRequest(requestUrl)

    def _clientPoll(self, jobId):
        """Client-side poll"""
        status = 'PENDING'
        while status == 'RUNNING' or status == 'PENDING':
            status = self._serviceGetStatus(jobId)
            # print >> sys.stderr, status
            if status == 'RUNNING' or status == 'PENDING':
                time.sleep(self._checkInterval)

    def _serviceGetStatus(self, jobId):
        """Get job status"""
        requestUrl = self._baseUrl + '/status/' + jobId
        return self._restRequest(requestUrl)

    def _restRequest(self, url):
        """Wrapper for a REST (HTTP GET) request"""
        try:
            gcontext = ssl._create_unverified_context()
            reqH = urlopen(url, context=gcontext)
            result = reqH.read()
            reqH.close()
            return result
        except Exception as exc:
            return ''

def main(argv):
    opts, args = getopt.getopt(argv, "", ["sequence=", "outfile="])
    sequence = None
    filename = None
    for opt, arg in opts:
        if opt == '--sequence':
            fh = open(arg, 'r')
            sequence = fh.read()
            fh.close()
        elif opt == '--outfile':
            filename = arg + '.xml'

    service = UnpBlastService(sequence)
    service.RunService()
    service.WriteResultFile(filename)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception as exc:
        print(exc)
        sys.exit(1)
