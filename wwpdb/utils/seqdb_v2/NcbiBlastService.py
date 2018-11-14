"""
File:    NcbiBlastService.py
Author:  Zukang Feng
Update:  24-May-2010
Version: 001  Initial version

"""

__author__  = "Zukang Feng"
__email__   = "zfeng@rcsb.rutgers.edu"
__version__ = "V0.001"

import re, os, sys, time, urllib, urllib2, ssl
import getopt

class NcbiBlastService:
    """ Utility class for running blastn service using nt database from NCBI site
    """
    def __init__(self, sequence):
        self._sequence = sequence;
        self._baseUrl = 'https://www.ncbi.nlm.nih.gov/blast/Blast.cgi'
        self._result = None
        self._checkInterval = 10

    def RunService(self):
        """Main function for running blastn service"""
        params = {}
        params['CMD'] = 'Put'
        params['DATABASE'] = 'nt'
        params['FILTER'] = 'L'
        params['QUERY'] = self._sequence
        params['PROGRAM'] = 'blastn'
        # Submit the job
        return_data = self._serviceSubmit(params)
        jobid = self._GetJobId(return_data)
        if jobid:
            # print >> sys.stderr, jobid
            time.sleep(10)
            self._result = self._getResultfromServer(jobid)

    def GetResult(self):
        return self._result

    def WriteResultFile(self, filename):
        fh = open(filename, 'w')
        if self._result:
            fh.write(self._result)
        fh.close()

    def _serviceSubmit(self, params):
        """Submit job"""
        # Get the data for the options
        try:
            requestData = urllib.urlencode(params)
            request = urllib2.Request(self._baseUrl, requestData)
            # We are ignoring the certificates
            gcontext = ssl._create_unverified_context()
            reqH = urllib2.urlopen(request, context=gcontext)
            return_data = reqH.read()
            reqH.close()    
            return return_data
        except Exception, exc:
            #print exc
            return ''

    def _GetJobId(self, data):
        """parse return html file and search for
              <!--QBlastInfoBegin
                  RID = X8SRH0TH012
                  RTOE = 11
              QBlastInfoEnd
              -->
           return job ID X8SRH0TH012
        """

        start = data.find('QBlastInfoBegin')
        end = data.find('QBlastInfoEnd')

        if start >= 0 and end > start + 15:
            list = data[start+15:end].split('\n')
            for e in list:
                if e.find('RID') >= 0:
                    list1 = e.strip().split(' ')
                    if len(list1) == 3:
                        return list1[2]

        return ""

    def _getResultfromServer(self, jobId):
        """Get result from server"""
        self._clientPoll(jobId)

        # status is READY, get the result xml file
        params = {}
        params['CMD'] = 'Get'
        params['RID'] = jobId
        params['NCBI_GI'] = 'yes'
        params['FORMAT_TYPE'] = 'XML'
        return self._serviceSubmit(params)

    def _clientPoll(self, jobId):
        """Client-side poll"""
        status = self._serviceGetStatus(jobId)
        while status != 'READY':
            # print >> sys.stderr, status
            time.sleep(self._checkInterval)
            status = self._serviceGetStatus(jobId)

    def _serviceGetStatus(self, jobId):
        """Get job status"""
        params = {}
        params['CMD'] = 'Get'
        params['RID'] = jobId
        # Submit the job
        return_data = self._serviceSubmit(params)
        return self._GetJobStatus(return_data)

    def _GetJobStatus(self, data):
        """parse return html file and search for
               QBlastInfoBegin
                       Status=READY
               QBlastInfoEnd
           return status READY
        """

        start = data.find('QBlastInfoBegin')
        end = data.find('QBlastInfoEnd')

        if start >= 0 and end > start + 15:
            list = data[start+15:end].split('\n')
            for e in list:
                if e.find('Status') >= 0:
                    list1 = e.strip().split('=')
                    if len(list1) == 2:
                        return list1[1]

        return ""

def main(argv):
    opts, args = getopt.getopt(argv, "s:o:", ["sequence=", "outfile="])
    sequence = None
    filename = None
    for opt, arg in opts:
        if opt in ('-s', '--sequence'):
            fh = open(arg, 'r')
            sequence = fh.read()
            fh.close()
        elif opt in ('-o', '--outfile'):
            filename = arg + '.xml'

    if sequence and filename:
        sequence = re.sub('[\t \n]', '', sequence)
        service = NcbiBlastService(sequence)
        service.RunService()
        service.WriteResultFile(filename)

if __name__ == "__main__":
    try:
        main(sys.argv[1:])
        sys.exit(0)
    except Exception, exc:
        print exc
        sys.exit(1)
