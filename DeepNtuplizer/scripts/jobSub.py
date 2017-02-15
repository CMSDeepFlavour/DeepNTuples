#!/usr/bin/env python2


import sys, os, time
import shutil
from pdb import set_trace
from glob import glob
from argparse import ArgumentParser
import re
#import URAnalysis.Utilities.prettyjson as prettyjson
#from URAnalysis.PlotTools.data_views import get_best_style, log #not for the same purpose, but the same mechanism
#import rootpy
#log.setLevel(log.CRITICAL)

swdir = os.path.realpath(os.environ['CMSSW_BASE'])
jobid = 'jobid' #os.environ['jobid']
inputdir = os.path.join(swdir, 'inputs')
inputdir = os.path.join(inputdir, jobid)

parser = ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('configfile')
parser.add_argument('jobdir')
parser.add_argument('--file',default='samples.cfg',help='file containing a sample list')
parser.add_argument('--nosubmit',default=False,help='no submission')

args = parser.parse_args()


if os.path.isdir(args.jobdir):
    print (args.jobdir), 'exists: EXIT'
    sys.exit(-1)
os.mkdir(args.jobdir)
configFile=os.path.abspath(args.configfile)
shutil.copy(configFile, args.jobdir)
configFile=os.path.abspath(os.path.join(args.jobdir, os.path.basename(configFile)))
print ('submitting jobs for '+configFile)

samplesdir='DeepNTuples.DeepNtuplizer.samples.'

#format: njobs  sample  output  args1 args2 ... (simple whitespace)
lines = [line.rstrip('\n') for line in open(args.file)]

for sampledescription in lines:
    
    if sampledescription.strip().startswith("#"):
        continue
    
    entries= [s.strip() for s in sampledescription.split('  ') if s]
    if len(entries) < 3:
        continue
        
    nJobs=entries[0]
    sample=samplesdir+entries[1]
    print('preparing '+sample)
    
    outputFile=entries[2]
    jobargs=''
    if len(entries) >3:
        jobargs=entries[3]
    jobpath = os.path.join(
          args.jobdir, 
          outputFile
          )
    jobpath=os.path.abspath(jobpath)
    os.mkdir(jobpath) 
    sheelscp=os.path.abspath(os.path.join(jobpath, 'batchscript.sh'))
    
    #create full output path on eos
    cernboxpath='/eos/user/'+os.environ['USER'][0]+'/'+os.environ['USER']+'/'
    #print (cernboxpath)
    ntupleOutDir=cernboxpath+'DeepNtuples/'+time.strftime('%a_%H%M%S')+'_'+args.jobdir+'/'+outputFile+'/'
    os.makedirs(ntupleOutDir)
    #print (ntupleOutDir)
    
    #link to ntupleOutDir
    os.symlink(ntupleOutDir,jobpath+'/output')
       
    condorfile ="""executable            = {batchscriptpath}
arguments             = {configfile} inputScript={sample} outputFile={ntupledir}{outputfile} nJobs={njobs} job=$(ProcId) {options}
output                = batch/con_out.$(ProcId).out
error                 = batch/con_out.$(ProcId).err
log                   = batch/con_out.$(ProcId).log
send_credential        = True
use_x509userproxy = True
queue {njobs}
""".format(
          batchscriptpath=sheelscp,
          configfile=configFile, 
          sample=sample,
          ntupledir=ntupleOutDir,
          outputfile=outputFile,
          njobs=nJobs, 
          options=jobargs
          )
    
    conf = open(os.path.join(jobpath, 'condor.sub'), 'w')
    conf.write(condorfile)
    conf.close()
    print("wrote condor file for "+outputFile)
    os.mkdir(jobpath+'/batch')
    
    #create individual condor files for resubmission
    for job in range(0,int(nJobs)):
         jobcondorfile ="""executable            = {batchscriptpath}
arguments             = {configfile} inputScript={sample} outputFile={ntupledir}{outputfile} nJobs={njobs} job={job} {options}
output                = batch/con_out.{job}.out
error                 = batch/con_out.{job}.err
log                   = batch/con_out.{job}.log
send_credential        = True
use_x509userproxy = True
queue 1
         """.format(
              batchscriptpath=sheelscp,
              configfile=configFile, 
              sample=sample,
              ntupledir=ntupleOutDir,
              outputfile=outputFile,
              njobs=nJobs, 
              options=jobargs,
              job=str(job)
          )
         jconf = open(os.path.join(jobpath+'/batch', 'condor_'+str(job)+'.sub'), 'w')
         jconf.write(jobcondorfile)
         jconf.close()
         
         
    #create script
    shellscript = """#!/bin/bash
echo "JOBSUB::RUN job running"
trap "echo JOBSUB::FAIL job killed" SIGTERM
cd {basedir}
eval `scramv1 runtime -sh`
cd {jobdir}
cmsRun "$@" 2>&1
exitstatus=$?
if [ $exitstatus != 0 ]
then
echo JOBSUB::FAIL job failed with status $exitstatus
else
echo JOBSUB::SUCC job ended sucessfully
fi
    """.format(
        basedir=swdir,
        jobdir=os.path.abspath(args.jobdir)
        )
    
    shellsc = open(sheelscp, 'w')
    shellsc.write(shellscript)
    shellsc.close()
    os.system('chmod +x '+sheelscp)
    
    #add a 'touch for the .out file to make the check realise it is there
    
    if not args.nosubmit:
         os.system('cd ' + jobpath + ' && condor_submit condor.sub') # gives back submitted to cluster XXX message - use




exit()

