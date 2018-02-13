#!/usr/bin/env python2


from __future__ import print_function
import sys, os, time
import shutil
from pdb import set_trace
from glob import glob
import re

########################## Parsing and environment ############################
import subprocess

from helpers import submitjob, createClusterInfo, resetJobOutput

def doSub():
    
    from argparse import ArgumentParser as argps

    swdir = os.path.realpath(os.environ['CMSSW_BASE'])
    issgesched = len(os.getenv('SGE_CELL','')) > 0
    jobid = 'jobid' #os.environ['jobid']
    inputdir = os.path.join(swdir, 'inputs')
    inputdir = os.path.join(inputdir, jobid)
    
    parser = argps('submit analyzer to the batch queues')
    parser.add_argument('configfile')
    parser.add_argument('jobdir')
    parser.add_argument('--file',default='samples.cfg',help='file containing a sample list')
    parser.add_argument('--nosubmit',default=False,help='no submission')
    parser.add_argument('--outpath',default='',help='set path to store the .root output')
    parser.add_argument('--walltime',default='21600',help='set job wall time in seconds')
    parser.add_argument('--maxsize',default='2000',help='set maximum allowed size of output ntuple')
    
    args = parser.parse_args()
    
    jobruntime=args.walltime
    maxSize=args.maxsize
    
    #eosGlobalOutDir='/eos/user/'+os.environ['USER'][0]+'/'+os.environ['USER']+'/DeepNtuples'
    eosGlobalOutDir='/eos/cms/store/cmst3/group/dehep/DeepJet/NTuples/'
    if len(args.outpath):
        eosGlobalOutDir=args.outpath
        if not os.path.isdir(eosGlobalOutDir):
            print('please specify a valid output path')
            sys.exit(-1)
    
    
    if os.path.isdir(args.jobdir):
        print (args.jobdir, 'exists: EXIT')
        sys.exit(-1)
    configFile=os.path.realpath(args.configfile)
    
    ###### check for grid certificate
    
    #check grid proxy with voms-proxy-info
    # grep for timeleft and require at least 3 hours 
    # and check wether the path lives in AFS (voms-proxy-info -path)
    #save this output for auto wget on file lists
    checkcert = subprocess.Popen(['voms-proxy-info','-path'], stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE)
    sout, serr = checkcert.communicate()
    certpath=sout
    checkcert = subprocess.Popen(['voms-proxy-info','-timeleft'], stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE)
    sout, serr = checkcert.communicate()
    certtime=sout
    
    if float(certtime) < 2*60*60:
        print('grid proxy loses validity in less than 2 hours, please renew before job submission.')
        exit()
        
    usercertfile=os.getenv('HOME')+'/.globus/usercert.pem'
    userkeyfile=os.getenv('HOME')+'/.globus/userkey.pem'
    
    nousercertsfound=False
    if not os.path.isfile(usercertfile):
        print('pleace locate your grid certificate file in ~/.globus/usercert.pem')
        nousercertsfound=True
    if not os.path.isfile(userkeyfile):
        print('pleace locate your grid key file in ~/.globus/userkey.pem')
        nousercertsfound=True
        
    if nousercertsfound:
        print('please follow the Twiki https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid')
        exit()
        
        
    #recreates samples directory (removing old one avoids possible errors in creating importsamples)
    samplescriptdir=os.getenv('HOME')+'/.deepntuples_scripts_tmp'
    if not os.path.isdir(samplescriptdir):
        os.mkdir(samplescriptdir)
    
    samplescriptdir+='/'
    
    #https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid
    
    #testurl='https://cmsweb.cern.ch/das/makepy?dataset=/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16MiniAOD-PhaseIFall16PUFlat20to50_PhaseIFall16_81X_upgrade2017_realistic_v26-v1/MINIAODSIM&instance=prod/global'
    
    #make a system call to wget.. urllib is not powerful enoguh apparently
    
    
    #checkcert = subprocess.Popen(['wget','--certificate='+certpath,'-o bla.py', testurl], stdout=subprocess.PIPE, 
    #                                 stderr=subprocess.PIPE)
    
    #print(checkcert.communicate())
    #exit()#testing
    
    ######## all checks done. From now on it just runs
    
    # create output dir
    
    os.mkdir(args.jobdir)
    shutil.copy(configFile, args.jobdir)
    configFile=os.path.realpath(os.path.join(args.jobdir, os.path.basename(configFile)))
    
    
    globalOutDir=eosGlobalOutDir+'/'+time.strftime('%a_%H%M%S')+'_'+args.jobdir
    globalOutDir=os.path.realpath(globalOutDir)
    
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
            
        
        #check if sufficient files
        
        
        ######check out from DAS
        samplename=entries[1]
        isdasname=False
        #do differently if it is a DAS name
        if '/' in samplename: #this is DAS name
            isdasname=True
            if '*' in samplename:
                print('no wildcards in sample names allowed')
                exit()
        
        print('preparing\n'+samplename)
        sample=""
        if not isdasname:
            sample=samplesdir+samplename
        else:
            import string
            chars = re.escape(string.punctuation)
            scriptfile=re.sub(r'['+chars+']', '', str(samplename))
            scriptfile=scriptfile
            if not os.path.isfile(samplescriptdir+scriptfile+'.py'):
                cmd = 'dasgoclient -query="file dataset=%s"' % (samplename)
                dasquery = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                sout = dasquery.communicate()[0]
                filelist = ['"%s",' % f for f in sout.strip().split('\n')]

                template_sample = os.path.join(os.environ['CMSSW_BASE'], 'src/DeepNTuples/DeepNtuplizer/python/samples/samples_template.py')
                dest_file = samplescriptdir+scriptfile+'.py'
                with open(template_sample) as temp:
                    s = temp.read().replace('_FILES_', '\n'.join(filelist))
                    with open(dest_file, 'w') as fout:
                        fout.write(s)

            sample=scriptfile
        
        
        
        sys.path.append(samplescriptdir)
        sys.path.append(swdir+'/src/DeepNTuples/DeepNtuplizer/python/samples')
        importsample=sample.split('.')[-1]
        
        #print(importsample)
        cmssource = __import__(importsample)
        #print(cmssource.source)
        
        nJobs=entries[0]
        totalfiles=len(cmssource.source.fileNames)+len(cmssource.source.secondaryFileNames)
        if int(nJobs)>totalfiles:
            print('reduced number of jobs to number of files (',totalfiles,') from ', nJobs)
            nJobs=str(totalfiles)
        
        outputFile=entries[2]
        jobargs=''
        if len(entries) >3:
            jobargs=entries[3]
        jobpath = os.path.join(
              args.jobdir, 
              outputFile
              )
        jobpath=os.path.realpath(jobpath)
        os.mkdir(jobpath) 
        sheelscp=os.path.realpath(os.path.join(jobpath, 'batchscript.sh'))
        
        #create full output path on eos
        #print (eosGlobalOutDir)
        ntupleOutDir=globalOutDir+'/'+outputFile+'/output/'
        os.makedirs(ntupleOutDir)
        #print (ntupleOutDir)
        
        #link to ntupleOutDir
        os.symlink(ntupleOutDir,jobpath+'/output')
        #ntupleOutDir=jobpath+'/output/'

        hostname=os.getenv('HOSTNAME')
        if not hostname:
            raise Exception("hostname could not be determined")

	#create directory in /tmp/ to host the .out and .log files
	#logDir=globalOutDir+'/'+outputFile+'/batch/' #'/tmp/'+os.environ['USER']+'/batch/'
        helperdir='/tmp/'+os.environ['USER']+'/'+ str(os.getpid()) +'/'+outputFile+'/batch/'
        logDir=jobpath+'/batch/'
        os.makedirs(logDir)
        os.makedirs(helperdir)
        os.symlink(helperdir,jobpath+'/helper')
        infostring='helper files in \n'+helperdir+' on machine\n'+hostname
        #print(infostring+'\n delete after all jobs are finished (check.py does it automatically)')
        
	    #create a txt file and save the eos path there
        address = open(jobpath+'/hostinfo.txt','w')
        #address.write(logDir)
        
            
        address.write(hostname)
        address.close()
        
        #print(ntupleOutDir)
#The maximum wall time of a condor job is defined in the MaxRuntime parameter in seconds.
# 3 hours (10800s) seems to be currently enough

        condorfile ="""executable            = {batchscriptpath}
arguments             = {configfile} inputScript={sample} nJobs={njobs} job=$(ProcId) {options}
log                   = {logdir}con_out.$(ProcId).log
getenv = True
environment = "NTUPLEOUTFILEPATH={ntupledir}{outputfile}_$(ProcId).root"
use_x509userproxy = True
+MaxRuntime = {maxruntime}
max_transfer_output_mb = {maxsize}
RequestCpus = 2
transfer_output_remaps = "stdout.txt={logdir}con_out.$(ProcId).out"
max_retries = 20
queue {njobs}
    """.format(
              batchscriptpath=sheelscp,
              configfile=configFile, 
              sample=sample,
              ntupledir=ntupleOutDir,
              outputfile=outputFile,
              njobs=nJobs, 
              options=jobargs,
              maxruntime=jobruntime,
	      maxsize=maxSize,
	      logdir=logDir
              )
        
        conf = open(os.path.join(jobpath, 'condor.sub'), 'w')
        conf.write(condorfile)
        conf.close()
        print("wrote condor file for "+outputFile)
        
        allsgescripts=[]
        
        #create individual condor files for resubmission
        for job in range(0,int(nJobs)):
            jobcondorfile ="""executable            = {batchscriptpath}
arguments = {configfile} inputScript={sample} nJobs={njobs} job={job} {options}
log   = {logdir}con_out.{job}.log
getenv = True
environment = "NTUPLEOUTFILEPATH={ntupledir}{outputfile}_{job}.root"
use_x509userproxy = True
+MaxRuntime= {maxruntime}
RequestCpus = 2
max_transfer_output_mb = {maxsize}
transfer_output_remaps = "stdout.txt={logdir}con_out.{job}.out"
max_retries = 3
queue 1
             """.format(
                  batchscriptpath=sheelscp,
                  configfile=configFile, 
                  sample=sample,
                  ntupledir=ntupleOutDir,
                  outputfile=outputFile,
                  njobs=nJobs, 
                  options=jobargs,
                  job=str(job),
                  maxruntime=jobruntime,
                  maxsize=maxSize,
		          logdir=logDir
                  )
            jconf = open(os.path.join(logDir,'condor_'+str(job)+'.sub'), 'w')
            jconf.write(jobcondorfile)
            jconf.close()
        
            jobsgefile="""
#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#$ -l site=hh
#$ -l distro=sld6
#
#(the cpu time for this job)
#$ -l h_rt=05:55:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=4096M
#$ -l cvmfs
#(stderr and stdout are merged together to stdout)
#$ -j y
#$ -m a
#$ -cwd -V
#( -l h_stack=1536M) #try with small stack
#$ -pe local 1 -R y
#$ -P af-cms

export LOGDIR={logdir}
export JOB={job}
export NTUPLEOUTFILEPATH={ntupledir}{outputfile}_{job}.root

{batchscriptpath} {configfile} inputScript={sample} nJobs={njobs} job={job} {options}
            """.format(
                  batchscriptpath=sheelscp,
                  configfile=configFile, 
                  sample=sample,
                  ntupledir=ntupleOutDir,
                  outputfile=outputFile,
                  njobs=nJobs, 
                  options=jobargs,
                  job=str(job),
                  #maxruntime=jobruntime,
                  #maxsize=maxSize,
                  logdir=logDir
                  )
            if issgesched or True:
                sgefile=os.path.join(logDir,'sge_'+str(job)+'.sh')
                jconf = open(sgefile, 'w')
                jconf.write(jobsgefile)
                jconf.close()
                allsgescripts.append(sgefile)
                os.system('chmod +x '+sgefile)
        
            resetJobOutput(jobpath,job)
             
             
        #create script
        shellscript = """#!/bin/bash
        
workdir=""
if [ $SGE_CELL ]
then
workdir=`mktemp -d -t DeepNTuplesXXXXXX`        
cd $workdir
workdir=$workdir"/"
fi

exec > "$PWD/stdout.txt" 2>&1
echo "JOBSUB::RUN job running"
trap "echo JOBSUB::FAIL job killed" SIGTERM
export OUTPUT=$workdir"{outputfile}"
cd {basedir}
eval `scramv1 runtime -sh`
export PYTHONPATH={sampleScriptdir}:$PYTHONPATH
which cmsRun
cd -
cmsRun -n 1 "$@" outputFile=$OUTPUT 
exitstatus=$?
if [ $exitstatus != 0 ]
then
   echo JOBSUB::FAIL job failed with status $exitstatus
else
   pwd
   ls -ltr $OUTPUT*.root
   if [ $SGE_CELL ]
   then
     cp $OUTPUT*.root $NTUPLEOUTFILEPATH
   else
     eoscp $OUTPUT*.root $NTUPLEOUTFILEPATH
   fi
   exitstatus=$?
   rm -f $OUTPUT*.root
   if [ $exitstatus != 0 ]
   then
     echo JOBSUB::FAIL job failed with status $exitstatus
   else
     echo JOBSUB::SUCC job ended sucessfully
   fi
fi
rm -f $OUTPUT*.root
if [ $workdir ]
then
# JOB is only defined for SGE submit
cp $workdir/stdout.txt $LOGDIR/con_out.$JOB.out
rm -rf $workdir
fi
exit $exitstatus
        """.format(
            sampleScriptdir=samplescriptdir,
            basedir=swdir,
            jobdir=os.path.realpath(args.jobdir),
	        outputfile=outputFile
            )
        
        shellsc = open(sheelscp, 'w')
        shellsc.write(shellscript)
        shellsc.close()
        os.system('chmod +x '+sheelscp)
        os.system('touch '+logDir+'/nJobs.'+str(nJobs))
        
        sgelist=open(jobpath+'/sge_sub.sh','w')
        for line in allsgescripts:
            sgelist.write('qsub '+line+'\n')
        os.system('chmod +x '+jobpath+'/sge_sub.sh')
        sgelist.close()
        
        #add a 'touch for the .out file to make the check realise it is there
        
        if not args.nosubmit and not issgesched:
            print('submitting '+outputFile)
             #os.system('cd ' + jobpath + ' && echo condor_submit condor.sub') # gives back submitted to cluster XXX message - use
            cluster=submitjob(jobpath,'condor.sub')
            #print(cluster)
            for job in range(0,int(nJobs)):
                createClusterInfo(jobpath,job,cluster,True)
    
        
        if not args.nosubmit and issgesched:
            pass
        
    exit()
    
    
doSub()
