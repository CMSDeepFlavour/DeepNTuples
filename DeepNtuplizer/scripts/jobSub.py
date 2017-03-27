#!/usr/bin/env python2


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
    jobid = 'jobid' #os.environ['jobid']
    inputdir = os.path.join(swdir, 'inputs')
    inputdir = os.path.join(inputdir, jobid)
    
    parser = argps('submit analyzer to the batch queues')
    parser.add_argument('configfile')
    parser.add_argument('jobdir')
    parser.add_argument('--file',default='samples.cfg',help='file containing a sample list')
    parser.add_argument('--nosubmit',default=False,help='no submission')
    
    args = parser.parse_args()
    
    
    if os.path.isdir(args.jobdir):
        print (args.jobdir), 'exists: EXIT'
        sys.exit(-1)
    configFile=os.path.abspath(args.configfile)
    
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
        
        
    #make samples directory
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
    configFile=os.path.abspath(os.path.join(args.jobdir, os.path.basename(configFile)))
    
    cernboxpath='/eos/user/'+os.environ['USER'][0]+'/'+os.environ['USER']+'/'
    globalOutDir=cernboxpath+'DeepNtuples/'+time.strftime('%a_%H%M%S')+'_'+args.jobdir
    
    
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
                sampleurl='https://cmsweb.cern.ch/das/makepy?dataset='+samplename+'&instance=prod/global'
                print(scriptfile)
                #check if already saved a file list from das?
                ##get from das etc, prepare query
                print('getting script file from DAS')
                dasquery = subprocess.Popen(['wget','--certificate',usercertfile,
                                             '--private-key',userkeyfile,
                                             '--no-check-certificate',
                                             '-O',samplescriptdir+scriptfile+'.py',
                                             sampleurl],
                                             stdout=subprocess.PIPE, 
                                             stderr=subprocess.PIPE)
                sout, serr = dasquery.communicate()
            sample=scriptfile
        
        
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
        #print (cernboxpath)
        ntupleOutDir=globalOutDir+'/'+outputFile+'/'
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
             resetJobOutput(jobpath,job)
             
             
        #create script
        shellscript = """#!/bin/bash
    echo "JOBSUB::RUN job running"
    trap "echo JOBSUB::FAIL job killed" SIGTERM
    export PYTHONPATH={sampleScriptdir}:$PYTHONPATH
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
            sampleScriptdir=samplescriptdir,
            basedir=swdir,
            jobdir=os.path.abspath(args.jobdir)
            )
        
        shellsc = open(sheelscp, 'w')
        shellsc.write(shellscript)
        shellsc.close()
        os.system('chmod +x '+sheelscp)
        os.system('touch '+jobpath+'/batch/nJobs.'+str(nJobs))
        
        #add a 'touch for the .out file to make the check realise it is there
        
        if not args.nosubmit:
            print('submitting '+outputFile)
             #os.system('cd ' + jobpath + ' && echo condor_submit condor.sub') # gives back submitted to cluster XXX message - use
            cluster=submitjob(jobpath,'condor.sub')
            #print(cluster)
            for job in range(0,int(nJobs)):
               createClusterInfo(jobpath,job,cluster,True)
    
    
    exit()
doSub()
