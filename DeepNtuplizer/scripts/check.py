#!/usr/bin/env python2

# read stdout from each job to get job status (cmsRun success or not)
# use log files to evaluate batch status
# get job ids from file list
# 
# print number of suvvessful jobs
#


# usage: check.py singledir


#for each job in 


import os
import glob
from argparse import ArgumentParser as ARgpsrs


from helpers import submitjob, createClusterInfo, getCondorStatus, resetJobOutput, bcolors

#R: running, X: removed

    

parse = ARgpsrs('check the status of the jobs on the batch queue, pass a directory and an action. Action can be: none, remove, resubmit. The action remove removes failed output. After that, the failed jobs can be resubmitted with resubmit')

parse.add_argument('dirs', metavar='N', nargs='+',
                    help='dirs to check')

parse.add_argument('--action',default='none',
                    help='create file lists even though jobs are not finished yet')

#parser.add_argument('dir')
#parser.add_argument('action')
args = parse.parse_args()
dirs=args.dirs

action=args.action

issgesched = len(os.getenv('SGE_CELL','')) > 0

if not  issgesched:
    print('getting condor batch status...\n')
    clustersandjobs,statuses=getCondorStatus()

#if len(clustersandjobs) < 1:
#    raise Exception("condor_q not available, cannot check job status")

#for i in range(len(clustersandjobs)):
#    print(clustersandjobs[i],statuses[i])



for dir in dirs:
    
    print('\nchecking dir '+dir)
    
    if len(dir.split('/'))>1:
        print ('please run this script directly in the parent directory of the job directory')
        exit()
    with open(dir+'/hostinfo.txt','r') as f:
        hostname=f.readline().strip()
        
    if not hostname == os.getenv('HOSTNAME'):
        raise Exception("check must be run on host "+hostname)
    
    #stdoutfiles=glob.glob(dir+"/batch/con_out*.out")
    clusterfiles=glob.glob(dir+"/helper/condorcluster_*")
    nJobsFile=glob.glob(dir+"/batch/nJobs*")
    
    
    nJobs=int(nJobsFile[0].split('.')[1])
    
    jobstatus_list=[]
    
    for i in range(nJobs):
        jobstatus_list.append('L')
        
    for clf in clusterfiles:
        jcl=clf.split('_')[-1]
        jidx=int(jcl.split('.')[-3])
        jcl=jcl.split('.')[1]+'.'+jcl.split('.')[2]
        
        try:
            idxofjob=clustersandjobs.index(jcl)
            jobstatus_list[jidx]=statuses[idxofjob]
        except:
            jobstatus_list[jidx]='L'
        
    
    failedjobs=[]
    idlejobs=[]
    succjobs=[]
    runjobs=[]
    lostjobs=[]
    holdjobs=[]
    #exit()
    
    #stati: idle, running, failed(file), lost(file)
    
    nsucc=0
    nhold=0
    nfail=0
    nlost=0
    nidle=0
    nrunning=0
    
    for i in range(nJobs):
        jobno=i
        filename = dir+"/batch/con_out."+str(jobno)+ ".out"
        if os.path.isfile(dir+"/helper/"+str(jobno)+'.succ'):
            jobstatus_list[jobno]='S'
        elif os.path.isfile(filename) and 'JOBSUB::SUCC' in open(filename).read():
            jobstatus_list[jobno]='S'
        elif os.path.isfile(filename) and  'JOBSUB::FAIL' in open(filename).read():
            jobstatus_list[jobno]='F'
            os.system('touch '+dir+"/batch/"+str(jobno)+'.succ')

           
        if jobstatus_list[jobno]=='S' and not os.path.isfile(dir+'/output/'+dir+'_'+str(jobno)+'.root'):
            jobstatus_list[jobno]='F'
            os.system('rm -f '+dir+"/helper/"+str(jobno)+'.succ')
        
        j=jobstatus_list[jobno]
        if j=='S':
            nsucc+=1
            succjobs.append(filename)
        elif j=='F':
            nfail+=1
            failedjobs.append(filename)
        elif j=='H':
            nhold+=1
            holdjobs.append(filename)
        elif j=='L' or j=='X':
            nlost+=1
            lostjobs.append(filename)
        elif j=='I':
            nidle+=1
            idlejobs.append(filename)
        elif j=='R':
            nrunning+=1
            runjobs.append(filename)
    

    ntupleOutDir=os.path.abspath(dir+'/output/')
    
###
    #get some numbers
    
    nJobs=nJobs
    if float(nidle)/float(nJobs)>0.8:
        print (bcolors.BOLD+'idle:    '  + str(nidle)+bcolors.ENDC)
    else:
        print ('idle:    '  + str(nidle))
    if float(nrunning)/float(nJobs)>0.8:
        print (bcolors.BOLD+'running: '  + str(nrunning)+bcolors.ENDC)
    else:
        print ('running: '  + str(nrunning))
    if float(nsucc)/float(nJobs)>0.95:
        print ('\x1b[6;30;42m'+bcolors.BOLD+'success: '  + str(nsucc)+'\x1b[0m')
    elif float(nsucc)/float(nJobs)>0.9:
        print (bcolors.BOLD+'success: '  + str(nsucc)+bcolors.ENDC)
    else:
        print ('success: '  + str(nsucc))
    if float(nfail)/float(nJobs)>0.9:
        print ('\x1b[5;32;31m'+bcolors.BOLD+'failed:  '  + str(nfail)+'\x1b[0m')
    else:
        print ('failed:  '  + str(nfail))
    if float(nlost)/float(nJobs)>0.8:
        print ('\x1b[5;32;31m'+bcolors.BOLD+'lost:    '  + str(nlost)+'\x1b[0m')
    else:
        print ('lost:    '  + str(nlost))
    if float(nhold)/float(nJobs)>0.8:
        print (bcolors.BOLD+'hold:    '  + str(nhold)+bcolors.ENDC)
    else:
        print ('hold:    '  + str(nhold))
        
        
        
    
    for f in failedjobs:
        jobno=os.path.basename(f).split('.')[1]
        
        # remove the failed output
        if action == 'remove':
            resetJobOutput(dir,jobno)
            print('removed output of job '+ str(jobno))
    
        elif action == 'resubmit':
            
            print('resubmitting '+dir+' ' +str(jobno)+ '...')
            if issgesched:
                resubfile=dir+'/batch/sge_'+ str(jobno)+'.sh'
                resetJobOutput(dir,jobno)
                os.system('qsub '+resubfile)
            else:
                resetJobOutput(dir,jobno)
                cluster=submitjob(dir,'batch/condor_'+ jobno+'.sub',jobno)
                createClusterInfo(dir,jobno,cluster,False)#single submission
            
        else:
            print('failed job, see: ' + f)
            
    if action == 'resubmit' and not issgesched:
        for f in lostjobs:
            jobno=os.path.basename(f).split('.')[1]
            print('resubmitting '+dir+' ' +str(jobno)+ '...')
            resetJobOutput(dir,jobno)
            cluster=submitjob(dir,'batch/condor_'+ jobno+'.sub',jobno)
            createClusterInfo(dir,jobno,cluster,False)
                    
            #os.system('cd ' + dir + ' && condor_submit batch/condor_'+ jobno+'.sub && touch batch/con_out.'+ jobno +'.out')
           
    
    if (action == 'merge' or nsucc==nJobs or action == 'filelist') and nsucc>0:
         
        succoutfile=[]
        trainvalfiles=[]
        testfiles=[]
        
        
        tocombine=' '
        combinedsize=0
        succfilellist=[]
        for i in range(nsucc):
            
            f=succjobs[i]
            jobno=os.path.basename(f).split('.')[1]
            outputFile=dir+'_'+jobno
            succfilellist.append(outputFile+'.root')
            jobfrac = float(i)/float(nsucc)
            if jobfrac < 0.9:
                trainvalfiles.append(outputFile+'.root')
            else:
                testfiles.append(outputFile+'.root')
            
            fulloutfile=os.path.join(ntupleOutDir,outputFile+'.root')

            combinedsize+=os.path.getsize(fulloutfile)
            if combinedsize < 4e9:
                tocombine = tocombine+' '+fulloutfile
            else:
                succoutfile.append(tocombine)
                tocombine=' '
                combinedsize=0
        # add remaining ones
        succoutfile.append(tocombine)
        fileListFile=open(ntupleOutDir+'/all_samples.txt','w')
        for l in succfilellist:
           fileListFile.write("%s\n" % l)
        fileListFile.close()
        fileListFile=open(ntupleOutDir+'/train_val_samples.txt','w')
        for l in trainvalfiles:
           fileListFile.write("%s\n" % l)
        fileListFile.close()
        fileListFile=open(ntupleOutDir+'/test_samples.txt','w')
        for l in testfiles:
           fileListFile.write("%s\n" % l)
        fileListFile.close()
        
        
        idx=0
        #processes=[]
        if action == 'merge':
           print('merging...')
           for out in succoutfile:
               outputroot=ntupleOutDir+'/'+dir+'_merged_'+str(idx)+'.root'
	       print (outputroot)
               idx+=1
               #can run in parallel
               os.system('hadd '+outputroot+ out) #use subprocess for parallelisatio later
           
           print ('merged to '+str(idx) +' files')
            #os.system('hadd '+outputroot+ succoutfile)
             
        
            
    if (action == 'clean'):         
        tmpdir=os.path.realpath(dir+"/helper")
        print('removing temp dir '+tmpdir)
        os.system("rm -rf "+tmpdir)
    else:
        print('if no jobs are running anymore, please run "check.py <dirs> --action clean" to clean temporary files')
              
          
         
         
     
