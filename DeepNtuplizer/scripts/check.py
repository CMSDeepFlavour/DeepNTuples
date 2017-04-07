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
    
    stdoutfiles=glob.glob(dir+"/batch/con_out*.out")
    clusterfiles=glob.glob(dir+"/batch/condorcluster_*")
    nJobsFile=glob.glob(dir+"/batch/nJobs*")
    
    #print (stdoutfiles)
    
    nJobs=int(nJobsFile[0].split('.')[1])
    
    failedjobs=[]
    idlejobs=[]
    jobstatus=[]
    
    for i in range(nJobs):
        jobstatus.append('L')
        
    for clf in clusterfiles:
        jcl=clf.split('_')[-1]
        jidx=int(jcl.split('.')[-3])
        jcl=jcl.split('.')[1]+'.'+jcl.split('.')[2]
        
        try:
            idxofjob=clustersandjobs.index(jcl)
            jobstatus[jidx]=statuses[idxofjob]
        except:
            jobstatus[jidx]='L'
        
    
    succjobs=[]
    runjobs=[]
    lostjobs=[]
    holdjobs=[]
    #exit()
    
    #stati: idle, running, failed(file), lost(file)
    
        
    
    for i in range(len(stdoutfiles)):
        filename = stdoutfiles[i]
        #print (filename)
        jobno=int(os.path.basename(filename).split('.')[1])
        jobstat='X'
        if 'JOBSUB::FAIL' in open(filename).read():
            failedjobs.append(filename)
        elif 'JOBSUB::SUCC' in open(filename).read():
            succjobs.append(filename)
        else:
            jobstat=jobstatus[jobno]
            if jobstat == 'I':
                idlejobs.append(filename)
            elif jobstat == 'R':
                runjobs.append(filename)
            elif jobstat == 'H':
                holdjobs.append(filename)
            else:
                lostjobs.append(filename)
    
    
    

    ntupleOutDir=os.path.abspath(dir+'/output/')
    
    
    
    
    for f in failedjobs:
        jobno=os.path.basename(f).split('.')[1]
        
        # remove the failed output
        if action == 'remove':
            resetJobOutput(dir,jobno)
            print('removed output of job '+ str(jobno))
    
        elif action == 'resubmit':
            
            print('resubmitting '+dir+' ' +str(jobno)+ '...')
            resetJobOutput(dir,jobno)
            cluster=submitjob(dir,'batch/condor_'+ jobno+'.sub')
            createClusterInfo(dir,jobno,cluster,False)#single submission
            
        else:
            print('failed job, see: ' + f)
            
            
            #os.system('cd ' + dir + ' && condor_submit batch/condor_'+ jobno+'.sub && touch batch/con_out.'+ jobno +'.out')
           
    allValidjobs=nJobs
    if float(len(idlejobs))/float(allValidjobs)>0.8:
        print (bcolors.BOLD+'idle:    '  + str(len(idlejobs))+bcolors.ENDC)
    else:
        print ('idle:    '  + str(len(idlejobs)))
    if float(len(runjobs))/float(allValidjobs)>0.8:
        print (bcolors.BOLD+'running: '  + str(len(runjobs))+bcolors.ENDC)
    else:
        print ('running: '  + str(len(runjobs)))
    if float(len(succjobs))/float(allValidjobs)>0.95:
        print ('\x1b[6;30;42m'+bcolors.BOLD+'success: '  + str(len(succjobs))+'\x1b[0m')
    elif float(len(succjobs))/float(allValidjobs)>0.9:
        print (bcolors.BOLD+'success: '  + str(len(succjobs))+bcolors.ENDC)
    else:
        print ('success: '  + str(len(succjobs)))
    if float(len(failedjobs))/float(allValidjobs)>0.9:
        print ('\x1b[5;32;31m'+bcolors.BOLD+'failed:  '  + str(len(failedjobs))+'\x1b[0m')
    else:
        print ('failed:  '  + str(len(failedjobs)))
    if float(len(lostjobs))/float(allValidjobs)>0.8:
        print ('\x1b[5;32;31m'+bcolors.BOLD+'lost:    '  + str(len(lostjobs))+'\x1b[0m')
    else:
        print ('lost:    '  + str(len(lostjobs)))
    if float(len(holdjobs))/float(allValidjobs)>0.8:
        print (bcolors.BOLD+'hold:    '  + str(len(holdjobs))+bcolors.ENDC)
    else:
        print ('hold:    '  + str(len(holdjobs)))
    
    if (action == 'merge' or len(succjobs)==nJobs or action == 'filelist') and len(succjobs)>0:
         
        succoutfile=[]
        trainvalfiles=[]
        testfiles=[]
        
        nsucc=len(succjobs)
        
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
               idx+=1
               #can run in parallel
               os.system('hadd '+outputroot+ out) #use subprocess for parallelisatio later
           
           print ('merged to '+str(idx) +' files')
            #os.system('hadd '+outputroot+ succoutfile)
             
             
             
          
          
          
         
         
     
