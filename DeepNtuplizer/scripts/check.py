#!/usr/bin/env python2

# read stdout from each job to get job status (cmsRun success or not)
# use log files to evaluate batch status
# get job ids from file list
# 
# print number of suvvessful jobs
#


# usage: check.py singledir


#for each job in 


import glob
from argparse import ArgumentParser

parser = ArgumentParser('check the status of the jobs on the batch queue')
parser.add_argument('dir')
args = parser.parse_args()
dir=args.dir

stdoutfiles=glob.glob(dir+"/batch/*.out")

#print (stdoutfiles)

failedjobs=[]
succjobs=[]
runjobs=[]
#exit()

for filename in stdoutfiles:
    if 'JOBSUB::FAIL' in open(filename).read():
        failedjobs.append(filename)
    elif 'JOBSUB::SUCC' in open(filename).read():
        succjobs.append(filename)
    elif 'JOBSUB::RUN' in open(filename).read():
        runjobs.append(filename)


print ('idle: ' + str(len(stdoutfiles)-len(failedjobs)-len(succjobs)-len(runjobs)))
print ('failed: ' + str(len(failedjobs)))
print ('success: ' + str(len(succjobs)))
print ('running: ' + str(len(runjobs)))

        