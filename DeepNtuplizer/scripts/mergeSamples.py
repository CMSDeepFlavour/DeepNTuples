#!/usr/bin/env python2


from argparse import ArgumentParser
import os
import subprocess

print("This script is still experimental and not fully completed")

parser = ArgumentParser('merge samples')
parser.add_argument('nsamples')
parser.add_argument('outdir')
parser.add_argument('infiles', metavar='N', nargs='+',
                    help='sample list files')

args = parser.parse_args()

if not os.path.isdir(args.outdir):

    allins=''
    for l in args.infiles:
        allins+=' '+l
        
    os.system('createMergeList '+str(args.nsamples)+' '+args.outdir+' '+allins)
    
    
#read number of jobs
file=open(args.outdir+'/nentries','r')
nJobs=file.read()

listtoberun=[]
listsucc=[]

for j in range(int(nJobs)):
    
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
        continue
    
    listtoberun.append(j)

print('successful: ',listsucc)

import multiprocessing as mp


def worker(j):
    print('starting '+str(j))
    os.system('merge '+args.outdir+'/mergeconfig '+str(j))


pool = mp.Pool(processes=mp.cpu_count(),) 
pool.map(worker, listtoberun)


for j in range(int(nJobs)):
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
    
if len(listsucc) == int(nJobs):
    print('merge successful, creating file list')
    file=open(args.outdir+'/samples.txt','w')
    for filenumber in listsucc:
        file.write('ntuple_merged_'+str(filenumber)+'.root\n')
    file.close()






