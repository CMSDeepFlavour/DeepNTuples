#!/usr/bin/env python2

from argparse import ArgumentParser
import os
import subprocess

def syscall(cmd):
  print 'Executing: %s' % cmd
  retval = os.system(cmd)
  if retval != 0:
    raise RuntimeError('Command failed!')

print("This script is still experimental and not fully completed")

parser = ArgumentParser('merge samples')
parser.add_argument('nsamples')
parser.add_argument('outdir')
parser.add_argument('infiles', metavar='N', nargs='+',
                    help='sample list files')
parser.add_argument('--batch', action='store_true')

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    allins=''
    for l in args.infiles:
        allins+=' '+l
        
    syscall('createMergeList '+str(args.nsamples)+' '+args.outdir+' '+allins)
    
    
#read number of jobs
file=open(args.outdir+'/nentries','r')
nJobs=int(file.read())

listtoberun=[]
listsucc=[]

for j in range(nJobs):
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
        continue
    listtoberun.append(j)

print('successful: ',listsucc)

if args.batch:
	script = '''#! /bin/bash

set -o errexit
set -o nounset

cd {cwd}
eval `scramv1 runtime -sh`
cd -

merge {outdir}/mergeconfig $1
'''.format(cwd=os.getcwd(), outdir=args.outdir)
	from time import time
	dname = '%s/merging_%d' % (os.getcwd(), int(time()))
	os.makedirs(dname)
	with open('%s/batch.sh' % dname, 'w') as batch:
		batch.write(script)
	
	header = '''
executable = {dname}/batch.sh
getenv = True
use_x509userproxy = True
+MaxRuntime = 21600
max_transfer_output_mb = 2000
max_retries = 4
'''.format(dname=dname)
	
	template = '''
arguments  = {args}
log        = con_out.{idx}.log
output     = con_out.{idx}.out
error      = con_out.{idx}.err
queue
'''
	jobs = [
		template.format(args=val, idx=idx) 
		for idx, val in enumerate(listtoberun)
		]
	with open('%s/condor.sub' % dname, 'w') as condor:
		condor.write(header)
		condor.write(''.join(jobs))
	os.chdir(dname)
	syscall('condor_submit condor.sub')
	print 'Once all the jobs are run please run again this command to ensure everything worked'
else:
	import multiprocessing as mp

	def worker(j):
		print('starting '+str(j))
		os.system('merge '+args.outdir+'/mergeconfig '+str(j))

	pool = mp.Pool(processes=mp.cpu_count(),) 
	pool.map(worker, listtoberun)


for j in range(nJobs):
    if os.path.exists(args.outdir+'/'+str(j)+'.succ'):
        listsucc.append(j)
    
if len(listsucc) == nJobs:
    print('merge successful, creating file list')
    file=open(args.outdir+'/samples.txt','w')
    for filenumber in listsucc:
        file.write('ntuple_merged_'+str(filenumber)+'.root\n')
    file.close()






