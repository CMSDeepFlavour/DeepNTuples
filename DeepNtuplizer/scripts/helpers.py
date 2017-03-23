

class bcolors:
    
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def resetJobOutput(parentdir,jobno):
    import os
    rootfile=parentdir+'/output/'+parentdir+'_'+str(jobno)+'.root'
    outs=parentdir+'/batch/con_*'+str(jobno)+'.*'
    clusterf=parentdir+"/batch/condorcluster_*"+str(jobno)
    
    os.system('rm -f '+rootfile+' '+outs+' '+clusterf)
    exit()
    os.system('touch '+parentdir+'/batch/con_out.'+ jobno +'.out')
    


def submitjob(path,condorfile):
    import subprocess
    print(condorfile)
    proc = subprocess.Popen(['cd ' + path + ' && condor_submit '+ condorfile], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print (out)
    cluster=out.split()[-1][0:-1]
    return cluster

def createClusterInfo(path,job,cluster):
    import os
    clfile= open(os.path.join(path+'/batch', 'condorcluster_'+str(cluster)+'.'+str(job)), 'w')
    clfile.write(str(cluster))
    clfile.close()

def getCondorStatus():
    
    import subprocess
    proc = subprocess.Popen(['condor_q -nobatch'], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    out=out.split('\n')
    clustersandjobs=[]
    statuses=[]
    isdata=False
    for line in out:
        if len(line)<10 : continue
        tline=line.split()
        if tline[0] == 'ID':
            isdata=True
        if tline[0] == 'ID' or tline[1] == 'jobs;' or not isdata:
            continue

        clustersandjobs.append((tline[0]))
        statuses.append(tline[5])

    return clustersandjobs,statuses

def readStatuses(cluster, nJobs, condorstatus):
    out=[]
    for i in range(nJobs):
        out.append('I')