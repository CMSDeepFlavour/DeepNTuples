

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
    clusterf=parentdir+"/helper/condorcluster_"+str(jobno)+'*'
    submitindicator=parentdir+"/helper/"+str(jobno)+".submitted"
    os.system('rm -f '+rootfile+' '+outs+' '+clusterf+' '+submitindicator)
    #os.system('touch '+parentdir+'/batch/con_out.'+ str(jobno) +'.out')
    


def submitjob(path,condorfile,jobno=-1):
    import subprocess
    import os
    print(condorfile)
    proc = subprocess.Popen(['cd ' + path + ' && condor_submit '+ condorfile], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    print (out)
    if err and len(err):
        print (err)
    elif jobno>-1:
        os.system('touch '+path +'/helper/'+str(jobno)+".submitted")
        #os.system('touch '+path+'/batch/con_out.'+ str(jobno) +'.out')
    cluster=0
    try:
        cluster=out.split()[-1][0:-1]
    except:
        0
    return cluster

def createClusterInfo(path,job,cluster,batchsub):
    import os
    addstring='0'
    if batchsub:
        addstring=str(job)
    clfile= open(os.path.join(path+'/helper', 'condorcluster_'+str(job)+'.'+str(cluster)+'.'+addstring), 'w')
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
