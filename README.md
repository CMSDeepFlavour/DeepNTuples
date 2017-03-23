# DeepNTuples
NTuple framework for DeepFlavour

Installation
============

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src/
cmsenv
git cms-init
git clone https://github.com/mverzett/DeepNTuples
# Add JetToolBox
cd DeepNTuples
git submodule init
git submodule update

# Add DeepFlavour -- To be updated once the 80X PR is done
cd -
git cms-merge-topic -u mverzett:DeepFlavour-from-CMSSW_8_0_21
mkdir RecoBTag/DeepFlavour/data/
cd RecoBTag/DeepFlavour/data/
wget http://home.fnal.gov/~verzetti//DeepFlavour/training/DeepFlavourNoSL.json
cd -
```

It is important to create your grid proxy in a location that is accessible by other nodes (there is no security issue, your full credentials are still needed for access). For this purpose, redirect the grid proxy location by adding the following to your login script:

```
export X509_USER_PROXY=${HOME}/.gridproxy.pem
```

Production
==========

The jobs can be submitted in the production directory using the following syntax
```
jobSub.py --file <sample file> DeepNtuplizer.py <batch directory>
```
For an example of sample files, please refer to the ones already in the directory.

The status of the jobs can be checked with
```
cd <batch directory>
check.py <sample subdirectories to be checked>
```

The check.py script provides additional options to resubmit failed jobs or to create sample lists in case a satisfying fraction of jobs ended successfully. In this case do:
```
check.py <sample subdirectories to be checked> --action filelist
```

This will create file lists that can be further processed by the DeepJet framework