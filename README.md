# DeepNTuples
NTuple framework for DeepFlavour

Installation (CMSSW 8_0_25)
============

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src/
cmsenv
git cms-init
git clone https://github.com/CMSDeepFlavour/DeepNTuples
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
#compile
scram b -j 4
```
Installation (CMSSW 8_1_X)
============

```
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src/
cmsenv
git cms-init
# Add DeepFlavour -- To be updated once the 80X PR is done
git cms-merge-topic -u cms-btv-pog:DeepFlavour-from-CMSSW_8_1_0
git clone https://github.com/CMSDeepFlavour/DeepNTuples
# Add JetToolBox
cd DeepNTuples
git submodule init
git submodule update

#compile
scram b -j 4
```

Installation (CMSSW 8_4_X and 9_0_X)
============

```
cmsrel CMSSW_8_0_25
cd CMSSW_8_0_25/src/
cmsenv
git cms-init
git clone https://github.com/CMSDeepFlavour/DeepNTuples
# Add JetToolBox
cd DeepNTuples
git submodule init
git submodule update

#DeepCSV is already in the release, but with different names, which will become the defaults in the close future
sed -i 's|deepFlavourJetTags|pfDeepCSVJetTags|g' DeepNTuples/DeepNtuplizer/production/DeepNtuplizer.py
#compile
scram b -j 4
```

Further settings
============

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
The large job output (root files) will NOT be stored in the batch directory. A directory in your cernbox eos space will be created for those output files. The batch directory will contain a symlink to this directory.

The status of the jobs can be checked with
```
cd <batch directory>
check.py <sample subdirectories to be checked>
```

The check.py script provides additional options to resubmit failed jobs or to create sample lists in case a satisfying fraction of jobs ended successfully. 
In this case do:
```
check.py <sample subdirectories to be checked> --action filelist
```
This will create file lists that can be further processed by the DeepJet framework
For resubmitting failed jobs, do:
```
check.py <sample subdirectories to be checked> --action resubmit
```

When the file lists are created, the part used for training of the ttbar and QCD samples (or in principle any other process) can be merged using the executable:
```
mergeSamples <no of jets per file> <output dir> <file lists 1> <file lists 2> <file lists 3> ...
```
Note that there should be at least two separate file lists. One for ttbar and one for QCD. 
Otherwise the jets are not randomized properly.

For example:
```
mergeSamples 400000 merged ntuple_ttbar*/train_val_samples.txt ntuple_qcd*/train_val_samples.txt

```
This will take a significant amount of time - likely more than the ntuple production itself. It is therefore recommended to run the command within 'screen'.

