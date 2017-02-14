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
