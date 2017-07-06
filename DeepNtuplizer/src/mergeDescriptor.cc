/*
 * mergeDescriptor.cc
 *
 *  Created on: 22 May 2017
 *      Author: jkiesele
 */




#include "../interface/mergeDescriptor.h"
#include <fstream>


#include <dirent.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include "TROOT.h"
#include <stdio.h>


#include "DeepNTuples/DeepNtuplizer/interface/ntuple_bTagVars.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_JetInfo.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_pfCands.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_SV.h"

static bool debug=true;

TString createTempName(){
    TString tin ="/tmp/mergeParallel_XXXXXXXXX";
    char t[tin.Length()];
    strcpy(t, tin.Data());
    int f=mkstemp(t);
    //std::cout << t << std::endl;
    close(f);
    TString n(t);
    return n;
}

TString prependXRootD(const TString& path){

    return path; //not used
    TString full_path = realpath(path, NULL);
    if(full_path.BeginsWith("/eos/cms/")){
        TString append="root://eoscms.cern.ch//";
        TString s_remove="/eos/cms/";
        TString newpath (full_path(s_remove.Length(),full_path.Length()));
        newpath=append+newpath;
        return newpath;
    }
    return path;
}

void setPreCache(TChain* tree){
    return ; //don't do anything for now
    tree->SetCacheSize(100e6);//100MB precache (eos is slow) - but increases CPU a lot...
}

bool FileExists (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;
    DIR *pDir;
    bool bExists = false;
    pDir = opendir (pzPath);
    if (pDir != NULL){
        bExists = true;
        (void) closedir (pDir);
    }
    return bExists;
}


void mergeDescriptor::writeToFile(std::string filename){
    std::ofstream file(filename);
    serializedWrite(whichchain_perfile,file);
    serializedWrite(infiles,file);
    serializedWrite(outpath,file);
    serializedWrite(fractions,file);
    serializedWrite(startentries,file);
    file.close();
}
void mergeDescriptor::readFromFile(std::string filename, int pickone){
    std::ifstream file(filename);
    if(pickone<0){
        serializedRead(whichchain_perfile,file);
        serializedRead(infiles,file);
        serializedRead(outpath,file);
        serializedRead(fractions,file);
        serializedRead(startentries,file);
    }
    else{
        whichchain_perfile=std::vector<std::vector<size_t> >(1,std::vector<size_t>());
        serializedReadFromVector(whichchain_perfile.at(0),file,(size_t)pickone);

        serializedRead(infiles,file);//not sorted per outfile

        serializedRead(outpath,file);
        serializedRead(fractions,file);

        startentries=std::vector<std::vector<size_t> >(1,std::vector<size_t> ());
        serializedReadFromVector(startentries.at(0),file,(size_t)pickone);
    }
    file.close();
}


std::vector<TChain* > mergeDescriptor::createChains(
        std::vector<size_t>& entriesperchain,
        size_t& totalentries, bool usexrootd){

    static int ntimescalled=0;

    if(debug){
        std::cout << "creating chains" <<std::endl;
    }

    if(branchinfos.size())
        for(auto& b: branchinfos)
            delete b;
    branchinfos.clear();
    entriesperchain=std::vector<size_t>(infiles.size(),0);

    branchinfos.push_back(new ntuple_JetInfo());
    branchinfos.push_back(new ntuple_SV());
    //branchinfos.push_back(new ntuple_SV("LooseIVF_"));
    branchinfos.push_back(new ntuple_bTagVars());
    branchinfos.push_back(new ntuple_pfCands());

    std::vector<TChain* > chains;
    for(size_t i=0;i<infiles.size();i++){
        TString chainname="";
        chainname+=i;
        chainname+="_";
        chainname+=ntimescalled;
        chains.push_back(new TChain(chainname,chainname)); //to get ahead of root background lsiting problems...
    }

    for(size_t i=0;i<infiles.size();i++){
        for(const auto& f:infiles.at(i)){
            TString xrootdedpath=f;
            if(usexrootd)
                xrootdedpath=prependXRootD(xrootdedpath);
            chains.at(i)->Add(xrootdedpath+"/deepntuplizer/tree");
        }
        for(auto& bi:branchinfos){
            bi->setIsRead(true);
            bi->initBranches(chains.at(i));
        }
        entriesperchain.at(i) = chains.at(i)->GetEntries();
        setPreCache(chains.at(i));
        totalentries+=entriesperchain.at(i);
    }
    ntimescalled++;
    return chains;
}
