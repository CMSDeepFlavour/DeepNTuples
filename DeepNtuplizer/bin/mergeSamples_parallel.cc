/*
 * mergeSamples.cc
 *
 *  Created on: 23 Feb 2017
 *      Author: jkiesele
 *
 *
 *  Executable to merge sample collections, e.g. QCD and ttbar
 *
 *  Notes:
 *    - should work on samples.txt outputs from check.py
 *    - other inputs: output directory, number of samples to be merged to (if not given, use first nsamples)
 *    - writes out root files and at the end a samples.txt file again
 *    - should work with an unlimited number of input collections
 *
 */

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include "TRandom3.h"
#include <libgen.h>
#include <string>
#include <fstream>
#include "TDirectory.h"
#include <libgen.h>
#include "TStyle.h"
#include "TLine.h"

//for control plots
#include "TH1D.h"
#include "TCanvas.h"
#include <ctime>

#include "DeepNTuples/DeepNtuplizer/interface/ntuple_bTagVars.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_JetInfo.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_pfCands.h"
#include "DeepNTuples/DeepNtuplizer/interface/ntuple_SV.h"

#include <dirent.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include "TROOT.h"

bool debug=false;


void setPreCache(TChain* tree){
    return ; //don't do anything for now
    tree->SetCacheSize(100e6);//100MB precache (eos is slow) - but increases CPU a lot...
}

inline bool FileExists (const std::string& name) {
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


std::vector<TString> readSampleFile(const TString& file, const TString& addpath){
    std::string line;
    std::vector<TString>  out;
    std::ifstream myfile (file.Data());
    if (myfile.is_open()){
        while ( getline (myfile,line) ){
            out.push_back(addpath+"/"+(TString)line);
        }
        myfile.close();
    }
    else
        throw std::runtime_error("could not read file");
    return out;
}

std::vector<TChain* > createChains(const std::vector<std::vector<TString> >& infiles,
        std::vector<ntuple_content*>& branchinfos,
        std::vector<size_t>& entriesperchain,
        size_t& totalentries){

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
            TString xrootdedpath=f;//prependXRootD(f);
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


void prepareSplitting(const std::vector<std::vector<TString> >& infiles,
        const size_t entriesperfile,const TString& outpath,
        std::vector<double>& fractions, std::vector<std::vector<size_t> >& whichchain_perfile,
        std::vector<std::vector<size_t> >&startentries){

    // input files read, create TChains
    std::vector<ntuple_content*> branchinfos;
    std::vector<size_t> entriesperchain;
    size_t totalentries=0;
    std::vector<TChain* > chains = createChains(infiles,
            branchinfos,
            entriesperchain,totalentries);

    //std::cout << "setting up trees... \nPlease ignore wrong branch type errors - they do no harm." <<std::endl;
    //read in files




    //create "binning"
    fractions=std::vector<double>(chains.size());
    double alloldfracs=0;
    for(size_t i=0;i<fractions.size();i++){
        fractions.at(i) = (double)entriesperchain.at(i)/(double)totalentries + alloldfracs;
        alloldfracs=fractions.at(i);
    }

    std::cout << "entries per chain: ";

    for(const auto& e:entriesperchain)
        std::cout << e << " ";
    std::cout <<std::endl;
    for(const auto& e:entriesperchain){
        if(!e)
            throw std::runtime_error("at least one sample with zero contribution");
    }
    std::cout << "entry fractions: " ;
    for(const auto& e:fractions)
        std::cout << e << " ";
    std::cout <<std::endl;


    TRandom3 rand;


    std::vector<double> histobins(1,0);
    histobins.insert(histobins.end(),fractions.begin(),fractions.end());
    TH1D hist("samplefraction","samplefraction",histobins.size()-1,&histobins.at(0));
    std::vector<TH1D> contr_hists;
    for(size_t i=0;i<chains.size();i++){
        TString basen="";
        TH1D hcontr(basen,basen,500,0,totalentries);
        contr_hists.push_back(hcontr);
    }


    std::cout << "looping for indices" <<std::endl;
    std::vector<size_t> entryindices(chains.size(),0);

    size_t alloutentries=0;
    whichchain_perfile.clear();
    std::vector<size_t> perfiletmp;
    startentries.push_back(entryindices);
    while(1){

        size_t treeindex=std::lower_bound(fractions.begin(), fractions.end(),
                rand.Uniform(0,1))-fractions.begin();

        perfiletmp.push_back(treeindex);

        hist.Fill(histobins.at(treeindex));
        contr_hists.at(treeindex).Fill(alloutentries);

        entryindices.at(treeindex)++;
        if(entryindices.at(treeindex)>=entriesperchain.at(treeindex)){
            whichchain_perfile.push_back(perfiletmp);
            perfiletmp.clear();
            break;
        }

        if(perfiletmp.size()>=entriesperfile){
            whichchain_perfile.push_back(perfiletmp);
            perfiletmp.clear();
            startentries.push_back(entryindices);
        }
        alloutentries++;
        if(alloutentries>=totalentries){
            throw std::runtime_error("serious problem");
        }
    }
    std::cout << "prepared randomized indices" <<std::endl;


    TH1D fracs=hist;
    fracs.Scale(1./(double)totalentries);
    TCanvas cv;
    gStyle->SetOptStat(0);
    fracs.Draw();
    cv.Print(outpath+"/samplefractions.pdf");

    hist.Scale(1,"width");
    hist.Scale(1./(double)totalentries);
    hist.Draw();
    cv.Print(outpath+"/samplefractions_binwidth.pdf");


    double min=totalentries,max=0;
    for(const auto& h:contr_hists){
        double tmax=h.GetMaximum();
        double tmin=h.GetMinimum();
        if(tmax>max)max=tmax;
        if(tmin<min)min=tmin;
    }

    contr_hists.at(0).GetYaxis()->SetRangeUser(min*0.95,max*1.05);
    contr_hists.at(0).Draw("AXIS");
    contr_hists.at(0).GetYaxis()->SetRangeUser(min*0.95,max*1.05);

    for(size_t i=0;i<contr_hists.size();i++){
        int color=i+1;
        if(i>8)color+=30;
        contr_hists.at(i).SetLineColor(color);
        contr_hists.at(i).Draw("same");
    }
    cv.SetCanvasSize(cv.GetWindowWidth()*5,cv.GetWindowHeight());
    cv.SetLeftMargin(0.02);
    cv.SetRightMargin(0.02);
    double count=0;
    for(size_t i=0;i<totalentries/entriesperfile+1;i++){
        if(count){
            TLine * l = new TLine(count,min*0.95,count,max*1.05);
            l->SetLineColorAlpha(kBlack,0.5);
            l->Draw("same");
        }
        count+=entriesperfile;
    }
    cv.Print(outpath+"/sample_contributions.pdf");

    //launch make the children
    for(auto& c:chains)
        delete c; //clear before forking

    gROOT->Reset();//really annoying root behaviour with all the background globals

}


#include <signal.h>
#include <execinfo.h>
std::string gl_currentfailfile="";
void kill_callback_handler(int signum){
    system(("touch "+gl_currentfailfile).data());
    exit(signum);
}
void segfault_callback_handler(int signum){
    system(("touch "+gl_currentfailfile).data());


    const int lines=100;
    void *array[lines];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, lines);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", signum);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(signum);
}

void fillTree(const std::vector<std::vector<size_t> > whichchain_perfile,
        std::vector<std::vector<TString> > infiles,
        const TString outpath,
        const std::vector<double>& fractions,
        const std::vector<std::vector<size_t> >&startentries,
        const size_t outfileindex

){
    //forked part

    TString idxstring=""; //for testing
    idxstring+=outfileindex;
    gl_currentfailfile=outpath+"/"+idxstring+".fail";
    //make sure the fail path is set before anything else

    TFile * fout= 0;
    TTree * outtree =0;


    std::vector<ntuple_content*> branchinfos;
    std::vector<size_t> entriesperchain;
    size_t totalentries=0;
    std::vector<TChain* > chains = createChains(infiles,
            branchinfos,
            entriesperchain,totalentries);




    size_t outtreeentry=0;
    size_t alloutentries=0;
    std::vector<size_t> entryindices=startentries.at(outfileindex);


    TString outfilename="ntuple_merged_"+idxstring+".root";
    //goes in the indexed loop later
    //std::cout << "new output file " <<outfilename <<std::endl;
    fout= new TFile(outpath+"/"+outfilename,"RECREATE");
    TDirectory* tdir=fout->mkdir("deepntuplizer","deepntuplizer");
    tdir->cd();
    outtree = chains.at(0)->CloneTree(0);//new TTree("tree","tree");//

    if(debug){
        std::cout << "looping " << outfileindex << " starting from: ";
        for(const auto& e:entryindices){
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }

    for(size_t i=0;i<whichchain_perfile.at(outfileindex).size();i++){

        if(debug){
            if(i<5)
                std::cout << whichchain_perfile.at(outfileindex).at(i) <<std::endl;
        }

        //associate the next entry randomly to one of the sources
        size_t treeindex=whichchain_perfile.at(outfileindex).at(i);

        chains.at(treeindex)->GetEntry(entryindices.at(treeindex));
        entryindices.at(treeindex)++;

        outtree->Fill();

        outtreeentry++;
        alloutentries++;

    }
    outtree->Write();
    //check file as ok
    fout->Close();
    if(debug){
        std::cout << idxstring << " file success"<<std::endl;
    }
    system(("touch "+outpath+"/"+idxstring+".succ").Data());
    delete fout;
    exit(0);

}










int main(int argc, char *argv[]){

    //simple opt parsing
    TString helpmessage="\n\
    run this program with:\n \
    <nEntriesPerFile (about 400k is a good choice)> <outputdir> <input sample.txt 0> <input sample.txt 1> <input sample.txt 2> ...\
             ";

    if(argc < 4 || atof(argv[1]) == 0){
        std::cout << helpmessage <<std::endl;
        return -1;
    }
    size_t entriesperfile=atof(argv[1]);
    const TString outpath=argv[2];

    //check if exists
    if(DirectoryExists(outpath.Data())){
        std::cout << "ERROR: output directory must not exists ("<<outpath << ")"<<std::endl;
        return -1;
    }

    std::cout << "WARNING: this executable needs to be monitored and killed by hand if needed. If one of the merging jobs fails, it might otherwise continue forever, since only a fraction of guards is implmented" << std::endl;


    system(("mkdir -p "+outpath).Data());

    std::vector<std::vector<TString> > infiles;

    // unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
    int numCPU = sysconf(_SC_NPROCESSORS_ONLN);

    size_t maxrunning=numCPU;

    for(int i=3;i<argc;i++){
        TString samplefile=argv[i];
        TString inpath=dirname(argv[i]);
        inpath+="/";
        infiles.push_back(readSampleFile(samplefile,inpath));
    }

    std::vector<double> fractions;
    std::vector<std::vector<size_t> > whichchain_perfile;
    std::vector<std::vector<size_t> >startentries;
    prepareSplitting(infiles,
            entriesperfile,outpath,
            fractions,  whichchain_perfile,startentries);


    std::vector<size_t> childpids(whichchain_perfile.size());
    size_t nrunning=0;

    //loop and increase at each fork




    time_t now;
    time_t started;
    time(&started);
    time(&now);

    //wait for them
    enum job_stat{js_succ,js_run,js_wait,js_fail};

    std::vector<job_stat> statuses(childpids.size(),js_wait);
    size_t nextspawn=0;
    size_t seccounter=0;
    while(1){
        for(size_t i=nextspawn;i<childpids.size();i++){

            if(nrunning>=maxrunning)
                break;
            if(debug){
                std::cout << "starting " <<i << std::endl;
            }
            childpids.at(i)=fork();
            if(childpids.at(i)==0){
                signal(SIGSEGV, segfault_callback_handler);
                signal(SIGKILL, kill_callback_handler);
                //child process
                fillTree(whichchain_perfile,
                        infiles,
                        outpath,
                        fractions,startentries,i);
                exit(0);
            }
            else{
                usleep(100000);
                statuses.at(i)=js_run;
            }
            nextspawn++;
            nrunning++;
        }
        size_t nsucc=0;
        size_t nfail=0;
        bool newfilefinished=false;
        for(size_t i=0;i<childpids.size();i++){
            TString idxstring="";
            idxstring+=i;
            if(debug){
                //std::cout << "checking " <<i << std::endl;
            }
            if(FileExists((outpath+"/"+idxstring+".succ").Data())){
                if(statuses.at(i)!=js_succ){
                    nrunning--;
                    newfilefinished=true;
                    if(debug){
                        std::cout << "check success " <<i << std::endl;
                    }
                }
                statuses.at(i)=js_succ;
                nsucc++;
            }
            if(FileExists((outpath+"/"+idxstring+".fail").Data())){
                if(statuses.at(i)!=js_fail){
                    nrunning--;
                    newfilefinished=true;
                    if(debug){
                        std::cout << "check success " <<i << std::endl;
                    }
                }
                statuses.at(i)=js_fail;
                nfail++;
            }

        }
        if(seccounter>10){//just for displaying info
            std::cout << "processed " << nsucc << " out of "<< childpids.size() << " files. ";
            if(nfail)
                std::cout << nfail << " failed. Aborting.";
            std::cout<<std::endl;
            seccounter=0;
            if(nfail)
                exit(-1);
        }
        seccounter++;
        if(nsucc == childpids.size()){
            //all good
            //this needs an additional exit condition in case of job fails
            break;
        }
        if(!newfilefinished)
            sleep(1);
    }

    std::cout << "merging successful, creating sample list.." <<std::endl;

    std::ofstream outlist((outpath+"/samples.txt").Data());
    //all successful
    //check for .succ files
    for(size_t i=0;i<childpids.size();i++){
        TString idxstring="";
        idxstring+=i;
        TString outfilename="ntuple_merged_"+idxstring+".root";
        outlist << outfilename<<"\n";
    }


    outlist << '\n';
    outlist.close();

    return 0;
}








