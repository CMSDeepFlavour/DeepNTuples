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


void setPreCache(TChain* tree){
	tree->SetCacheSize(100e6);//100MB precache (eos is slow)
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


int main(int argc, char *argv[]){

	std::vector<ntuple_content*> branchinfos;
	branchinfos.push_back(new ntuple_JetInfo());
	branchinfos.push_back(new ntuple_SV());
	branchinfos.push_back(new ntuple_bTagVars());
	branchinfos.push_back(new ntuple_pfCands());



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
	system(("mkdir -p "+outpath).Data());

	std::vector<std::vector<TString> > infiles;

	for(int i=3;i<argc;i++){
		TString samplefile=argv[i];
		TString inpath=dirname(argv[i]);
		inpath+="/";
		infiles.push_back(readSampleFile(samplefile,inpath));
	}

	// input files read, create TChains

	std::vector<TChain* > chains;
	for(size_t i=0;i<infiles.size();i++)
		chains.push_back(new TChain());

	//std::cout << "setting up trees... \nPlease ignore wrong branch type errors - they do no harm." <<std::endl;
	//read in files
	std::vector<size_t> entriesperchain(chains.size());
	size_t totalentries=0;
	for(size_t i=0;i<infiles.size();i++){
		for(const auto& f:infiles.at(i))
			chains.at(i)->Add(f+"/deepntuplizer/tree");
		for(auto& bi:branchinfos){
			bi->setIsRead(true);
			bi->initBranches(chains.at(i));
		}
		entriesperchain.at(i) = chains.at(i)->GetEntries();
		setPreCache(chains.at(i));
		totalentries+=entriesperchain.at(i);
	}

	//create "binning"
	std::vector<double> fractions(chains.size());
	double alloldfracs=0;
	for(size_t i=0;i<fractions.size();i++){
		fractions.at(i) = (double)entriesperchain.at(i)/(double)totalentries + alloldfracs;
		alloldfracs=fractions.at(i);
	}

	std::cout << "entry fractions: " ;
	for(const auto& e:fractions)
		std::cout << e << " ";
	std::cout <<std::endl;


	TRandom3 rand;

	bool lastwrite=false;
	size_t outfileindex=0;
	size_t outtreeentry=0;
	size_t alloutentries=0;
	std::vector<size_t> entryindices(chains.size(),0);
	TFile * fout= 0;
	TTree * outtree =0;
	bool newfile=true;
	std::ofstream outlist((outpath+"/samples.txt").Data());


	//control histogram
	std::vector<double> histobins(1,0);
	histobins.insert(histobins.end(),fractions.begin(),fractions.end());
	TH1D hist("samplefraction","samplefraction",histobins.size()-1,&histobins.at(0));
	std::vector<TH1D> contr_hists;
	for(size_t i=0;i<chains.size();i++){
		TString basen="";
		TH1D hcontr(basen,basen,500,0,totalentries);
		contr_hists.push_back(hcontr);
	}

	std::cout << "looping" <<std::endl;
	time_t now;
	time_t started;
	time(&started);
	time(&now);

	while(1){

		if(newfile || lastwrite){
			if(fout){
				outtree->Write();
				fout->Close();
			}
			if(lastwrite)
				break;
			TString idxstring=""; //for testing
			idxstring+=outfileindex;
			TString outfilename="ntuple_merged_"+idxstring+".root";
			//goes in the indexed loop later
			outlist << outfilename<<"\n";
			std::cout << "new output file " <<outfilename <<std::endl;

			time(&now);
			float runningsecond=difftime(now,started);
			float processedfraction=(double)alloutentries/(double)totalentries;
			float expectedtotal=runningsecond/processedfraction;
			int eta=expectedtotal-runningsecond;
			int min=0;
			if(eta>120)
				min=eta/60;

			if(outfileindex){
				if(min)
					std::cout << "total progress: " << (int)(processedfraction * 100) << "%, ETA " << min << " min"<<std::endl;
				else
					std::cout << "total progress: " << (int)(processedfraction * 100) << "%, ETA " << eta << " sec"<<std::endl;

			}
			fout= new TFile(outpath+"/"+outfilename,"RECREATE");
			TDirectory* tdir=fout->mkdir("deepntuplizer","deepntuplizer");
			tdir->cd();
			outtree = chains.at(0)->CloneTree(0);//new TTree("tree","tree");//
			outtreeentry=0;
			outfileindex++;
		}

		//associate the next entry randomly to one of the sources
		size_t treeindex=std::lower_bound(fractions.begin(), fractions.end(), rand.Uniform(0,1))-fractions.begin();

		chains.at(treeindex)->GetEntry(entryindices.at(treeindex));
		entryindices.at(treeindex)++;

		outtree->Fill();

		contr_hists.at(treeindex).Fill(alloutentries);

		outtreeentry++;
		alloutentries++;


		newfile= outtreeentry>=entriesperfile;

		hist.Fill(histobins.at(treeindex));

		if(entryindices.at(treeindex)>=entriesperchain.at(treeindex))
			lastwrite=true;

		continue;
		//testing
		if(alloutentries>=100000)
			lastwrite=true;

	}

	std::cout << "total progress: 100%"<<std::endl;
	std::cout << "merging done. Discarded " << totalentries-alloutentries << " of " << totalentries << "entries" <<std::endl;

	for(size_t i=0;i<entriesperchain.size();i++)
		std::cout << "used " << entryindices.at(i) << " out of " << entriesperchain.at(i) << " entries for chain " <<i <<std::endl;

	outlist << '\n';
	outlist.close();

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

	return 0;
}


















