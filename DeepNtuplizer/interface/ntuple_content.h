/*
 * ntuple_content.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TTree.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <iostream>
#include <math.h>
#include <iostream>
#include "TString.h"
/**
 * Base class for modules to inherit from.
 */
class ntuple_content{
public:
    ntuple_content():ntuple_content(0.4) {}
    ntuple_content(double jetR):vertices_(0),jetR_(jetR),pupInfo_(0),rhoInfo_(0),read_(false){}
    virtual ~ntuple_content();

    virtual void getInput(const edm::ParameterSet& iConfig){}
    virtual void initBranches(TTree* )=0;
    virtual void readEvent(const edm::Event& iEvent)=0;
    virtual void readSetup(const edm::EventSetup& iSetup){}
    //use either of these functions

    virtual bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0)=0;

    void setPrimaryVertices(const reco::VertexCollection* v){
        vertices_=v;
    }
    void setPuInfo(const std::vector<PileupSummaryInfo> *v){
	pupInfo_ =v;
    }
    void setRhoInfo(const double *v){
        rhoInfo_ =v;
    }

    void setIsRead(bool isread){read_=isread;}

    std::vector<TString> getListOfBranches(){
        if(allbranches_.size())
            return allbranches_;
        else{
            TTree *t=new TTree();
            initBranches(t);
            return allbranches_;
        }
    }

protected:
    const reco::VertexCollection * vertices()const;
    const double* rhoInfo()const;
    const std::vector<PileupSummaryInfo> * pupInfo()const;


    template <class T>
    void addBranch(TTree* t, const char* name,  T*, const char* leaflist=0);

    double jetR() const { return jetR_; }

    static inline const float& catchInfs(const float& in,const float& replace_value){
        if(in==in){
            if(std::isinf(in))
                return replace_value;
            else if(in < -1e32 || in > 1e32)
                return replace_value;
            return in;
        }
        return replace_value;
    }

    static inline float catchInfsAndBound(const float& in,const float& replace_value, const float& lowerbound, const float& upperbound){
        float withoutinfs=catchInfs(in,replace_value);
        if(withoutinfs<lowerbound) return lowerbound;
        if(withoutinfs>upperbound) return upperbound;
        return withoutinfs;
    }

private:
    const reco::VertexCollection* vertices_;
    double jetR_;
    const std::vector<PileupSummaryInfo> * pupInfo_;
    const double* rhoInfo_;
    bool read_;
    std::vector<TString> allbranches_;
};

template <class T>
void ntuple_content::addBranch(TTree* t, const char* name,  T* address, const char* leaflist){

    if(read_ ){
        t->SetBranchAddress(name,address);
    }
    else{
        if(leaflist)
            t->Branch(name  ,address  ,leaflist );
        else
            t->Branch(name  ,address);
    }
    allbranches_.push_back((TString)name);

}




#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_ */
