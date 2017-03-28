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
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include <iostream>
#include <math.h>
#include <iostream>
/**
 * Base class for modules to inherit from.
 */
class ntuple_content{
public:
	ntuple_content():vertices_(0),read_(false){}
	virtual ~ntuple_content();

	virtual void getInput(const edm::ParameterSet& iConfig){}
	virtual void initBranches(TTree* )=0;
	virtual void readEvent(const edm::Event& iEvent)=0;

	//use either of these functions

	virtual bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0)=0;

	void setPrimaryVertices(const reco::VertexCollection* v){
		vertices_=v;
	}

	void setIsRead(bool isread){read_=isread;}

protected:
	const reco::VertexCollection * vertices()const;


	template <class T>
	void addBranch(TTree* t, const char* name,  T*, const char* leaflist=0);


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
	bool read_;
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

}




#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_ */
