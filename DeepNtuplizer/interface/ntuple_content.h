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

/**
 * Base class for modules to inherit from.
 */
class ntuple_content{
public:
	ntuple_content():vertices_(0){}
	virtual ~ntuple_content();

	virtual void getInput(const edm::ParameterSet& iConfig){}
	virtual void initBranches(TTree* )=0;
	virtual void readEvent(const edm::Event& iEvent)=0;

	//use either of these functions

	virtual bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0)=0;

	void setPrimaryVertices(const reco::VertexCollection* v){
		vertices_=v;
	}


protected:
	const reco::VertexCollection * vertices()const;

private:
	const reco::VertexCollection* vertices_;
};




#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_CONTENT_H_ */
