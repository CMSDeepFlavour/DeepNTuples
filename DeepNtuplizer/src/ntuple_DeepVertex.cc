/*
 * ntuple_DeepVertex.cc
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat
 */


#include "../interface/ntuple_DeepVertex.h"

#include "DataFormats/GeometrySurface/interface/Line.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"





ntuple_DeepVertex::ntuple_DeepVertex(double jetR):ntuple_content(jetR){
}
ntuple_DeepVertex::~ntuple_DeepVertex(){}


void ntuple_DeepVertex::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_DeepVertex::initBranches(TTree* tree){
    
    addBranch(tree,"n_seeds",&n_seeds, "n_seeds/i");
    addBranch(tree,"nSeeds",&nSeeds, "nSeeds/f");

    addBranch(tree,"seed_pt",&seed_pt, "seed_pt[n_seeds]/f");
    addBranch(tree,"seed_eta",&seed_eta, "seed_eta[n_seeds]/f");
    addBranch(tree,"seed_phi",&seed_phi, "seed_phi[n_seeds]/f");
    addBranch(tree,"seed_mass",&seed_mass, "seed_mass[n_seeds]/f");
    
    addBranch(tree,"seed_dz", &seed_dz, "seed_dz[n_seeds]/f");
    addBranch(tree,"seed_dxy", &seed_dxy, "seed_dxy[n_seeds]/f");
    addBranch(tree,"seed_3D_ip", &seed_3D_ip, "seed_3D_ip[n_seeds]/f");
    addBranch(tree,"seed_3D_sip", &seed_3D_sip, "seed_3D_sip[n_seeds]/f");
    addBranch(tree,"seed_2D_ip", &seed_2D_ip, "seed_2D_ip[n_seeds]/f");
    addBranch(tree,"seed_2D_sip", &seed_2D_sip, "seed_2D_sip[n_seeds]/f");
    
    addBranch(tree,"seed_3D_signedIp", &seed_3D_signedIp, "seed_3D_signedIp[n_seeds]/f");
    addBranch(tree,"seed_3D_signedSip", &seed_3D_signedSip, "seed_3D_signedSip[n_seeds]/f");
    addBranch(tree,"seed_2D_signedIp", &seed_2D_signedIp, "seed_2D_signedIp[n_seeds]/f");
    addBranch(tree,"seed_2D_signedSip", &seed_2D_signedSip, "seed_2D_signedSip[n_seeds]/f");
    
    addBranch(tree,"seed_chi2reduced",&seed_chi2reduced, "seed_chi2reduced[n_seeds]/f");
    addBranch(tree,"seed_nPixelHits",&seed_nPixelHits, "seed_nPixelHits[n_seeds]/f");
    addBranch(tree,"seed_nHits",&seed_nHits, "seed_nHits[n_seeds]/f");
    addBranch(tree,"seed_jetAxisDistance",&seed_jetAxisDistance, "seed_jetAxisDistance[n_seeds]/f");
    addBranch(tree,"seed_jetAxisDlength",&seed_jetAxisDlength, "seed_jetAxisDlength[n_seeds]/f");
    
    addBranch(tree,"seed_n_NearTracks",&seed_n_NearTracks, "seed_n_NearTracks[n_seeds]/i");
    addBranch(tree,"seed_nNearTracks",&seed_nNearTracks, "seed_nNearTracks[n_seeds]/f");
    
    
    
    // near Tracks
    
    addBranch(tree,"n_NearTracksTotal",&n_NearTracksTotal, "n_NearTracksTotal/i");
    addBranch(tree,"nNearTracksTotal",&nNearTracksTotal, "nNearTracksTotal/f");
    
    addBranch(tree,"nearTracks_pt", &nearTracks_pt, "nearTracks_pt[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_eta", &nearTracks_eta, "nearTracks_eta[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_phi", &nearTracks_phi, "nearTracks_phi[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_mass", &nearTracks_mass, "nearTracks_mass[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dz", &nearTracks_dz, "nearTracks_dz[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dxy", &nearTracks_dxy, "nearTracks_dxy[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_3D_ip", &nearTracks_3D_ip, "nearTracks_3D_ip[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_3D_sip", &nearTracks_3D_sip, "nearTracks_3D_sip[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_2D_ip", &nearTracks_2D_ip, "nearTracks_2D_ip[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_2D_sip", &nearTracks_2D_sip, "nearTracks_2D_sip[n_NearTracksTotal]/f");

    addBranch(tree,"nearTracks_PCAdist", &nearTracks_PCAdist, "nearTracks_PCAdist[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAdsig", &nearTracks_PCAdsig, "nearTracks_PCAdsig[n_NearTracksTotal]/f");
    
    addBranch(tree,"nearTracks_PCAonSeed_x", &nearTracks_PCAonSeed_x, "nearTracks_PCAonSeed_x[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonSeed_y", &nearTracks_PCAonSeed_y, "nearTracks_PCAonSeed_y[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonSeed_z", &nearTracks_PCAonSeed_z, "nearTracks_PCAonSeed_z[n_NearTracksTotal]/f");

    addBranch(tree,"nearTracks_PCAonSeed_xerr", &nearTracks_PCAonSeed_xerr, "nearTracks_PCAonSeed_xerr[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonSeed_yerr", &nearTracks_PCAonSeed_yerr, "nearTracks_PCAonSeed_yerr[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonSeed_zerr", &nearTracks_PCAonSeed_zerr, "nearTracks_PCAonSeed_zerr[n_NearTracksTotal]/f");

    addBranch(tree,"nearTracks_PCAonTrack_x", &nearTracks_PCAonTrack_x, "nearTracks_PCAonTrack_x[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonTrack_y", &nearTracks_PCAonTrack_y, "nearTracks_PCAonTrack_y[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonTrack_z", &nearTracks_PCAonTrack_z, "nearTracks_PCAonTrack_z[n_NearTracksTotal]/f");

    addBranch(tree,"nearTracks_PCAonTrack_xerr", &nearTracks_PCAonTrack_xerr, "nearTracks_PCAonTrack_xerr[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonTrack_yerr", &nearTracks_PCAonTrack_yerr, "nearTracks_PCAonTrack_yerr[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonTrack_zerr", &nearTracks_PCAonTrack_zerr, "nearTracks_PCAonTrack_zerr[n_NearTracksTotal]/f"); 

    addBranch(tree,"nearTracks_dotprodTrack", &nearTracks_dotprodTrack, "nearTracks_dotprodTrack[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dotprodSeed", &nearTracks_dotprodSeed, "nearTracks_dotprodSeed[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dotprodTrackSeed2D", &nearTracks_dotprodTrackSeed2D, "nearTracks_dotprodTrackSeed2D[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dotprodTrackSeed3D", &nearTracks_dotprodTrackSeed3D, "nearTracks_dotprodTrackSeed3D[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dotprodTrackSeedVectors2D", &nearTracks_dotprodTrackSeedVectors2D, "nearTracks_dotprodTrackSeedVectors2D[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_dotprodTrackSeedVectors3D", &nearTracks_dotprodTrackSeedVectors3D, "nearTracks_dotprodTrackSeedVectors3D[n_NearTracksTotal]/f");
    
    addBranch(tree,"nearTracks_PCAonSeed_pvd", &nearTracks_PCAonSeed_pvd, "nearTracks_PCAonSeed_pvd[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAonTrack_pvd", &nearTracks_PCAonTrack_pvd, "nearTracks_PCAonTrack_pvd[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAjetAxis_dist",&nearTracks_PCAjetAxis_dist,"nearTracks_PCAjetAxis_dist[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAjetMomenta_dotprod",&nearTracks_PCAjetMomenta_dotprod,"nearTracks_PCAjetMomenta_dotprod[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAjetDirs_DEta",&nearTracks_PCAjetDirs_DEta,"nearTracks_PCAjetDirs_DEta[n_NearTracksTotal]/f");
    addBranch(tree,"nearTracks_PCAjetDirs_DPhi",&nearTracks_PCAjetDirs_DPhi,"nearTracks_PCAjetDirs_DPhi[n_NearTracksTotal]/f");


}


void ntuple_DeepVertex::readEvent(const edm::Event& iEvent){

    iEvent.getByToken(CandidateToken, tracks);

}


void ntuple_DeepVertex::readSetup(const edm::EventSetup& iSetup){

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}




bool ntuple_DeepVertex::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    // pv info
    const reco::Vertex &pv = vertices()->at(0);
    GlobalPoint pvp(pv.x(),pv.y(),pv.z());

    
    std::vector<reco::TransientTrack> selectedTracks;
    std::vector<float> masses;
    
    
   for(size_t k = 0; k<tracks->size(); ++k) {
        if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5 && std::fabs(pvp.z()-builder->build(tracks->ptrAt(k)).track().vz())<0.5) {
            selectedTracks.push_back(builder->build(tracks->ptrAt(k)));
            masses.push_back(tracks->ptrAt(k)->mass());
        }
    }
    
    double jet_radius = jetR();
    GlobalVector direction(jet.px(), jet.py(), jet.pz());
    
    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){

        //is the track in the jet cone?
        float angular_distance=std::sqrt(std::pow(jet.eta()-it->track().eta(),2) + std::pow(jet.phi()-it->track().phi(),2) );
        if (angular_distance>jet_radius) { continue; }
        
        // is it a seed track?
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*it, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*it, pv);
		std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
		float length=999;
        if (closest.isValid()) length=(closest.globalPosition() - pvp).mag();
        
        // shouldn't it be like this, including the minimal 3DIP cuts? more conform with IVF! https://github.com/cms-sw/cmssw/blob/09c3fce6626f70fd04223e7dacebf0b485f73f54/RecoVertex/AdaptiveVertexFinder/src/TracksClusteringFromDisplacedSeed.cc#L96 
        // bool is_seed_candidate = (ip.first && ip.second.value() >= min3DIPValue && ip.second.significance() >= min3DIPSignificance &&
//             ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance &&
//             it->track().normalizedChi2()<5. && std::fabs(it->track().dxy(pv.position())) < 2 &&
//             std::fabs(it->track().dz(pv.position())) < 17  && jet_dist.second.value()<0.07 && length<5. );

        // this is what is in the DeepVertex code: https://github.com/leonardogiannini/Analyzer/blob/master/plugins/AnalyzerSignedIP_MINIAOD_wArr.cc#L886-L889
        bool is_seed_candidate = (ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 &&
            ip.second.value() <= max3DIPValue && ip.second.significance() <= max3DIPSignificance &&
            it->track().normalizedChi2()<5. && std::fabs(it->track().dxy(pv.position())) < 2 &&
            std::fabs(it->track().dz(pv.position())) < 17  && jet_dist.second.value()<0.07 && length<5. );
        
        if (!is_seed_candidate){continue;}
        
        std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*it,direction, pv); 
        //n_seeds++;
        
        nearTracks.clear();
        //now that we found a seed, loop over all other tracks and look for neighbours
        for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin();tt!=selectedTracks.end(); ++tt ) {
            VertexDistance3D distanceComputer;
            TwoTrackMinimumDistance dist;
            if(*tt==*it) continue;
            if(std::fabs(pvp.z()-tt->track().vz())>0.1) continue;
            if(dist.calculate(tt->impactPointState(),it->impactPointState())) {
                GlobalPoint ttPoint          = dist.points().first;
                GlobalError ttPointErr       = tt->impactPointState().cartesianError().position();
                GlobalPoint seedPosition     = dist.points().second;
                GlobalError seedPositionErr  = it->impactPointState().cartesianError().position();
                Measurement1D m = distanceComputer.distance(VertexState(seedPosition,seedPositionErr), VertexState(ttPoint, ttPointErr));
                GlobalPoint cp(dist.crossingPoint()); 
                
                GlobalVector PairMomentum(it->track().px()+tt->track().px(), it->track().py()+tt->track().py(), it->track().pz()+tt->track().pz());
                GlobalVector  PCA_pv(cp-pvp);

                float PCAseedFromPV =  (dist.points().second-pvp).mag();
                float PCAtrackFromPV =  (dist.points().first-pvp).mag();               
                float distance = dist.distance();

                GlobalVector trackDir2D(tt->impactPointState().globalDirection().x(),tt->impactPointState().globalDirection().y(),0.); 
                GlobalVector seedDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.); 
                GlobalVector trackPCADir2D(dist.points().first.x()-pvp.x(),dist.points().first.y()-pvp.y(),0.); 
                GlobalVector seedPCADir2D(dist.points().second.x()-pvp.x(),dist.points().second.y()-pvp.y(),0.); 
                
                float dotprodTrack = (dist.points().first-pvp).unit().dot(tt->impactPointState().globalDirection().unit());
                float dotprodSeed = (dist.points().second-pvp).unit().dot(it->impactPointState().globalDirection().unit());                    
                float dotprodTrackSeed2D = trackDir2D.unit().dot(seedDir2D.unit());
                float dotprodTrackSeed3D = it->impactPointState().globalDirection().unit().dot(tt->impactPointState().globalDirection().unit());
                float dotprodTrackSeed2DV = trackPCADir2D.unit().dot(seedPCADir2D.unit());
                float dotprodTrackSeed3DV = (dist.points().second-pvp).unit().dot((dist.points().first-pvp).unit());

                std::pair<bool,Measurement1D> t_ip = IPTools::absoluteImpactParameter3D(*tt,pv);        
                std::pair<bool,Measurement1D> t_ip2d = IPTools::absoluteTransverseImpactParameter(*tt,pv);

                myTrack.set_values(tt->track().pt(), tt->track().eta(), tt->track().phi(),  tt->track().dz(pv.position()), tt->track().dxy(pv.position()), distance,  m.significance(), seedPosition.x(), seedPosition.y(), seedPosition.z(), seedPositionErr.cxx(), seedPositionErr.cyy(), seedPositionErr.czz(),  ttPoint.x(),  ttPoint.y(),  ttPoint.z(),  ttPointErr.cxx(),  ttPointErr.cyy(),  ttPointErr.czz(), dotprodTrack, dotprodSeed );
                myTrack.set_index(-1);
                myTrack.set_distances(PCAseedFromPV, PCAtrackFromPV);
                myTrack.set_vars(masses[tt-selectedTracks.begin()],t_ip2d.second.value() , t_ip2d.second.significance(), t_ip.second.value() , t_ip.second.significance(), dotprodTrackSeed2D, dotprodTrackSeed3D, dotprodTrackSeed2DV, dotprodTrackSeed3DV ); 
                
                Line::PositionType pos(pvp);
                Line::DirectionType dir(direction);
                Line::DirectionType pairMomentumDir(PairMomentum);
                Line jetLine(pos,dir);   
                Line PCAMomentumLine(cp,pairMomentumDir);
                float PCA_JetAxis_dist=jetLine.distance(cp).mag();
                float dotprodMomenta=PairMomentum.unit().dot(direction.unit());
                float dEta=std::fabs(PCA_pv.eta()-jet.eta());
                float dPhi=std::fabs(PCA_pv.phi()-jet.phi());

                myTrack.setSeedMass(masses[it-selectedTracks.begin()]);                    
                myTrack.set_JetAxisVars(PCA_JetAxis_dist,dotprodMomenta,dEta,dPhi);
                nearTracks.push_back(myTrack);
            
            }
        }            
         
        std::sort (nearTracks.begin(), nearTracks.end(), sortfunction2());
        if (nearTracks.size() > 20){nearTracks.resize(20);}
        SortedSeedsMap.insert(std::make_pair(-ipSigned.second.significance(), std::make_pair(&(*it), nearTracks)));
            
    }
    
       
    unsigned int seeds_max_counter=0;
    unsigned int neartracks_max_counter=0;
    for(std::multimap<double,std::pair<const reco::TransientTrack*,const std::vector<trackVars2> > >::const_iterator im = SortedSeedsMap.begin(); im != SortedSeedsMap.end(); im++){
        
        if(seeds_max_counter>=10) {
            //n_seeds = 10;
            break;
        }
        
        std::pair<bool,Measurement1D> ipSigned = IPTools::signedImpactParameter3D(*im->second.first,direction, pv);        
        std::pair<bool,Measurement1D> ip2dSigned = IPTools::signedTransverseImpactParameter(*im->second.first,direction, pv);  
        std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*im->second.first, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::absoluteTransverseImpactParameter(*im->second.first, pv);	
        
        seed_pt[seeds_max_counter]=im->second.first->track().pt();
        seed_eta[seeds_max_counter]=im->second.first->track().eta();
        seed_phi[seeds_max_counter]=im->second.first->track().phi();
        seed_mass[seeds_max_counter]=im->second.second.at(0).seedMass;
        seed_dz[seeds_max_counter]=im->second.first->track().dz(pv.position());
        seed_dxy[seeds_max_counter]=im->second.first->track().dxy(pv.position());
        seed_3D_ip[seeds_max_counter]=ip.second.value();
        seed_3D_sip[seeds_max_counter]=ip.second.significance();
        seed_2D_ip[seeds_max_counter]=ip2d.second.value();
        seed_2D_sip[seeds_max_counter]=ip2d.second.significance();
        seed_3D_signedIp[seeds_max_counter]=ipSigned.second.value();
        seed_3D_signedSip[seeds_max_counter]=ipSigned.second.significance();
        seed_2D_signedIp[seeds_max_counter]=ip2dSigned.second.value();
        seed_2D_signedSip[seeds_max_counter]=ip2dSigned.second.significance();		
        seed_chi2reduced[seeds_max_counter]=im->second.first->track().normalizedChi2();
        seed_nPixelHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidPixelHits();
        seed_nHits[seeds_max_counter]=im->second.first->track().hitPattern().numberOfValidHits();

        std::pair<double, Measurement1D> jet_distance =IPTools::jetTrackDistance(*im->second.first, direction, pv);
        seed_jetAxisDistance[seeds_max_counter]=std::fabs(jet_distance.second.value());

        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(im->second.first->impactPointState(),pv, direction,im->second.first->field());
        if (closest.isValid()) seed_jetAxisDlength[seeds_max_counter]=(closest.globalPosition() - pvp).mag(); 
        else seed_jetAxisDlength[seeds_max_counter]= -99;
        
        seed_n_NearTracks[seeds_max_counter]=im->second.second.size();
        seed_nNearTracks[seeds_max_counter]=im->second.second.size();
        
        // FILL NEAREAST VARIABLES
        for(unsigned int i=0; i< im->second.second.size(); i++) {

            if((neartracks_max_counter+i)>=200) break;
            
			nearTracks_pt[neartracks_max_counter+i]=im->second.second.at(i).pt;
            nearTracks_eta[neartracks_max_counter+i]=im->second.second.at(i).eta;
            nearTracks_phi[neartracks_max_counter+i]=im->second.second.at(i).phi;
            nearTracks_dz[neartracks_max_counter+i]=im->second.second.at(i).dz;
            nearTracks_dxy[neartracks_max_counter+i]=im->second.second.at(i).dxy;
            nearTracks_mass[neartracks_max_counter+i]=im->second.second.at(i).mass;
            nearTracks_3D_ip[neartracks_max_counter+i]=im->second.second.at(i).t3Dip;
            nearTracks_3D_sip[neartracks_max_counter+i]=im->second.second.at(i).t3Dsip;
            nearTracks_2D_ip[neartracks_max_counter+i]=im->second.second.at(i).t2Dip;
            nearTracks_2D_sip[neartracks_max_counter+i]=im->second.second.at(i).t2Dsip;
            nearTracks_PCAdist[neartracks_max_counter+i]=im->second.second.at(i).dist;
            nearTracks_PCAdsig[neartracks_max_counter+i]=im->second.second.at(i).dsig;
            nearTracks_PCAonSeed_x[neartracks_max_counter+i]=im->second.second.at(i).PCA_sx;
            nearTracks_PCAonSeed_y[neartracks_max_counter+i]=im->second.second.at(i).PCA_sy;
            nearTracks_PCAonSeed_z[neartracks_max_counter+i]=im->second.second.at(i).PCA_sz;
            nearTracks_PCAonSeed_xerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_sxerr;
            nearTracks_PCAonSeed_yerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_syerr;
            nearTracks_PCAonSeed_zerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_szerr;
            nearTracks_PCAonTrack_x[neartracks_max_counter+i]=im->second.second.at(i).PCA_tx;
            nearTracks_PCAonTrack_y[neartracks_max_counter+i]=im->second.second.at(i).PCA_ty;
            nearTracks_PCAonTrack_z[neartracks_max_counter+i]=im->second.second.at(i).PCA_tz;
            nearTracks_PCAonTrack_xerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_txerr;
            nearTracks_PCAonTrack_yerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_tyerr;
            nearTracks_PCAonTrack_zerr[neartracks_max_counter+i]=im->second.second.at(i).PCA_tzerr;
            nearTracks_dotprodTrack[neartracks_max_counter+i]=im->second.second.at(i).dotprodTrack;
            nearTracks_dotprodSeed[neartracks_max_counter+i]=im->second.second.at(i).dotprodSeed;
            nearTracks_dotprodTrackSeed2D[neartracks_max_counter+i]=im->second.second.at(i).dotprodTrackSeed2D;
            nearTracks_dotprodTrackSeed3D[neartracks_max_counter+i]=im->second.second.at(i).dotprodTrackSeed3D;
            nearTracks_dotprodTrackSeedVectors2D[neartracks_max_counter+i]=im->second.second.at(i).dotprodTrackSeedVectors2D;
            nearTracks_dotprodTrackSeedVectors3D[neartracks_max_counter+i]=im->second.second.at(i).dotprodTrackSeedVectors3D;

            nearTracks_PCAonSeed_pvd[neartracks_max_counter+i]=im->second.second.at(i).seedPCA_pv;
            nearTracks_PCAonTrack_pvd[neartracks_max_counter+i]=im->second.second.at(i).trackPCA_pv;

            nearTracks_PCAjetAxis_dist[neartracks_max_counter+i]=im->second.second.at(i).PCA_JetAxis_distance;
            nearTracks_PCAjetMomenta_dotprod[neartracks_max_counter+i]=im->second.second.at(i).PCAPair_Jet_dotprod;

            nearTracks_PCAjetDirs_DEta[neartracks_max_counter+i]=im->second.second.at(i).PCAAxis_JetAxis_DEta;
            nearTracks_PCAjetDirs_DPhi[neartracks_max_counter+i]=im->second.second.at(i).PCAAxis_JetAxis_DPhi;

        }
        
        // *********
        neartracks_max_counter += im->second.second.size();
        seeds_max_counter++; 
    }
    n_NearTracksTotal = neartracks_max_counter;
    nNearTracksTotal = neartracks_max_counter;
    n_seeds = seeds_max_counter;
    nSeeds = seeds_max_counter;
    
    SortedSeedsMap.clear();
    nearTracks.clear();
    masses.clear();

    return true;
}

