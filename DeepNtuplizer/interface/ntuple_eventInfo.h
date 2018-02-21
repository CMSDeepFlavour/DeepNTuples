/*
 * ntuple_eventInfo.h
 *
 *  Created on: 20 Feb 2018
 *      Author: dwalter
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_

#include "ntuple_content.h"

#include <vector>
#include <TObjString.h>


#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

/*
 * For MC weights such as pileup, lhe, ... later: lepton scalefactors
 */
class ntuple_eventInfo: public ntuple_content{
public:
    ntuple_eventInfo():ntuple_content(),
    isData_(false)
{}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);


    void setLHEToken(edm::EDGetTokenT<LHEEventProduct> lheToken) {
        lheToken_ = lheToken;
    }


private:

    edm::EDGetTokenT<LHEEventProduct> lheToken_;

    edm::Handle<LHEEventProduct> lheInfo;

    bool isData_;
    bool useLHEWeights_;

    // global variables

    float ntrueInt_;


    std::vector<double> pupWeights = {0., 0.9059208563049067, 0.4689946658078284, 0.4646384395209148, 0.35969952879169353, 0.20581098696265324, 0.1272975292012222, 0.11484070280016173, 0.18997631713506571, 0.2134442228038988, 0.2643260644534947, 0.38290419689123334, 0.5233231080095909, 0.6531355038202962, 0.7413472652417477, 0.8226986246368242, 0.909166618129074, 0.9754246515040906, 1.0116680456990728, 1.0419003275306946, 1.0762295787102802, 1.1034517536607453, 1.1240057591005752, 1.1262609366258989, 1.1198925605086503, 1.1078563848872773, 1.083911370342871, 1.0496801331833463, 1.0169474777740053, 0.9813073284102061, 0.9426647112682524, 0.9233559974737714, 0.894377685151131, 0.884518640484559, 0.8804996845439541, 0.888865624721422, 0.9107720191678217, 0.953824079957025, 1.0112115792549745, 1.0947266350670766, 1.1906625100966544, 1.3370686699103784, 1.510754710039529, 1.7098870624229212, 1.9391935863130325, 2.1529480239752585, 2.396292845558282, 2.7067102383965915, 2.9081731387455574, 3.126536829527805, 3.263323444233829, 3.3766746819738014, 3.3825183331846596, 3.332122539734733, 3.309580249670872, 2.962530063365449, 2.609107759179973, 2.4216997456337266, 2.228050917243899, 1.963846134469284, 1.844961576497657, 1.431817556288543, 1.5742387193813208, 1.731514171818786, 2.21989485741883, 2.8198818386713915, 3.8909891214054464, 5.4911613943458395, 8.610687700958811, 12.634799153629881, 20.32733145948569};

    /////////branches
    float event_weight_;


};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_EVENTINFO_H_ */
