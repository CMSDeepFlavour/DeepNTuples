
//#if !defined(__CINT__) && !defined(__MAKECINT__)
//#include "DataFormats/FWLite/interface/Handle.h"
//#include "DataFormats/FWLite/interface/Event.h"
//Headers for the data items
//#include <iostream>
//#include <TFile.h>
//#endif

void make_plots(){
   #include "DataFormats/FWLite/interface/Handle.h"
   std::cout<<" start make_plots "<<std::endl;

   TFile file("/afs/desy.de/user/d/dwalter/data/TT/output1_0_1.root");

   /*fwlite::Event ev(&file);

   for( ev.toBegin(); ! ev.atEnd(); ++ev) {
       fwlite::Handle<vector<pat::Jet> > h_jets;
       h_jets.getByLabel(ev,"slimmedJets");
       // now can access data
       std::cout <<" size "<<h_jets.ptr()->size()<<std::endl;
   }*/
}