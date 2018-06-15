/*
 * treeReader.h
 *
 *  Created on: 17 May 2018
 *      Author: dwalter
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_TRIGGERINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_TRIGGERINFO_H_


#include <map>
#include <vector>
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Framework/interface/Event.h"

class TriggerInfo {

public:
  TriggerInfo(const edm::Event& iEvent,
	      const edm::EDGetTokenT<edm::TriggerResults>& triggerBitsToken);
  bool IsTriggered(std::string triggername) const ;
  bool Exists(std::string triggername) const ;
  bool IsAnyTriggered(std::vector< std::string > triggers) const ;
  bool IsAnyTriggered(std::string triggerstring) const;
  void ListTriggers() ;
  std::map<std::string, bool> GetTriggers() const;
private:
  std::map<std::string, bool> triggers;
  std::vector<std::string> triggernames;
};

#endif
