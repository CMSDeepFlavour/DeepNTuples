#include "../interface/triggerInfo.h"

using namespace std;
TriggerInfo::TriggerInfo(const edm::Event& iEvent,
			 const edm::EDGetTokenT<edm::TriggerResults>& triggerBitsToken
			 ){

  edm::Handle<edm::TriggerResults> h_triggerBits;
  iEvent.getByToken(triggerBitsToken, h_triggerBits);
  // todo : use trigger objects
  //    edm::Handle<pat::TriggerObjectStandAloneCollection> h_triggerObjects;
  //    iEvent.getByToken(triggerObjectsToken, h_triggerObjects);

  const edm::TriggerNames &names = iEvent.triggerNames(*h_triggerBits);

  for (unsigned int i = 0; i < h_triggerBits->size(); ++i) {
    string name=names.triggerName(i);
//        cout << "Trigger #" << i << "  " << name << endl;
    triggernames.push_back(name);
    triggers[name]=h_triggerBits->accept(i);
  }
}

bool TriggerInfo::Exists(std::string triggername) const {
  if(triggers.count(triggername)==0){
    return false;
  }
  return true;
}


bool TriggerInfo::IsTriggered(std::string triggername) const {
  if(triggername=="any"||triggername=="Any"||triggername=="none"||triggername=="None") return true;
  int asterix = triggername.find("*");
  if(asterix==int(std::string::npos)){
    return Exists(triggername)&&triggers.at(triggername);
  }
  else{
    if(asterix==int(triggername.length()-1)){
      triggername.pop_back();
      auto it=triggers.lower_bound(triggername);
      while(it!=triggers.end()&&it->first.find(triggername) == 0){
	if(it->second) return true;
	it++;
      }
      return false;
    }
    else{
      cerr << "* wildcard only allowed at the end" << endl;
      return false;
    }
  }
  return false;

}

bool TriggerInfo::IsAnyTriggered(std::vector< std::string > triggerNames) const {
  for(auto name=triggerNames.begin(); name!=triggerNames.end();name++){
    if(IsTriggered(*name)) return true;
  }
  return false;
}

bool TriggerInfo::IsAnyTriggered(std::string triggerstring) const{
    std::vector<std::string> triggernames;

    int posOr;
    while(true){
        posOr = triggerstring.find("||");
        if(posOr==int(std::string::npos)){
            triggernames.push_back(triggerstring);
            return IsAnyTriggered(triggernames);
        }
        triggernames.push_back(triggerstring.substr(0,posOr));
        triggerstring=triggerstring.substr(posOr+2,triggerstring.length());
        //std::cout<<"triggerstring = "<<triggerstring<<std::endl;
    }
}

std::map<std::string, bool> TriggerInfo::GetTriggers() const{
  return triggers;
}


void TriggerInfo::ListTriggers() {
  int counter=0;
  for(unsigned int i=0; i<triggernames.size();i++){
    cout << "Trigger #" << counter << " : " << triggernames[i] << endl;
    counter++;
  }
}
