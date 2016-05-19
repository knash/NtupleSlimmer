#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>


class  SlimUserData_TriggerUserData : public edm::EDProducer {
public:
  SlimUserData_TriggerUserData( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
 };


SlimUserData_TriggerUserData::SlimUserData_TriggerUserData(const edm::ParameterSet& iConfig)
 {   
   produces<std::vector<float>>("triggerBitTree");
   produces<std::vector<int>>("triggerPrescaleTree");  
   produces<std::vector<std::string>>("triggerNameTree");  
 }


void SlimUserData_TriggerUserData::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {


  std::auto_ptr<std::vector<float>> triggerBitTree(new std::vector<float>()); 
  edm::Handle<std::vector<float>> triggerBitTreeHandle;

  std::auto_ptr<std::vector<std::string>> triggerNameTree(new std::vector<std::string>()); 
  edm::Handle<std::vector<std::string>> triggerNameTreeHandle;

  std::auto_ptr<std::vector<int>> triggerPrescaleTree(new std::vector<int>()); 
  edm::Handle<std::vector<int>> triggerPrescaleTreeHandle;

  iEvent.getByLabel("TriggerUserData" , "triggerBitTree", triggerBitTreeHandle);
  iEvent.getByLabel("TriggerUserData" , "triggerPrescaleTree", triggerPrescaleTreeHandle);
  iEvent.getByLabel("TriggerUserData" , "triggerNameTree", triggerNameTreeHandle);


  for( size_t i=0; i<triggerNameTreeHandle->size(); i++ ) 
	{
	std::string tname = triggerNameTreeHandle->at(i);
	if (tname.find("HLT_PFHT800") || tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") ||  tname.find("HLT_PFHT475"))
		{
  		triggerNameTree->push_back(triggerNameTreeHandle->at(i));
  		triggerPrescaleTree->push_back(triggerPrescaleTreeHandle->at(i));
  		triggerBitTree->push_back(triggerBitTreeHandle->at(i));
		}
	}
  iEvent.put(triggerNameTree       ,"triggerNameTree");
  iEvent.put(triggerPrescaleTree   ,"triggerPrescaleTree");
  iEvent.put(triggerBitTree        ,"triggerBitTree");


 }
// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_TriggerUserData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_TriggerUserData::endJob() {
}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_TriggerUserData);
