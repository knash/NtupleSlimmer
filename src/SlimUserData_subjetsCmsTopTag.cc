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


class  SlimUserData_subjetsCmsTopTag : public edm::EDProducer {
public:
  SlimUserData_subjetsCmsTopTag( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
 };


SlimUserData_subjetsCmsTopTag::SlimUserData_subjetsCmsTopTag(const edm::ParameterSet& iConfig)
 {   

   produces<std::vector<float>>("subjetCmsTopTagCSV");
   produces<std::vector<float>>("subjetAK8CSV");  

 }


void SlimUserData_subjetsCmsTopTag::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<float>> subjetCmsTopTagCSV(new std::vector<float>()); 
  edm::Handle<std::vector<float>> subjetCmsTopTagCSVHandle;

  iEvent.getByLabel("jetsAK8" , "subjetCmsTopTagCSV", subjetCmsTopTagCSVHandle);


  for( size_t i=0; i<subjetCmsTopTagCSV->size(); i++ ) 
	{
	subjetCmsTopTagCSV->push_back(subjetCmsTopTagCSVHandle->at(i));      
	}
  iEvent.put(subjetCmsTopTagCSV        ,"subjetCmsTopTagCSV");


  std::auto_ptr<std::vector<float>> subjetSoftDropCSV(new std::vector<float>()); 
  edm::Handle<std::vector<float>> subjetSoftDropCSVHandle;

  iEvent.getByLabel("jetsAK8" , "subjetSoftDropCSV", subjetSoftDropCSVHandle);


  for( size_t i=0; i<subjetSoftDropCSV->size(); i++ ) 
	{
	subjetSoftDropCSV->push_back(subjetSoftDropCSVHandle->at(i));      
	}
  iEvent.put(subjetSoftDropCSV        ,"subjetSoftDropCSV");


 }
// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_subjetsCmsTopTag::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_subjetsCmsTopTag::endJob() {
}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_subjetsCmsTopTag);
