#include <memory>
#include <cmath>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"


// dR and dPhi
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// Muons
//#include "ttbarDM/TopPlusDMAna/interface/Muons.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// trigger
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" // gives access to the (release cycle dependent) trigger object codes

#include <TFile.h>
#include <TH1F.h>
#include <TGraphAsymmErrors.h>
#include<vector>

using namespace reco;
using namespace edm;
using namespace std;
using namespace trigger;

class  SlimUserData_counter : public edm::one::EDProducer<edm::BeginRunProducer,edm::EndRunProducer> {
public:
  SlimUserData_counter( const edm::ParameterSet & );   

private:
  int nev;
  void produce( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
  void beginRunProduce(edm::Run& , const edm::EventSetup &);

  void endRunProduce(edm::Run&, const edm::EventSetup &);

 };


SlimUserData_counter::SlimUserData_counter(const edm::ParameterSet& iConfig)
 {   
  produces<std::vector<int>, edm::InRun>("nevr");
  nev = 0;
 }

void 
SlimUserData_counter::beginRunProduce(edm::Run& iRun, edm::EventSetup const& iSetup)
{
}


void SlimUserData_counter::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
nev+=1;
}

void SlimUserData_counter::endRunProduce(edm::Run& run, edm::EventSetup const& es)
{
   //	std::auto_ptr<int> nevr(new int); 
        std::auto_ptr<std::vector<int>> nevr(new std::vector<int>()); 
	nevr->push_back(nev);
	std::cout<<"Total Events Processed: "<<nevr->at(0)<<std::endl;
	run.put(nevr,"nevr");
}


// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_counter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_counter::endJob() {
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SlimUserData_counter);
