
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include <TH1.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"


using namespace reco;
class  SlimUserData_GENFilter : public edm::EDFilter {
public:

  explicit SlimUserData_GENFilter( const edm::ParameterSet & );   

 // ~SlimUserData_GENFilter();
private:
  bool ISDATA_ ;
  std::vector<std::string> TrigNames;
  std::vector<int> TrigIndices;
  std::vector<std::string> FiltNames;
  std::vector<int> FiltIndices;
  bool filter( edm::Event &, const edm::EventSetup & );
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob() ;
 };


SlimUserData_GENFilter::SlimUserData_GENFilter(const edm::ParameterSet& iConfig):
   ISDATA_ (iConfig.getUntrackedParameter<bool>("ISDATA", false))
 {   
   edm::EDGetTokenT<std::vector<reco::GenParticle>>(consumes<std::vector<reco::GenParticle>>(edm::InputTag("filteredPrunedGenParticles" , ""))); 
 }


bool SlimUserData_GENFilter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup) {
 

  edm::Handle<std::vector<reco::GenParticle>  > GPHandle;
  iEvent.getByLabel("filteredPrunedGenParticles" , "", GPHandle);
  
  //std::cout<<GPHandle->at(0).pdgId()<<std::endl;
  int nfound=0;
  for( uint32_t igp=0 ; igp<GPHandle->size() ; igp++ ) 
	{
	if (std::abs(GPHandle->at(igp).pdgId()) == 1000024) nfound+=1;	
	if (nfound==2)break;
	}

  if (nfound!=2) return 0;

  else 
	{
	return 1;
	}
 }




void 
SlimUserData_GENFilter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{		
}



// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_GENFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_GENFilter::endJob() 
{

}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_GENFilter);








