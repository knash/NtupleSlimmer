

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>


class  SlimUserData_Filter : public edm::EDFilter {
public:
  explicit SlimUserData_Filter( const edm::ParameterSet & );   
 // ~SlimUserData_Filter();
private:
  bool filter( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
 };


SlimUserData_Filter::SlimUserData_Filter(const edm::ParameterSet& iConfig)
 {   
   produces<std::vector<bool>>("HT800bit"); 
   produces<std::vector<bool>>("DijetBit"); 
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSPt")));

   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("METUserData" , "triggerBitTree")));
   edm::EDGetTokenT<std::vector<std::string>>(consumes<std::vector<std::string>>(edm::InputTag("METUserData" ,"triggerNameTree"))); 

   edm::EDGetTokenT<bool>(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer" , "HBHENoiseFilterResult"))); 


   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("TriggerUserData" , "triggerBitTree"))); 
   edm::EDGetTokenT<std::vector<std::string>>(consumes<std::vector<std::string>>(edm::InputTag("TriggerUserData" , "triggerNameTree"))); 


 }


bool SlimUserData_Filter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<std::vector<bool>> HT800bit(new std::vector<bool>());    
  std::auto_ptr<std::vector<bool>> DijetBit(new std::vector<bool>()); 
             
  std::auto_ptr<std::vector<float>> jetAK8Pt(new std::vector<float>()); 
  edm::Handle<std::vector<float>> jetAK8PtHandle;
  iEvent.getByLabel("jetsAK8CHS" , "jetAK8CHSPt", jetAK8PtHandle);
  edm::Handle<std::vector<float>> metBitTreeHandle;
  edm::Handle<std::vector<std::string>> metNameTreeHandle;
  iEvent.getByLabel("METUserData" , "triggerBitTree", metBitTreeHandle);
  iEvent.getByLabel("METUserData" ,"triggerNameTree", metNameTreeHandle);




  std::auto_ptr<std::vector<std::string>> filters(new std::vector<std::string>({"Flag_goodVertices","Flag_CSCTightHaloFilter","Flag_eeBadScFilter"})) ;
  edm::Handle<bool> HBHENoiseFilterHandle;
  iEvent.getByLabel("HBHENoiseFilterResultProducer" , "HBHENoiseFilterResult", HBHENoiseFilterHandle);
  //std::cout<<"NewFilter"<<std::endl;
  if (not *HBHENoiseFilterHandle) 
	{
	//std::cout<<*HBHENoiseFilterHandle<<std::endl; 
	return 0;
	}
  //filters->push_back();
  for( size_t i=0; i<metNameTreeHandle->size(); i++ ) 
		{
  		for( size_t j=0; j<filters->size(); j++ ) 
			{
			if (filters->at(j)==metNameTreeHandle->at(i))

				{
				//std::cout<<filters->at(j)<<std::endl;
			//	std::cout<<metBitTreeHandle->at(i)<<std::endl;
				if (not metBitTreeHandle->at(i))
					{
					//std::cout<<"Throw Away"<<std::endl;
					return 0;
					}
				//std::cout<<std::endl;
				}
			}
			
		}
	






  if (iEvent.eventAuxiliary().isRealData())
	{


  	edm::Handle<std::vector<float>> triggerBitTreeHandle;
  	edm::Handle<std::vector<std::string>> triggerNameTreeHandle;


  	iEvent.getByLabel("TriggerUserData" , "triggerBitTree", triggerBitTreeHandle);
  	iEvent.getByLabel("TriggerUserData" , "triggerNameTree", triggerNameTreeHandle);



	bool trigpass = 0;
  	for( size_t i=0; i<triggerNameTreeHandle->size(); i++ ) 
		{

		std::string tname = triggerNameTreeHandle->at(i);

		if (tname.find("HLT_PFHT800") != std::string::npos || tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") != std::string::npos)
			{
  			if (triggerBitTreeHandle->at(i)) trigpass=1;
			if (tname.find("HLT_PFHT800") != std::string::npos)
				{
  				HT800bit->push_back(triggerBitTreeHandle->at(i));    
				}
			if (tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45") != std::string::npos)
				{
  				DijetBit->push_back(triggerBitTreeHandle->at(i)); 
				}

			}
		}
  	if (not trigpass) return 0;
	}
  if (jetAK8PtHandle->size()<2 ) return 0;
  if (jetAK8PtHandle->at(0)<300. || jetAK8PtHandle->at(1)<300. ) return 0;
  iEvent.put(DijetBit,"DijetBit");
  iEvent.put(HT800bit,"HT800bit");
  return 1;


 }
// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_Filter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_Filter::endJob() {
}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_Filter);








