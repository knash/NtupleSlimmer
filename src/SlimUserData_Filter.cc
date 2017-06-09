
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include <TH1.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"



class  SlimUserData_Filter : public edm::EDFilter {
public:

  explicit SlimUserData_Filter( const edm::ParameterSet & );   

 // ~SlimUserData_Filter();
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


SlimUserData_Filter::SlimUserData_Filter(const edm::ParameterSet& iConfig):
   ISDATA_ (iConfig.getUntrackedParameter<bool>("ISDATA", false))
 {   
   produces<std::vector<bool>>("filtersbit"); 
   produces<std::vector<bool>>("HT800bit"); 
   produces<std::vector<bool>>("HT900bit"); 
   produces<std::vector<bool>>("JET450bit"); 
   produces<std::vector<bool>>("JET260bit"); 
   produces<std::vector<bool>>("HT475bit"); 
   produces<std::vector<bool>>("DijetBit"); 
   edm::EDGetTokenT<std::vector<std::string>>(mayConsume<std::vector<std::string>,edm::InRun>(edm::InputTag("TriggerUserData" , "triggerNameTree")));
   edm::EDGetTokenT<std::vector<std::string>>(mayConsume<std::vector<std::string>,edm::InRun>(edm::InputTag("METUserData","triggerNameTree")));
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSPt")));
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiPt")));
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("METUserData" , "triggerBitTree")));
   edm::EDGetTokenT<bool>(consumes<bool>(edm::InputTag("HBHENoiseFilterResultProducer" , "HBHENoiseFilterResult"))); 
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("TriggerUserData" , "triggerBitTree"))); 
   edm::EDGetTokenT<bool>(consumes<bool>(edm::InputTag("BadChargedCandidateFilter" , ""))); 
   edm::EDGetTokenT<bool>(consumes<bool>(edm::InputTag("BadPFMuonFilter" , ""))); 

 }


bool SlimUserData_Filter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup) {
  std::auto_ptr<std::vector<bool>> filtersbit(new std::vector<bool>()); 
  std::auto_ptr<std::vector<bool>> HT800bit(new std::vector<bool>()); 
  std::auto_ptr<std::vector<bool>> HT900bit(new std::vector<bool>()); 
  std::auto_ptr<std::vector<bool>> HT475bit(new std::vector<bool>());       
  std::auto_ptr<std::vector<bool>> DijetBit(new std::vector<bool>()); 
  std::auto_ptr<std::vector<bool>> JET450bit(new std::vector<bool>());       
  std::auto_ptr<std::vector<bool>> JET260bit(new std::vector<bool>());             


  std::auto_ptr<std::vector<float>> jetAK8CHSPt(new std::vector<float>()); 
  edm::Handle<std::vector<float>> jetAK8CHSPtHandle;
  iEvent.getByLabel("jetsAK8CHS" , "jetAK8CHSPt", jetAK8CHSPtHandle);

  std::auto_ptr<std::vector<float>> jetAK8PuppiPt(new std::vector<float>()); 
  edm::Handle<std::vector<float>> jetAK8PuppiPtHandle;
  iEvent.getByLabel("jetsAK8Puppi" , "jetAK8PuppiPt", jetAK8PuppiPtHandle);

  edm::Handle<std::vector<float>> metBitTreeHandle;
  iEvent.getByLabel("METUserData" , "triggerBitTree", metBitTreeHandle);

  edm::Handle<bool> BadChargedCandidateFilterHandle;
  edm::Handle<bool> BadPFMuonFilterHandle;
  iEvent.getByLabel("BadPFMuonFilter" , "", BadPFMuonFilterHandle);
  iEvent.getByLabel("BadChargedCandidateFilter" , "", BadChargedCandidateFilterHandle);
  bool BPFmu = BadPFMuonFilterHandle.product();
  bool BCCFmu = BadChargedCandidateFilterHandle.product();
  for( size_t i=0; i<FiltNames.size(); i++ ) 
	{

	std::string filtname = FiltNames.at(i);
	int filtind = FiltIndices.at(i);
	filtersbit->push_back(metBitTreeHandle->at(filtind));
	}
  filtersbit->push_back(BPFmu);
  filtersbit->push_back(BCCFmu);
  if (ISDATA_)
  {
 

  	edm::Handle<std::vector<float>> triggerBitTreeHandle;
  	iEvent.getByLabel("TriggerUserData" , "triggerBitTree", triggerBitTreeHandle);
  	//edm::Handle<std::vector<float>> triggerPrescaleHandle;
  	//iEvent.getByLabel("TriggerUserData" , "triggerPrescaleTree", triggerPrescaleHandle);

	bool trigpass = 0;
  	for( size_t i=0; i<TrigNames.size(); i++ ) 
		{

		std::string tname = TrigNames.at(i);
		int tind = TrigIndices.at(i);
		
		//std::cout<<"tname "<< tname << std::endl;
		//std::cout<<"ind "<< tind << std::endl;
		//std::cout<<triggerBitTreeHandle->at(tind)<<std::endl<< std::endl;
  		if (triggerBitTreeHandle->at(tind)) trigpass=1;
		if (tname.find("HLT_PFHT800") != std::string::npos)
			{
  			HT800bit->push_back(triggerBitTreeHandle->at(tind));
		//	std::cout<<"tname "<< tname << std::endl;
		//	std::cout<<"ind "<< tind << std::endl;
		//	std::cout<<triggerBitTreeHandle->at(tind)<< std::endl;
		    
			}
		if (tname.find("HLT_PFHT900") != std::string::npos)
			{
  			HT900bit->push_back(triggerBitTreeHandle->at(tind));  
		//	std::cout<<"tname "<< tname << std::endl;
		//	std::cout<<"ind "<< tind << std::endl;
		//	std::cout<<triggerBitTreeHandle->at(tind)<< std::endl;  
		//	std::cout<<HT900bit->at(0)<< std::endl;  
			}
		if (tname.find("HLT_PFHT475") != std::string::npos)
			{
  			HT475bit->push_back(triggerBitTreeHandle->at(tind)); 
		//	std::cout<<"tname "<< tname << std::endl;
		//	std::cout<<"ind "<< tind << std::endl;
		//	std::cout<<triggerBitTreeHandle->at(tind)<< std::endl;
  		//	std::cout<<HT475bit->at(0)<< std::endl;
			}
		if (tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20") != std::string::npos)
			{
  			DijetBit->push_back(triggerBitTreeHandle->at(tind)); 
			}

		if (tname.find("HLT_AK8PFJet450") != std::string::npos)
			{
  			JET450bit->push_back(triggerBitTreeHandle->at(tind));    
			}

		if (tname.find("HLT_AK8PFJet260") != std::string::npos)
			{
  			JET260bit->push_back(triggerBitTreeHandle->at(tind));    
			}

		}

  	if (not trigpass) return 0;

  	iEvent.put(DijetBit,"DijetBit");

  	iEvent.put(HT800bit,"HT800bit");
  	iEvent.put(HT900bit,"HT900bit");
  	iEvent.put(HT475bit,"HT475bit");
  	iEvent.put(JET450bit,"JET450bit");
  	iEvent.put(JET260bit,"JET260bit");
	
  }

  float ptcutoff = 200.0;
 
  iEvent.put(filtersbit,"filtersbit");
  if (jetAK8CHSPtHandle->size()<2 && jetAK8PuppiPtHandle->size()<2) return 0;

  else if (jetAK8CHSPtHandle->size()<2)
	{
	if (jetAK8PuppiPtHandle->at(0)<ptcutoff || jetAK8PuppiPtHandle->at(1)<ptcutoff) return 0;
	else return 1;
	}

  else if (jetAK8PuppiPtHandle->size()<2)
	{
	if (jetAK8CHSPtHandle->at(0)<ptcutoff || jetAK8CHSPtHandle->at(1)<ptcutoff) return 0;
	else return 1;
	}

  else if ((jetAK8CHSPtHandle->at(0)<ptcutoff || jetAK8CHSPtHandle->at(1)<ptcutoff) && (jetAK8PuppiPtHandle->at(0)<ptcutoff || jetAK8PuppiPtHandle->at(1)<ptcutoff)) return 0;

  else return 1;


 }




void 
SlimUserData_Filter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
 	edm::Handle<std::vector<std::string>> metNameTreeHandle;

  	TrigNames.clear();
  	TrigIndices.clear();
  	FiltNames.clear();
  	FiltIndices.clear();
	if (ISDATA_)
		{

  			edm::Handle<std::vector<std::string>> triggerNameTreeHandle;

  			iRun.getByLabel("TriggerUserData" , "triggerNameTree", triggerNameTreeHandle);

		  	for( size_t i=0; i<triggerNameTreeHandle->size(); i++ ) 
				{
				std::string tname = triggerNameTreeHandle->at(i);
				if (tname.find("HLT_PFHT800") != std::string::npos || tname.find("HLT_PFHT900") != std::string::npos   ||    tname.find("HLT_AK8PFJet450") != std::string::npos  || tname.find("HLT_AK8PFJet260") != std::string::npos  || tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20") != std::string::npos  || tname.find("HLT_PFHT475") != std::string::npos)
					{		
					std::cout<<"Found at index "<< i << std::endl;
					TrigNames.push_back(tname);
					TrigIndices.push_back(i);
					}
				}
		}








    		iRun.getByLabel("METUserData","triggerNameTree", metNameTreeHandle);
  		std::vector<std::string> filters ={"Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_globalTightHalo2016Filter","Flag_goodVertices"};
		if (ISDATA_) filters.push_back("Flag_eeBadScFilter");
  		for( size_t i=0; i<metNameTreeHandle->size(); i++ ) 
			{
			std::cout<<metNameTreeHandle->at(i)<<std::endl;
			std::string fname = metNameTreeHandle->at(i);
  			for( size_t j=0; j<filters.size(); j++ ) 
				{

				if (filters[j]==metNameTreeHandle->at(i))
					{
						FiltNames.push_back(fname);
						FiltIndices.push_back(i);
						std::cout<<"Found at index "<< i << std::endl;
						//std::cout<<metBitTreeHandle->at(i)<<std::endl;
					}
				}
			}
			


		
}



// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_Filter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_Filter::endJob() 
{

}

#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_Filter);








