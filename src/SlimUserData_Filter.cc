
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
  bool filter( edm::Event &, const edm::EventSetup & );
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void beginJob() ;
  void endJob() ;
 };


SlimUserData_Filter::SlimUserData_Filter(const edm::ParameterSet& iConfig):
   ISDATA_ (iConfig.getUntrackedParameter<bool>("ISDATA", false))
 {   
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


 }


bool SlimUserData_Filter::filter( edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  if (ISDATA_)
  {
 /*
  edm::Handle<std::vector<float>> metBitTreeHandle;

  std::auto_ptr<std::vector<std::string>> filters(new std::vector<std::string>({"Flag_goodVertices","Flag_CSCTightHaloFilter","Flag_eeBadScFilter"})) ;
  for( size_t i=0; i<metNameTreeHandle->size(); i++ ) 
		{
  		for( size_t j=0; j<filters->size(); j++ ) 
			{
			if (filters->at(j)==metNameTreeHandle->at(i))

				{
				std::cout<<filters->at(j)<<std::endl;
				std::cout<<metBitTreeHandle->at(i)<<std::endl;
				if (not metBitTreeHandle->at(i))
					{
					//std::cout<<"Throw Away"<<std::endl;
					return 0;
					}
				//std::cout<<std::endl;
				}
			}
			
		}
	








	*/
  	edm::Handle<std::vector<float>> triggerBitTreeHandle;
  	iEvent.getByLabel("TriggerUserData" , "triggerBitTree", triggerBitTreeHandle);

	bool trigpass = 0;
  	for( size_t i=0; i<TrigNames.size(); i++ ) 
		{

		std::string tname = TrigNames.at(i);
		int tind = TrigIndices.at(i);
		
		//std::cout<<"tname "<< tname << std::endl;
		//std::cout<<"ind "<< tind << std::endl<< std::endl;
		
  		if (triggerBitTreeHandle->at(tind)) trigpass=1;
		if (tname.find("HLT_PFHT800") != std::string::npos)
			{
  			HT800bit->push_back(triggerBitTreeHandle->at(tind));    
			}
		if (tname.find("HLT_PFHT900") != std::string::npos)
			{
  			HT900bit->push_back(triggerBitTreeHandle->at(tind));    
			}
		if (tname.find("HLT_PFHT475") != std::string::npos)
			{
  			HT475bit->push_back(triggerBitTreeHandle->at(tind));    
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
 

  if (jetAK8CHSPtHandle->size()<2 && jetAK8PuppiPtHandle->size()<2) return 0;

  else if (jetAK8CHSPtHandle->size()<2)
	{
	if (jetAK8PuppiPtHandle->at(0)<250. || jetAK8PuppiPtHandle->at(1)<250.) return 0;
	else return 1;
	}

  else if (jetAK8PuppiPtHandle->size()<2)
	{
	if (jetAK8CHSPtHandle->at(0)<250. || jetAK8CHSPtHandle->at(1)<250.) return 0;
	else return 1;
	}

  else if ((jetAK8CHSPtHandle->at(0)<250. || jetAK8CHSPtHandle->at(1)<250.) && (jetAK8PuppiPtHandle->at(0)<250. || jetAK8PuppiPtHandle->at(1)<250.)) return 0;

  else return 1;


 }




void 
SlimUserData_Filter::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{


	if (ISDATA_)
		{

  			edm::Handle<std::vector<std::string>> triggerNameTreeHandle;
 			edm::Handle<std::vector<std::string>> metNameTreeHandle;
  			iRun.getByLabel("TriggerUserData" , "triggerNameTree", triggerNameTreeHandle);

		  	for( size_t i=0; i<triggerNameTreeHandle->size(); i++ ) 
				{
				std::string tname = triggerNameTreeHandle->at(i);
				if (tname.find("HLT_PFHT800") != std::string::npos || tname.find("HLT_PFHT900") != std::string::npos   ||    tname.find("HLT_AK8PFJet450") != std::string::npos  || tname.find("HLT_AK8PFJet260") != std::string::npos  || tname.find("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20") != std::string::npos  || tname.find("HLT_PFHT475") != std::string::npos)
					{		
					TrigNames.push_back(tname);
					TrigIndices.push_back(i);
					}
				}

    			/*iRun.getByLabel("METUserData","triggerNameTree", metNameTreeHandle);
  			std::auto_ptr<std::vector<std::string>> filters(new std::vector<std::string>({"Flag_HBHENoiseIsoFilter","Flag_HBHENoiseFilter","Flag_goodVertices","Flag_CSCTightHaloFilter","Flag_eeBadScFilter"})) ;
  			for( size_t i=0; i<metNameTreeHandle->size(); i++ ) 
				{
  				for( size_t j=0; j<filters->size(); j++ ) 
					{
					std::cout<<metNameTreeHandle->at(i)<<std::endl;
					if (filters->at(j)==metNameTreeHandle->at(i))

						{
							std::cout<<"Found at index "<< i << std::endl;
							//std::cout<<metBitTreeHandle->at(i)<<std::endl;
						}
					}
				}
			*/


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








