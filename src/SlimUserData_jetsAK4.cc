
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"


//#include <errno.h>
#include <Math/VectorUtil.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TF1.h>
#include "SlimUserData_JecCorr.h"







class  SlimUserData_jetsAK4 : public edm::EDProducer {
public:
  SlimUserData_jetsAK4( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
  bool reapplyjec_,reapplyjer_;
  std::string jes_,jer_,era_;

 };


SlimUserData_jetsAK4::SlimUserData_jetsAK4(const edm::ParameterSet& iConfig):
   reapplyjec_ (iConfig.getParameter<bool>("reapplyjec")),
   reapplyjer_ (iConfig.getParameter<bool>("reapplyjer")),
   jes_ (iConfig.getParameter<std::string>("jes")),
   jer_ (iConfig.getParameter<std::string>("jer")),
   era_ (iConfig.getParameter<std::string>("era"))


 {   

   produces<std::vector<float>>("jetAK4PuppiPt");
   produces<std::vector<float>>("jetAK4PuppiPhi"); 
   produces<std::vector<float>>("jetAK4PuppiEta");  
   produces<std::vector<float>>("jetAK4PuppiMass");

   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi","jetAK4PuppiSmearedPt")));



  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi","jetAK4PuppiPtResolution")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi","jetAK4PuppiJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi","jetAK4PuppiJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi","jetAK4PuppiJERSFDown")));  


  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppijecFactor0"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiPt"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiPhi"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiEta"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiE"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppijetArea"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiGenJetPt" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiGenJetEta" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiGenJetPhi" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiGenJetE" )));

  edm::EDGetTokenT<int>(consumes<int>(edm::InputTag("vertexInfo", "npv"    )));
  edm::EDGetTokenT<double>(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll", ""    )));





  if (jes_=="nominal"&&jer_=="nominal")
	{



  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiCMVAv2"     )));  

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiPartonFlavour"  )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppijecFactor0"  )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppichargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppineutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppineutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppiNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppichargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppineutralMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK4Puppi", "jetAK4PuppichargedEmEnergyFrac" )));

  		produces<std::vector<float>>("jetAK4PuppiTight"); 
   		produces<std::vector<float>>("jetAK4PuppiLoose"); 
 

   		produces<std::vector<float>>("jetAK4PuppiCSV"); 
   		produces<std::vector<float>>("jetAK4PuppiCMVAv2"); 



   		produces<std::vector<float>>("jetAK4PuppiPartonFlavour"); 
   		produces<std::vector<float>>("jetAK4PuppijecFactor0");

		}

 }



void SlimUserData_jetsAK4::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {



  bool ISDATA = iEvent.eventAuxiliary().isRealData();
  unsigned int runnum = iEvent.eventAuxiliary().run();


  std::auto_ptr<int> npv(new int()); 
  std::auto_ptr<double> Rho(new double()); 

  edm::Handle<int> npvHandle;
  edm::Handle<double> RhoHandle;


  std::auto_ptr<std::vector<float>> jetAK4PuppiPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiSmearedPt(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK4PuppiPhi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK4PuppiEta(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK4PuppiMass(new std::vector<float>()); 
      

  std::auto_ptr<std::vector<float>> jetAK4PuppiTight(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiLoose(new std::vector<float>()); 


  std::auto_ptr<std::vector<float>> jetAK4PuppijetArea(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiE(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppijecFactor0(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiGenJetPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiGenJetEta(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK4PuppiGenJetPhi(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK4PuppiCSV(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK4PuppiCMVAv2(new std::vector<float>());
              
  std::auto_ptr<std::vector<float>> jetAK4PuppiPartonFlavour(new std::vector<float>());  
     

  edm::Handle<std::vector<float>> jetAK4PuppiEHandle;
  edm::Handle<std::vector<float>> jetAK4PuppijecFactor0Handle;
  edm::Handle<std::vector<float>> jetAK4PuppijetAreaHandle;
  edm::Handle<std::vector<float>> jetAK4PuppiGenJetPtHandle;
  edm::Handle<std::vector<float>> jetAK4PuppiGenJetEtaHandle;
  edm::Handle<std::vector<float>> jetAK4PuppiGenJetPhiHandle;

  edm::Handle<std::vector<float>> jetAK4PuppiPtHandle; 

  edm::Handle<std::vector<float>> jetAK4PuppiSmearedPtHandle;             
  edm::Handle<std::vector<float>> jetAK4PuppiPhiHandle;        
  edm::Handle<std::vector<float>> jetAK4PuppiEtaHandle;        
  edm::Handle<std::vector<float>> jetAK4PuppiMassHandle;   


  edm::Handle<std::vector<float>> jetAK4PuppiJERSFHandle;
  edm::Handle<std::vector<float>> jetAK4PuppiJERSFUpHandle;
  edm::Handle<std::vector<float>> jetAK4PuppiJERSFDownHandle;

  edm::Handle<std::vector<float>> jetAK4PuppiRESHandle;

  edm::Handle<std::vector<float>> jetAK4PuppiCSVHandle;     
  edm::Handle<std::vector<float>> jetAK4PuppiCMVAv2Handle;
           
  edm::Handle<std::vector<float>> jetAK4PuppiPartonFlavourHandle;


  edm::Handle<std::vector<float>>  jetAK4PuppichargedHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppineutralEmEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppineutralHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppiNumConstituentsHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppichargedMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppineutralMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK4PuppichargedEmEnergyFracHandle;  


  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiPt"        ,jetAK4PuppiPtHandle);
  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiSmearedPt"        ,jetAK4PuppiSmearedPtHandle);
  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiPhi"       ,jetAK4PuppiPhiHandle);  
  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiEta"       ,jetAK4PuppiEtaHandle);  
  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiE"        ,jetAK4PuppiEHandle);

  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppijecFactor0"        ,jetAK4PuppijecFactor0Handle);
  iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppijetArea"        ,jetAK4PuppijetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiGenJetPt"        ,jetAK4PuppiGenJetPtHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiGenJetEta"        ,jetAK4PuppiGenJetEtaHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiGenJetPhi"        ,jetAK4PuppiGenJetPhiHandle);


  iEvent.getByLabel("vertexInfo", "npv"        ,npvHandle);
  iEvent.getByLabel("fixedGridRhoFastjetAll", ""        ,RhoHandle);

  iEvent.getByLabel( "jetsAK4Puppi", "jetAK4PuppiPtResolution",  jetAK4PuppiRESHandle);

  iEvent.getByLabel( "jetsAK4Puppi", "jetAK4PuppiJERSF",  jetAK4PuppiJERSFHandle);
  iEvent.getByLabel( "jetsAK4Puppi", "jetAK4PuppiJERSFUp",  jetAK4PuppiJERSFUpHandle);
  iEvent.getByLabel( "jetsAK4Puppi", "jetAK4PuppiJERSFDown",  jetAK4PuppiJERSFDownHandle);



  if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppichargedHadronEnergyFrac"       ,jetAK4PuppichargedHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppineutralEmEnergyFrac"       ,jetAK4PuppineutralEmEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppineutralHadronEnergyFrac"       ,jetAK4PuppineutralHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiNumConstituents"       ,jetAK4PuppiNumConstituentsHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppichargedMultiplicity"       ,jetAK4PuppichargedMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppineutralMultiplicity"       ,jetAK4PuppineutralMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppichargedEmEnergyFrac"       ,jetAK4PuppichargedEmEnergyFracHandle);  








  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiCSVv2"       ,jetAK4PuppiCSVHandle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiCMVAv2"       ,jetAK4PuppiCMVAv2Handle);  
  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppiPartonFlavour"   ,jetAK4PuppiPartonFlavourHandle);  

  	iEvent.getByLabel("jetsAK4Puppi", "jetAK4PuppijecFactor0"   ,jetAK4PuppijecFactor0Handle);  
	}

JME::JetResolutionScaleFactor res_sfP;
JME::JetResolution resoP;

std::string JERFileP_ = "jecfiles/Spring16_25nsV10_MC_SF_AK4PFPuppi.txt";
std::string RESFileP_ = "jecfiles/Spring16_25nsV10_MC_PtResolution_AK4PFPuppi.txt";
res_sfP = JME::JetResolutionScaleFactor(JERFileP_);
resoP = JME::JetResolution(RESFileP_);


boost::shared_ptr<FactorizedJetCorrector> jecAK4PUPPI_;
std::vector<std::string> jecAK4PUPPIPayloadNames_;
std::string runtxt_;

if(ISDATA) 
		{

		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV4";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV4";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV4"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV4"; //IOV H:[280919,Infinity] f



  		jecAK4PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2Relative_AK4PFPuppi.txt");
  		jecAK4PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L3Absolute_AK4PFPuppi.txt");
  		jecAK4PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2L3Residual_AK4PFPuppi.txt");
		}
else
		{

		runtxt_ = "V4";
  		jecAK4PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L2Relative_AK4PFPuppi.txt");
  		jecAK4PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L3Absolute_AK4PFPuppi.txt");
		}

std::vector<JetCorrectorParameters> vParPUPPI;

for ( std::vector<std::string>::const_iterator payloadBegin = jecAK4PUPPIPayloadNames_.begin(), payloadEnd = jecAK4PUPPIPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) 
		{
  	  	JetCorrectorParameters pars(*ipayload);
 	   	vParPUPPI.push_back(pars);
 		}
jecAK4PUPPI_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vParPUPPI) );




JetCorrectionUncertainty *jecUncPUPPI = new JetCorrectionUncertainty("jecfiles/"+era_+"V4_MC_Uncertainty_AK4PFPuppi.txt");

  //Puppi
  for( size_t i=0; i<jetAK4PuppiPtHandle->size(); i++ ) 
	{

	float shift = 1.0;
	float JECcorr = 1.0;

        TLorentzVector v1;
	v1.SetPtEtaPhiE(jetAK4PuppiPtHandle->at(i),jetAK4PuppiEtaHandle->at(i),jetAK4PuppiPhiHandle->at(i),jetAK4PuppiEHandle->at(i));
	float calcmass =  v1.M();




	float uncorrpt = fmaxf(1.0,jetAK4PuppijecFactor0Handle->at(i)*jetAK4PuppiPtHandle->at(i));	
	//float uncorrE = fmaxf(1.0,jetAK4PuppijecFactor0Handle->at(i)*jetAK4PuppiEHandle->at(i));	

   	if (reapplyjec_)
		{
		JECcorr = fmaxf(0.0,jetAK4PuppijecFactor0Handle->at(i)*JEC_Corr(jecAK4PUPPI_,uncorrpt,jetAK4PuppiEtaHandle->at(i),ISDATA,jetAK4PuppijetAreaHandle->at(i),*RhoHandle.product()));
		}


	float puppicorrpt = jetAK4PuppiPtHandle->at(i)*JECcorr;
   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert(jecUncPUPPI,puppicorrpt,jetAK4PuppiEtaHandle->at(i),jes_));
		}
	//std::cout<<"JECUPDATE "<<JECcorr*shift<<std::endl;

	if (not ISDATA) 
		{
			float SFapp=1.0;
			float RESapp=1.0;

			if (reapplyjer_)
				{	
					
				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(puppicorrpt).setJetEta(jetAK4PuppiEtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") 
						{
						SFapp = res_sfP.getScaleFactor(jetParam);
						}
					else if (jer_ == "up")
						{ 
						SFapp = res_sfP.getScaleFactor(jetParam, Variation::UP);	
						}
 					else if (jer_ == "down") 
						{
						SFapp = res_sfP.getScaleFactor(jetParam, Variation::DOWN);
						}

				
					RESapp = resoP.getResolution(jetParam)*puppicorrpt;

				}
			else
				{

					if (jer_ == "nominal") 
						{
						SFapp = jetAK4PuppiJERSFHandle->at(i);
						}
					else if (jer_ == "up")
						{
						SFapp = jetAK4PuppiJERSFUpHandle->at(i);
						}
 					else if (jer_ == "down")
						{
						SFapp = jetAK4PuppiJERSFDownHandle->at(i);
						}
					RESapp = jetAK4PuppiRESHandle->at(i)*puppicorrpt;

				}	


			float DRval =  deltaR(jetAK4PuppiEtaHandle->at(i), jetAK4PuppiPhiHandle->at(i), jetAK4PuppiGenJetEtaHandle->at(i), jetAK4PuppiGenJetPhiHandle->at(i)) ;
			float ptval = std::fabs(puppicorrpt-jetAK4PuppiGenJetPtHandle->at(i));


			bool DRmatch = DRval<0.2 ;
			bool ptmatch = ptval<3.0*RESapp;



			float ptcorrunsmear = puppicorrpt;
			float ptgen =jetAK4PuppiGenJetPtHandle->at(i);

			float smearcorr = 1.0+(SFapp-1.0)*((ptcorrunsmear-ptgen)/ptcorrunsmear);



			if (DRmatch and ptmatch)
				{
					shift = shift*fmaxf(0.0,smearcorr);	
				}
			else
				{

					/*
					std::cout<<"Hybrid: "<<std::endl;
					std::cout<<"DeltaR match: "<<DRval<<std::endl;
					std::cout<<"DeltaPt match: "<<ptval <<std::endl;
					std::cout<<"RES: "<<RESapp<<std::endl;
					std::cout<<"RES threesig: "<<3.0*RESapp<<std::endl;
				        std::cout<<"SF: "<<SFapp<<std::endl;
					*/



					Double_t ptresfac =std::sqrt(SFapp*SFapp-1.0 );
					TRandom3 *r1 = new TRandom3(0); 
					Double_t ptsmearify = r1->Gaus(0.0,RESapp/puppicorrpt)*ptresfac; 
					shift = shift*fmaxf(0.0,(1.0+ptsmearify));
				        //std::cout<<"smearify: "<<ptsmearify<<std::endl;
					delete r1;

				}

		/*
		std::cout<<"JER = "<<jer_<<std::endl;
		std::cout<<"JES = "<<jes_<<std::endl;

		std::cout<<"PUPPTRAW = "<<puppicorrpt<<std::endl;
		std::cout<<"SHIFT = "<<shift<<std::endl;
		std::cout<<"JECCORR = "<<JECcorr<<std::endl;
		std::cout<<"PTcor = "<<puppicorrpt*shift<<std::endl;
		*/
		}		


	jetAK4PuppiPt->push_back(puppicorrpt*shift);      
	jetAK4PuppiPhi->push_back(jetAK4PuppiPhiHandle->at(i));       
	jetAK4PuppiEta->push_back(jetAK4PuppiEtaHandle->at(i));       


	jetAK4PuppiMass->push_back(calcmass*shift*JECcorr);  

   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK4PuppiCSV->push_back(jetAK4PuppiCSVHandle->at(i));
		jetAK4PuppiCMVAv2->push_back(jetAK4PuppiCMVAv2Handle->at(i));              
		jetAK4PuppiPartonFlavour->push_back(jetAK4PuppiPartonFlavourHandle->at(i)); 
		jetAK4PuppijecFactor0->push_back(jetAK4PuppijecFactor0Handle->at(i));   
   

	  	float CHF = jetAK4PuppichargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK4PuppineutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK4PuppineutralHadronEnergyFracHandle->at(i);
	  	float CM = jetAK4PuppichargedMultiplicityHandle->at(i);  
	  	float CEMF = jetAK4PuppichargedEmEnergyFracHandle->at(i);

		float NC = jetAK4PuppineutralMultiplicityHandle->at(i) + jetAK4PuppichargedMultiplicityHandle->at(i); 

		float TJ = 0.0;
		if ((NHF<0.9) and (NEMF<0.9) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) TJ = 1.0;

		float LJ = 0.0;
		if ((NHF<0.99) and (NEMF<0.99) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) LJ = 1.0;


		jetAK4PuppiTight->push_back(TJ); 
		jetAK4PuppiLoose->push_back(LJ);




		}

	}


  delete jecUncPUPPI;
 

  iEvent.put(jetAK4PuppiPt,"jetAK4PuppiPt");
  iEvent.put(jetAK4PuppiPhi,"jetAK4PuppiPhi"); 
  iEvent.put(jetAK4PuppiEta,"jetAK4PuppiEta");  
  iEvent.put(jetAK4PuppiMass,"jetAK4PuppiMass");


   if (jes_=="nominal"&&jer_=="nominal")
	{
	
  	iEvent.put(jetAK4PuppiCSV,"jetAK4PuppiCSV"); 
  	iEvent.put(jetAK4PuppiCMVAv2,"jetAK4PuppiCMVAv2");
 
  	iEvent.put(jetAK4PuppiPartonFlavour,"jetAK4PuppiPartonFlavour"); 

  	iEvent.put(jetAK4PuppiTight,"jetAK4PuppiTight"); 
  	iEvent.put(jetAK4PuppiLoose,"jetAK4PuppiLoose");

  	iEvent.put(jetAK4PuppijecFactor0,"jetAK4PuppijecFactor0");

	}	

 }
// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_jetsAK4::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_jetsAK4::endJob() {
}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_jetsAK4);
