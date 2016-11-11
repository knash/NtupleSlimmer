
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include <errno.h>
#include <Math/VectorUtil.h>
//#include <TLorentzVector.h>






float JES_Uncert(float pt,float eta,std::string val,std::string era_)
	{
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(era_+"_MC_Uncertainty_AK8PFchs.txt");
	int sign = 0;
	if (val=="up") sign = 1;
	else if (val=="down") sign = -1;


  	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float shift = (1+sign*jecUnc->getUncertainty(true));

	delete jecUnc;
	return  shift;
	}

float JEC_Corr(float  pt,float eta  , bool isdata,float Area,const double Rho,std::string era_)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;

  	if(isdata) 
		{
 		jecAK8PayloadNames_.push_back(era_+"_DATA_L1FastJet_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L2L3Residual_AK8PFchs.txt");
		}
  	else
		{
  		jecAK8PayloadNames_.push_back(era_+"_MC_L1FastJet_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_MC_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_MC_L3Absolute_AK8PFchs.txt");
		}
  	std::vector<JetCorrectorParameters> vPar;

  	for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) 
		{
  	  	JetCorrectorParameters pars(*ipayload);
 	   	vPar.push_back(pars);
 		}
  	jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );


        jecAK8_->setJetEta( eta );
        jecAK8_->setJetPt ( pt );
        //jecAK8_->setJetE  ( energy);
        jecAK8_->setJetA  ( Area );
        jecAK8_->setRho   ( Rho );
        //jecAK8_->setNPV   ( NPV );


        float corr = jecAK8_->getCorrection();

	return corr;



	}





float Mass_Corr(float  pt,float eta,float energy  , bool isdata,float Area,const double Rho,const int NPV,std::string era_)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;

  	if(isdata) 
		{
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_DATA_L2L3Residual_AK8PFchs.txt");
		}
  	else
		{
  		jecAK8PayloadNames_.push_back(era_+"_MC_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+"_MC_L3Absolute_AK8PFchs.txt");
		}
  	std::vector<JetCorrectorParameters> vPar;

  	for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) 
		{
  	  	JetCorrectorParameters pars(*ipayload);
 	   	vPar.push_back(pars);
 		}
  	jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );


        jecAK8_->setJetEta( eta );
        jecAK8_->setJetPt ( pt );
        jecAK8_->setJetE  ( energy);
        jecAK8_->setJetA  ( Area );
        jecAK8_->setRho   ( Rho );
        jecAK8_->setNPV   ( NPV );
	//std::cout<<" "<<std::endl;
	//std::cout<<"eta "<<eta<<std::endl;
	//std::cout<<"pt "<<pt<<std::endl;
	//std::cout<<"energy "<<energy<<std::endl;
	//std::cout<<"area "<<Area<<std::endl;
	//std::cout<<"rho "<<Rho<<std::endl;
	//std::cout<<"npv "<<NPV<<std::endl;

        float corr = jecAK8_->getCorrection();
	//std::cout<<"corr "<<corr<<std::endl;


	return corr;



	}













class  SlimUserData_jetsAK8 : public edm::EDProducer {
public:
  SlimUserData_jetsAK8( const edm::ParameterSet & );   

private:
  void produce( edm::Event &, const edm::EventSetup & );
  void beginJob() ;
  void endJob() ;
  bool reapplyjec_,reapplyjer_;
  std::string jes_,jer_,era_;

 };


SlimUserData_jetsAK8::SlimUserData_jetsAK8(const edm::ParameterSet& iConfig):
   reapplyjec_ (iConfig.getParameter<bool>("reapplyjec")),
   reapplyjer_ (iConfig.getParameter<bool>("reapplyjer")),
   jes_ (iConfig.getParameter<std::string>("jes")),
   jer_ (iConfig.getParameter<std::string>("jer")),
   era_ (iConfig.getParameter<std::string>("era"))


 {   

   produces<std::vector<float>>("jetAK8Pt");
   produces<std::vector<float>>("jetAK8Phi"); 
   produces<std::vector<float>>("jetAK8Eta");  
   produces<std::vector<float>>("jetAK8Mass");

   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSSmearedPt")));



  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFDown")));  


  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPt"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPhi"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSEta"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSMass"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSE"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjecFactor0"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjetArea"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetPt" )));

  edm::EDGetTokenT<int>(consumes<int>(edm::InputTag("eventUserData", "npv"    )));
  edm::EDGetTokenT<double>(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll", ""    )));


  produces<std::vector<float>>("jetAK8prunedMass");
  produces<std::vector<float>>("jetAK8softDropMass");
  produces<std::vector<float>>("jetAK8softDropMassuncorr");
  produces<std::vector<float>>("jetAK8prunedMassuncorr");

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSprunedMass"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMass"  )));

  if (jes_=="nominal"&&jer_=="nominal")
	{

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCMVA"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCMVAv2"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPartonFlavour"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSfilteredMass")));

  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopMass"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStrimmedMass"   )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjecFactor0"  )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSminmass"   )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSnSubJets" )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau1"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau2"    )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau3"    )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"  ))); 
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex1" )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac" )));





   		produces<std::vector<float>>("jetAK8Tight"); 
   		produces<std::vector<float>>("jetAK8Loose"); 
 


   		produces<std::vector<float>>("jetAK8CSV"); 
   		produces<std::vector<float>>("jetAK8CMVA"); 
   		produces<std::vector<float>>("jetAK8CMVAv2"); 
   		produces<std::vector<float>>("jetAK8PartonFlavour"); 
   		produces<std::vector<float>>("jetAK8filteredMass"); 

   		//produces<std::vector<float>>("jetAK8topMass");
   		produces<std::vector<float>>("jetAK8trimmedMass");
   		produces<std::vector<float>>("jetAK8jecFactor0");
   		//produces<std::vector<float>>("jetAK8minmass"); 
   		//produces<std::vector<float>>("jetAK8nSubJets"); 
   		produces<std::vector<float>>("jetAK8tau1");
   		produces<std::vector<float>>("jetAK8tau2");
   		produces<std::vector<float>>("jetAK8tau3");
   		produces<std::vector<float>>("jetAK8vSubjetIndex0");
   		produces<std::vector<float>>("jetAK8vSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8topSubjetIndex0");
   		//produces<std::vector<float>>("jetAK8topSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8topSubjetIndex2");
   		//produces<std::vector<float>>("jetAK8topSubjetIndex3");
		}

 }


void SlimUserData_jetsAK8::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {



  bool ISDATA = iEvent.eventAuxiliary().isRealData();







  std::auto_ptr<std::vector<float>> jetAK8Pt(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8SmearedPt(new std::vector<float>());         
        
  std::auto_ptr<std::vector<float>> jetAK8Phi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8Eta(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8Mass(new std::vector<float>()); 
      

  std::auto_ptr<std::vector<float>> jetAK8Tight(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8Loose(new std::vector<float>()); 


  std::auto_ptr<std::vector<float>> jetAK8jetArea(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8E(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8jecFactor0(new std::vector<float>()); 
  std::auto_ptr<int> npv(new int()); 
  std::auto_ptr<double> Rho(new double()); 
  std::auto_ptr<std::vector<float>> jetAK8GenJetPt(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8CSV(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8CMVA(new std::vector<float>());  
  std::auto_ptr<std::vector<float>> jetAK8CMVAv2(new std::vector<float>());              
  std::auto_ptr<std::vector<float>> jetAK8PartonFlavour(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8filteredMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8prunedMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8softDropMass(new std::vector<float>());     
  //std::auto_ptr<std::vector<float>> jetAK8topMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8trimmedMass(new std::vector<float>());       
  std::auto_ptr<std::vector<float>> jetAK8softDropMassuncorr(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8prunedMassuncorr(new std::vector<float>());

 // std::auto_ptr<std::vector<float>> jetAK8CHSJERSF(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8CHSJERSFUp(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8CHSJERSFDown(new std::vector<float>());




  //std::auto_ptr<std::vector<float>> jetAK8minmass(new std::vector<float>());      
//  std::auto_ptr<std::vector<float>> jetAK8nSubJets(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8tau1(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8tau2(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8tau3(new std::vector<float>());    
  std::auto_ptr<std::vector<float>> jetAK8vSubjetIndex0(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8vSubjetIndex1(new std::vector<float>());      
       
 // std::auto_ptr<std::vector<float>> jetAK8topSubjetIndex0(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8topSubjetIndex1(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8topSubjetIndex2(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8topSubjetIndex3(new std::vector<float>());      

  edm::Handle<std::vector<float>> jetAK8EHandle;
  edm::Handle<std::vector<float>> jetAK8jecFactor0Handle;
  edm::Handle<std::vector<float>> jetAK8jetAreaHandle;
  edm::Handle<std::vector<float>> jetAK8GenJetPtHandle;

  edm::Handle<int> npvHandle;
  edm::Handle<double> RhoHandle;

  edm::Handle<std::vector<float>> jetAK8PtHandle; 

  edm::Handle<std::vector<float>> jetAK8SmearedPtHandle;             
  edm::Handle<std::vector<float>> jetAK8PhiHandle;        
  edm::Handle<std::vector<float>> jetAK8EtaHandle;        
  edm::Handle<std::vector<float>> jetAK8MassHandle;   


  edm::Handle<std::vector<float>> jetAK8CHSJERSFHandle;
  edm::Handle<std::vector<float>> jetAK8CHSJERSFUpHandle;
  edm::Handle<std::vector<float>> jetAK8CHSJERSFDownHandle;


  edm::Handle<std::vector<float>> jetAK8CSVHandle;     
  edm::Handle<std::vector<float>> jetAK8CMVAHandle;        
  edm::Handle<std::vector<float>> jetAK8CMVAv2Handle;           
  edm::Handle<std::vector<float>> jetAK8PartonFlavourHandle;    
  edm::Handle<std::vector<float>> jetAK8filteredMassHandle;    
  edm::Handle<std::vector<float>> jetAK8prunedMassHandle;    
  edm::Handle<std::vector<float>> jetAK8softDropMassHandle;   
 // edm::Handle<std::vector<float>> jetAK8topMassHandle;    
  edm::Handle<std::vector<float>> jetAK8trimmedMassHandle;     

  //edm::Handle<std::vector<float>> jetAK8minmassHandle;    
 // edm::Handle<std::vector<float>> jetAK8nSubJetsHandle;    
  edm::Handle<std::vector<float>> jetAK8tau1Handle;       
  edm::Handle<std::vector<float>> jetAK8tau2Handle;       
  edm::Handle<std::vector<float>> jetAK8tau3Handle;  


  edm::Handle<std::vector<float>> jetAK8vSubjetIndex0Handle;    
  edm::Handle<std::vector<float>> jetAK8vSubjetIndex1Handle;  
    

  edm::Handle<std::vector<float>>  jetAK8CHSchargedHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSneutralEmEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSneutralHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSNumConstituentsHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSchargedMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSchargedEmEnergyFracHandle;  





 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex0Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex1Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex2Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex3Handle;  

  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPt"        ,jetAK8PtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSSmearedPt"        ,jetAK8SmearedPtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPhi"       ,jetAK8PhiHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSEta"       ,jetAK8EtaHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSMass"      ,jetAK8MassHandle);




  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSE"        ,jetAK8EHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"        ,jetAK8jecFactor0Handle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjetArea"        ,jetAK8jetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPt"        ,jetAK8GenJetPtHandle);

  iEvent.getByLabel("eventUserData", "npv"        ,npvHandle);
  iEvent.getByLabel("fixedGridRhoFastjetAll", ""        ,RhoHandle);

  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSF",  jetAK8CHSJERSFHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFUp",  jetAK8CHSJERSFUpHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFDown",  jetAK8CHSJERSFDownHandle);


  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSprunedMass"   ,jetAK8prunedMassHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSsoftDropMass"   ,jetAK8softDropMassHandle); 
  if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedHadronEnergyFrac"       ,jetAK8CHSchargedHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSneutralEmEnergyFrac"       ,jetAK8CHSneutralEmEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSneutralHadronEnergyFrac"       ,jetAK8CHSneutralHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSNumConstituents"       ,jetAK8CHSNumConstituentsHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedMultiplicity"       ,jetAK8CHSchargedMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac"       ,jetAK8CHSchargedEmEnergyFracHandle);  






  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCSVv2"       ,jetAK8CSVHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCMVA"       ,jetAK8CMVAHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCMVAv2"       ,jetAK8CMVAv2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPartonFlavour"   ,jetAK8PartonFlavourHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSfilteredMass"   ,jetAK8filteredMassHandle);  

  	///iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopMass"   ,jetAK8topMassHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStrimmedMass"   ,jetAK8trimmedMassHandle);   
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"   ,jetAK8jecFactor0Handle);  
  	//iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSminmass"   ,jetAK8minmassHandle);  
  	//iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSnSubJets"   ,jetAK8nSubJetsHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau1"      ,jetAK8tau1Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau2"      ,jetAK8tau2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau3"      ,jetAK8tau3Handle); 

  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"   ,jetAK8vSubjetIndex0Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex1"   ,jetAK8vSubjetIndex1Handle);

 
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex0"   ,jetAK8topSubjetIndex0Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex1"   ,jetAK8topSubjetIndex1Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex2"   ,jetAK8topSubjetIndex2Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex3"   ,jetAK8topSubjetIndex3Handle); 
 
	}








  for( size_t i=0; i<jetAK8PtHandle->size(); i++ ) 
	{
	float shift = 1.0;
	float JECcorr = 1.0;







	float uncorrpt = fmaxf(1.0,jetAK8jecFactor0Handle->at(i)*jetAK8PtHandle->at(i));	
	float uncorrE = fmaxf(1.0,jetAK8jecFactor0Handle->at(i)*jetAK8EHandle->at(i));	

   	if (reapplyjec_)
		{
		JECcorr = fmaxf(0.0,jetAK8jecFactor0Handle->at(i)*JEC_Corr(uncorrpt,jetAK8EtaHandle->at(i),ISDATA,jetAK8jetAreaHandle->at(i),*RhoHandle.product(),era_));

		}

	//float prunedjes = 1.0;
	//float prunedjer = 1.0;


	//float prunedjerunc = 1.0;
	
	float prunedshift = 1.0;


   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert(jetAK8PtHandle->at(i),jetAK8EtaHandle->at(i),jes_,era_));

		float sign = (shift-1.0)/std::fabs(shift-1.0);
		float prunedjesunc = 1.+sign*(std::sqrt( (1.-shift)*(1.-shift) + (0.002)*(0.002)));
		

		prunedshift=prunedshift*prunedjesunc;		

		}


        
	if (not ISDATA) 
		{
			float SFapp=1.0;
			if (reapplyjer_)
				{	
					JME::JetResolutionScaleFactor res_sf;
					std::string JERFile_ = era_+"_MC_SF_AK8PFchs.txt";
					res_sf = JME::JetResolutionScaleFactor(JERFile_);
				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(jetAK8PtHandle->at(i)).setJetEta(jetAK8EtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") SFapp = res_sf.getScaleFactor(jetParam);
					else if (jer_ == "up") SFapp = res_sf.getScaleFactor(jetParam, Variation::UP);
 					else if (jer_ == "down") SFapp = res_sf.getScaleFactor(jetParam, Variation::DOWN);
  
				}
			else
				{

					if (jer_ == "nominal") SFapp = jetAK8CHSJERSFHandle->at(i);
					else if (jer_ == "up") SFapp = jetAK8CHSJERSFUpHandle->at(i);
 					else if (jer_ == "down") SFapp = jetAK8CHSJERSFDownHandle->at(i);
				
				}	

			float ptJERCor = jetAK8GenJetPtHandle->at(i)+SFapp*(jetAK8PtHandle->at(i)-jetAK8GenJetPtHandle->at(i));
			shift = shift*fmaxf(0.0,ptJERCor/jetAK8GenJetPtHandle->at(i));
			prunedshift = prunedshift*fmaxf(0.0,ptJERCor/jetAK8GenJetPtHandle->at(i));


		}		

	jetAK8Pt->push_back(jetAK8PtHandle->at(i)*shift*JECcorr);      
	jetAK8Phi->push_back(jetAK8PhiHandle->at(i));       
	jetAK8Eta->push_back(jetAK8EtaHandle->at(i));       
	jetAK8Mass->push_back(jetAK8MassHandle->at(i)*shift*JECcorr);  



	float corrsdmass = Mass_Corr(uncorrpt,jetAK8EtaHandle->at(i),uncorrE,ISDATA,jetAK8jetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product(),era_);

	jetAK8softDropMass->push_back(corrsdmass*jetAK8softDropMassHandle->at(i)*shift);  
	jetAK8softDropMassuncorr->push_back(jetAK8softDropMassHandle->at(i)*shift);  
	jetAK8prunedMass->push_back(corrsdmass*0.99*jetAK8prunedMassHandle->at(i)*prunedshift);   
	jetAK8prunedMassuncorr->push_back(jetAK8prunedMassHandle->at(i)*prunedshift);   


 
   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8CSV->push_back(jetAK8CSVHandle->at(i));
		jetAK8CMVA->push_back(jetAK8CMVAHandle->at(i));       
		jetAK8CMVAv2->push_back(jetAK8CMVAv2Handle->at(i));              
		jetAK8PartonFlavour->push_back(jetAK8PartonFlavourHandle->at(i)); 
		  

		//TLV.SetPtEtaPhiM(ptcorr1,jetAK8EtaHandle->at(i),jetAK8PhiHandle->at(i),jetAK8MassHandle->at(i));
		


		jetAK8filteredMass->push_back(jetAK8filteredMassHandle->at(i));   
		//jetAK8topMass->push_back(jetAK8topMassHandle->at(i));   
		jetAK8trimmedMass->push_back(jetAK8trimmedMassHandle->at(i));    
		jetAK8jecFactor0->push_back(jetAK8jecFactor0Handle->at(i));   
		//jetAK8minmass->push_back(jetAK8minmassHandle->at(i));   
		//jetAK8nSubJets->push_back(jetAK8nSubJetsHandle->at(i));   
		jetAK8tau1->push_back(jetAK8tau1Handle->at(i));      
		jetAK8tau2->push_back(jetAK8tau2Handle->at(i));      
		jetAK8tau3->push_back(jetAK8tau3Handle->at(i));    

		jetAK8vSubjetIndex0->push_back(jetAK8vSubjetIndex0Handle->at(i));   
		jetAK8vSubjetIndex1->push_back(jetAK8vSubjetIndex1Handle->at(i)); 
  
		//jetAK8topSubjetIndex0->push_back(jetAK8topSubjetIndex0Handle->at(i));   
		//jetAK8topSubjetIndex1->push_back(jetAK8topSubjetIndex1Handle->at(i));   
		//jetAK8topSubjetIndex2->push_back(jetAK8topSubjetIndex2Handle->at(i));   
		//jetAK8topSubjetIndex3->push_back(jetAK8topSubjetIndex3Handle->at(i)); 







	  	float CHF = jetAK8CHSchargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK8CHSneutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK8CHSneutralHadronEnergyFracHandle->at(i);
	  	float NC = jetAK8CHSNumConstituentsHandle->at(i);
	  	float CM = jetAK8CHSchargedMultiplicityHandle->at(i);  
	  	float CEMF = jetAK8CHSchargedEmEnergyFracHandle->at(i);

		float TJ = 0.0;
		if ((NHF<0.9) and (NEMF<0.9) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) TJ = 1.0;

		float LJ = 0.0;
		if ((NHF<0.99) and (NEMF<0.99) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) LJ = 1.0;


		jetAK8Tight->push_back(TJ); 
		jetAK8Loose->push_back(LJ);




		}

	}


  iEvent.put(jetAK8Pt,"jetAK8Pt");
  iEvent.put(jetAK8Phi,"jetAK8Phi"); 
  iEvent.put(jetAK8Eta,"jetAK8Eta");  
  iEvent.put(jetAK8Mass,"jetAK8Mass");


  iEvent.put(jetAK8prunedMass,"jetAK8prunedMass");
  iEvent.put(jetAK8softDropMass,"jetAK8softDropMass");
  iEvent.put(jetAK8softDropMassuncorr,"jetAK8softDropMassuncorr");
  iEvent.put(jetAK8prunedMassuncorr,"jetAK8prunedMassuncorr");

   if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.put(jetAK8CSV,"jetAK8CSV"); 
  	iEvent.put(jetAK8CMVA,"jetAK8CMVA"); 
  	iEvent.put(jetAK8CMVAv2,"jetAK8CMVAv2"); 
  	iEvent.put(jetAK8PartonFlavour,"jetAK8PartonFlavour"); 
  	iEvent.put(jetAK8filteredMass,"jetAK8filteredMass"); 

  	iEvent.put(jetAK8Tight,"jetAK8Tight"); 
  	iEvent.put(jetAK8Loose,"jetAK8Loose");


  	//iEvent.put(jetAK8topMass,"jetAK8topMass");
  	iEvent.put(jetAK8trimmedMass,"jetAK8trimmedMass");
  	iEvent.put(jetAK8jecFactor0,"jetAK8jecFactor0");
  	///iEvent.put(jetAK8minmass,"jetAK8minmass"); 
  	//iEvent.put(jetAK8nSubJets,"jetAK8nSubJets"); 
  	iEvent.put(jetAK8tau1,"jetAK8tau1");
  	iEvent.put(jetAK8tau2,"jetAK8tau2");
  	iEvent.put(jetAK8tau3,"jetAK8tau3");

  	iEvent.put(jetAK8vSubjetIndex0,"jetAK8vSubjetIndex0");
  	iEvent.put(jetAK8vSubjetIndex1,"jetAK8vSubjetIndex1");

  	//iEvent.put(jetAK8topSubjetIndex0,"jetAK8topSubjetIndex0");
  	//iEvent.put(jetAK8topSubjetIndex1,"jetAK8topSubjetIndex1");
  	//iEvent.put(jetAK8topSubjetIndex2,"jetAK8topSubjetIndex2");
  	//iEvent.put(jetAK8topSubjetIndex3,"jetAK8topSubjetIndex3");

	}	
 }
// ------------ method called once each job just before starting event loop  ------------
void 
SlimUserData_jetsAK8::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_jetsAK8::endJob() {
}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_jetsAK8);
