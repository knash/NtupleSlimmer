#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include <errno.h>
#include <Math/VectorUtil.h>
//#include <TLorentzVector.h>






float JER_Uncert(float pt,float eta,std::string val,float AK8GENPT)
	{


	float pTgen = AK8GENPT;

	float regs[5][2] = {{0.0,0.8},{0.8,1.3},{1.3,1.9},{1.9,2.5},{2.5,3.0}};
	float SF[5][2] = {{1.061,0.023},{1.088,0.029},{1.106,0.030},{1.126,0.094},{1.343,0.123}};
	float SFmatch[2];
        for( size_t i=0; i<sizeof(regs); i++ ) 
		{
		if (regs[i][0]<=fabs(eta) && fabs(eta)<regs[i][1])
			{
			SFmatch[0]=SF[i][0];
			SFmatch[1]=SF[i][1];
			break;
			}
		}
		
	float SFapp;
	if (val == "nominal") SFapp = SFmatch[0];
	else if (val == "up") SFapp = SFmatch[0]+SFmatch[1];
 	else if (val == "down") SFapp = SFmatch[0]-SFmatch[1];
	else SFapp = 0;

	float ptJERCor = pTgen+SFapp*(pt-pTgen);
	//std::cout<<"eta "<<fabs(eta)<<std::endl;
	//std::cout<<"SF "<<SFapp<<std::endl;
	//std::cout<<"ptcorr "<<ptJERCor<<std::endl;
	return fmaxf(1.0,ptJERCor);

	}




float JES_Uncert(float pt,float eta,std::string val)
	{
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt");
	int sign = 0;
	if (val=="up") sign = 1;
	else if (val=="down") sign = -1;


  	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float ptCor = pt*(1+sign*jecUnc->getUncertainty(true));

	delete jecUnc;
	return  fmaxf(1.0,ptCor);
	}


float Mass_Corr(float  pt,float eta,float energy  , bool isdata,float Area,const double Rho,const int NPV)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;

  	if(isdata) 
		{
  		jecAK8PayloadNames_.push_back("Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back("Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back("Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt");
		}
  	else
		{
  		jecAK8PayloadNames_.push_back("Fall15_25nsV2_MC_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back("Fall15_25nsV2_MC_L3Absolute_AK8PFchs.txt");
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
  std::string jes_,jer_;

 };


SlimUserData_jetsAK8::SlimUserData_jetsAK8(const edm::ParameterSet& iConfig):

   jes_ (iConfig.getParameter<std::string>("jes")),
   jer_ (iConfig.getParameter<std::string>("jer"))


 {   

   produces<std::vector<float>>("jetAK8Pt");
   produces<std::vector<float>>("jetAK8Phi"); 
   produces<std::vector<float>>("jetAK8Eta");  
   produces<std::vector<float>>("jetAK8Mass");

   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSPt")));






  
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




  if (jes_=="nominal"&&jer_=="nominal")
	{

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCMVA"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCMVAv2"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPartonFlavour"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSfilteredMass")));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSprunedMass"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMass"  )));
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

 	
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopSubjetIndex0"  )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopSubjetIndex1"  )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopSubjetIndex2"   )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopSubjetIndex3" )));
 
   		produces<std::vector<float>>("jetAK8CSV"); 
   		produces<std::vector<float>>("jetAK8CMVA"); 
   		produces<std::vector<float>>("jetAK8CMVAv2"); 
   		produces<std::vector<float>>("jetAK8PartonFlavour"); 
   		produces<std::vector<float>>("jetAK8filteredMass"); 
   		produces<std::vector<float>>("jetAK8prunedMass");
   		produces<std::vector<float>>("jetAK8softDropMass");
   		produces<std::vector<float>>("jetAK8softDropMassuncorr");
   		produces<std::vector<float>>("jetAK8prunedMassuncorr");
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
  std::auto_ptr<std::vector<float>> jetAK8Phi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8Eta(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8Mass(new std::vector<float>()); 
      
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
  edm::Handle<std::vector<float>> jetAK8PhiHandle;        
  edm::Handle<std::vector<float>> jetAK8EtaHandle;        
  edm::Handle<std::vector<float>> jetAK8MassHandle;   
  
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
     
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex0Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex1Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex2Handle;    
 // edm::Handle<std::vector<float>> jetAK8topSubjetIndex3Handle;  

  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPt"        ,jetAK8PtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPhi"       ,jetAK8PhiHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSEta"       ,jetAK8EtaHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSMass"      ,jetAK8MassHandle);

  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSE"        ,jetAK8EHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"        ,jetAK8jecFactor0Handle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjetArea"        ,jetAK8jetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPt"        ,jetAK8GenJetPtHandle);

  iEvent.getByLabel("eventUserData", "npv"        ,npvHandle);
  iEvent.getByLabel("fixedGridRhoFastjetAll", ""        ,RhoHandle);




  if (jes_=="nominal"&&jer_=="nominal")
	{
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCSVv2"       ,jetAK8CSVHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCMVA"       ,jetAK8CMVAHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCMVAv2"       ,jetAK8CMVAv2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPartonFlavour"   ,jetAK8PartonFlavourHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSfilteredMass"   ,jetAK8filteredMassHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSprunedMass"   ,jetAK8prunedMassHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSsoftDropMass"   ,jetAK8softDropMassHandle); 
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
	float ptcorr = jetAK8PtHandle->at(i);
   	if (jes_!="nominal")
		{
		ptcorr = JES_Uncert(jetAK8PtHandle->at(i),jetAK8EtaHandle->at(i),jes_);
		}


	float ptcorr1 = ptcorr;



	if (not ISDATA) ptcorr1 = JER_Uncert(ptcorr,jetAK8EtaHandle->at(i),jer_,jetAK8GenJetPtHandle->at(i));
		
	jetAK8Pt->push_back(ptcorr1);      
	jetAK8Phi->push_back(jetAK8PhiHandle->at(i));       
	jetAK8Eta->push_back(jetAK8EtaHandle->at(i));       
	jetAK8Mass->push_back(jetAK8MassHandle->at(i));   
   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8CSV->push_back(jetAK8CSVHandle->at(i));
		jetAK8CMVA->push_back(jetAK8CMVAHandle->at(i));       
		jetAK8CMVAv2->push_back(jetAK8CMVAv2Handle->at(i));              
		jetAK8PartonFlavour->push_back(jetAK8PartonFlavourHandle->at(i)); 
		  

		//TLV.SetPtEtaPhiM(ptcorr1,jetAK8EtaHandle->at(i),jetAK8PhiHandle->at(i),jetAK8MassHandle->at(i));
		
		float uncorrpt = fmaxf(1.0,jetAK8jecFactor0Handle->at(i)*ptcorr1);	
		float uncorrE = fmaxf(1.0,jetAK8jecFactor0Handle->at(i)*jetAK8EHandle->at(i));	
		float corrsdmass = Mass_Corr(uncorrpt,jetAK8EtaHandle->at(i),uncorrE,ISDATA,jetAK8jetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product());
		jetAK8softDropMass->push_back(corrsdmass*jetAK8softDropMassHandle->at(i));  
		jetAK8softDropMassuncorr->push_back(jetAK8softDropMassHandle->at(i));  
		jetAK8prunedMass->push_back(corrsdmass*jetAK8prunedMassHandle->at(i));   
		jetAK8prunedMassuncorr->push_back(jetAK8prunedMassHandle->at(i));   
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
		}

	}


  iEvent.put(jetAK8Pt,"jetAK8Pt");
  iEvent.put(jetAK8Phi,"jetAK8Phi"); 
  iEvent.put(jetAK8Eta,"jetAK8Eta");  
  iEvent.put(jetAK8Mass,"jetAK8Mass");


   if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.put(jetAK8CSV,"jetAK8CSV"); 
  	iEvent.put(jetAK8CMVA,"jetAK8CMVA"); 
  	iEvent.put(jetAK8CMVAv2,"jetAK8CMVAv2"); 
  	iEvent.put(jetAK8PartonFlavour,"jetAK8PartonFlavour"); 
  	iEvent.put(jetAK8filteredMass,"jetAK8filteredMass"); 
  	iEvent.put(jetAK8prunedMass,"jetAK8prunedMass");
  	iEvent.put(jetAK8softDropMass,"jetAK8softDropMass");
  	iEvent.put(jetAK8softDropMassuncorr,"jetAK8softDropMassuncorr");
  	iEvent.put(jetAK8prunedMassuncorr,"jetAK8prunedMassuncorr");

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
