
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

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
//#include <errno.h>
#include <Math/VectorUtil.h>
#include <TRandom3.h>

//#include <TLorentzVector.h>




float JES_Uncert_CHS(float pt,float eta,std::string val,std::string era_, unsigned int runnum)
	{
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(era_+"V2_MC_Uncertainty_AK8PFchs.txt");
	int sign = 0;
	if (val=="up") sign = 1;
	else if (val=="down") sign = -1;


  	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float shift = (1+sign*jecUnc->getUncertainty(true));

	delete jecUnc;
	return  shift;
	}

float JES_Uncert_Puppi(float pt,float eta,std::string val,std::string era_, unsigned int runnum)
	{
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(era_+"V2_MC_Uncertainty_AK8PFPuppi.txt");
	int sign = 0;
	if (val=="up") sign = 1;
	else if (val=="down") sign = -1;


  	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float shift = (1+sign*jecUnc->getUncertainty(true));

	delete jecUnc;
	return  shift;
	}

float JEC_Corr_CHS(float  pt,float eta  , bool isdata,float Area,const double Rho,std::string era_,unsigned int runnum)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;
	std::string runtxt_;
  	if(isdata) 
		{

		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV2";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV2";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV2"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV2"; //IOV H:[280919,Infinity] f



 		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L1FastJet_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2L3Residual_AK8PFchs.txt");
		}
  	else
		{

		runtxt_ = "V2";
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L1FastJet_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L3Absolute_AK8PFchs.txt");
		}
	//std::cout<<runtxt_<<std::endl;
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

float JEC_Corr_Puppi(float  pt,float eta  , bool isdata,float Area,const double Rho,std::string era_,unsigned int runnum)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;
	std::string runtxt_;
  	if(isdata) 
		{

		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV2";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV2";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV2"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV2"; //IOV H:[280919,Infinity] f



  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2Relative_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L3Absolute_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2L3Residual_AK8PFPuppi.txt");
		}
  	else
		{

		runtxt_ = "V2";
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L2Relative_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L3Absolute_AK8PFPuppi.txt");
		}
	//std::cout<<runtxt_<<std::endl;
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



float Mass_Corr_CHS(float  pt,float eta,float energy  , bool isdata,float Area,const double Rho,const int NPV,std::string era_, unsigned int runnum)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;
	std::string runtxt_;
  	if(isdata) 
		{



		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV2";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV2";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV2"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV2"; //IOV H:[280919,Infinity] f




  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2L3Residual_AK8PFchs.txt");
		}
  	else
		{		
		runtxt_ = "V2";
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L2Relative_AK8PFchs.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L3Absolute_AK8PFchs.txt");
		}
  	std::vector<JetCorrectorParameters> vPar;
	//std::cout<<runtxt_<<std::endl;
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


float Mass_Corr_Puppi(float  pt,float eta,float energy  , bool isdata,float Area,const double Rho,const int NPV,std::string era_, unsigned int runnum)
	{
        boost::shared_ptr<FactorizedJetCorrector> jecAK8_;
  	std::vector<std::string> jecAK8PayloadNames_;
	std::string runtxt_;
  	if(isdata) 
		{



		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV2";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV2";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV2"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV2"; //IOV H:[280919,Infinity] f




  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2Relative_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L3Absolute_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_DATA_L2L3Residual_AK8PFPuppi.txt");
		}
  	else
		{		
		runtxt_ = "V2";
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L2Relative_AK8PFPuppi.txt");
  		jecAK8PayloadNames_.push_back(era_+runtxt_+"_MC_L3Absolute_AK8PFPuppi.txt");
		}
  	std::vector<JetCorrectorParameters> vPar;
	//std::cout<<runtxt_<<std::endl;
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

   produces<std::vector<float>>("jetAK8CHSPt");
   produces<std::vector<float>>("jetAK8CHSPhi"); 
   produces<std::vector<float>>("jetAK8CHSEta");  
   produces<std::vector<float>>("jetAK8CHSMass");

   produces<std::vector<float>>("jetAK8PuppiPt");
   produces<std::vector<float>>("jetAK8PuppiPhi"); 
   produces<std::vector<float>>("jetAK8PuppiEta");  
   produces<std::vector<float>>("jetAK8PuppiMass");

   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSSmearedPt")));
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiSmearedPt")));



  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSPtResolution")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFDown")));  


  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPt"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPhi"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSEta"    )));
//  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSMass"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSE"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjecFactor0"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjetArea"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetPt" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetEta" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetPhi" )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiPtResolution")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSFDown")));  


  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiPt"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiPhi"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiEta"    )));
//  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiMass"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiE"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppijecFactor0"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppijetArea"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiGenJetPt" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiGenJetEta" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiGenJetPhi" )));

  edm::EDGetTokenT<int>(consumes<int>(edm::InputTag("vertexInfo", "npv"    )));
  edm::EDGetTokenT<double>(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll", ""    )));


  produces<std::vector<float>>("jetAK8CHSprunedMass");
  produces<std::vector<float>>("jetAK8CHSsoftDropMass");
  produces<std::vector<float>>("jetAK8CHSsoftDropMassuncorr");
  produces<std::vector<float>>("jetAK8CHSprunedMassuncorr");

  produces<std::vector<float>>("jetAK8PuppiprunedMass");
  produces<std::vector<float>>("jetAK8PuppisoftDropMass");
  produces<std::vector<float>>("jetAK8PuppisoftDropMassuncorr");
  produces<std::vector<float>>("jetAK8PuppiprunedMassuncorr");

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSprunedMassCHS"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS"  )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiprunedMass"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppisoftDropMass"  )));


  if (jes_=="nominal"&&jer_=="nominal")
	{

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSCMVAv2"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPartonFlavour"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSfilteredMassCHS")));

  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStopMass"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStrimmedMassCHS"   )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjecFactor0"  )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSminmass"   )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSnSubJets" )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau1CHS"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau2CHS"    )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHStau3CHS"    )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"  ))); 
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex1" )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac" )));




  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiCMVAv2"     )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiPartonFlavour"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppifilteredMass")));

  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppitopMass"  )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppitrimmedMass"   )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppijecFactor0"  )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8Puppiminmass"   )));
  	   	//edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppinSubJets" )));  
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8Puppitau1"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8Puppitau2"    )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8Puppitau3"    )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppivSubjetIndex0"  ))); 
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppivSubjetIndex1" )));

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedEmEnergyFrac" )));





   		produces<std::vector<float>>("jetAK8CHSTight"); 
   		produces<std::vector<float>>("jetAK8CHSLoose"); 

  		produces<std::vector<float>>("jetAK8PuppiTight"); 
   		produces<std::vector<float>>("jetAK8PuppiLoose"); 
 


   		produces<std::vector<float>>("jetAK8CHSCSV"); 
   		produces<std::vector<float>>("jetAK8CHSCMVAv2"); 
   		produces<std::vector<float>>("jetAK8CHSPartonFlavour"); 
   		produces<std::vector<float>>("jetAK8CHSfilteredMass"); 

   		//produces<std::vector<float>>("jetAK8CHStopMass");
   		produces<std::vector<float>>("jetAK8CHStrimmedMass");
   		produces<std::vector<float>>("jetAK8CHSjecFactor0");
   		//produces<std::vector<float>>("jetAK8CHSminmass"); 
   		//produces<std::vector<float>>("jetAK8CHSnSubJets"); 
   		produces<std::vector<float>>("jetAK8CHStau1");
   		produces<std::vector<float>>("jetAK8CHStau2");
   		produces<std::vector<float>>("jetAK8CHStau3");
   		produces<std::vector<float>>("jetAK8CHSvSubjetIndex0");
   		produces<std::vector<float>>("jetAK8CHSvSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex0");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex2");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex3");


   		produces<std::vector<float>>("jetAK8PuppiCSV"); 
   		produces<std::vector<float>>("jetAK8PuppiCMVAv2"); 
   		produces<std::vector<float>>("jetAK8PuppiPartonFlavour"); 
   		produces<std::vector<float>>("jetAK8PuppifilteredMass"); 

   		//produces<std::vector<float>>("jetAK8PuppitopMass");
   		produces<std::vector<float>>("jetAK8PuppitrimmedMass");
   		produces<std::vector<float>>("jetAK8PuppijecFactor0");
   		//produces<std::vector<float>>("jetAK8Puppiminmass"); 
   		//produces<std::vector<float>>("jetAK8PuppinSubJets"); 
   		produces<std::vector<float>>("jetAK8Puppitau1");
   		produces<std::vector<float>>("jetAK8Puppitau2");
   		produces<std::vector<float>>("jetAK8Puppitau3");
   		produces<std::vector<float>>("jetAK8PuppivSubjetIndex0");
   		produces<std::vector<float>>("jetAK8PuppivSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8PuppitopSubjetIndex0");
   		//produces<std::vector<float>>("jetAK8PuppitopSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8PuppitopSubjetIndex2");
   		//produces<std::vector<float>>("jetAK8PuppitopSubjetIndex3");



		}

 }



void SlimUserData_jetsAK8::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {



  bool ISDATA = iEvent.eventAuxiliary().isRealData();
  unsigned int runnum = iEvent.eventAuxiliary().run();


  //std::cout<<runnum<<std::endl;



  std::auto_ptr<std::vector<float>> jetAK8CHSPt(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8CHSSmearedPt(new std::vector<float>());         
        
  std::auto_ptr<std::vector<float>> jetAK8CHSPhi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8CHSEta(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8CHSMass(new std::vector<float>()); 
      

  std::auto_ptr<std::vector<float>> jetAK8CHSTight(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSLoose(new std::vector<float>()); 


  std::auto_ptr<std::vector<float>> jetAK8CHSjetArea(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSE(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSjecFactor0(new std::vector<float>()); 
  std::auto_ptr<int> npv(new int()); 
  std::auto_ptr<double> Rho(new double()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetEta(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetPhi(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8CHSCSV(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8CHSCMVAv2(new std::vector<float>());              
  std::auto_ptr<std::vector<float>> jetAK8CHSPartonFlavour(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSfilteredMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSprunedMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSsoftDropMass(new std::vector<float>());     
  //std::auto_ptr<std::vector<float>> jetAK8CHStopMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHStrimmedMass(new std::vector<float>());       
  std::auto_ptr<std::vector<float>> jetAK8CHSsoftDropMassuncorr(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8CHSprunedMassuncorr(new std::vector<float>());

 // std::auto_ptr<std::vector<float>> jetAK8CHSJERSF(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8CHSJERSFUp(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8CHSJERSFDown(new std::vector<float>());




  //std::auto_ptr<std::vector<float>> jetAK8CHSminmass(new std::vector<float>());      
//  std::auto_ptr<std::vector<float>> jetAK8CHSnSubJets(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHStau1(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8CHStau2(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8CHStau3(new std::vector<float>());    
  std::auto_ptr<std::vector<float>> jetAK8CHSvSubjetIndex0(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSvSubjetIndex1(new std::vector<float>());      
       
 // std::auto_ptr<std::vector<float>> jetAK8CHStopSubjetIndex0(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8CHStopSubjetIndex1(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8CHStopSubjetIndex2(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8CHStopSubjetIndex3(new std::vector<float>());      

  edm::Handle<std::vector<float>> jetAK8CHSEHandle;
  edm::Handle<std::vector<float>> jetAK8CHSjecFactor0Handle;
  edm::Handle<std::vector<float>> jetAK8CHSjetAreaHandle;
  edm::Handle<std::vector<float>> jetAK8CHSGenJetPtHandle;
  edm::Handle<std::vector<float>> jetAK8CHSGenJetEtaHandle;
  edm::Handle<std::vector<float>> jetAK8CHSGenJetPhiHandle;

  edm::Handle<int> npvHandle;
  edm::Handle<double> RhoHandle;

  edm::Handle<std::vector<float>> jetAK8CHSPtHandle; 

  edm::Handle<std::vector<float>> jetAK8CHSSmearedPtHandle;             
  edm::Handle<std::vector<float>> jetAK8CHSPhiHandle;        
  edm::Handle<std::vector<float>> jetAK8CHSEtaHandle;        
  edm::Handle<std::vector<float>> jetAK8CHSMassHandle;   


  edm::Handle<std::vector<float>> jetAK8CHSJERSFHandle;
  edm::Handle<std::vector<float>> jetAK8CHSJERSFUpHandle;
  edm::Handle<std::vector<float>> jetAK8CHSJERSFDownHandle;

  edm::Handle<std::vector<float>> jetAK8CHSRESHandle;

  edm::Handle<std::vector<float>> jetAK8CHSCSVHandle;     
  edm::Handle<std::vector<float>> jetAK8CHSCMVAv2Handle;           
  edm::Handle<std::vector<float>> jetAK8CHSPartonFlavourHandle;    
  edm::Handle<std::vector<float>> jetAK8CHSfilteredMassHandle;    
  edm::Handle<std::vector<float>> jetAK8CHSprunedMassHandle;    
  edm::Handle<std::vector<float>> jetAK8CHSsoftDropMassHandle;   
 // edm::Handle<std::vector<float>> jetAK8CHStopMassHandle;    
  edm::Handle<std::vector<float>> jetAK8CHStrimmedMassHandle;     

  //edm::Handle<std::vector<float>> jetAK8CHSminmassHandle;    
 // edm::Handle<std::vector<float>> jetAK8CHSnSubJetsHandle;    
  edm::Handle<std::vector<float>> jetAK8CHStau1Handle;       
  edm::Handle<std::vector<float>> jetAK8CHStau2Handle;       
  edm::Handle<std::vector<float>> jetAK8CHStau3Handle;  


  edm::Handle<std::vector<float>> jetAK8CHSvSubjetIndex0Handle;    
  edm::Handle<std::vector<float>> jetAK8CHSvSubjetIndex1Handle;  
    

  edm::Handle<std::vector<float>>  jetAK8CHSchargedHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSneutralEmEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSneutralHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSNumConstituentsHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSchargedMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSneutralMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK8CHSchargedEmEnergyFracHandle;  





 // edm::Handle<std::vector<float>> jetAK8CHStopSubjetIndex0Handle;    
 // edm::Handle<std::vector<float>> jetAK8CHStopSubjetIndex1Handle;    
 // edm::Handle<std::vector<float>> jetAK8CHStopSubjetIndex2Handle;    
 // edm::Handle<std::vector<float>> jetAK8CHStopSubjetIndex3Handle;  





  std::auto_ptr<std::vector<float>> jetAK8PuppiPt(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8PuppiSmearedPt(new std::vector<float>());         
        
  std::auto_ptr<std::vector<float>> jetAK8PuppiPhi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8PuppiEta(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8PuppiMass(new std::vector<float>()); 
      

  std::auto_ptr<std::vector<float>> jetAK8PuppiTight(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppiLoose(new std::vector<float>()); 


  std::auto_ptr<std::vector<float>> jetAK8PuppijetArea(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppiE(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppijecFactor0(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppiGenJetPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppiGenJetEta(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8PuppiGenJetPhi(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8PuppiCSV(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8PuppiCMVAv2(new std::vector<float>());              
  std::auto_ptr<std::vector<float>> jetAK8PuppiPartonFlavour(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppifilteredMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppiprunedMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMass(new std::vector<float>());     
  //std::auto_ptr<std::vector<float>> jetAK8PuppitopMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppitrimmedMass(new std::vector<float>());       
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMassuncorr(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8PuppiprunedMassuncorr(new std::vector<float>());

 // std::auto_ptr<std::vector<float>> jetAK8PuppiJERSF(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8PuppiJERSFUp(new std::vector<float>());
  //std::auto_ptr<std::vector<float>> jetAK8PuppiJERSFDown(new std::vector<float>());




  //std::auto_ptr<std::vector<float>> jetAK8Puppiminmass(new std::vector<float>());      
//  std::auto_ptr<std::vector<float>> jetAK8PuppinSubJets(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8Puppitau1(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8Puppitau2(new std::vector<float>());         
  std::auto_ptr<std::vector<float>> jetAK8Puppitau3(new std::vector<float>());    
  std::auto_ptr<std::vector<float>> jetAK8PuppivSubjetIndex0(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppivSubjetIndex1(new std::vector<float>());      
       
 // std::auto_ptr<std::vector<float>> jetAK8PuppitopSubjetIndex0(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8PuppitopSubjetIndex1(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8PuppitopSubjetIndex2(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8PuppitopSubjetIndex3(new std::vector<float>());      

  edm::Handle<std::vector<float>> jetAK8PuppiEHandle;
  edm::Handle<std::vector<float>> jetAK8PuppijecFactor0Handle;
  edm::Handle<std::vector<float>> jetAK8PuppijetAreaHandle;
  edm::Handle<std::vector<float>> jetAK8PuppiGenJetPtHandle;
  edm::Handle<std::vector<float>> jetAK8PuppiGenJetEtaHandle;
  edm::Handle<std::vector<float>> jetAK8PuppiGenJetPhiHandle;

  edm::Handle<std::vector<float>> jetAK8PuppiPtHandle; 

  edm::Handle<std::vector<float>> jetAK8PuppiSmearedPtHandle;             
  edm::Handle<std::vector<float>> jetAK8PuppiPhiHandle;        
  edm::Handle<std::vector<float>> jetAK8PuppiEtaHandle;        
  edm::Handle<std::vector<float>> jetAK8PuppiMassHandle;   


  edm::Handle<std::vector<float>> jetAK8PuppiJERSFHandle;
  edm::Handle<std::vector<float>> jetAK8PuppiJERSFUpHandle;
  edm::Handle<std::vector<float>> jetAK8PuppiJERSFDownHandle;

  edm::Handle<std::vector<float>> jetAK8PuppiRESHandle;

  edm::Handle<std::vector<float>> jetAK8PuppiCSVHandle;     
  edm::Handle<std::vector<float>> jetAK8PuppiCMVAv2Handle;           
  edm::Handle<std::vector<float>> jetAK8PuppiPartonFlavourHandle;    
  edm::Handle<std::vector<float>> jetAK8PuppifilteredMassHandle;    
  edm::Handle<std::vector<float>> jetAK8PuppiprunedMassHandle;    
  edm::Handle<std::vector<float>> jetAK8PuppisoftDropMassHandle;   
 // edm::Handle<std::vector<float>> jetAK8PuppitopMassHandle;    
  edm::Handle<std::vector<float>> jetAK8PuppitrimmedMassHandle;     

  //edm::Handle<std::vector<float>> jetAK8PuppiminmassHandle;    
 // edm::Handle<std::vector<float>> jetAK8PuppinSubJetsHandle;    
  edm::Handle<std::vector<float>> jetAK8Puppitau1Handle;       
  edm::Handle<std::vector<float>> jetAK8Puppitau2Handle;       
  edm::Handle<std::vector<float>> jetAK8Puppitau3Handle;  


  edm::Handle<std::vector<float>> jetAK8PuppivSubjetIndex0Handle;    
  edm::Handle<std::vector<float>> jetAK8PuppivSubjetIndex1Handle;  
    

  edm::Handle<std::vector<float>>  jetAK8PuppichargedHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppineutralEmEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppineutralHadronEnergyFracHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppiNumConstituentsHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppichargedMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppineutralMultiplicityHandle;  
  edm::Handle<std::vector<float>>  jetAK8PuppichargedEmEnergyFracHandle;  





 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex0Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex1Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex2Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex3Handle;  


  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPt"        ,jetAK8CHSPtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSSmearedPt"        ,jetAK8CHSSmearedPtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPhi"       ,jetAK8CHSPhiHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSEta"       ,jetAK8CHSEtaHandle);  
 // iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSMass"      ,jetAK8MassHandle);




  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSE"        ,jetAK8CHSEHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"        ,jetAK8CHSjecFactor0Handle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjetArea"        ,jetAK8CHSjetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPt"        ,jetAK8CHSGenJetPtHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetEta"        ,jetAK8CHSGenJetEtaHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPhi"        ,jetAK8CHSGenJetPhiHandle);


  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiPt"        ,jetAK8PuppiPtHandle);
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiSmearedPt"        ,jetAK8PuppiSmearedPtHandle);
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiPhi"       ,jetAK8PuppiPhiHandle);  
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiEta"       ,jetAK8PuppiEtaHandle);  
 // iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiMass"      ,jetAK8PuppiMassHandle);




  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiE"        ,jetAK8PuppiEHandle);
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppijecFactor0"        ,jetAK8PuppijecFactor0Handle);
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppijetArea"        ,jetAK8PuppijetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiGenJetPt"        ,jetAK8PuppiGenJetPtHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiGenJetEta"        ,jetAK8PuppiGenJetEtaHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiGenJetPhi"        ,jetAK8PuppiGenJetPhiHandle);


  iEvent.getByLabel("vertexInfo", "npv"        ,npvHandle);
  iEvent.getByLabel("fixedGridRhoFastjetAll", ""        ,RhoHandle);


  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSPtResolution",  jetAK8CHSRESHandle);

  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSF",  jetAK8CHSJERSFHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFUp",  jetAK8CHSJERSFUpHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFDown",  jetAK8CHSJERSFDownHandle);


  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSprunedMassCHS"   ,jetAK8CHSprunedMassHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS"   ,jetAK8CHSsoftDropMassHandle); 


  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiPtResolution",  jetAK8PuppiRESHandle);

  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSF",  jetAK8PuppiJERSFHandle);
  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSFUp",  jetAK8PuppiJERSFUpHandle);
  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSFDown",  jetAK8PuppiJERSFDownHandle);


  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiprunedMass"   ,jetAK8PuppiprunedMassHandle);  
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppisoftDropMass"   ,jetAK8PuppisoftDropMassHandle); 


  if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedHadronEnergyFrac"       ,jetAK8CHSchargedHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSneutralEmEnergyFrac"       ,jetAK8CHSneutralEmEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSneutralHadronEnergyFrac"       ,jetAK8CHSneutralHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSNumConstituents"       ,jetAK8CHSNumConstituentsHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedMultiplicity"       ,jetAK8CHSchargedMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSneutralMultiplicity"       ,jetAK8CHSneutralMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac"       ,jetAK8CHSchargedEmEnergyFracHandle);  



  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCSVv2"       ,jetAK8CHSCSVHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSCMVAv2"       ,jetAK8CHSCMVAv2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPartonFlavour"   ,jetAK8CHSPartonFlavourHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSfilteredMassCHS"   ,jetAK8CHSfilteredMassHandle);  

  	///iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopMass"   ,jetAK8CHStopMassHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStrimmedMassCHS"   ,jetAK8CHStrimmedMassHandle);   
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"   ,jetAK8CHSjecFactor0Handle);  
  	//iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSminmass"   ,jetAK8CHSminmassHandle);  
  	//iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSnSubJets"   ,jetAK8CHSnSubJetsHandle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau1CHS"      ,jetAK8CHStau1Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau2CHS"      ,jetAK8CHStau2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau3CHS"      ,jetAK8CHStau3Handle); 

  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"   ,jetAK8CHSvSubjetIndex0Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex1"   ,jetAK8CHSvSubjetIndex1Handle);

 
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex0"   ,jetAK8CHStopSubjetIndex0Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex1"   ,jetAK8CHStopSubjetIndex1Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex2"   ,jetAK8CHStopSubjetIndex2Handle);  
  //	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStopSubjetIndex3"   ,jetAK8CHStopSubjetIndex3Handle); 



  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedHadronEnergyFrac"       ,jetAK8PuppichargedHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralEmEnergyFrac"       ,jetAK8PuppineutralEmEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralHadronEnergyFrac"       ,jetAK8PuppineutralHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiNumConstituents"       ,jetAK8PuppiNumConstituentsHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedMultiplicity"       ,jetAK8PuppichargedMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralMultiplicity"       ,jetAK8PuppineutralMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedEmEnergyFrac"       ,jetAK8PuppichargedEmEnergyFracHandle);  



  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiCSVv2"       ,jetAK8PuppiCSVHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiCMVAv2"       ,jetAK8PuppiCMVAv2Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiPartonFlavour"   ,jetAK8PuppiPartonFlavourHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppifilteredMass"   ,jetAK8PuppifilteredMassHandle);  

  	///iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitopMass"   ,jetAK8PuppitopMassHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitrimmedMass"   ,jetAK8PuppitrimmedMassHandle);   
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppijecFactor0"   ,jetAK8PuppijecFactor0Handle);  
  	//iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppiminmass"   ,jetAK8PuppiminmassHandle);  
  	//iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppinSubJets"   ,jetAK8PuppinSubJetsHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau1"      ,jetAK8Puppitau1Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau2"      ,jetAK8Puppitau2Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau3"      ,jetAK8Puppitau3Handle); 

  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppivSubjetIndex0"   ,jetAK8PuppivSubjetIndex0Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppivSubjetIndex1"   ,jetAK8PuppivSubjetIndex1Handle);

 
  //	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitopSubjetIndex0"   ,jetAK8PuppitopSubjetIndex0Handle);  
  //	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitopSubjetIndex1"   ,jetAK8PuppitopSubjetIndex1Handle);  
  //	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitopSubjetIndex2"   ,jetAK8PuppitopSubjetIndex2Handle);  
  //	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitopSubjetIndex3"   ,jetAK8PuppitopSubjetIndex3Handle); 
 
	}





  //CHS
  for( size_t i=0; i<jetAK8CHSPtHandle->size(); i++ ) 
	{
	float shift = 1.0;
	float JECcorr = 1.0;

        TLorentzVector v1;
	v1.SetPtEtaPhiE(jetAK8CHSPtHandle->at(i),jetAK8CHSEtaHandle->at(i),jetAK8CHSPhiHandle->at(i),jetAK8CHSEHandle->at(i));
	float calcmass =  v1.M();




	float uncorrpt = fmaxf(1.0,jetAK8CHSjecFactor0Handle->at(i)*jetAK8CHSPtHandle->at(i));	
	float uncorrE = fmaxf(1.0,jetAK8CHSjecFactor0Handle->at(i)*jetAK8CHSEHandle->at(i));	

   	if (reapplyjec_)
		{
		JECcorr = fmaxf(0.0,jetAK8CHSjecFactor0Handle->at(i)*JEC_Corr_CHS(uncorrpt,jetAK8CHSEtaHandle->at(i),ISDATA,jetAK8CHSjetAreaHandle->at(i),*RhoHandle.product(),era_, runnum));
		}

	//float prunedjes = 1.0;
	//float prunedjer = 1.0;


	//float prunedjerunc = 1.0;
	
	float prunedshift = 1.0;


   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert_CHS(jetAK8CHSPtHandle->at(i),jetAK8CHSEtaHandle->at(i),jes_,era_,runnum));

		float sign = (shift-1.0)/std::fabs(shift-1.0);
		float prunedjesunc = 1.+sign*(std::sqrt( (1.-shift)*(1.-shift) + (0.002)*(0.002)));
		

		prunedshift=prunedshift*prunedjesunc;		

		}


        
	if (not ISDATA) 
		{

			float SFapp=1.0;
			float RESapp=1.0;
			if (reapplyjer_)
				{	
					//JME::JetResolutionScaleFactor res_sf;
					//JME::JetResolution reso;
					//std::string JERFile_ = era_+"_MC_SF_AK8PFchs.txt";
					



/////////////////////////////////////////////////////TOUPDATE!	
					//Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt
					//Spring16_25nsV10_MC_SF_AK8PFchs.txt			
					std::string JERFile_ = "Spring16_25nsV10_MC_SF_AK8PFchs.txt";
					std::string RESFile_ = "Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";
/////////////////////////////////////////////////////TOUPDATE!	


					JME::JetResolutionScaleFactor res_sf = JME::JetResolutionScaleFactor(JERFile_);

					JME::JetResolution reso = JME::JetResolution(RESFile_);

				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(jetAK8CHSPtHandle->at(i)).setJetEta(jetAK8CHSEtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") SFapp = res_sf.getScaleFactor(jetParam);
					else if (jer_ == "up") SFapp = res_sf.getScaleFactor(jetParam, Variation::UP);
 					else if (jer_ == "down") SFapp = res_sf.getScaleFactor(jetParam, Variation::DOWN);

					RESapp = reso.getResolution(jetParam)*jetAK8CHSPtHandle->at(i);

				}
			else
				{

					if (jer_ == "nominal") SFapp = jetAK8CHSJERSFHandle->at(i);
					else if (jer_ == "up") SFapp = jetAK8CHSJERSFUpHandle->at(i);
 					else if (jer_ == "down") SFapp = jetAK8CHSJERSFDownHandle->at(i);
					RESapp = jetAK8CHSRESHandle->at(i)*jetAK8CHSPtHandle->at(i);

				}	

			bool DRmatch =  deltaR(jetAK8CHSEtaHandle->at(i), jetAK8CHSPhiHandle->at(i), jetAK8CHSGenJetEtaHandle->at(i), jetAK8CHSGenJetPhiHandle->at(i))<0.4/2.0 ;
			bool ptmatch = std::fabs(jetAK8CHSPtHandle->at(i)-jetAK8CHSGenJetPtHandle->at(i))<3.0*RESapp;
			float ptJERCor = jetAK8CHSGenJetPtHandle->at(i)+SFapp*(jetAK8CHSPtHandle->at(i)-jetAK8CHSGenJetPtHandle->at(i));

			if (DRmatch and ptmatch)
				{
					shift = shift*fmaxf(0.0,ptJERCor/jetAK8CHSGenJetPtHandle->at(i));
					prunedshift = prunedshift*fmaxf(0.0,ptJERCor/jetAK8CHSGenJetPtHandle->at(i));
				}
			else
				{
					//std::cout<<"Hybrid: "<<std::endl;
					//std::cout<<"DeltaR match: "<<deltaR(jetAK8EtaHandle->at(i), jetAK8PhiHandle->at(i), jetAK8GenJetEtaHandle->at(i), jetAK8GenJetPhiHandle->at(i)) <<std::endl;
					//std::cout<<"DeltaPt match: "<<std::fabs(jetAK8PtHandle->at(i)-jetAK8GenJetPtHandle->at(i)) <<std::endl;
					//std::cout<<"RES: "<<RESapp<<std::endl;
					//std::cout<<"RES threesig: "<<3.0*RESapp<<std::endl;
					//std::cout<<"SF: "<<SFapp<<std::endl;
					Double_t sigma =std::sqrt( SFapp*SFapp-1 ) * RESapp;// √(SF^2-1) * sigma_MC_PT.;
					TRandom3 *r = new TRandom3(0); 
					Double_t smearify = r->Gaus(0.0,sigma); 

					shift = shift*fmaxf(0.0,(1.0+smearify/jetAK8CHSPtHandle->at(i)));
					prunedshift = prunedshift*fmaxf(0.0,(1.0+smearify/jetAK8CHSPtHandle->at(i)));
					
				}


		}		
	
	jetAK8CHSPt->push_back(jetAK8CHSPtHandle->at(i)*shift*JECcorr);      
	jetAK8CHSPhi->push_back(jetAK8CHSPhiHandle->at(i));       
	jetAK8CHSEta->push_back(jetAK8CHSEtaHandle->at(i));       


	jetAK8CHSMass->push_back(calcmass*shift*JECcorr);  



	float corrsdmass = Mass_Corr_CHS(uncorrpt,jetAK8CHSEtaHandle->at(i),uncorrE,ISDATA,jetAK8CHSjetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product(),era_,runnum);

	jetAK8CHSsoftDropMass->push_back(corrsdmass*jetAK8CHSsoftDropMassHandle->at(i)*shift);  
	jetAK8CHSsoftDropMassuncorr->push_back(jetAK8CHSsoftDropMassHandle->at(i)*shift);  
	jetAK8CHSprunedMass->push_back(corrsdmass*0.99*jetAK8CHSprunedMassHandle->at(i)*prunedshift);   
	jetAK8CHSprunedMassuncorr->push_back(jetAK8CHSprunedMassHandle->at(i)*prunedshift);   


 
   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8CHSCSV->push_back(jetAK8CHSCSVHandle->at(i));
		jetAK8CHSCMVAv2->push_back(jetAK8CHSCMVAv2Handle->at(i));              
		jetAK8CHSPartonFlavour->push_back(jetAK8CHSPartonFlavourHandle->at(i)); 
		  

		//TLV.SetPtEtaPhiM(ptcorr1,jetAK8EtaHandle->at(i),jetAK8PhiHandle->at(i),jetAK8MassHandle->at(i));
		


		jetAK8CHSfilteredMass->push_back(jetAK8CHSfilteredMassHandle->at(i));   
		//jetAK8CHStopMass->push_back(jetAK8CHStopMassHandle->at(i));   
		jetAK8CHStrimmedMass->push_back(jetAK8CHStrimmedMassHandle->at(i));    
		jetAK8CHSjecFactor0->push_back(jetAK8CHSjecFactor0Handle->at(i));   
		//jetAK8CHSminmass->push_back(jetAK8CHSminmassHandle->at(i));   
		//jetAK8CHSnSubJets->push_back(jetAK8CHSnSubJetsHandle->at(i));   
		jetAK8CHStau1->push_back(jetAK8CHStau1Handle->at(i));      
		jetAK8CHStau2->push_back(jetAK8CHStau2Handle->at(i));      
		jetAK8CHStau3->push_back(jetAK8CHStau3Handle->at(i));    

		jetAK8CHSvSubjetIndex0->push_back(jetAK8CHSvSubjetIndex0Handle->at(i));   
		jetAK8CHSvSubjetIndex1->push_back(jetAK8CHSvSubjetIndex1Handle->at(i)); 
  
		//jetAK8CHStopSubjetIndex0->push_back(jetAK8CHStopSubjetIndex0Handle->at(i));   
		//jetAK8CHStopSubjetIndex1->push_back(jetAK8CHStopSubjetIndex1Handle->at(i));   
		//jetAK8CHStopSubjetIndex2->push_back(jetAK8CHStopSubjetIndex2Handle->at(i));   
		//jetAK8CHStopSubjetIndex3->push_back(jetAK8CHStopSubjetIndex3Handle->at(i)); 







	  	float CHF = jetAK8CHSchargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK8CHSneutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK8CHSneutralHadronEnergyFracHandle->at(i);
	  	//float NC = jetAK8CHSNumConstituentsHandle->at(i);
	  	float CM = jetAK8CHSchargedMultiplicityHandle->at(i);  
	  	float CEMF = jetAK8CHSchargedEmEnergyFracHandle->at(i);

		float NC = jetAK8CHSneutralMultiplicityHandle->at(i) + jetAK8CHSchargedMultiplicityHandle->at(i); 

		float TJ = 0.0;
		if ((NHF<0.9) and (NEMF<0.9) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) TJ = 1.0;

		float LJ = 0.0;
		if ((NHF<0.99) and (NEMF<0.99) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) LJ = 1.0;


		jetAK8CHSTight->push_back(TJ); 
		jetAK8CHSLoose->push_back(LJ);




		}

	}

  //Puppi
  for( size_t i=0; i<jetAK8PuppiPtHandle->size(); i++ ) 
	{
	float shift = 1.0;
	float JECcorr = 1.0;

        TLorentzVector v1;
	v1.SetPtEtaPhiE(jetAK8PuppiPtHandle->at(i),jetAK8PuppiEtaHandle->at(i),jetAK8PuppiPhiHandle->at(i),jetAK8PuppiEHandle->at(i));
	float calcmass =  v1.M();




	float uncorrpt = fmaxf(1.0,jetAK8PuppijecFactor0Handle->at(i)*jetAK8PuppiPtHandle->at(i));	
	float uncorrE = fmaxf(1.0,jetAK8PuppijecFactor0Handle->at(i)*jetAK8PuppiEHandle->at(i));	

   	if (reapplyjec_)
		{
		JECcorr = fmaxf(0.0,jetAK8PuppijecFactor0Handle->at(i)*JEC_Corr_Puppi(uncorrpt,jetAK8PuppiEtaHandle->at(i),ISDATA,jetAK8PuppijetAreaHandle->at(i),*RhoHandle.product(),era_, runnum));
		}

	//float prunedjes = 1.0;
	//float prunedjer = 1.0;


	//float prunedjerunc = 1.0;
	
	float prunedshift = 1.0;


   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert_Puppi(jetAK8PuppiPtHandle->at(i),jetAK8PuppiEtaHandle->at(i),jes_,era_,runnum));

		float sign = (shift-1.0)/std::fabs(shift-1.0);
		float prunedjesunc = 1.+sign*(std::sqrt( (1.-shift)*(1.-shift) + (0.002)*(0.002)));
		

		prunedshift=prunedshift*prunedjesunc;		

		}


        
	if (not ISDATA) 
		{

			float SFapp=1.0;
			float RESapp=1.0;
			if (reapplyjer_)
				{	
					//JME::JetResolutionScaleFactor res_sf;
					//JME::JetResolution reso;
					//std::string JERFile_ = era_+"_MC_SF_AK8PFPuppi.txt";
					



/////////////////////////////////////////////////////TOUPDATE!	
					//Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt
					//Spring16_25nsV10_MC_SF_AK8PFPuppi.txt			
					std::string JERFile_ = "Spring16_25nsV10_MC_SF_AK8PFPuppi.txt";
					std::string RESFile_ = "Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt";
/////////////////////////////////////////////////////TOUPDATE!	


					JME::JetResolutionScaleFactor res_sf = JME::JetResolutionScaleFactor(JERFile_);

					JME::JetResolution reso = JME::JetResolution(RESFile_);

				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(jetAK8PuppiPtHandle->at(i)).setJetEta(jetAK8PuppiEtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") SFapp = res_sf.getScaleFactor(jetParam);
					else if (jer_ == "up") SFapp = res_sf.getScaleFactor(jetParam, Variation::UP);
 					else if (jer_ == "down") SFapp = res_sf.getScaleFactor(jetParam, Variation::DOWN);

					RESapp = reso.getResolution(jetParam)*jetAK8PuppiPtHandle->at(i);

				}
			else
				{

					if (jer_ == "nominal") SFapp = jetAK8PuppiJERSFHandle->at(i);
					else if (jer_ == "up") SFapp = jetAK8PuppiJERSFUpHandle->at(i);
 					else if (jer_ == "down") SFapp = jetAK8PuppiJERSFDownHandle->at(i);
					RESapp = jetAK8PuppiRESHandle->at(i)*jetAK8PuppiPtHandle->at(i);

				}	

			bool DRmatch =  deltaR(jetAK8PuppiEtaHandle->at(i), jetAK8PuppiPhiHandle->at(i), jetAK8PuppiGenJetEtaHandle->at(i), jetAK8PuppiGenJetPhiHandle->at(i))<0.4/2.0 ;
			bool ptmatch = std::fabs(jetAK8PuppiPtHandle->at(i)-jetAK8PuppiGenJetPtHandle->at(i))<3.0*RESapp;
			float ptJERCor = jetAK8PuppiGenJetPtHandle->at(i)+SFapp*(jetAK8PuppiPtHandle->at(i)-jetAK8PuppiGenJetPtHandle->at(i));

			if (DRmatch and ptmatch)
				{
					shift = shift*fmaxf(0.0,ptJERCor/jetAK8PuppiGenJetPtHandle->at(i));
					prunedshift = prunedshift*fmaxf(0.0,ptJERCor/jetAK8PuppiGenJetPtHandle->at(i));
				}
			else
				{
					//std::cout<<"Hybrid: "<<std::endl;
					//std::cout<<"DeltaR match: "<<deltaR(jetAK8EtaHandle->at(i), jetAK8PhiHandle->at(i), jetAK8GenJetEtaHandle->at(i), jetAK8GenJetPhiHandle->at(i)) <<std::endl;
					//std::cout<<"DeltaPt match: "<<std::fabs(jetAK8PtHandle->at(i)-jetAK8GenJetPtHandle->at(i)) <<std::endl;
					//std::cout<<"RES: "<<RESapp<<std::endl;
					//std::cout<<"RES threesig: "<<3.0*RESapp<<std::endl;
					//std::cout<<"SF: "<<SFapp<<std::endl;
					Double_t sigma =std::sqrt( SFapp*SFapp-1 ) * RESapp;// √(SF^2-1) * sigma_MC_PT.;
					TRandom3 *r = new TRandom3(0); 
					Double_t smearify = r->Gaus(0.0,sigma); 

					shift = shift*fmaxf(0.0,(1.0+smearify/jetAK8PuppiPtHandle->at(i)));
					prunedshift = prunedshift*fmaxf(0.0,(1.0+smearify/jetAK8PuppiPtHandle->at(i)));
					
				}


		}		
	
	jetAK8PuppiPt->push_back(jetAK8PuppiPtHandle->at(i)*shift*JECcorr);      
	jetAK8PuppiPhi->push_back(jetAK8PuppiPhiHandle->at(i));       
	jetAK8PuppiEta->push_back(jetAK8PuppiEtaHandle->at(i));       


	jetAK8PuppiMass->push_back(calcmass*shift*JECcorr);  



	float corrsdmass = Mass_Corr_Puppi(uncorrpt,jetAK8PuppiEtaHandle->at(i),uncorrE,ISDATA,jetAK8PuppijetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product(),era_,runnum);

	jetAK8PuppisoftDropMass->push_back(corrsdmass*jetAK8PuppisoftDropMassHandle->at(i)*shift);  
	jetAK8PuppisoftDropMassuncorr->push_back(jetAK8PuppisoftDropMassHandle->at(i)*shift);  
	jetAK8PuppiprunedMass->push_back(corrsdmass*0.99*jetAK8PuppiprunedMassHandle->at(i)*prunedshift);   
	jetAK8PuppiprunedMassuncorr->push_back(jetAK8PuppiprunedMassHandle->at(i)*prunedshift);   


 
   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8PuppiCSV->push_back(jetAK8PuppiCSVHandle->at(i));
		jetAK8PuppiCMVAv2->push_back(jetAK8PuppiCMVAv2Handle->at(i));              
		jetAK8PuppiPartonFlavour->push_back(jetAK8PuppiPartonFlavourHandle->at(i)); 
		  

		//TLV.SetPtEtaPhiM(ptcorr1,jetAK8EtaHandle->at(i),jetAK8PhiHandle->at(i),jetAK8MassHandle->at(i));
		


		jetAK8PuppifilteredMass->push_back(jetAK8PuppifilteredMassHandle->at(i));   
		//jetAK8PuppitopMass->push_back(jetAK8PuppitopMassHandle->at(i));   
		jetAK8PuppitrimmedMass->push_back(jetAK8PuppitrimmedMassHandle->at(i));    
		jetAK8PuppijecFactor0->push_back(jetAK8PuppijecFactor0Handle->at(i));   
		//jetAK8Puppiminmass->push_back(jetAK8PuppiminmassHandle->at(i));   
		//jetAK8PuppinSubJets->push_back(jetAK8PuppinSubJetsHandle->at(i));   
		jetAK8Puppitau1->push_back(jetAK8Puppitau1Handle->at(i));      
		jetAK8Puppitau2->push_back(jetAK8Puppitau2Handle->at(i));      
		jetAK8Puppitau3->push_back(jetAK8Puppitau3Handle->at(i));    

		jetAK8PuppivSubjetIndex0->push_back(jetAK8PuppivSubjetIndex0Handle->at(i));   
		jetAK8PuppivSubjetIndex1->push_back(jetAK8PuppivSubjetIndex1Handle->at(i)); 
  
		//jetAK8PuppitopSubjetIndex0->push_back(jetAK8PuppitopSubjetIndex0Handle->at(i));   
		//jetAK8PuppitopSubjetIndex1->push_back(jetAK8PuppitopSubjetIndex1Handle->at(i));   
		//jetAK8PuppitopSubjetIndex2->push_back(jetAK8PuppitopSubjetIndex2Handle->at(i));   
		//jetAK8PuppitopSubjetIndex3->push_back(jetAK8PuppitopSubjetIndex3Handle->at(i)); 







	  	float CHF = jetAK8PuppichargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK8PuppineutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK8PuppineutralHadronEnergyFracHandle->at(i);
	  	//float NC = jetAK8PuppiNumConstituentsHandle->at(i);
	  	float CM = jetAK8PuppichargedMultiplicityHandle->at(i);  
	  	float CEMF = jetAK8PuppichargedEmEnergyFracHandle->at(i);

		float NC = jetAK8PuppineutralMultiplicityHandle->at(i) + jetAK8PuppichargedMultiplicityHandle->at(i); 

		float TJ = 0.0;
		if ((NHF<0.9) and (NEMF<0.9) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) TJ = 1.0;

		float LJ = 0.0;
		if ((NHF<0.99) and (NEMF<0.99) and (NC>1) and (CHF>0.) and (CM>0) and (CEMF<0.99)) LJ = 1.0;


		jetAK8PuppiTight->push_back(TJ); 
		jetAK8PuppiLoose->push_back(LJ);




		}

	}




  iEvent.put(jetAK8CHSPt,"jetAK8CHSPt");
  iEvent.put(jetAK8CHSPhi,"jetAK8CHSPhi"); 
  iEvent.put(jetAK8CHSEta,"jetAK8CHSEta");  
  iEvent.put(jetAK8CHSMass,"jetAK8CHSMass");


  iEvent.put(jetAK8CHSprunedMass,"jetAK8CHSprunedMass");
  iEvent.put(jetAK8CHSsoftDropMass,"jetAK8CHSsoftDropMass");
  iEvent.put(jetAK8CHSsoftDropMassuncorr,"jetAK8CHSsoftDropMassuncorr");
  iEvent.put(jetAK8CHSprunedMassuncorr,"jetAK8CHSprunedMassuncorr");


  iEvent.put(jetAK8PuppiPt,"jetAK8PuppiPt");
  iEvent.put(jetAK8PuppiPhi,"jetAK8PuppiPhi"); 
  iEvent.put(jetAK8PuppiEta,"jetAK8PuppiEta");  
  iEvent.put(jetAK8PuppiMass,"jetAK8PuppiMass");


  iEvent.put(jetAK8PuppiprunedMass,"jetAK8PuppiprunedMass");
  iEvent.put(jetAK8PuppisoftDropMass,"jetAK8PuppisoftDropMass");
  iEvent.put(jetAK8PuppisoftDropMassuncorr,"jetAK8PuppisoftDropMassuncorr");
  iEvent.put(jetAK8PuppiprunedMassuncorr,"jetAK8PuppiprunedMassuncorr");




   if (jes_=="nominal"&&jer_=="nominal")
	{

  	iEvent.put(jetAK8CHSCSV,"jetAK8CHSCSV"); 
  	iEvent.put(jetAK8CHSCMVAv2,"jetAK8CHSCMVAv2"); 
  	iEvent.put(jetAK8CHSPartonFlavour,"jetAK8CHSPartonFlavour"); 
  	iEvent.put(jetAK8CHSfilteredMass,"jetAK8CHSfilteredMass"); 

  	iEvent.put(jetAK8CHSTight,"jetAK8CHSTight"); 
  	iEvent.put(jetAK8CHSLoose,"jetAK8CHSLoose");


  	//iEvent.put(jetAK8CHStopMass,"jetAK8CHStopMass");
  	iEvent.put(jetAK8CHStrimmedMass,"jetAK8CHStrimmedMass");
  	iEvent.put(jetAK8CHSjecFactor0,"jetAK8CHSjecFactor0");
  	///iEvent.put(jetAK8CHSminmass,"jetAK8CHSminmass"); 
  	//iEvent.put(jetAK8CHSnSubJets,"jetAK8CHSnSubJets"); 
  	iEvent.put(jetAK8CHStau1,"jetAK8CHStau1");
  	iEvent.put(jetAK8CHStau2,"jetAK8CHStau2");
  	iEvent.put(jetAK8CHStau3,"jetAK8CHStau3");

  	iEvent.put(jetAK8CHSvSubjetIndex0,"jetAK8CHSvSubjetIndex0");
  	iEvent.put(jetAK8CHSvSubjetIndex1,"jetAK8CHSvSubjetIndex1");

  	//iEvent.put(jetAK8CHStopSubjetIndex0,"jetAK8CHStopSubjetIndex0");
  	//iEvent.put(jetAK8CHStopSubjetIndex1,"jetAK8CHStopSubjetIndex1");
  	//iEvent.put(jetAK8CHStopSubjetIndex2,"jetAK8CHStopSubjetIndex2");
  	//iEvent.put(jetAK8CHStopSubjetIndex3,"jetAK8CHStopSubjetIndex3");

  	iEvent.put(jetAK8PuppiCSV,"jetAK8PuppiCSV"); 
  	iEvent.put(jetAK8PuppiCMVAv2,"jetAK8PuppiCMVAv2"); 
  	iEvent.put(jetAK8PuppiPartonFlavour,"jetAK8PuppiPartonFlavour"); 
  	iEvent.put(jetAK8PuppifilteredMass,"jetAK8PuppifilteredMass"); 

  	iEvent.put(jetAK8PuppiTight,"jetAK8PuppiTight"); 
  	iEvent.put(jetAK8PuppiLoose,"jetAK8PuppiLoose");


  	//iEvent.put(jetAK8PuppitopMass,"jetAK8PuppitopMass");
  	iEvent.put(jetAK8PuppitrimmedMass,"jetAK8PuppitrimmedMass");
  	iEvent.put(jetAK8PuppijecFactor0,"jetAK8PuppijecFactor0");
  	///iEvent.put(jetAK8Puppiminmass,"jetAK8Puppiminmass"); 
  	//iEvent.put(jetAK8PuppinSubJets,"jetAK8PuppinSubJets"); 
  	iEvent.put(jetAK8Puppitau1,"jetAK8Puppitau1");
  	iEvent.put(jetAK8Puppitau2,"jetAK8Puppitau2");
  	iEvent.put(jetAK8Puppitau3,"jetAK8Puppitau3");

  	iEvent.put(jetAK8PuppivSubjetIndex0,"jetAK8PuppivSubjetIndex0");
  	iEvent.put(jetAK8PuppivSubjetIndex1,"jetAK8PuppivSubjetIndex1");

  	//iEvent.put(jetAK8PuppitopSubjetIndex0,"jetAK8PuppitopSubjetIndex0");
  	//iEvent.put(jetAK8PuppitopSubjetIndex1,"jetAK8PuppitopSubjetIndex1");
  	//iEvent.put(jetAK8PuppitopSubjetIndex2,"jetAK8PuppitopSubjetIndex2");
  	//iEvent.put(jetAK8PuppitopSubjetIndex3,"jetAK8PuppitopSubjetIndex3");

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
