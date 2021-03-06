
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
//#include <TLorentzVector.h>




float getPUPPIweight(float puppipt, float puppieta ){

 TFile* file = TFile::Open( "puppiCorr.root","READ");

 TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
 TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
 TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");


  float genCorr  = 1.;
  float recoCorr = 1.;
  float totalWeight = 1.;
  genCorr =  puppisd_corrGEN->Eval( puppipt );
  if( fabs(puppieta)  <= 1.3 ){
    recoCorr = puppisd_corrRECO_cen->Eval( puppipt );
  }
  else{
    recoCorr = puppisd_corrRECO_for->Eval( puppipt );
  }

  totalWeight = genCorr * recoCorr;
  delete file;
  delete puppisd_corrGEN;
  delete puppisd_corrRECO_cen;
  delete puppisd_corrRECO_for;
  return totalWeight;

}



float Mass_Corr(boost::shared_ptr<FactorizedJetCorrector> jecAK8_,float  pt,float eta,float energy  , bool isdata,float Area,const double Rho,const int NPV)
	{

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
/*
   produces<std::vector<float>>("jetAK8CHSPt");
   produces<std::vector<float>>("jetAK8CHSPhi"); 
   produces<std::vector<float>>("jetAK8CHSEta");  
   produces<std::vector<float>>("jetAK8CHSMass");
*/
  // produces<std::vector<float>>("jetAK8CHSPtPuppi");
  // produces<std::vector<float>>("jetAK8CHSEtaPuppi");  

   produces<std::vector<float>>("jetAK8PuppiPt");
   produces<std::vector<float>>("jetAK8PuppiPhi"); 
   produces<std::vector<float>>("jetAK8PuppiEta");  
   produces<std::vector<float>>("jetAK8PuppiMass");

   //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSSmearedPt")));
   edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiSmearedPt")));



  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8CHS", "subjetAK8CHSjecFactor0")));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8CHS", "subjetAK8CHSPt")));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8CHS", "subjetAK8CHSPhi")));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8CHS", "subjetAK8CHSEta")));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8CHS", "subjetAK8CHSE")));



/*
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSPtResolution")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS","jetAK8CHSJERSFDown")));  


  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPt"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPhi"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSEta"    )));
//  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSMass"    )));


 // edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSPtPuppi"  )));
 // edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSEtaPuppi"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSE"    )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjecFactor0"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSjetArea"    )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetPt" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetEta" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetPhi" )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSGenJetE" )));
*/
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiPtResolution")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSF")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSFUp")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi","jetAK8PuppiJERSFDown")));  



  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppijecFactor0")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiPt")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiPhi")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiEta")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiE")));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiGenJetPt")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiGenJetPhi")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiGenJetEta")));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("subjetsAK8Puppi", "subjetAK8PuppiGenJetE")));



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
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiGenJetE" )));

  edm::EDGetTokenT<int>(consumes<int>(edm::InputTag("vertexInfo", "npv"    )));
  edm::EDGetTokenT<double>(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll", ""    )));

/*
  produces<std::vector<float>>("jetAK8CHSprunedMass");
  produces<std::vector<float>>("jetAK8CHSsoftDropMass");
 // produces<std::vector<float>>("jetAK8CHSsoftDropMassPuppi");

  produces<std::vector<float>>("jetAK8CHSsoftDropMassuncorr");
  produces<std::vector<float>>("jetAK8CHSprunedMassuncorr");
*/
  produces<std::vector<float>>("jetAK8PuppiprunedMass");
  produces<std::vector<float>>("jetAK8PuppisoftDropMass");
  produces<std::vector<float>>("jetAK8PuppiCorrectedsoftDropMass");
  //produces<std::vector<float>>("jetAK8CHSCorrectedsoftDropMassPuppi");

  produces<std::vector<float>>("jetAK8PuppiCorrectedsoftDropMassUnsmear");
 // produces<std::vector<float>>("jetAK8CHSCorrectedsoftDropMassUnsmearPuppi");

  produces<std::vector<float>>("jetAK8PuppisoftDropMassuncorr");
  produces<std::vector<float>>("jetAK8PuppisoftDropMassForTopPUPPIAK8JEC");
  produces<std::vector<float>>("jetAK8PuppiprunedMassuncorr");

  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSprunedMassCHS"  )));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS"  )));
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSsoftDropMassPuppi"  )));

  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiprunedMass"  )));
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppisoftDropMass"  )));

  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"  ))); 
  //edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetIndex1" )));

 // edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetPuppiIndex0"  ))); 
 // edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSvSubjetPuppiIndex1" )));



  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppivSubjetIndex0"  ))); 
  edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppivSubjetIndex1" )));


 

  produces<std::vector<float>>("jetAK8PuppivSubjetIndex0");
  produces<std::vector<float>>("jetAK8PuppivSubjetIndex1");




  if (jes_=="nominal"&&jer_=="nominal")
	{
/*
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

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSneutralMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8CHS", "jetAK8CHSchargedEmEnergyFrac" )));

*/



  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiCSVv2"     )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiCMVAv2"     )));  

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiDoubleBAK8"  ))); 
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

  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralEmEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralHadronEnergyFrac" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppiNumConstituents" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppineutralMultiplicity" )));
  	   	edm::EDGetTokenT<std::vector<float>>(consumes<std::vector<float>>(edm::InputTag("jetsAK8Puppi", "jetAK8PuppichargedEmEnergyFrac" )));





   		//produces<std::vector<float>>("jetAK8CHSTight"); 
   		//produces<std::vector<float>>("jetAK8CHSLoose"); 

  		produces<std::vector<float>>("jetAK8PuppiTight"); 
   		produces<std::vector<float>>("jetAK8PuppiLoose"); 
 

/*
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

   		//produces<std::vector<float>>("jetAK8CHSvSubjetPuppiIndex0");
   		//produces<std::vector<float>>("jetAK8CHSvSubjetPuppiIndex1");

   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex0");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex1");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex2");
   		//produces<std::vector<float>>("jetAK8CHStopSubjetIndex3");

*/
   		produces<std::vector<float>>("jetAK8PuppiCSV"); 
   		produces<std::vector<float>>("jetAK8PuppiCMVAv2"); 

   		produces<std::vector<float>>("jetAK8PuppiDoubleBAK8"); 

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



  		produces<std::vector<float>>("jetAK8PuppisoftDropMassForTopRAW" );
  		produces<std::vector<float>>("jetAK8PuppisoftDropMassForTopuncorr" );

  		produces<std::vector<float>>("sM0" );
  		produces<std::vector<float>>("sM1" );
  		produces<std::vector<float>>("sSDM0" );


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


/*
  std::auto_ptr<std::vector<float>> jetAK8CHSPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSPtPuppi(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8CHSSmearedPt(new std::vector<float>());         
        
  std::auto_ptr<std::vector<float>> jetAK8CHSPhi(new std::vector<float>());          
  std::auto_ptr<std::vector<float>> jetAK8CHSEtaPuppi(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSEta(new std::vector<float>());          
    
  std::auto_ptr<std::vector<float>> jetAK8CHSMass(new std::vector<float>()); 
      

  std::auto_ptr<std::vector<float>> jetAK8CHSTight(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSLoose(new std::vector<float>()); 


  std::auto_ptr<std::vector<float>> jetAK8CHSjetArea(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSE(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSjecFactor0(new std::vector<float>()); 
*/
  std::auto_ptr<int> npv(new int()); 
  std::auto_ptr<double> Rho(new double()); 

/*
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetPt(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetEta(new std::vector<float>()); 
  std::auto_ptr<std::vector<float>> jetAK8CHSGenJetPhi(new std::vector<float>()); 

  std::auto_ptr<std::vector<float>> jetAK8CHSCSV(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8CHSCMVAv2(new std::vector<float>());              
  std::auto_ptr<std::vector<float>> jetAK8CHSPartonFlavour(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSfilteredMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSprunedMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8CHSsoftDropMass(new std::vector<float>());     
  //std::auto_ptr<std::vector<float>> jetAK8CHSsoftDropMassPuppi(new std::vector<float>());     
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
       
 // std::auto_ptr<std::vector<float>> jetAK8CHSvSubjetPuppiIndex0(new std::vector<float>());      
 // std::auto_ptr<std::vector<float>> jetAK8CHSvSubjetPuppiIndex1(new std::vector<float>());      
       

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
*/
  edm::Handle<int> npvHandle;
  edm::Handle<double> RhoHandle;
/*
  edm::Handle<std::vector<float>> jetAK8CHSPtHandle; 
  edm::Handle<std::vector<float>> jetAK8CHSPtPuppiHandle; 

  edm::Handle<std::vector<float>> jetAK8CHSSmearedPtHandle;             
  edm::Handle<std::vector<float>> jetAK8CHSPhiHandle;        
  edm::Handle<std::vector<float>> jetAK8CHSEtaHandle;     
  edm::Handle<std::vector<float>> jetAK8CHSEtaPuppiHandle;        
   
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
 // edm::Handle<std::vector<float>> jetAK8CHSsoftDropMassPuppiHandle;   
 // edm::Handle<std::vector<float>> jetAK8CHStopMassHandle;    
  edm::Handle<std::vector<float>> jetAK8CHStrimmedMassHandle;     

  //edm::Handle<std::vector<float>> jetAK8CHSminmassHandle;    
 // edm::Handle<std::vector<float>> jetAK8CHSnSubJetsHandle;    
  edm::Handle<std::vector<float>> jetAK8CHStau1Handle;       
  edm::Handle<std::vector<float>> jetAK8CHStau2Handle;       
  edm::Handle<std::vector<float>> jetAK8CHStau3Handle;  


  edm::Handle<std::vector<float>> jetAK8CHSvSubjetIndex0Handle;    
  edm::Handle<std::vector<float>> jetAK8CHSvSubjetIndex1Handle;  

  //edm::Handle<std::vector<float>> jetAK8CHSvSubjetPuppiIndex0Handle;    
//  edm::Handle<std::vector<float>> jetAK8CHSvSubjetPuppiIndex1Handle;  
    
    
  //edm::Handle<std::vector<float>> subjetAK8CHSjecFactor0Handle;    
  //edm::Handle<std::vector<float>> subjetAK8CHSPtHandle;  
  //edm::Handle<std::vector<float>> subjetAK8CHSPhiHandle;   
  //edm::Handle<std::vector<float>> subjetAK8CHSEtaHandle;         
  //edm::Handle<std::vector<float>> subjetAK8CHSEHandle;         

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



*/

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

  std::auto_ptr<std::vector<float>> jetAK8PuppiDoubleBAK8(new std::vector<float>());
              
  std::auto_ptr<std::vector<float>> jetAK8PuppiPartonFlavour(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppifilteredMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppiprunedMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMass(new std::vector<float>());  
  std::auto_ptr<std::vector<float>> jetAK8PuppiCorrectedsoftDropMass(new std::vector<float>()); 
  //std::auto_ptr<std::vector<float>> jetAK8CHSCorrectedsoftDropMassPuppi(new std::vector<float>());  

  std::auto_ptr<std::vector<float>> jetAK8PuppiCorrectedsoftDropMassUnsmear(new std::vector<float>()); 
 // std::auto_ptr<std::vector<float>> jetAK8CHSCorrectedsoftDropMassUnsmearPuppi(new std::vector<float>()); 
 
  //std::auto_ptr<std::vector<float>> jetAK8PuppitopMass(new std::vector<float>());      
  std::auto_ptr<std::vector<float>> jetAK8PuppitrimmedMass(new std::vector<float>());       
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMassuncorr(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMassForTopPUPPIAK8JEC(new std::vector<float>());
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



  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMassForTopRAW(new std::vector<float>());
  std::auto_ptr<std::vector<float>> jetAK8PuppisoftDropMassForTopuncorr(new std::vector<float>());


  std::auto_ptr<std::vector<float>> sM0(new std::vector<float>());
  std::auto_ptr<std::vector<float>> sM1(new std::vector<float>());
  std::auto_ptr<std::vector<float>> sSDM0(new std::vector<float>());


       
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

  edm::Handle<std::vector<float>> jetAK8PuppiDoubleBAK8Handle;
           
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



  edm::Handle<std::vector<float>>  subjetAK8PuppijecFactor0Handle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiPtHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiPhiHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiEtaHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiEHandle;

  edm::Handle<std::vector<float>>  subjetAK8PuppiGenJetPtHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiGenJetPhiHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiGenJetEtaHandle;
  edm::Handle<std::vector<float>>  subjetAK8PuppiGenJetEHandle;







 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex0Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex1Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex2Handle;    
 // edm::Handle<std::vector<float>> jetAK8PuppitopSubjetIndex3Handle;  


/*
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPt"        ,jetAK8CHSPtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSSmearedPt"        ,jetAK8CHSSmearedPtHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPhi"       ,jetAK8CHSPhiHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSEta"       ,jetAK8CHSEtaHandle);  
 // iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSMass"      ,jetAK8MassHandle);

 // iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSPtPuppi"        ,jetAK8CHSPtPuppiHandle);
 // iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSEtaPuppi"       ,jetAK8CHSEtaPuppiHandle);  


  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSE"        ,jetAK8CHSEHandle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"        ,jetAK8CHSjecFactor0Handle);
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjetArea"        ,jetAK8CHSjetAreaHandle);

  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPt"        ,jetAK8CHSGenJetPtHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetEta"        ,jetAK8CHSGenJetEtaHandle);
  if (not ISDATA) iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSGenJetPhi"        ,jetAK8CHSGenJetPhiHandle);

*/
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

/*
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSPtResolution",  jetAK8CHSRESHandle);

  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSF",  jetAK8CHSJERSFHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFUp",  jetAK8CHSJERSFUpHandle);
  iEvent.getByLabel( "jetsAK8CHS", "jetAK8CHSJERSFDown",  jetAK8CHSJERSFDownHandle);


  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSprunedMassCHS"   ,jetAK8CHSprunedMassHandle);  
  iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSsoftDropMassCHS"   ,jetAK8CHSsoftDropMassHandle); 
  //iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSsoftDropMassPuppi"   ,jetAK8CHSsoftDropMassPuppiHandle); 

*/
  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiPtResolution",  jetAK8PuppiRESHandle);

  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSF",  jetAK8PuppiJERSFHandle);
  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSFUp",  jetAK8PuppiJERSFUpHandle);
  iEvent.getByLabel( "jetsAK8Puppi", "jetAK8PuppiJERSFDown",  jetAK8PuppiJERSFDownHandle);


  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiprunedMass"   ,jetAK8PuppiprunedMassHandle);  
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppisoftDropMass"   ,jetAK8PuppisoftDropMassHandle); 

  //iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex0"   ,jetAK8CHSvSubjetIndex0Handle);  
  //iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetIndex1"   ,jetAK8CHSvSubjetIndex1Handle);


  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppivSubjetIndex0"   ,jetAK8PuppivSubjetIndex0Handle);  
  iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppivSubjetIndex1"   ,jetAK8PuppivSubjetIndex1Handle);




  //iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetPuppiIndex0"   ,jetAK8CHSvSubjetPuppiIndex0Handle);  
 // iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSvSubjetPuppiIndex1"   ,jetAK8CHSvSubjetPuppiIndex1Handle);


  iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppijecFactor0",subjetAK8PuppijecFactor0Handle);
  iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiPt",subjetAK8PuppiPtHandle);
  iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiPhi",subjetAK8PuppiPhiHandle);
  iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiEta",subjetAK8PuppiEtaHandle);
  iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiE",subjetAK8PuppiEHandle);

  if (not ISDATA) iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiGenJetPt",subjetAK8PuppiGenJetPtHandle);
  if (not ISDATA) iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiGenJetPhi",subjetAK8PuppiGenJetPhiHandle);
  if (not ISDATA) iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiGenJetEta",subjetAK8PuppiGenJetEtaHandle);
  if (not ISDATA) iEvent.getByLabel("subjetsAK8Puppi", "subjetAK8PuppiGenJetE",subjetAK8PuppiGenJetEHandle);


  if (jes_=="nominal"&&jer_=="nominal")
	{
/*
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


  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStrimmedMassCHS"   ,jetAK8CHStrimmedMassHandle);   
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHSjecFactor0"   ,jetAK8CHSjecFactor0Handle);  
 
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau1CHS"      ,jetAK8CHStau1Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau2CHS"      ,jetAK8CHStau2Handle);  
  	iEvent.getByLabel("jetsAK8CHS", "jetAK8CHStau3CHS"      ,jetAK8CHStau3Handle); 

*/
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedHadronEnergyFrac"       ,jetAK8PuppichargedHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralEmEnergyFrac"       ,jetAK8PuppineutralEmEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralHadronEnergyFrac"       ,jetAK8PuppineutralHadronEnergyFracHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiNumConstituents"       ,jetAK8PuppiNumConstituentsHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedMultiplicity"       ,jetAK8PuppichargedMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppineutralMultiplicity"       ,jetAK8PuppineutralMultiplicityHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppichargedEmEnergyFrac"       ,jetAK8PuppichargedEmEnergyFracHandle);  








  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiCSVv2"       ,jetAK8PuppiCSVHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiDoubleBAK8"       ,jetAK8PuppiDoubleBAK8Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiCMVAv2"       ,jetAK8PuppiCMVAv2Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppiPartonFlavour"   ,jetAK8PuppiPartonFlavourHandle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppifilteredMass"   ,jetAK8PuppifilteredMassHandle);  

  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppitrimmedMass"   ,jetAK8PuppitrimmedMassHandle);   
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8PuppijecFactor0"   ,jetAK8PuppijecFactor0Handle);  
  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau1"      ,jetAK8Puppitau1Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau2"      ,jetAK8Puppitau2Handle);  
  	iEvent.getByLabel("jetsAK8Puppi", "jetAK8Puppitau3"      ,jetAK8Puppitau3Handle); 


	}

  //float JMR_PRUNED = 1.07;
  //float JMS_PRUNED = 1.0;

  //float SIGMA_JMR_PRUNED_EX = 0.103;
  //float SIGMA_JMS_PRUNED_EX = 0.02;

  float JMR_PUPPISD = 1.0;
  float JMS_PUPPISD = 1.0;

  float SIGMA_JMR_PUPPISD = 0.2;
  float SIGMA_JMS_PUPPISD = 0.0094;

////TO UPDATE!

JME::JetResolutionScaleFactor res_sfC;
JME::JetResolution resoC;
JME::JetResolutionScaleFactor res_sfP;
JME::JetResolution resoP;


std::string JERFileC_ = "jecfiles/Spring16_25nsV10_MC_SF_AK8PFchs.txt";
std::string RESFileC_ = "jecfiles/Spring16_25nsV10_MC_PtResolution_AK8PFchs.txt";
res_sfC = JME::JetResolutionScaleFactor(JERFileC_);
resoC = JME::JetResolution(RESFileC_);
//resoC.dump();



std::string JERFileP_ = "jecfiles/Spring16_25nsV10_MC_SF_AK8PFPuppi.txt";
std::string RESFileP_ = "jecfiles/Spring16_25nsV10_MC_PtResolution_AK8PFPuppi.txt";
res_sfP = JME::JetResolutionScaleFactor(JERFileP_);
resoP = JME::JetResolution(RESFileP_);









boost::shared_ptr<FactorizedJetCorrector> jecAK8CHS_(new FactorizedJetCorrector);
std::vector<std::string> jecAK8CHSPayloadNames_;
std::string runtxt_;

if(ISDATA) 
		{
		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV4";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV4";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV4"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV4"; //IOV H:[280919,Infinity] f

  		jecAK8CHSPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2Relative_AK8PFchs.txt");
  		jecAK8CHSPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L3Absolute_AK8PFchs.txt");
  		jecAK8CHSPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2L3Residual_AK8PFchs.txt");
		}
else
		{	
		runtxt_ = "V4";
  		jecAK8CHSPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L2Relative_AK8PFchs.txt");
  		jecAK8CHSPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L3Absolute_AK8PFchs.txt");
		}

std::vector<JetCorrectorParameters> vParCHS;
for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8CHSPayloadNames_.begin(), payloadEnd = jecAK8CHSPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) 
		{

  	  	JetCorrectorParameters parsCHS(*ipayload);
 	   	vParCHS.push_back(parsCHS);
 		}

jecAK8CHS_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vParCHS) );








boost::shared_ptr<FactorizedJetCorrector> jecAK8PUPPI_;
std::vector<std::string> jecAK8PUPPIPayloadNames_;

if(ISDATA) 
		{

		if (runnum>=1 and runnum<=276811) runtxt_ = "BCDV4";  // IOV BCD:[1,276811]  (For Runs B/C/D)
		if (runnum>=276831 and runnum<=278801) runtxt_ = "EFV4";  //IOV EF:[276831,278801]  (For Runs E/early F)
		if (runnum>=278802 and runnum<=280385) runtxt_ = "GV4"; //IOV G:[278802,280385] (For Runs late F/G)
		if (runnum>=280919) runtxt_ = "HV4"; //IOV H:[280919,Infinity] f



  		jecAK8PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2Relative_AK8PFPuppi.txt");
  		jecAK8PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L3Absolute_AK8PFPuppi.txt");
  		jecAK8PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_DATA_L2L3Residual_AK8PFPuppi.txt");
		}
else
		{

		runtxt_ = "V4";
  		jecAK8PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L2Relative_AK8PFPuppi.txt");
  		jecAK8PUPPIPayloadNames_.push_back("jecfiles/"+era_+runtxt_+"_MC_L3Absolute_AK8PFPuppi.txt");
		}

std::vector<JetCorrectorParameters> vParPUPPI;

for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PUPPIPayloadNames_.begin(), payloadEnd = jecAK8PUPPIPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) 
		{
  	  	JetCorrectorParameters pars(*ipayload);
 	   	vParPUPPI.push_back(pars);
 		}
jecAK8PUPPI_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vParPUPPI) );




JetCorrectionUncertainty *jecUncPUPPI = new JetCorrectionUncertainty("jecfiles/"+era_+"V4_MC_Uncertainty_AK8PFPuppi.txt");
JetCorrectionUncertainty *jecUncCHS = new JetCorrectionUncertainty("jecfiles/"+era_+"V4_MC_Uncertainty_AK8PFchs.txt");

/*
std::cout<<"- "<<std::endl;
std::cout<<"- "<<std::endl;
std::cout<<"- "<<std::endl;
std::cout<<"newevent "<<std::endl;
std::cout<<"JES "<<jes_<<std::endl;
std::cout<<"JER "<<jer_<<std::endl;
*/



  /*CHS
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
		JECcorr = fmaxf(0.0,jetAK8CHSjecFactor0Handle->at(i)*JEC_Corr(jecAK8CHS_,uncorrpt,jetAK8CHSEtaHandle->at(i),ISDATA,jetAK8CHSjetAreaHandle->at(i),*RhoHandle.product()));
		}

	
	float chscorrpt = jetAK8CHSPtHandle->at(i)*JECcorr;

	float prunedshift = JMS_PRUNED;


   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert(jecUncCHS,chscorrpt,jetAK8CHSEtaHandle->at(i),jes_));

		float sign = (shift-1.0)/std::fabs(shift-1.0);
		float prunedjesunc = 1.+sign*(std::sqrt( (1.-shift)*(1.-shift) + (SIGMA_JMS_PRUNED_EX)*(SIGMA_JMS_PRUNED_EX)));


		prunedshift=prunedshift*prunedjesunc;		

		}
	

	if (not ISDATA) 
		{


			float SFapp=1.0;
			float RESapp=1.0;






			if (reapplyjer_)
				{	
					
				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(chscorrpt).setJetEta(jetAK8CHSEtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") 
						{
						SFapp = res_sfC.getScaleFactor(jetParam);

						}
					else if (jer_ == "up")
						{ 
						SFapp = res_sfC.getScaleFactor(jetParam, Variation::UP);	

						}
 					else if (jer_ == "down") 
						{
						SFapp = res_sfC.getScaleFactor(jetParam, Variation::DOWN);

						}

					RESapp = resoC.getResolution(jetParam)*chscorrpt;

				}
			else
				{




					if (jer_ == "nominal") 
						{
						SFapp = jetAK8CHSJERSFHandle->at(i);
						}
					else if (jer_ == "up")
						{
						SFapp = jetAK8CHSJERSFUpHandle->at(i);
						}
 					else if (jer_ == "down")
						{
						SFapp = jetAK8CHSJERSFDownHandle->at(i);
						}
					RESapp = jetAK8CHSRESHandle->at(i)*chscorrpt;




				}	

			bool DRmatch =  deltaR(jetAK8CHSEtaHandle->at(i), jetAK8CHSPhiHandle->at(i), jetAK8CHSGenJetEtaHandle->at(i), jetAK8CHSGenJetPhiHandle->at(i))<0.4 ;
			bool ptmatch = std::fabs(chscorrpt-jetAK8CHSGenJetPtHandle->at(i))<3.0*RESapp;


			float ptcorrunsmear = chscorrpt;
			float ptgen =jetAK8CHSGenJetPtHandle->at(i);

			float smearcorr = 1.0+(SFapp-1.0)*((ptcorrunsmear-ptgen)/ptcorrunsmear);


			
			if (DRmatch and ptmatch)
				{
					shift = shift*fmaxf(0.0,smearcorr);
					
					prunedshift = prunedshift*fmaxf(0.0,smearcorr);
				}
			else
				{
					//std::cout<<"Hybrid: "<<std::endl;
					//std::cout<<"DeltaR match: "<<deltaR(jetAK8EtaHandle->at(i), jetAK8PhiHandle->at(i), jetAK8GenJetEtaHandle->at(i), jetAK8GenJetPhiHandle->at(i)) <<std::endl;
					//std::cout<<"DeltaPt match: "<<std::fabs(jetAK8PtHandle->at(i)-jetAK8GenJetPtHandle->at(i)) <<std::endl;
					//std::cout<<"RES: "<<RESapp<<std::endl;
					//std::cout<<"RES threesig: "<<3.0*RESapp<<std::endl;
					//std::cout<<"SF: "<<SFapp<<std::endl;
					Double_t sigma =std::sqrt( SFapp*SFapp-1.0 ) * RESapp;// √(SF^2-1) * sigma_MC_PT.;
					TRandom3 *r = new TRandom3(0); 
					Double_t smearify = r->Gaus(0.0,sigma); 

					shift = shift*fmaxf(0.0,(1.0+smearify/chscorrpt));
					prunedshift = prunedshift*fmaxf(0.0,(1.0+smearify/chscorrpt));
				
					delete r;
				}




		}		
	

	jetAK8CHSPt->push_back(chscorrpt*shift);      
	jetAK8CHSPhi->push_back(jetAK8CHSPhiHandle->at(i));       
	jetAK8CHSEta->push_back(jetAK8CHSEtaHandle->at(i));       


	jetAK8CHSMass->push_back(calcmass*shift*JECcorr);  


	//Here we recalculate the SD mass and eliminate all corrections
	 

	float corrsdmass = Mass_Corr(jecAK8CHS_,uncorrpt,jetAK8CHSEtaHandle->at(i),uncorrE,ISDATA,jetAK8CHSjetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product());


	jetAK8CHSsoftDropMass->push_back(corrsdmass*jetAK8CHSsoftDropMassHandle->at(i)*shift); 

	jetAK8CHSsoftDropMassuncorr->push_back(jetAK8CHSsoftDropMassHandle->at(i)*shift);  
	jetAK8CHSprunedMass->push_back(corrsdmass*jetAK8CHSprunedMassHandle->at(i)*prunedshift);   
	jetAK8CHSprunedMassuncorr->push_back(jetAK8CHSprunedMassHandle->at(i)*prunedshift);   


 
   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8CHSCSV->push_back(jetAK8CHSCSVHandle->at(i));
		jetAK8CHSCMVAv2->push_back(jetAK8CHSCMVAv2Handle->at(i));              
		jetAK8CHSPartonFlavour->push_back(jetAK8CHSPartonFlavourHandle->at(i)); 
		  

		jetAK8CHSfilteredMass->push_back(jetAK8CHSfilteredMassHandle->at(i));   

		jetAK8CHStrimmedMass->push_back(jetAK8CHStrimmedMassHandle->at(i));    
		jetAK8CHSjecFactor0->push_back(jetAK8CHSjecFactor0Handle->at(i));   
	 
		jetAK8CHStau1->push_back(jetAK8CHStau1Handle->at(i));      
		jetAK8CHStau2->push_back(jetAK8CHStau2Handle->at(i));      
		jetAK8CHStau3->push_back(jetAK8CHStau3Handle->at(i));    

		jetAK8CHSvSubjetIndex0->push_back(jetAK8CHSvSubjetIndex0Handle->at(i));   
		jetAK8CHSvSubjetIndex1->push_back(jetAK8CHSvSubjetIndex1Handle->at(i)); 
 




	  	float CHF = jetAK8CHSchargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK8CHSneutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK8CHSneutralHadronEnergyFracHandle->at(i);
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
  */
  //Puppi
  for( size_t i=0; i<jetAK8PuppiPtHandle->size(); i++ ) 
	{

	float shift = 1.0;
	float JECcorr = 1.0;
	float SDshift = 1.0;
        TLorentzVector v1;
	v1.SetPtEtaPhiE(jetAK8PuppiPtHandle->at(i),jetAK8PuppiEtaHandle->at(i),jetAK8PuppiPhiHandle->at(i),jetAK8PuppiEHandle->at(i));
	float calcmass =  v1.M();




	float uncorrpt = fmaxf(1.0,jetAK8PuppijecFactor0Handle->at(i)*jetAK8PuppiPtHandle->at(i));	
	float uncorrE = fmaxf(1.0,jetAK8PuppijecFactor0Handle->at(i)*jetAK8PuppiEHandle->at(i));	

   	if (reapplyjec_)
		{
		JECcorr = fmaxf(0.0,jetAK8PuppijecFactor0Handle->at(i)*JEC_Corr(jecAK8PUPPI_,uncorrpt,jetAK8PuppiEtaHandle->at(i),ISDATA,jetAK8PuppijetAreaHandle->at(i),*RhoHandle.product()));
		}


	float puppicorrpt = jetAK8PuppiPtHandle->at(i)*JECcorr;
   	if (jes_!="nominal")
		{
		shift = fmaxf(0.0,JES_Uncert(jecUncPUPPI,puppicorrpt,jetAK8PuppiEtaHandle->at(i),jes_));


		//float SDjesunc = 1.0;
		if (jes_=="up") JMS_PUPPISD += SIGMA_JMS_PUPPISD;
		else if (jes_=="down") JMS_PUPPISD -= SIGMA_JMS_PUPPISD;
		SDshift=SDshift*JMS_PUPPISD;		
		//std::cout<<"sdjes "<<SDjesunc<<std::endl;


		}
	//std::cout<<"JECUPDATE "<<JECcorr*shift<<std::endl;


	float puppiEta = jetAK8PuppiEtaHandle->at(i);

	float puppiCorr = getPUPPIweight(puppicorrpt,puppiEta);

 
	float ptsmfacM0 = 0.0;
	float ptsmfacM1 = 0.0;
	float smfac = 1.0;
	if (not ISDATA) 
		{


			float SFSDMASS = 1.0;
			float SFapp=1.0;
			float RESapp=1.0;
			if (reapplyjer_)
				{	
					
				  	JME::JetParameters jetParam;
				    	jetParam.setJetPt(puppicorrpt).setJetEta(jetAK8PuppiEtaHandle->at(i)).setRho(*RhoHandle.product());

					if (jer_ == "nominal") 
						{
						SFapp = res_sfP.getScaleFactor(jetParam);
						SFSDMASS = JMR_PUPPISD;
						}
					else if (jer_ == "up")
						{ 
						SFapp = res_sfP.getScaleFactor(jetParam, Variation::UP);	
						SFSDMASS = JMR_PUPPISD+SIGMA_JMR_PUPPISD;
						}
 					else if (jer_ == "down") 
						{
						SFapp = res_sfP.getScaleFactor(jetParam, Variation::DOWN);
						SFSDMASS = JMR_PUPPISD-SIGMA_JMR_PUPPISD;
						}

				
					RESapp = resoP.getResolution(jetParam)*puppicorrpt;

				}
			else
				{

					if (jer_ == "nominal") 
						{
						SFapp = jetAK8PuppiJERSFHandle->at(i);
						SFSDMASS = JMR_PUPPISD;
						}
					else if (jer_ == "up")
						{
						SFapp = jetAK8PuppiJERSFUpHandle->at(i);
						SFSDMASS = JMR_PUPPISD+SIGMA_JMR_PUPPISD;
						}
 					else if (jer_ == "down")
						{
						SFapp = jetAK8PuppiJERSFDownHandle->at(i);
						SFSDMASS = JMR_PUPPISD-SIGMA_JMR_PUPPISD;
						}
					RESapp = jetAK8PuppiRESHandle->at(i)*puppicorrpt;

				}	


			float DRval =  deltaR(jetAK8PuppiEtaHandle->at(i), jetAK8PuppiPhiHandle->at(i), jetAK8PuppiGenJetEtaHandle->at(i), jetAK8PuppiGenJetPhiHandle->at(i)) ;
			float ptval = std::fabs(puppicorrpt-jetAK8PuppiGenJetPtHandle->at(i));


			bool DRmatch = DRval<0.4 ;
			bool ptmatch = ptval<3.0*RESapp;




			//std::cout<<"GENSDMASS "<<GENSDMASS<<std::endl;
			
			//std::cout<<"PT SMEAR FACTOR "<<SFapp<<std::endl;
			//std::cout<<"PRE "<<puppicorrpt<<std::endl;
			//float ptJERCor = jetAK8PuppiGenJetPtHandle->at(i)+SFapp*(puppicorrpt-jetAK8PuppiGenJetPtHandle->at(i));
			//float SDmassJERCor = GENSDMASS+SFSDMASS*(puppicorrsdmass-GENSDMASS);

			float ptcorrunsmear = puppicorrpt;
			float ptgen =jetAK8PuppiGenJetPtHandle->at(i);

			float smearcorr = 1.0+(SFapp-1.0)*((ptcorrunsmear-ptgen)/ptcorrunsmear);

			//Currently debugging
			

			Double_t massresfac =std::sqrt(fmaxf(0.0,(SFSDMASS*SFSDMASS-1.0)));// √(SF^2-1) * sigma_MC_PT.;
			TRandom3 *massr = new TRandom3(0); 
			float SFSDRESO =  10.1/80.4; 
			
			Double_t masssmearify = massr->Gaus(0.0,SFSDRESO)*massresfac; 
			//std::cout<<"massresfac "<<massresfac<<std::endl;
			//std::cout<<"masssmearify "<<masssmearify<<std::endl;
			smfac = smfac*fmaxf(0.0,(1.0+masssmearify));
			//std::cout<<"puppicorrsdmass "<<SDshift*jetAK8PuppisoftDropMassHandle->at(i)*puppiCorr<<std::endl;
			//std::cout<<"sdsmear "<<smfac<<std::endl;
			SDshift = SDshift*smfac;

			delete massr;


			if (DRmatch and ptmatch)
				{
					shift = shift*fmaxf(0.0,smearcorr);
					
						
				}
			else
				{


					//std::cout<<"Hybrid: "<<std::endl;
					//std::cout<<"DeltaR match: "<<DRval<<std::endl;
					//std::cout<<"DeltaPt match: "<<ptval <<std::endl;
					//std::cout<<"RES: "<<RESapp<<std::endl;
					//std::cout<<"RES threesig: "<<3.0*RESapp<<std::endl;
				        //std::cout<<"SF: "<<SFapp<<std::endl;
					//Double_t sigma =std::sqrt( SFapp*SFapp-1.0 ) * RESapp;// √(SF^2-1) * sigma_MC_PT.;
					//TRandom3 *r = new TRandom3(0); 
					//Double_t smearify = r->Gaus(0.0,sigma); 



					
					//delete r;




					Double_t ptresfac =std::sqrt(SFapp*SFapp-1.0 );// √(SF^2-1) * sigma_MC_PT.;
					TRandom3 *r1 = new TRandom3(0); 
					Double_t ptsmearify = r1->Gaus(0.0,RESapp/puppicorrpt)*ptresfac; 
					shift = shift*fmaxf(0.0,(1.0+ptsmearify));
					//std::cout<<"Mehtod0 "<<fmaxf(0.0,(1.0+ptsmearify))<<std::endl;
					//std::cout<<"Mehtod1 "<<fmaxf(0.0,(1.0+smearify/puppicorrpt))<<std::endl;
					ptsmfacM0 = fmaxf(0.0,(1.0+ptsmearify));

					delete r1;

				}

	/*
	std::cout<<"JER = "<<jer_<<std::endl;
	std::cout<<"JES = "<<jes_<<std::endl;

	std::cout<<"PUPPTRAW = "<<puppicorrpt<<std::endl;
	std::cout<<"SHIFT = "<<shift<<std::endl;
	std::cout<<"JECCORR = "<<JECcorr<<std::endl;
	std::cout<<"PTcor = "<<puppicorrpt*shift*JECcorr<<std::endl;

	std::cout<<"PUPSDMASSRAW = "<<jetAK8PuppisoftDropMassHandle->at(i)<<std::endl;
	std::cout<<"PUPSDsdcorr = "<<puppiCorr<<std::endl;
	std::cout<<"PUPSDsdcorrected = "<<puppicorrsdmass<<std::endl;
	std::cout<<"PUPSDsmear = "<<SDsmearcorr<<std::endl;
	std::cout<<"PUPSDUNsmear = "<<puppicorrsdmass*SDshift*(1.0/smfac)<<std::endl;
	//std::cout<<"PUPSDsmearPT = "<<SDsmearcorrPT<<std::endl;
	std::cout<<"SDshift = "<<SDshift<<std::endl;

	std::cout<<"Final = "<<puppicorrsdmass*SDshift<<std::endl<<std::endl;
	*/



		}		



	float SJ0index = jetAK8PuppivSubjetIndex0Handle->at(i);
	float SJ1index = jetAK8PuppivSubjetIndex1Handle->at(i);

	float SDFromCALC = -1.0;
	float SDFromCALCuncorr = -1.0;
	float SDFromCALCPUPPIAK8JEC  = -1.0;


	if (SJ0index>=0)
		{

		TLorentzVector vSJ;
		vSJ.SetPtEtaPhiE(subjetAK8PuppiPtHandle->at(SJ0index),subjetAK8PuppiEtaHandle->at(SJ0index),subjetAK8PuppiPhiHandle->at(SJ0index),subjetAK8PuppiEHandle->at(SJ0index));


		TLorentzVector vSJuncorr;
		vSJuncorr = vSJ*(subjetAK8PuppijecFactor0Handle->at(SJ0index));



		if (SJ1index>=0)
			{
				TLorentzVector vSJ1;
				vSJ1.SetPtEtaPhiE(subjetAK8PuppiPtHandle->at(SJ1index),subjetAK8PuppiEtaHandle->at(SJ1index),subjetAK8PuppiPhiHandle->at(SJ1index),subjetAK8PuppiEHandle->at(SJ1index));

				vSJuncorr += vSJ1*(subjetAK8PuppijecFactor0Handle->at(SJ1index));
				vSJ+=vSJ1;


			}
		//GENSDMASS = (vgen0+vgen1).M();

		SDFromCALC = vSJ.M();
		SDFromCALCuncorr = vSJuncorr.M();

		float JECcorrFull = JECcorr/jetAK8PuppijecFactor0Handle->at(i);
		
		//Here, use the L2L3 just as the fat jet.  Need to eliminate the subjet corrections first.  Currently only apply and vary JES
		SDFromCALCPUPPIAK8JEC = vSJuncorr.M()*JECcorrFull*shift;	


		}

	//std::cout<<"SDFromCALC "<<SDFromCALC<<std::endl;
	//std::cout<<"SDFromCALCuncorr "<<SDFromCALCuncorr<<std::endl;
	//std::cout<<"SDFromCALCPUPPIAK8JEC "<<SDFromCALCPUPPIAK8JEC<<std::endl;



	float puppicorrsdmass = SDshift*jetAK8PuppisoftDropMassHandle->at(i)*puppiCorr;
	jetAK8PuppiPt->push_back(puppicorrpt*shift);      
	jetAK8PuppiPhi->push_back(jetAK8PuppiPhiHandle->at(i));       
	jetAK8PuppiEta->push_back(jetAK8PuppiEtaHandle->at(i));       


	jetAK8PuppiMass->push_back(calcmass*shift*JECcorr);  



	float corrsdmass = Mass_Corr(jecAK8PUPPI_,uncorrpt,jetAK8PuppiEtaHandle->at(i),uncorrE,ISDATA,jetAK8PuppijetAreaHandle->at(i),*RhoHandle.product(),*npvHandle.product());




	jetAK8PuppisoftDropMass->push_back(corrsdmass*jetAK8PuppisoftDropMassHandle->at(i)*shift);

	//std::cout<<"sdshift "<<SDshift<<std::endl;

	jetAK8PuppisoftDropMassForTopPUPPIAK8JEC->push_back(SDFromCALCPUPPIAK8JEC);

	jetAK8PuppiCorrectedsoftDropMass->push_back(puppicorrsdmass);    
	jetAK8PuppiCorrectedsoftDropMassUnsmear->push_back(puppicorrsdmass*(1.0/smfac));    
	jetAK8PuppisoftDropMassuncorr->push_back(jetAK8PuppisoftDropMassHandle->at(i)*shift);  

   	if (jes_=="nominal"&&jer_=="nominal")
		{
		jetAK8PuppiCSV->push_back(jetAK8PuppiCSVHandle->at(i));
		jetAK8PuppiDoubleBAK8->push_back(jetAK8PuppiDoubleBAK8Handle->at(i));



		jetAK8PuppiCMVAv2->push_back(jetAK8PuppiCMVAv2Handle->at(i));              
		jetAK8PuppiPartonFlavour->push_back(jetAK8PuppiPartonFlavourHandle->at(i)); 
		  
		


		jetAK8PuppifilteredMass->push_back(jetAK8PuppifilteredMassHandle->at(i));   

		jetAK8PuppitrimmedMass->push_back(jetAK8PuppitrimmedMassHandle->at(i));    
		jetAK8PuppijecFactor0->push_back(jetAK8PuppijecFactor0Handle->at(i));   
   
		jetAK8Puppitau1->push_back(jetAK8Puppitau1Handle->at(i));      
		jetAK8Puppitau2->push_back(jetAK8Puppitau2Handle->at(i));      
		jetAK8Puppitau3->push_back(jetAK8Puppitau3Handle->at(i));    

		jetAK8PuppivSubjetIndex0->push_back(jetAK8PuppivSubjetIndex0Handle->at(i));   
		jetAK8PuppivSubjetIndex1->push_back(jetAK8PuppivSubjetIndex1Handle->at(i)); 
  



		jetAK8PuppisoftDropMassForTopRAW->push_back(SDFromCALC);
		jetAK8PuppisoftDropMassForTopuncorr->push_back(SDFromCALCuncorr);


		sM0->push_back(ptsmfacM0);
		sM1->push_back(ptsmfacM1);
		sSDM0->push_back(smfac);

	  	float CHF = jetAK8PuppichargedHadronEnergyFracHandle->at(i);
	  	float NEMF = jetAK8PuppineutralEmEnergyFracHandle->at(i);
	  	float NHF = jetAK8PuppineutralHadronEnergyFracHandle->at(i);
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


  delete jecUncPUPPI;
  delete jecUncCHS;
  /*
  iEvent.put(jetAK8CHSPt,"jetAK8CHSPt");
  iEvent.put(jetAK8CHSPhi,"jetAK8CHSPhi"); 
  iEvent.put(jetAK8CHSEta,"jetAK8CHSEta");  
  iEvent.put(jetAK8CHSMass,"jetAK8CHSMass");


  iEvent.put(jetAK8CHSprunedMass,"jetAK8CHSprunedMass");
  iEvent.put(jetAK8CHSsoftDropMass,"jetAK8CHSsoftDropMass");
  iEvent.put(jetAK8CHSsoftDropMassuncorr,"jetAK8CHSsoftDropMassuncorr");
  iEvent.put(jetAK8CHSprunedMassuncorr,"jetAK8CHSprunedMassuncorr");
  */

  iEvent.put(jetAK8PuppiPt,"jetAK8PuppiPt");
  iEvent.put(jetAK8PuppiPhi,"jetAK8PuppiPhi"); 
  iEvent.put(jetAK8PuppiEta,"jetAK8PuppiEta");  
  iEvent.put(jetAK8PuppiMass,"jetAK8PuppiMass");


 // iEvent.put(jetAK8PuppiprunedMass,"jetAK8PuppiprunedMass");
  iEvent.put(jetAK8PuppisoftDropMass,"jetAK8PuppisoftDropMass");
  iEvent.put(jetAK8PuppiCorrectedsoftDropMass,"jetAK8PuppiCorrectedsoftDropMass");
  //iEvent.put(jetAK8CHSCorrectedsoftDropMassPuppi,"jetAK8CHSCorrectedsoftDropMassPuppi");

  iEvent.put(jetAK8PuppiCorrectedsoftDropMassUnsmear,"jetAK8PuppiCorrectedsoftDropMassUnsmear");
 // iEvent.put(jetAK8CHSCorrectedsoftDropMassUnsmearPuppi,"jetAK8CHSCorrectedsoftDropMassUnsmearPuppi");



  iEvent.put(jetAK8PuppisoftDropMassuncorr,"jetAK8PuppisoftDropMassuncorr");
  //iEvent.put(jetAK8PuppiprunedMassuncorr,"jetAK8PuppiprunedMassuncorr");
  iEvent.put(jetAK8PuppisoftDropMassForTopPUPPIAK8JEC,"jetAK8PuppisoftDropMassForTopPUPPIAK8JEC");

   if (jes_=="nominal"&&jer_=="nominal")
	{
	/*
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
	*/
  	iEvent.put(jetAK8PuppiCSV,"jetAK8PuppiCSV"); 
  	iEvent.put(jetAK8PuppiCMVAv2,"jetAK8PuppiCMVAv2");

  	iEvent.put(jetAK8PuppiDoubleBAK8,"jetAK8PuppiDoubleBAK8"); 
 
  	iEvent.put(jetAK8PuppiPartonFlavour,"jetAK8PuppiPartonFlavour"); 
  	//iEvent.put(jetAK8PuppifilteredMass,"jetAK8PuppifilteredMass"); 

  	iEvent.put(jetAK8PuppiTight,"jetAK8PuppiTight"); 
  	iEvent.put(jetAK8PuppiLoose,"jetAK8PuppiLoose");


  	//iEvent.put(jetAK8PuppitopMass,"jetAK8PuppitopMass");
  	//iEvent.put(jetAK8PuppitrimmedMass,"jetAK8PuppitrimmedMass");
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

	//delete jecAK8CHS_
	//delete jecAK8PUPPI_

	iEvent.put(jetAK8PuppisoftDropMassForTopRAW,"jetAK8PuppisoftDropMassForTopRAW");
	iEvent.put(jetAK8PuppisoftDropMassForTopuncorr,"jetAK8PuppisoftDropMassForTopuncorr");


	iEvent.put(sM0,"sM0");
	iEvent.put(sM1,"sM1");
	iEvent.put(sSDM0,"sSDM0");

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
