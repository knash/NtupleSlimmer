#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"


float JES_Uncert(JetCorrectionUncertainty *jecUnc,float pt,float eta,std::string val);
float JEC_Corr(boost::shared_ptr<FactorizedJetCorrector> jecAK8_,float  pt,float eta  , bool isdata,float Area,const double Rho);

