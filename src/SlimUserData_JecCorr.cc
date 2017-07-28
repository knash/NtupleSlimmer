#include "SlimUserData_JecCorr.h"

float JES_Uncert(JetCorrectionUncertainty *jecUnc,float pt,float eta,std::string val)
	{
	int sign = 0;
	if (val=="up") sign = 1;
	else if (val=="down") sign = -1;


  	jecUnc->setJetEta(eta);
	jecUnc->setJetPt(pt);
	float shift = (1+sign*jecUnc->getUncertainty(true));
	return  shift;
	}

float JEC_Corr(boost::shared_ptr<FactorizedJetCorrector> jecAK8_,float  pt,float eta  , bool isdata,float Area,const double Rho)
	{


        jecAK8_->setJetEta( eta );
        jecAK8_->setJetPt ( pt );
        jecAK8_->setJetA  ( Area );
        jecAK8_->setRho   ( Rho );

        float corr = jecAK8_->getCorrection();

	return corr;
	}
