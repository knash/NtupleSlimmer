#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include <Math/VectorUtil.h>

#include <vector>


class  SlimUserData_weights : public edm::EDProducer {
public:
  SlimUserData_weights( const edm::ParameterSet & );   

private:
  int lha_pdf_id_;
  void produce( edm::Event &, const edm::EventSetup & );
  void beginRun(edm::Run const&, edm::EventSetup const&);
  void endJob() ;
  bool ISDATA_ ;
  std::string lhe_label_;
  // gen_token;
  edm::EDGetTokenT<GenEventInfoProduct> gen_token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
 };

namespace LHAPDF {
      void initPDFSet(int nset, const std::string& filename, int member=0);
      int numberPDF(int nset);
      void usePDFMember(int nset, int member);
      double xfx(int nset, double x, double Q, int fl);
      double getXmin(int nset, int member);
      double getXmax(int nset, int member);
      double getQ2min(int nset, int member);
      double getQ2max(int nset, int member);
      void extrapolate(bool extrapolate=true);
}
SlimUserData_weights::SlimUserData_weights(const edm::ParameterSet& iConfig):
   ISDATA_ (iConfig.getUntrackedParameter<bool>("ISDATA", false)),
   lhe_label_ (iConfig.getParameter<std::string>("lhe_label"))
 {   
    edm::EDGetTokenT<LHERunInfoProduct>(mayConsume<LHERunInfoProduct, edm::InRun>(edm::InputTag(lhe_label_, "")));
    edm::EDGetTokenT<GenEventInfoProduct>(consumes<GenEventInfoProduct>(edm::InputTag("generator")));
    edm::EDGetTokenT<LHEEventProduct>(consumes<LHEEventProduct>(edm::InputTag(lhe_label_, "")));   	
    produces<std::vector<float>>("pdfWeights");
    produces<std::vector<float>>("pdfWeightsNNPDF");
  // produces<std::vector<float>>("pdfWeightsCT10");
   //produces<std::vector<float>>("pdfWeightsMSTW2008");
  // produces<std::vector<float>>("pdfId");
 }

void SlimUserData_weights::produce( edm::Event& iEvent, const edm::EventSetup& iSetup) {
	
  	std::auto_ptr<std::vector<float>> pdfWeights(new std::vector<float>());   
  	std::auto_ptr<std::vector<float>> pdfWeights_temp(new std::vector<float>());    
  	//std::auto_ptr<std::vector<float>> pdfWeights_NNPDF(new std::vector<float>());   
  	//std::auto_ptr<std::vector<float>> pdfWeights_CT10(new std::vector<float>());  
  	//std::auto_ptr<std::vector<float>> pdfWeights_MSTW2008(new std::vector<float>());   

  	//int pdfId; 
  	std::auto_ptr<std::vector<std::string>> names(new std::vector<std::string>({"pdfWeightsNNPDF"}));//,"pdfWeightsCT10","pdfWeightsMSTW2008"})) ;


	//std::cout<<"PDFID"<<std::endl;
	//std::cout<<pdfId<<std::endl;
  	if (!ISDATA_) 
		{
    		edm::Handle<LHEEventProduct> lheEvtInfo;
    		iEvent.getByLabel(lhe_label_, lheEvtInfo);

    		//lha_pdf_id_ = lheRunInfo->heprup().PDFSUP.first;
    		//std::cout<<  lheEvtInfo.isValid() << std::endl;
		//std::cout<< lha_pdf_id_ << std::endl << std::endl;

    		if (lheEvtInfo.isValid()) 
			{

			//pdfId = lha_pdf_id_;
      			double lheOrigWeight = lheEvtInfo->originalXWGTUP();
      			size_t first = 9;
      			if (lha_pdf_id_ == 263000) first = 10;
      			if (lha_pdf_id_ == 263400) first = 111;
      			if (lheEvtInfo->weights().size()>=first+100) for (size_t i=first; i<first+100; ++i)
				{
        			pdfWeights->push_back(lheEvtInfo->weights()[i].wgt/lheOrigWeight);
     				}
    			}
		}


      		edm::Handle<GenEventInfoProduct> pdfstuff;
      		//gen_token = consumes<GenEventInfoProduct>(edm::InputTag("generator"));

      		//gen_token = iC.consumes<GenEventInfoProduct>(edm::InputTag("generator"));

      		iEvent.getByToken(gen_token, pdfstuff);



      		float q = pdfstuff->pdf()->scalePDF;

      		int id1 = pdfstuff->pdf()->id.first;
      		double x1 = pdfstuff->pdf()->x.first;
      	//	double pdf1 = pdfstuff->pdf()->xPDF.first;

      		int id2 = pdfstuff->pdf()->id.second;
      		double x2 = pdfstuff->pdf()->x.second;
      		//double pdf2 = pdfstuff->pdf()->xPDF.second;

		for(int ipdf=1; ipdf <=1; ++ipdf)
		{			
			
			




			LHAPDF::usePDFMember(ipdf,0);
			double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
			double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
			double w0 = xpdf1 * xpdf2;

			
			for(int i=1; i <=100; ++i)
				{

				LHAPDF::usePDFMember(ipdf,i);
		   		double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
		   		double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
		   		double weight = xpdf1_new * xpdf2_new / w0;
									
			
				pdfWeights_temp->push_back(weight);
				}

  			iEvent.put(pdfWeights_temp,names->at(ipdf-1));
			pdfWeights_temp.reset();
		}

		
  		iEvent.put(pdfWeights,"pdfWeights");
  		//iEvent.put(pdfId*,"pdfId");
}





 

void 
SlimUserData_weights::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  if (!ISDATA_) {
    LHAPDF::initPDFSet(1, "NNPDF30_lo_as_0130.LHgrid");
   // LHAPDF::initPDFSet(2, "CT10.LHgrid");
    //LHAPDF::initPDFSet(2, "MSTW2008nnlo68cl.LHgrid");

    edm::Handle<LHERunInfoProduct> lheRunInfo;
    iRun.getByLabel(lhe_label_, lheRunInfo);
    
    if (lheRunInfo.isValid()) {


      lha_pdf_id_ = lheRunInfo->heprup().PDFSUP.first;
      std::cout<<"LHE: LHA PDF ID = "<<lha_pdf_id_<<std::endl;
      std::cout<<"LHE:   --> For more info about the sets, check: https://lhapdf.hepforge.org/pdfsets.html"<<std::endl;
      
      // Check headers
      std::cout<<"LHE: Weight info in header:"<<std::endl;
      LHERunInfoProduct lheRunInfoProduct = *(lheRunInfo.product());
      typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
      size_t iHead = 0;
      for (headers_const_iterator header=lheRunInfoProduct.headers_begin(); header!=lheRunInfoProduct.headers_end(); header++){
        if (header->tag()=="initrwgt") {
          std::cout<<"LHE: "<<iHead<<" "<<header->tag()<<std::endl;
          for (auto line : header->lines()) {
	    std::cout<<"LHE: "<<line;
	    // Fix buggy powheg samples
	    if (lha_pdf_id_==-1 && line.find("weight id=\"2001\"")!=std::string::npos) {
	      if (line.find("PDF set = 260001")!=std::string::npos) lha_pdf_id_ = 260000;
	      else if (line.find("PDF set = 260401")!=std::string::npos) lha_pdf_id_ = 260400;
	    }
	  }
        }
        iHead++;
      }
      
    }
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
SlimUserData_weights::endJob() {
}


#include "FWCore/Framework/interface/MakerMacros.h"


DEFINE_FWK_MODULE(SlimUserData_weights);
