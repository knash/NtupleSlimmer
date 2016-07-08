
import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy

options = opts.VarParsing ('analysis')

options.register('sample',
                 #'/store/user/knash/JetHT/crab_JetHT_Run2015D-PromptReco-v4_B2GAnaFW_V8p4_25ns_v1/151222_180124/0000/B2GEDMNtuple_59.root',
              #   '/store/user/knash/JetHT/crab_JetHT_Run2015D-PromptReco-v4_B2GAnaFW_V8p4_25ns_v1/151222_180124/0000/B2GEDMNtuple_126.root',
		# '/store/user/knash/JetHT/crab_JetHT_Run2015D-05Oct2015-v1_B2GAnaFW_V8p4_25ns/151210_174906/0000/B2GEDMNtuple_9.root',
		#'/store/user/knash/WprimeToTB_TToHad_M-1200_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V8p4_M1200_RH_25ns/151113_172838/0000/B2GEDMNtuple_1.root',
		#'/store/user/knash/JetHT/crab_JetHT_Run2015D-PromptReco-v4_B2GAnaFW_V8p4_25ns_JECv7_v2/160324_125554/0000/B2GEDMNtuple_482.root',
		#'/store/user/lcorcodi/BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/crab_BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/160318_162851/0000/B2GEDMNtuple_2.root',#
		'/store/group/lpctlbsm/B2GAnaFW_80X_V1p0/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0/160609_173944/0001/B2GEDMNtuple_1141.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/CMSSW_7_6_3_patch2/src/Analysis/B2GAnaFW/test/B2GEDMNtuple.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')


options.register('type',
                 'NONE',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample type (DATA or MC)')


options.register('outputlabel',
                 '',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

### Events to process: 'maxEvents' is already registered by the framework

options.parseArguments()
process = cms.Process("slimntuple")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        options.sample
        )
)

#process.TriggerUserData = cms.EDProducer('SlimUserData_TriggerUserData')



process.jetsAK8 = cms.EDProducer(
    	'SlimUserData_jetsAK8',
   	reapplyjec = cms.bool(True),
   	reapplyjer = cms.bool(True),
	jes  = cms.string("nominal"),
	jer  = cms.string("nominal"),
	era  = cms.string("Spring16_25nsV6")
    )


process.jetsAK8jesup = process.jetsAK8.clone(
	jes  = cms.string("up"),
      )

process.jetsAK8jesdown = process.jetsAK8.clone(
	jes  = cms.string("down"),
      )


process.jetsAK8jerup = process.jetsAK8.clone(
	jer  = cms.string("up"),
      )

process.jetsAK8jerdown = process.jetsAK8.clone(
	jer  = cms.string("down"),
      )



process.Filter = cms.EDFilter('SlimUserData_Filter')
if options.type=='MC':

	process.weights = cms.EDProducer(
    		'SlimUserData_weights',
		ISDATA  = cms.untracked.bool(False),
		lhe_label = cms.string("externalLHEProducer"),
    		)

    	#process.pdfWeights = cms.EDProducer("PdfWeightProducer",
    	#	#  FixPOWHEG = cms.untracked.string(""),
      	#	GenTag = cms.untracked.InputTag("genPart"),
      	#	PdfInfoTag = cms.untracked.InputTag("generator"),
      	#	PdfSetNames = cms.untracked.vstring( 	   "CT10.LHgrid" 
	#						  ,"MSTW2008nnlo68cl.LHgrid"
        #		                                  ,"NNPDF23_nnlo_as_0118.LHgrid"
        #		                           )
	#	)

	process.p = cms.Path(
		process.weights
		*process.Filter
		*process.jetsAK8
		*process.jetsAK8jesup
		*process.jetsAK8jesdown
		*process.jetsAK8jerup
		*process.jetsAK8jerdown	

	    	)
elif options.type=='DATA':
	process.p = cms.Path(
		process.Filter
		*process.jetsAK8
	    	)
else:
	sys.exit("Must be either DATA or MC for type")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.f1 = cms.Path(process.Filter)
process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('B2GEDMNtuple_slim'+options.outputlabel+'.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents =  cms.vstring('f1')),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_*_*_slimntuple",
    "drop edmTriggerResults_*_*_slimntuple",
    "keep *_subjetsCmsTopTag_subjetCmsTopTagCSV_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSCSV*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSCMVA*_*",
    "keep *_eventUserData_*_*",
    "drop *_eventUserData_v*_*",
    "keep *_eventInfo_*_*",
   # "keep *_genPart_*Eta_*",
   # "keep *_genPart_*Pt_*",
   # "keep *_genPart_*Phi_*",
   # "keep *_genPart_*Mass_*",
   # "keep *_genPart_*ID_*",
   # "keep *_genPart_*Status_*",
    "keep *_met_*_*",
    "drop *_met_Px_*",    
    "drop *_met_Py_*",
  #  "keep *_genJetsAK8_*_*",
    "keep *_pdfWeights*_*_*",
    ),
    dropMetaData = cms.untracked.string('ALL'),
    )
#process.Timing = cms.Service("Timing", summaryOnly=cms.untracked.bool(True))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.endPath = cms.EndPath(process.edmNtuplesOut)


