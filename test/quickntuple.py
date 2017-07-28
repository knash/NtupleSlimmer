
import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy
import sys
import subprocess
options = opts.VarParsing ('analysis')

options.register('sample',
		#'/store/user/knash/WprimeToTB_TToHad_M-1200_RH_TuneCUETP8M1_13TeV-comphep-pythia8/crab_WPrime13TeV_B2GAnaFW_V8p4_M1200_RH_25ns/151113_172838/0000/B2GEDMNtuple_1.root',
		#'/store/user/knash/JetHT/crab_JetHT_Run2015D-PromptReco-v4_B2GAnaFW_V8p4_25ns_JECv7_v2/160324_125554/0000/B2GEDMNtuple_482.root',
		#'/store/user/lcorcodi/BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/crab_BstarToTW_M-1400_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/160318_162851/0000/B2GEDMNtuple_2.root',#
		#'/store/group/lpcrutgers/knash/WprimeToTB_TToHad_M-1500_RH_TuneCUETP8M1_13TeV-comphep-pythia8/RunIISpring16MiniAODv2_80X_reHLT_B2GAnaFW_80X_V2p1/161109_215328/0000/B2GEDMNtuple_47.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/CMSSW_7_6_3_patch2/src/Analysis/B2GAnaFW/test/B2GEDMNtuple.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/CMSSW_8_0_24_patch1/src/Analysis/NtupleSlimmer/test/B2GEDMNtuple_MC.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/WithPUPPI/CMSSW_8_0_24_patch1/src/Analysis/NtupleSlimmer/test/B2GEDMNtuple_74.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/WithPUPPI/CMSSW_8_0_24_patch1/src/Analysis/NtupleSlimmer/test/B2GEDMNtuple_40.root',
		#'/store/group/phys_b2g/B2GAnaFW_80X_V2p3/JetHT/Run2016F/JetHT/Run2016F-23Sep2016-v1_B2GAnaFW_80X_V2p3/161216_222052/0000/B2GEDMNtuple_103.root',
		#'/store/group/lpcrutgers/knash/SMS-T7WgStealth_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2_B2GAnaFW_80x_V2p4/170221_164255/0001/B2GEDMNtuple_1399.root',
		#'/store/user/lcorcodi/BstarToTW_M-2000_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2_PUMoriond17_B2GAnaFW_80X_V2p4/170212_044129/0000/B2GEDMNtuple_20.root',
		#'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/WithPUPPI/CMSSW_8_0_24_patch1/src/Analysis/NtupleSlimmer/test/B2GEDMNtuple_30.root',
		'file:///uscms_data/d3/knash/WPrime13TeV/B2GAnaFW/SlimNtuples_test/WithPUPPI/CMSSW_8_0_24_patch1/src/Analysis/NtupleSlimmer/test/B2GEDMNtuple_100.root',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample to analyze')


options.register('type',
                 'NONE',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Sample type (DATA or MC)')

options.register('genfilter',
                 'OFF',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'cut out non chargino')


options.register('outputlabel',
                 '',
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 'Output label')

### Events to process: 'maxEvents' is already registered by the framework

options.parseArguments()


process = cms.Process("slimntuple")







process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
        options.sample
        ),
)

#process.TriggerUserData = cms.EDProducer('SlimUserData_TriggerUserData')



process.jetsAK8 = cms.EDProducer(
    	'SlimUserData_jetsAK8',
   	reapplyjec = cms.bool(True),
   	reapplyjer = cms.bool(True),
	jes  = cms.string("nominal"),
	jer  = cms.string("nominal"),
	era  = cms.string("Summer16_23Sep2016")
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









process.jetsAK4 = cms.EDProducer(
    	'SlimUserData_jetsAK4',
   	reapplyjec = cms.bool(True),
   	reapplyjer = cms.bool(True),
	jes  = cms.string("nominal"),
	jer  = cms.string("nominal"),
	era  = cms.string("Summer16_23Sep2016")
    )


process.jetsAK4jesup = process.jetsAK4.clone(
	jes  = cms.string("up"),
      )

process.jetsAK4jesdown = process.jetsAK4.clone(
	jes  = cms.string("down"),
      )


process.jetsAK4jerup = process.jetsAK4.clone(
	jer  = cms.string("up"),
      )

process.jetsAK4jerdown = process.jetsAK4.clone(
	jer  = cms.string("down"),
      )







#process.EventCounter = cms.EDAnalyzer("EventCounter")

process.counter = cms.EDProducer('SlimUserData_counter')

process.Filter = cms.EDFilter('SlimUserData_Filter')
if options.type=='MC':
	process.p = cms.Path(process.counter)

	if options.genfilter == 'ON':
		print "WARNING"
		print "GENFILTER ON - Only for use with SMS signal!"
		process.GENFilter = cms.EDFilter(
			'SlimUserData_GENFilter',
			ISDATA  = cms.untracked.bool(False)
			)


		process.p*=process.GENFilter


	process.weights = cms.EDProducer(
    		'SlimUserData_weights',
		ISDATA  = cms.untracked.bool(False),
		lhe_label = cms.string("externalLHEProducer"),
    		)

	process.Filter = cms.EDFilter(
		'SlimUserData_Filter',
		ISDATA  = cms.untracked.bool(False)
		)
	process.p*=process.Filter*process.weights*process.jetsAK8*process.jetsAK8jesup*process.jetsAK8jesdown*process.jetsAK8jerup*process.jetsAK8jerdown*process.jetsAK4*process.jetsAK4jesup*process.jetsAK4jesdown*process.jetsAK4jerup*process.jetsAK4jerdown		

	
elif options.type=='DATA':
	process.Filter = cms.EDFilter(
		'SlimUserData_Filter',
		ISDATA  = cms.untracked.bool(True)
		)
	process.p = cms.Path(
        	#process.EventCounter
		process.counter
		*process.Filter
		*process.jetsAK8
		*process.jetsAK4
	    	)
else:
	sys.exit("Must be either DATA or MC for type")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.f1 = cms.Path(process.Filter)
if options.genfilter == 'ON':
	process.f1 *= process.GENFilter
process.edmNtuplesOut = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('B2GEDMNtuple_slim'+options.outputlabel+'.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents =  cms.vstring('f1')),
    outputCommands = cms.untracked.vstring(
    "drop *",
    "keep *_*_*_slimntuple",
    "drop edmTriggerResults_*_*_slimntuple",
   # "keep *_subjetsCmsTopTag_subjetCmsTopTagCSV_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSCSV*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSCMVA*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSPartonFlavour*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSEta*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSE*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSPhi*_*",
    "keep *_subjetsAK8CHS_subjetAK8CHSPt*_*",
    #"keep *_jetsAK4Puppi_jetAK4CHSCSVv2*_*",
    #"keep *_jetsAK4Puppi_jetAK4CHSEta*_*",
    #"keep *_jetsAK4Puppi_jetAK4CHSE_*",
    #"keep *_jetsAK4Puppi_jetAK4CHSPhi*_*",
    #"keep *_jetsAK4Puppi_jetAK4CHSPt_*", 
    "keep *_subjetsAK8Puppi_subjetAK8PuppiCSV*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiCMVA*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiPartonFlavour*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiEta*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiE*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiPhi*_*",
    "keep *_subjetsAK8Puppi_subjetAK8PuppiPt*_*",
    #"keep *_jetsAK4Puppi_jetAK4PuppiCSVv2*_*",
    #"keep *_jetsAK4Puppi_jetAK4PuppiEta*_*",
    #"keep *_jetsAK4Puppi_jetAK4PuppiMass*_*",
    ##"keep *_jetsAK4Puppi_jetAK4PuppiPhi*_*",
    #"keep *_jetsAK4Puppi_jetAK4PuppiPt*_*",
    #"keep *_genPart_genPartDau0ID*_*",  
    #"keep *_genPart_genPartDau1ID*_*",  
    "keep *_genPart_genPartEta*_*",  
    "keep *_genPart_genPartID*_*",  
    "keep *_genPart_genPartE*_*",  
    #"keep *_genPart_genPartMom0ID*_*",  
    #"keep *_genPart_genPartMom1ID*_*",  
    "keep *_genPart_genPartPhi*_*",  
    "keep *_genPart_genPartPt*_*",  
    "keep *_genPart_genPartStatus*_*",  
    "keep *_eventUserData_*_*",
    "drop *_eventUserData_v*_*",
    "keep *_generator_*_*",
    "keep *_eventInfo_*_*",
    "keep *_electrons_elEta_*",      
    "keep *_electrons_elE_*",      
    "keep *_electrons_elPhi_*",    
    "keep *_electrons_elPt_*",      
    #"keep *_electrons_elisLoose_*",      
    #"keep *_electrons_elisMedium_*",      
    #"keep *_electrons_elisTight_*",   
   # "keep *_electrons_elvidLoose_*",   
    "keep *_electrons_elvidMedium_*",   
   # "keep *_electrons_elvidTight_*",   
    "keep *_muons_muEta_*",   
    "keep *_muons_muE_*",   
    "keep *_muons_muPhi_*",    
    "keep *_muons_muPt_*",   
    #"keep *_muons_muIsLooseMuon_*",    
    "keep *_muons_muIsMediumMuon_*",   
    #"keep *_muons_muIsTightMuon_*",   
    #"keep *_metFull_metFullPhi_*",   
    #"keep *_metFull_metFullPt_*",   
    "keep *_metNoHF_metNoHFPhi_*",   
    "keep *_metNoHF_metNoHFPt_*",  
    "keep *_pdfWeights*_*_*",
    "keep *_filteredPrunedGenParticles_*_*"
    ),
    dropMetaData = cms.untracked.string('ALL'),
    )

   


#process.Timing = cms.Service("Timing", summaryOnly=cms.untracked.bool(True))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",ignoreTotal = cms.untracked.int32(1) )
process.endPath = cms.EndPath(process.edmNtuplesOut)


