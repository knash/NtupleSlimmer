from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_V1p1_76x_Slim_V12'
config.General.workArea = 'SlimTuples'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'quickntuple.py'
config.JobType.pyCfgParams = ['type=MC']
config.JobType.inputFiles = 	[
				"Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt",
				"Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt",
				"Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt",
				"Fall15_25nsV2_MC_L2Relative_AK8PFchs.txt",
				"Fall15_25nsV2_MC_L3Absolute_AK8PFchs.txt",
				"Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt"
				]

config.section_("Data")
config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/jkarancs-B2GAnaFW_76X_V1p1_RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1-bf3ef703e3bdb5dcc5320cf3ff6ce74d/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 3
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

