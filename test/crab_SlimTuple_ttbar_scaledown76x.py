from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8_B2GAnaFW_V2p1_76x_Slim_V12'
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
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8/lcorcodi-crab_TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8_B2GAnaFW_76X_V2p1-ce108855676c8f127cd37fa14e406ad2/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

