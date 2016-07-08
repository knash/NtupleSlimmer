from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TT_TuneCUETP8M1_13TeV-powheg-pythia8_B2GAnaFW_V1p1_80x_Slim_V1'
config.General.workArea = 'SlimTuples'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'quickntuple.py'
config.JobType.pyCfgParams = ['type=MC']
config.JobType.inputFiles = 	[
				"Spring16_25nsV6_DATA_L1FastJet_AK8PFchs.txt",
				"Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt",
				"Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt",
				"Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt",
				"Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt",
				"Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt",
				"Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt",
				"Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt",
				"Spring16_25nsV6_MC_SF_AK8PFchs.txt"
				]

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/srappocc-RunIISpring16MiniAODv2_B2GAnaFW_80x_V1p0-4e74e3854bbd13b3866f4a57304f402f/USER'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 8000
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

