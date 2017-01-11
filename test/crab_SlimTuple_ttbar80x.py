from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TT_TuneCUETP8M1_13TeV-powheg-pythia8_B2GAnaFW_V2p1_80x_Slim_V6_v1'
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
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/grauco-B2GAnaFW_80X_V2p1-edbed0685401a5848e7d61871b3a63d8/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 16
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/group/lpcrutgers/knash'
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.blacklist = ['T3_US_UMiss']
