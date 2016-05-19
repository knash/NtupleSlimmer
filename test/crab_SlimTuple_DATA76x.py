from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'JetHT_Run2015D-16Dec2015_B2GAnaFW_V1p3_76x_Slim_V12_v1'
config.General.workArea = 'SlimTuples'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'quickntuple.py'
config.JobType.pyCfgParams = ['type=DATA']
config.JobType.inputFiles = 	[
				"Fall15_25nsV2_DATA_L2Relative_AK8PFchs.txt",
				"Fall15_25nsV2_DATA_L3Absolute_AK8PFchs.txt",
				"Fall15_25nsV2_DATA_L2L3Residual_AK8PFchs.txt",
				"Fall15_25nsV2_MC_L2Relative_AK8PFchs.txt",
				"Fall15_25nsV2_MC_L3Absolute_AK8PFchs.txt",
				"Fall15_25nsV2_MC_Uncertainty_AK8PFchs.txt"
				]

config.section_("Data")
config.Data.inputDataset = '/JetHT/devdatta-B2GAnaFW_76X_V1p3-c4c3f6dc16a584b6ec5462c242b1cbe4/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 40
config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True
config.Data.publication = False

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'

