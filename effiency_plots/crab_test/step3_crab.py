import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DeepCore_11923_s3_1k_T23_0214_thr1'
config.General.workArea = 'workflow_crab_step3'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
#config.JobType.pluginName = 'PrivateMC'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'cmsDriver_step3.py'
#config.JobType.psetName = 'cmsDriver_step3_jetcore.py'
#config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores=8
config.JobType.maxMemoryMB=20000

config.section_("Data")
config.Data.inputDataset = '/DeepCore_11923_step12_1k/hboucham-DeepCore_11923_step12_1k-9d61c57e09a1759f7fec6afa547af408/USER'
#config.Data.inputDataset = '/DeepCore_test_30/hboucham-DeepCore_test_30-9d61c57e09a1759f7fec6afa547af408/USER' 
#config.Data.outputPrimaryDataset = 'DeepCore_12323_1k_step3'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'EventBased'
config.Data.inputDBS = 'phys03'
config.Data.unitsPerJob = 1
#NJOBS = 10 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
#config.Data.publication = True
config.Data.publication = False
config.Data.outputDatasetTag = 'DeepCore_11923_s3_1k_T23_0214_thr1'

config.section_("Site")

#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T3_US_FNALLPC'

config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

