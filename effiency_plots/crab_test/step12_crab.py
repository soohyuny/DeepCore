import CRABClient
from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'DeepCore_11923_step12_1k'
config.General.workArea = 'workflows_crab'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'cmsDriver_step12.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores=8
config.JobType.maxMemoryMB=15000

config.section_("Data")
config.Data.outputPrimaryDataset = 'DeepCore_11923_step12_1k'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 100
NJOBS = 10 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'DeepCore_11923_step12_1k'

config.section_("Site")

#config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.storageSite = 'T3_US_FNALLPC'
