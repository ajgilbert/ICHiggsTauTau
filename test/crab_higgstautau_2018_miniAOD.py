from CRABClient.UserUtilities import config
from CRABClient.UserUtilities import getUsernameFromCRIC
from multiprocessing import Process

config = config()

config.General.transferOutputs = True
config.General.workArea='Sep18_MC_102X_2018'

config.JobType.psetName = 'higgstautau_cfg_102X_Aug19_2018.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['EventTree.root']
#config.JobType.maxMemoryMB = 2500
cfgParams = ['isData=0', 'globalTag=102X_upgrade2018_realistic_v20','doHT=1']
config.JobType.allowUndistributedCMSSW = True
#config.Data.unitsPerJob = 1
#config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.outLFNDirBase='/store/user/{}/{}/'.format(getUsernameFromCRIC(), config.General.workArea)
#config.Data.outLFNDirBase='/store/user/dwinterb/{}/'.format(config.General.workArea)
config.Data.allowNonValidInputDataset = True
# config.Data.inputDBS = 'phys03'
config.Data.ignoreLocality = True


config.Site.whitelist   = ['T2_*','T1_*','T3_*']
config.Site.storageSite = 'T2_UK_London_IC'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    from CRABClient.ClientExceptions import ClientException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException, hte:
            print hte.headers
        except ClientException as cle:
            print cle

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    tasks=list()

    tasks.append(('DYJetsToLL', '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('DYJetsToLL_M-10-50-LO', '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('DY1JetsToLL-LO', '/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('DY2JetsToLL-LO', '/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('DY3JetsToLL-LO', '/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('DY4JetsToLL-LO', '/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('DYJetsToLL-LO', '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('EWKZ2Jets', '/EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM')) 
    tasks.append(('EWKWMinus2Jets', '/EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('EWKWPlus2Jets', '/EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('W1JetsToLNu-LO', '/W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('W2JetsToLNu-LO', '/W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('W3JetsToLNu-LO', '/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('W4JetsToLNu-LO', '/W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WJetsToLNu-LO', '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM')) 
    tasks.append(('TTTo2L2Nu', '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('TTToHadronic', '/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('TTToSemiLeptonic', '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WZTo1L1Nu2Q', '/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WZTo2L2Q', '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WZTo3LNu', '/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WZTo3LNu-ext1', '/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'))
    tasks.append(('WZTo1L3Nu', '/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('WWTo1L1Nu2Q', '/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WWToLNuQQ', '/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WWTo2L2Nu', '/WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('ZZTo2L2Nu-ext1', '/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'))
    tasks.append(('ZZTo2L2Nu-ext2', '/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM'))
    tasks.append(('ZZTo2L2Q', '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ZZTo4L', '/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v2/MINIAODSIM'))
    tasks.append(('ZZTo4L-ext', '/ZZTo4L_TuneCP5_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM'))

    tasks.append(('WGToLNuG','/WGToLNuG_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('Tbar-tW-ext1', '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'))
    tasks.append(('T-tW-ext1', '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'))
    tasks.append(('Tbar-t', '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('T-t', '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('GluGluHToTauTau_M-125', '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHToTauTau_M-125-ext1', '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'))
    tasks.append(('ZHToTauTau_M-125', '/ZHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WminusHToTauTau_M-125', '/WminusHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WplusHToTauTau_M-125', '/WplusHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))

    tasks.append(('GluGluHToTauTauUncorrelatedDecay_Filtered','/GluGluHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('GluGluHToTauTauUncorrelatedDecay','/GluGluHToTauTauUncorrelatedDecay_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('VBFHToTauTauUncorrelatedDecay_Filtered','/VBFHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('VBFHToTauTauUncorrelatedDecay','/VBFHToTauTauUncorrelatedDecay_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WminusHToTauTauUncorrelatedDecay_Filtered','/WminusHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WminusHToTauTauUncorrelatedDecay','/WminusHToTauTauUncorrelatedDecay_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WplusHToTauTauUncorrelatedDecay_Filtered','/WplusHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('WplusHToTauTauUncorrelatedDecay','/WplusHToTauTauUncorrelatedDecay_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ZHToTauTauUncorrelatedDecay_Filtered','/ZHToTauTauUncorrelatedDecay_Filtered_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('JJH0MToTauTauPlusOneJets_Filtered','/JJH0MToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0MToTauTauPlusOneJets','/JJH0MToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0MToTauTauPlusTwoJets_Filtered','/JJH0MToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0MToTauTauPlusTwoJets','/JJH0MToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0MToTauTauPlusZeroJets_Filtered','/JJH0MToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0MToTauTauPlusZeroJets','/JJH0MToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusOneJets_Filtered','/JJH0Mf05ph0ToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusOneJets','/JJH0Mf05ph0ToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered','/JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusTwoJets','/JJH0Mf05ph0ToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered','/JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0Mf05ph0ToTauTauPlusZeroJets','/JJH0Mf05ph0ToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusOneJets_Filtered','/JJH0PMToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusOneJets','/JJH0PMToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusTwoJets_Filtered','/JJH0PMToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusTwoJets','/JJH0PMToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusZeroJets_Filtered','/JJH0PMToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('JJH0PMToTauTauPlusZeroJets','/JJH0PMToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('VBFHiggs0L1ToTauTau','/VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0L1ZgToTauTau','/VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0L1Zgf05ph0ToTauTau','/VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0L1f05ph0ToTauTau','/VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0MToTauTau','/VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0Mf05ph0ToTauTau','/VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0PHToTauTau','/VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0PHf05ph0ToTauTau','/VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('VBFHiggs0PMToTauTau','/VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0L1ToTauTau','/WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0L1f05ph0ToTauTau','/WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0MToTauTau','/WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0Mf05ph0ToTauTau','/WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0PHToTauTau','/WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0PHf05ph0ToTauTau','/WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('WHiggs0PMToTauTau','/WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0L1ToTauTau','/ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0L1ZgToTauTau','/ZHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0L1Zgf05ph0ToTauTau','/ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0L1f05ph0ToTauTau','/ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0MToTauTau','/ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0Mf05ph0ToTauTau','/ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0PHToTauTau','/ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0PHf05ph0ToTauTau','/ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ZHiggs0PMToTauTau','/ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ttHiggs0MToTauTau','/ttHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ttHiggs0Mf05ph0ToTauTau','/ttHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))
    tasks.append(('ttHiggs0PMToTauTau','/ttHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'))

    tasks.append(('GluGluHToWWTo2L2Nu_M-125', '/GluGluHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('VBFHToWWTo2L2Nu_M-125', '/VBFHToWWTo2L2Nu_M125_13TeV_powheg2_JHUGenV714_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    tasks.append(('VVTo2L2Nu', '/VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ggZH_HToTauTau_ZToLL', '/ggZH_HToTauTau_ZToLL_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ggZH_HToTauTau_ZToNuNu', '/ggZH_HToTauTau_ZToNuNu_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ggZH_HToTauTau_ZToQQ', '/ggZH_HToTauTau_ZToQQ_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('ttHToTauTau', '/ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('HWminusJ_HToWW', '/HWminusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('HWplusJ_HToWW', '/HWplusJ_HToWW_M125_13TeV_powheg_jhugen724_pythia8_TuneCP5/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('HZJ_HToWW', '/HZJ_HToWW_M125_13TeV_powheg_jhugen714_pythia8_TuneCP5/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))
    tasks.append(('GluGluZH_HToWW', '/GluGluZH_HToWW_M125_13TeV_powheg_pythia8_TuneCP5_PSweights/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'))

    for task in tasks:
        print task[0]
        config.General.requestName = task[0]
        config.Data.inputDataset = task[1]

        if "DY4JetsToLL-LO" in task[0]:
            # talk to DBS to get list of files in this dataset
            from dbs.apis.dbsClient import DbsApi
            dbs = DbsApi('https://cmsweb.cern.ch/dbs/prod/global/DBSReader')

            dataset = task[1]
            fileDictList=dbs.listFiles(dataset=dataset)

            print ("dataset %s has %d files" % (dataset, len(fileDictList)))

            # DBS client returns a list of dictionaries, but we want a list of Logical File Names
            lfnList = [ dic['logical_file_name'] for dic in fileDictList ]
            config.Data.userInputFiles = lfnList
            config.Data.splitting = 'FileBased'
            config.Data.unitsPerJob = 1
            config.Data.inputDataset = None
        else:
            config.Data.inputDataset = task[1]

            #config.Data.unitsPerJob = 100000
            #config.Data.splitting = 'EventAwareLumiBased'
            config.Data.splitting = 'FileBased'
            config.Data.unitsPerJob = 1

            config.Data.userInputFiles = None

        if ("HToTauTau" in task[0] or 'JJH' in task[0]) and 'JHUGen' not in task[1]:
            if 'mcatnloFXFX' in task[0]:
                config.JobType.pyCfgParams = cfgParams + ['LHEWeights=True','includenpNLO=True','includeHTXS=True']
            else:
                config.JobType.pyCfgParams = cfgParams + ['LHEWeights=True','tauSpinner=True','includeHTXS=True']
        else:
            config.JobType.pyCfgParams = cfgParams

        print(config)

        p = Process(target=submit, args=(config,))
        p.start()
        p.join()

