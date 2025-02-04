import FWCore.ParameterSet.Config as cms
process = cms.Process("MAIN")
import sys

################################################################
# Read Options
################################################################
import FWCore.ParameterSet.VarParsing as parser
opts = parser.VarParsing ('analysis')

opts.register('file',
# 'root://xrootd.unl.edu//store/user/jbechtel/gc_storage/TauTau_data_2017_CMSSW944/TauEmbedding_TauTau_data_2017_CMSSW944_Run2017B/1/merged_0.root_'
# 'root://xrootd.unl.edu//store/user/sbrommer/gc_storage/ElTau_data_2017_CMSSW944/TauEmbedding_ElTau_data_2017_CMSSW944_Run2017B/25/merged_1524.root_',
'root://xrootd.unl.edu//store/mc/RunIIFall17MiniAODv2/VBFHToTauTau_M125_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/280000/E2C688CA-1F1C-E911-9847-0CC47AD98BC2.root',
# 'root://xrootd.unl.edu//store/mc/RunIIFall17MiniAODv2/SUSYGluGluToHToTauTau_M-120_TuneCP5_13TeV-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/C4C4D050-DE41-E811-A2A6-0025905B85B6.root'
#'root://xrootd.unl.edu///store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/80000/248C8431-B838-E811-B418-0025905B85D2.root',
parser.VarParsing.multiplicity.singleton,
parser.VarParsing.varType.string, "input file")
opts.register('globalTag', '102X_dataRun2_v12', parser.VarParsing.multiplicity.singleton,
# opts.register('globalTag', '102X_mc2017_realistic_v6', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "global tag")
opts.register('isData', 1, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as data?")
opts.register('isEmbed', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Process as embedded?")
opts.register('LHEWeights', False, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.bool, "Produce LHE weights for sample")
opts.register('LHETag', 'externalLHEProducer', parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.string, "Input tag for LHE weights")
opts.register('doHT', 0, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.int, "Store HT and number of outgoing partons?")
opts.register('includenpNLO', False, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.bool, "Store npNLO for sample (number of partons for NLO sample)")
opts.register('tauSpinner', False, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.bool, "Compute weights using tauspinner")
opts.register('includeHTXS', False, parser.VarParsing.multiplicity.singleton,
    parser.VarParsing.varType.bool, "Compute HTXS Gen level vars")

opts.parseArguments()
infile      = opts.file
if not infile: infile = "file:/tmp/file.root"
isData      = opts.isData
isEmbed      = opts.isEmbed
if isEmbed: isData = 0
tag         = opts.globalTag
doLHEWeights = opts.LHEWeights
if not isData: doHT     = opts.doHT
else: doHT = 0
includenpNLO = opts.includenpNLO

print 'isData      : '+str(isData)
print 'isEmbed      : '+str(isEmbed)
print 'globalTag   : '+str(tag)
print 'doHT        : '+str(doHT)


################################################################
# Standard setup
################################################################
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("EventTree.root"),
    closeFileFast = cms.untracked.bool(True)
)

################################################################
# Message Logging, summary, and number of events
################################################################
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000) 
)

process.MessageLogger.cerr.FwkReport.reportEvery = 50


################################################################
# Input files and global tags
################################################################
process.load("CondCore.CondDB.CondDB_cfi")
from CondCore.CondDB.CondDB_cfi import *

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(infile))
process.GlobalTag.globaltag = cms.string(tag)

process.options   = cms.untracked.PSet(
    FailPath=cms.untracked.vstring("FileReadError"),
    wantSummary = cms.untracked.bool(True),
    # numberOfThreads = cms.untracked.uint32(4)
)

import UserCode.ICHiggsTauTau.default_producers_cfi as producers
import UserCode.ICHiggsTauTau.default_selectors_cfi as selectors

process.load("PhysicsTools.PatAlgos.patSequences_cff")

################################################################
# Re-do PFTau reconstruction
################################################################
process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")

################################################################
# Object Selection
################################################################

# embed new tauID trainings

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig

tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, debug = False,
                    updatedTauName = updatedTauName,
                    toKeep = [
                            "2017v2", "dR0p32017v2", "newDM2017v2", #classic MVAIso tau-Ids
                            "deepTau2017v2p1",
                            "MVADM_2016_v1","MVADM_2017_v1", 
                               ])
tauIdEmbedder.runTauID()

process.selectedElectrons = cms.EDFilter("PATElectronRefSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 9.5 & abs(eta) < 2.6")
    )
process.selectedMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("slimmedMuons"),
    cut = cms.string("pt > 3 & abs(eta) < 2.6")
    )
process.selectedTaus = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag("slimmedTausNewID"),
    cut = cms.string('pt > 18.0 & abs(eta) < 2.6 & tauID("decayModeFindingNewDMs") > 0.5')
    )

process.icSelectionSequence = cms.Sequence()

process.icSelectionSequence += cms.Sequence(
    process.selectedElectrons+
    process.selectedMuons+
    process.rerunMvaIsolationSequence+
    getattr(process,updatedTauName)+
    process.selectedTaus
)

################################################################
# PF sequence for lepton isolation
################################################################
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")

process.pfPileUp = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("fromPV <= 1")
    )
process.pfNoPileUp = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string("fromPV > 1")
    )
process.pfAllNeutralHadrons = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfNoPileUp"),
    cut = cms.string("abs(pdgId) = 111 | abs(pdgId) = 130 | " \
        "abs(pdgId) = 310 | abs(pdgId) = 2112")
    )
process.pfAllChargedHadrons= cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfNoPileUp"),
    cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
        "abs(pdgId) = 999211 | abs(pdgId) = 2212")
    )
process.pfAllPhotons= cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfNoPileUp"),
    cut = cms.string("abs(pdgId) = 22")
    )
process.pfAllChargedParticles= cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfNoPileUp"),
    cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
        "abs(pdgId) = 999211 | abs(pdgId) = 2212 | " \
        "abs(pdgId) = 11 | abs(pdgId) = 13")
    )
process.pfPileUpAllChargedParticles= cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfPileUp"),
    cut = cms.string("abs(pdgId) = 211 | abs(pdgId) = 321 | " \
        "abs(pdgId) = 999211 | abs(pdgId) = 2212 | " \
        "abs(pdgId) = 11 | abs(pdgId) = 13")
    )
process.pfAllNeutralHadronsAndPhotons = cms.EDFilter("CandPtrSelector",
    src = cms.InputTag("pfNoPileUp"),
    cut = cms.string("abs(pdgId) = 111 | abs(pdgId) = 130 | " \
        "abs(pdgId) = 310 | abs(pdgId) = 2112 | abs(pdgId) = 22")
    )
process.pfParticleSelectionSequence = cms.Sequence(
    process.pfPileUp+
    process.pfNoPileUp+
    process.pfAllNeutralHadrons+
    process.pfAllChargedHadrons+
    process.pfAllPhotons+
    process.pfAllChargedParticles+
    process.pfPileUpAllChargedParticles+
    process.pfAllNeutralHadronsAndPhotons
    )

vtxLabel = cms.InputTag("offlineSlimmedPrimaryVertices")

################################################################
# Vertices
################################################################
process.icVertexProducer = producers.icVertexProducer.clone(
    branch  = cms.string("vertices"),
    input = vtxLabel,
    firstVertexOnly = cms.bool(True)
)

process.icGenVertexProducer = producers.icGenVertexProducer.clone(
    input = cms.InputTag("prunedGenParticles")
)

### refit PV excluding tau decay products tracks

import VertexRefit.TauRefit.AdvancedRefitVertexProducer_cfi as vertexrefit
import VertexRefit.TauRefit.MiniAODRefitVertexProducer_cfi as vertexfit

process.refitOfflineSlimmedPrimaryVerticesBSWithAllTracks = vertexfit.MiniAODRefitVertexBSProducer.clone()
process.icBSVertexProducer = producers.icVertexProducer.clone(
    branch  = cms.string("verticesBS"),
    input = cms.InputTag("refitOfflineSlimmedPrimaryVerticesBSWithAllTracks"),
    firstVertexOnly = cms.bool(True)
)

process.refitOfflineSlimmedPrimaryVertices = vertexrefit.AdvancedRefitVertexNoBSProducer.clone()
process.refitOfflineSlimmedPrimaryVertices.storeAsMap = cms.bool(True)
process.refitOfflineSlimmedPrimaryVertices.srcLeptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", updatedTauName)

process.refitOfflineSlimmedPrimaryVerticesBS = vertexrefit.AdvancedRefitVertexBSProducer.clone()
process.refitOfflineSlimmedPrimaryVerticesBS.storeAsMap = cms.bool(True)
process.refitOfflineSlimmedPrimaryVerticesBS.srcLeptons = cms.VInputTag("slimmedElectrons", "slimmedMuons", updatedTauName)


process.icRefitVertexProducer = producers.icRefitVertexProducer.clone(
  branch  = cms.string("refittedVertices"),
  input = "refitOfflineSlimmedPrimaryVertices",
  firstVertexOnly = cms.bool(False)
)

process.icRefitVertexProducerBS = producers.icRefitVertexProducer.clone(
  branch  = cms.string("refittedVerticesBS"),
  input = "refitOfflineSlimmedPrimaryVerticesBS",
  firstVertexOnly = cms.bool(False)
)

process.icVertexSequence = cms.Sequence(
  process.icVertexProducer+
  process.icGenVertexProducer+
  process.refitOfflineSlimmedPrimaryVerticesBSWithAllTracks+
  process.icBSVertexProducer+
  process.refitOfflineSlimmedPrimaryVertices+
  process.refitOfflineSlimmedPrimaryVerticesBS+
  process.icRefitVertexProducer+
  process.icRefitVertexProducerBS
)

if isData :
    process.icVertexSequence.remove(process.icGenVertexProducer)

################################################################
# PFCandidates
################################################################
process.icPFFromPackedProducer = cms.EDProducer('ICPFFromPackedProducer',
    branch  = cms.string("pfCandidates"),
    input   = cms.InputTag("icTauProducer", "requestedPFCandidates"),
    requestTracks       = cms.bool(False),
    requestGsfTracks    = cms.bool(False),
    inputUnpackedTracks = cms.InputTag("unpackedTracksAndVertices")
    )

process.icPFSequence = cms.Sequence()
process.icPFSequence += process.icPFFromPackedProducer

################################################################
# Electrons
################################################################
electronLabel = cms.InputTag("slimmedElectrons")

process.icElectronSequence = cms.Sequence()

# electron smear and scale
from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       # applyEnergyCorrections=True,
                       # applyVIDOnCorrectedEgamma=True,
                       runVID=False, #saves CPU time by not needlessly re-running VID
                       era='2017-Nov17ReReco')
process.icElectronSequence += cms.Sequence(
    process.egammaPostRecoSeq
    )

process.icElectronConversionCalculator = cms.EDProducer('ICElectronConversionCalculator',
    input       = electronLabel,
    beamspot    = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("reducedEgamma:reducedConversions")
)

process.load("RecoEgamma.ElectronIdentification.ElectronMVAValueMapProducer_cfi")
process.electronMVAValueMapProducer.src = electronLabel
process.electronMVAValueMapProducer.srcMiniAOD = electronLabel

# trying to add things for ID Fall17v2
process.electronMVAVariableHelper.src = electronLabel
process.electronMVAVariableHelper.srcMiniAOD = electronLabel
process.electronMVAVariableHelper.beamSpot = cms.InputTag("offlineBeamSpot")
process.electronMVAVariableHelper.beamSpotMiniAOD = cms.InputTag("offlineBeamSpot")
process.electronMVAVariableHelper.vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.electronMVAVariableHelper.vertexCollectionMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
process.electronMVAVariableHelper.conversions = cms.InputTag("reducedEgamma:reducedConversions")
process.electronMVAVariableHelper.conversionsMiniAOD = cms.InputTag("reducedEgamma:reducedConversions")

process.icElectronSequence+=cms.Sequence(
    process.electronMVAVariableHelper+
    process.electronMVAValueMapProducer
)

#Electron PF iso sequence:
process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff")
process.elPFIsoValueChargedAll03PFIdPFIso = process.elPFIsoValueChargedAll03PFId.clone()
process.elPFIsoDepositChargedAll.src  = electronLabel

process.elPFIsoValueCharged03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("charged_iso")
)
process.elPFIsoValueGamma03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("photon_iso")
)
process.elPFIsoValueNeutral03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("neutral_iso")
)
process.elPFIsoValuePU03PFIdPFIso = cms.EDProducer('ICRecoElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("pu_iso")
)

process.electronPFIsolationValuesSequence = cms.Sequence(
    process.elPFIsoValueCharged03PFIdPFIso+
    process.elPFIsoValueChargedAll03PFIdPFIso+
    process.elPFIsoValueGamma03PFIdPFIso+
    process.elPFIsoValueNeutral03PFIdPFIso+
    process.elPFIsoValuePU03PFIdPFIso
)


process.elPFIsoDepositGamma.ExtractorPSet.ComponentName = cms.string("CandViewExtractor")
process.icElectronSequence += cms.Sequence(
    process.elPFIsoDepositChargedAll+
    process.electronPFIsolationValuesSequence
)

process.elPFIsoValueChargedAll03PFIdPFIso.deposits[0].vetos = (
    cms.vstring('EcalEndcaps:ConeVeto(0.015)','EcalBarrel:ConeVeto(0.01)'))

process.elPFIsoValueCharged04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("charged_iso")
)
process.elPFIsoValueGamma04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("photon_iso")
)
process.elPFIsoValueNeutral04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("neutral_iso")
)
process.elPFIsoValuePU04PFIdPFIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("pu_iso")
)

process.elEcalPFClusterIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type     = cms.string("ecal_pf_cluster_iso")
)

process.elHcalPFClusterIso = cms.EDProducer('ICElectronIsolation',
    input        = electronLabel,
    deltaR       = cms.double(0.3),
    iso_type     = cms.string("hcal_pf_cluster_iso")
)

process.elPFIsoValueChargedAll04PFIdPFIso = process.elPFIsoValueChargedAll04PFId.clone()
process.elPFIsoValueChargedAll04PFIdPFIso.deposits[0].vetos = (
    cms.vstring('EcalEndcaps:ConeVeto(0.015)','EcalBarrel:ConeVeto(0.01)'))

process.electronPFIsolationValuesSequence +=cms.Sequence(
    process.elPFIsoValueCharged04PFIdPFIso+
    process.elPFIsoValueGamma04PFIdPFIso+
    process.elPFIsoValuePU04PFIdPFIso+
    process.elPFIsoValueNeutral04PFIdPFIso+
    process.elPFIsoValueChargedAll04PFIdPFIso+
    process.elEcalPFClusterIso+
    process.elHcalPFClusterIso
)

import RecoEgamma.EgammaTools.calibratedEgammas_cff
calibratedPatElectrons94Xv1 = RecoEgamma.EgammaTools.calibratedEgammas_cff.calibratedPatElectrons.clone(
    produceCalibratedObjs = False
    )

process.icElectronProducer = producers.icElectronProducer.clone(
    branch                    = cms.string("electrons"),
    input                     = cms.InputTag("selectedElectrons"),
    includeConversionMatches  = cms.bool(True),
    inputConversionMatches    = cms.InputTag("icElectronConversionCalculator"),
    includeVertexIP           = cms.bool(True),
    inputVertices             = vtxLabel,
    includeBeamspotIP         = cms.bool(True),
    inputBeamspot             = cms.InputTag("offlineBeamSpot"),
    inputPostCorr             = cms.string("ecalTrkEnergyPostCorr"),
    inputPreCorr              = cms.string("ecalTrkEnergyPreCorr"),
    inputErrPostCorr          = cms.string("ecalTrkEnergyErrPostCorr"),
    inputErrPreCorr           = cms.string("ecalTrkEnergyErrPreCorr"),
    inputScaleUp              = cms.string("energyScaleUp"),
    inputScaleDown            = cms.string("energyScaleDown"),
    inputSigmaUp              = cms.string("energySigmaUp"),
    inputSigmaDown            = cms.string("energySigmaDown"),
    includeFloats = cms.PSet(
       generalPurposeMVASpring16  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
       mvaRun2Fall17  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
       mvaRun2IsoFall17  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
       mvaRun2Fall17V2  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues"), # need raw values for v2
       mvaRun2IsoFall17V2  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2RawValues"),
    ),
    includeClusterIso        = cms.bool(True),
    includePFIso03           = cms.bool(True),
    includePFIso04           = cms.bool(True),
)


process.icElectronSequence += cms.Sequence(
    process.icElectronConversionCalculator+
    process.icElectronProducer
)

################################################################
# Muons
################################################################
process.icMuonSequence = cms.Sequence()
muons = cms.InputTag("slimmedMuons")

process.load("CommonTools.ParticleFlow.Isolation.pfMuonIsolation_cff")
process.load("CommonTools.ParticleFlow.deltaBetaWeights_cff")
process.muPFIsoDepositChargedAll.src  = muons
process.muPFIsoDepositNeutral.src     = muons
process.muPFIsoDepositNeutral.ExtractorPSet.inputCandView = cms.InputTag("pfWeightedNeutralHadrons")
process.muPFIsoDepositGamma.src       = muons
process.muPFIsoDepositGamma.ExtractorPSet.inputCandView = cms.InputTag("pfWeightedPhotons")
process.muPFIsoValueChargedAll03PFIso = process.muPFIsoValueChargedAll03.clone()
process.muPFIsoValueNeutral04PFWeights = process.muPFIsoValueNeutral04.clone()
process.muPFIsoValueGamma04PFWeights = process.muPFIsoValueGamma04.clone()
process.muPFIsoValueNeutral03PFWeights = process.muPFIsoValueNeutral03.clone()
process.muPFIsoValueGamma03PFWeights = process.muPFIsoValueGamma03.clone()

process.muPFIsoValueCharged03PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("charged_iso")
)
process.muPFIsoValueGamma03PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("photon_iso")
)
process.muPFIsoValueNeutral03PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("neutral_iso")
)
process.muPFIsoValuePU03PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.3),
    iso_type = cms.string("pu_iso")
)

process.muonPFIsolationValuesSequence = cms.Sequence(
    process.muPFIsoValueCharged03PFIso+
    process.muPFIsoValueChargedAll03PFIso+
    process.muPFIsoValueGamma03PFIso+
    process.muPFIsoValueNeutral03PFIso+
    process.muPFIsoValuePU03PFIso+
    process.muPFIsoValueNeutral04PFWeights+
    process.muPFIsoValueGamma04PFWeights+
    process.muPFIsoValueNeutral03PFWeights+
    process.muPFIsoValueGamma03PFWeights
)

process.muPFIsoValueCharged04PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("charged_iso")
)
process.muPFIsoValueGamma04PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("photon_iso")
)
process.muPFIsoValueNeutral04PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("neutral_iso")
)
process.muPFIsoValuePU04PFIso = cms.EDProducer('ICMuonIsolation',
    input        = muons,
    deltaR       = cms.double(0.4),
    iso_type = cms.string("pu_iso")
)
process.muPFIsoValueChargedAll04PFIso = process.muPFIsoValueChargedAll04.clone()
process.muonPFIsolationValuesSequence +=cms.Sequence(
    process.muPFIsoValueCharged04PFIso+
    process.muPFIsoValueGamma04PFIso+
    process.muPFIsoValuePU04PFIso+
    process.muPFIsoValueNeutral04PFIso+
    process.muPFIsoValueChargedAll04PFIso
)

process.icMuonSequence += cms.Sequence(
    process.pfDeltaBetaWeightingSequence+
    process.muPFIsoDepositChargedAll+
    process.muPFIsoDepositNeutral+
    process.muPFIsoDepositGamma+
    process.muonPFIsolationValuesSequence
)

process.icMuonProducer = producers.icMuonProducer.clone(
    branch                    = cms.string("muons"),
    input                     = cms.InputTag("selectedMuons"),
    isPF                      = cms.bool(False),
    includeVertexIP           = cms.bool(True),
    inputVertices             = vtxLabel,
    includeBeamspotIP         = cms.bool(True),
    inputBeamspot             = cms.InputTag("offlineBeamSpot"),
    includeDoubles = cms.PSet(
        neutral_pfweighted_iso_03 = cms.InputTag("muPFIsoValueNeutral03PFWeights"),
        neutral_pfweighted_iso_04 = cms.InputTag("muPFIsoValueNeutral04PFWeights"),
        gamma_pfweighted_iso_03 = cms.InputTag("muPFIsoValueGamma03PFWeights"),
        gamma_pfweighted_iso_04 = cms.InputTag("muPFIsoValueGamma04PFWeights")
    ),
    requestTracks           = cms.bool(False),
    includePFIso03           = cms.bool(True),
    includePFIso04           = cms.bool(True),
)
process.icMuonProducer.isPF = cms.bool(False)

process.icMuonSequence += cms.Sequence(
    process.icMuonProducer
)

################################################################
# Taus
################################################################

process.icTauProducer = cms.EDProducer("ICPFTauFromPatProducer",
    branch                  = cms.string("taus"),
    input                   = cms.InputTag("selectedTaus"),
    inputVertices           = vtxLabel,
    includeVertexIP         = cms.bool(True),
    requestTracks           = cms.bool(False),
    includeTotalCharged     = cms.bool(False),
    totalChargedLabel       = cms.string('totalCharged'),
    requestPFCandidates     = cms.bool(True),
    inputPFCandidates       = cms.InputTag("packedPFCandidates"),
    isSlimmed               = cms.bool(True),
    tauIDs = cms.PSet(),
    doSVRefit             = cms.bool(True)
)


process.icTauSequence = cms.Sequence(
    process.icTauProducer
)

################################################################
# photons, clusters, and superclusters
################################################################

process.icPi0SuperClusterProducer = cms.EDProducer('ICPi0SuperClusterProducer',
  branch  = cms.string("superClusters"),
  input   = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
  input_EERecHits = cms.InputTag("reducedEgamma", "reducedEERecHits"),
  input_EBRecHits = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
  input_taus = cms.InputTag("slimmedTaus")
)

process.icPhotonProducer = producers.icPhotonProducer.clone(
    branch                  = cms.string("photons"),
    input                   = cms.InputTag("slimmedPhotons"),
    includeHadTowOverEm     = cms.bool(True),
    includePFIso03           = cms.bool(True),
    pfIso03 = cms.PSet(
      chargedAll  = cms.InputTag("phPFIsoValueChargedAll03PFIso"),
      charged     = cms.InputTag("phPFIsoValueCharged03PFIso"),
      neutral     = cms.InputTag("phPFIsoValueNeutral03PFIso"),
      gamma       = cms.InputTag("phPFIsoValueGamma03PFIso"),
      pu          = cms.InputTag("phPFIsoValuePU03PFIso")
    ),
    includeIsoFromPat = cms.bool(True)
)


process.icPhotonSequence = cms.Sequence(
  process.icPhotonProducer+
  process.icPi0SuperClusterProducer
)


# ################################################################
# # Jets
# ################################################################
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection


#rebuild ak4 chs jets as in  https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#Jets
process.load('PhysicsTools.PatAlgos.slimming.unpackedTracksAndVertices_cfi')

#Reapply JECs:
if not (isData or isEmbed):
  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJets"),
    labelName = "UpdatedJEC",
    jetCorrections = ("AK4PFchs", cms.vstring(['L1FastJet','L2Relative','L3Absolute']), 'None'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = [
      'pfDeepFlavourJetTags:probb',
      'pfDeepFlavourJetTags:probbb',
      'pfDeepFlavourJetTags:problepb',
      'pfDeepFlavourJetTags:probc',
      'pfDeepFlavourJetTags:probuds',
      'pfDeepFlavourJetTags:probg'
    ],
  )
else:
 updateJetCollection(
   process,
   jetSource = cms.InputTag("slimmedJets"),
   labelName = "UpdatedJEC",
   jetCorrections = ("AK4PFchs", cms.vstring(['L1FastJet','L2Relative','L3Absolute','L2L3Residual']), 'None'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   btagDiscriminators = [
     'pfDeepFlavourJetTags:probb',
     'pfDeepFlavourJetTags:probbb',
     'pfDeepFlavourJetTags:problepb',
     'pfDeepFlavourJetTags:probc',
     'pfDeepFlavourJetTags:probuds',
     'pfDeepFlavourJetTags:probg'
   ],
  )

process.selectedSlimmedJetsAK4 = cms.EDFilter("PATJetRefSelector",
    src = cms.InputTag("selectedUpdatedPatJetsUpdatedJEC"),
    cut = cms.string("pt > 15")
    )

 # Jet energy corrections
 # ----------------------
process.ak4PFL1FastjetCHS = cms.EDProducer("L1FastjetCorrectorProducer",
    srcRho = cms.InputTag("fixedGridRhoFastjetAll"),
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet')
)
process.ak4PFL2RelativeCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)
process.ak4PFL3AbsoluteCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)
process.ak4PFResidualCHS = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)

#Corrections applied to miniaod slimmedJets
pfJECS = cms.PSet(
    L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
    L2Relative = cms.string("ak4PFL2RelativeCHS"),
    L3Absolute = cms.string("ak4PFL3AbsoluteCHS")
)
if isData or isEmbed: pfJECS = cms.PSet(
    L1FastJet  = cms.string("ak4PFL1FastjetCHS"),
    L2Relative = cms.string("ak4PFL2RelativeCHS"),
    L3Absolute = cms.string("ak4PFL3AbsoluteCHS"),
    L2L3Residual = cms.string("ak4PFResidualCHS")
)

 # b-tagging
 # ---------
process.load("RecoJets.JetAssociationProducers.ak4JTA_cff")
from RecoJets.JetAssociationProducers.ak4JTA_cff import ak4JetTracksAssociatorAtVertex
process.load("RecoBTag.Configuration.RecoBTag_cff")
import RecoBTag.Configuration.RecoBTag_cff as btag

process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfos.candidates = cms.InputTag("packedPFCandidates")

 # Pileup ID
 # ---------
 # Recalculated puJetId isn't the same as miniaod stored - should investigate
#stdalgos = cms.VPSet()
#from RecoJets.JetProducers.PileupJetIDParams_cfi import *
#stdalgos = cms.VPSet(full_5x_chs,cutbased)
process.load('RecoJets.JetProducers.PileupJetID_cfi')

process.pileupJetIdCalculator.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pileupJetIdCalculator.jets = cms.InputTag("ak4PFJetsCHS")
process.pileupJetIdCalculator.rho = cms.InputTag("fixedGridRhoFastjetAll")
process.pileupJetIdEvaluator.jets = cms.InputTag("ak4PFJetsCHS")
process.pileupJetIdEvaluator.rho = cms.InputTag("fixedGridRhoFastjetAll")

 # JER
 # --------
process.load('Configuration.StandardSequences.Services_cff')
process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")

process.slimmedJetsSmeared = cms.EDProducer('SmearedPATJetProducer',
        src = cms.InputTag('selectedUpdatedPatJetsUpdatedJEC'),
        enabled = cms.bool(True),
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        algo = cms.string('AK4PFchs'),
        algopt = cms.string('AK4PFchs_pt'),
        #resolutionFile = cms.FileInPath('Autumn18_V7_MC_PtResolution_AK4PFchs.txt'),
        #scaleFactorFile = cms.FileInPath('combined_SFs_uncertSources.txt'),

        genJets = cms.InputTag('slimmedGenJets'),
        dRMax = cms.double(0.2),
        dPtMaxFactor = cms.double(3),

        debug = cms.untracked.bool(False),
        # Systematic variation
        # 0: Nominal
        # -1: -1 sigma (down variation)
        # 1: +1 sigma (up variation)
        variation = cms.int32(0),  # If not specified, default to 0
)

process.slimmedJetsSmearedDown = process.slimmedJetsSmeared.clone(variation=cms.int32(-1))
process.slimmedJetsSmearedUp = process.slimmedJetsSmeared.clone(variation=cms.int32(1))

# process.icPFJetProducerFromPat = producers.icPFJetFromPatProducer.clone(
#     branch                    = cms.string("ak4PFJetsCHS"),
#     input                     = cms.InputTag("selectedSlimmedJetsAK4"),
#     srcConfig = cms.PSet(
#         isSlimmed               = cms.bool(True),
#         slimmedPileupIDLabel    = cms.string('pileupJetId:fullDiscriminant'),
#         includeJetFlavour       = cms.bool(True),
#         includeJECs             = cms.bool(True),
#         inputSVInfo             = cms.InputTag(""),
#         requestSVInfo           = cms.bool(False)
#     ),
#    destConfig = cms.PSet(
#         includePileupID         = cms.bool(True),
#         inputPileupID           = cms.InputTag("puJetMva", "fullDiscriminant"),
#         includeTrackBasedVars   = cms.bool(False),
#         inputTracks             = cms.InputTag("unpackedTracksAndVertices"),
#         inputVertices           = cms.InputTag("unpackedTracksAndVertices"),
#         requestTracks           = cms.bool(False)
#     )
# )

process.icPFJetProducerFromPatNew = producers.icPFJetFromPatNewProducer.clone(
    branch                    = cms.string("ak4PFJetsCHS"),
    input                     = cms.InputTag("selectedSlimmedJetsAK4"),
    inputSmear                = cms.InputTag("slimmedJetsSmeared"),
    inputSmearDown            = cms.InputTag("slimmedJetsSmearedUp"),
    inputSmearUp              = cms.InputTag("slimmedJetsSmearedDown"),
    srcConfig = cms.PSet(
      isSlimmed               = cms.bool(True),
      slimmedPileupIDLabel    = cms.string('pileupJetId:fullDiscriminant'),
      includeJetFlavour       = cms.bool(True),
      includeJECs             = cms.bool(True),
      inputSVInfo             = cms.InputTag(""),
      requestSVInfo           = cms.bool(False)
    ),
   destConfig = cms.PSet(
     includePileupID         = cms.bool(True),
     inputPileupID           = cms.InputTag("puJetMva", "fullDiscriminant"),
     includeTrackBasedVars   = cms.bool(False),
     inputTracks             = cms.InputTag("unpackedTracksAndVertices"),
     inputVertices           = cms.InputTag("unpackedTracksAndVertices"),
     requestTracks           = cms.bool(False)
    )
)
if isData or isEmbed:
    process.icPFJetProducerFromPatNew.doSmear = cms.bool(False)
else:
    process.icPFJetProducerFromPatNew.doSmear = cms.bool(True)

process.icPFJetSequence = cms.Sequence()

process.icPFJetSequence += cms.Sequence(
    process.patJetCorrFactorsUpdatedJEC+
    process.updatedPatJetsUpdatedJEC+
    process.pfImpactParameterTagInfosUpdatedJEC+
    process.pfInclusiveSecondaryVertexFinderTagInfosUpdatedJEC+
    process.pfDeepCSVTagInfosUpdatedJEC+
    process.pfDeepFlavourTagInfosUpdatedJEC+
    process.pfDeepFlavourJetTagsUpdatedJEC+
    process.patJetCorrFactorsTransientCorrectedUpdatedJEC+
    process.updatedPatJetsTransientCorrectedUpdatedJEC+
    process.selectedUpdatedPatJetsUpdatedJEC+
    process.selectedSlimmedJetsAK4+
    # process.icPFJetProducerFromPat +
    process.slimmedJetsSmeared+
    process.slimmedJetsSmearedDown+
    process.slimmedJetsSmearedUp+
    process.icPFJetProducerFromPatNew
    )

################################################################
# PF MET
################################################################
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load("RecoJets.JetProducers.ak4PFJets_cfi")

from RecoMET.METProducers.PFMET_cfi import pfMet
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
        process,
        isData = (bool(isData) or bool(isEmbed)),
        isEmbeddedSample=bool(isEmbed),
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
)

from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
makePuppiesFromMiniAOD( process, True );
runMetCorAndUncFromMiniAOD(process,
                           isData=(bool(isData) or bool(isEmbed)),
                           isEmbeddedSample=bool(isEmbed),
                           metType="Puppi",
                           postfix="PuppiModifiedMET",
                           jetFlavor="AK4PFPuppi",
                           fixEE2017 = False, 
                           )

process.icPfMetProducer = producers.icMetFromPatProducer.clone(
                         branch = cms.string("pfMetFromSlimmed"),
                         input = cms.InputTag("slimmedMETsModifiedMET"), # for 2017 apply re-correction
                         getUncorrectedMet=cms.bool(False),
                         includeMetUncertainties=cms.bool(True)
                         )

################################################################
# Puppi MET
################################################################

process.icPuppiMetProducer = producers.icMetFromPatProducer.clone(
                         input=cms.InputTag("slimmedMETsPuppiModifiedMET"),
                         branch = cms.string("puppiMet"),
                         getUncorrectedMet=cms.bool(False),
                         includeMetUncertainties=cms.bool(True)
                         )

process.icMetSequence = cms.Sequence(
  process.puppiMETSequence *
  process.fullPatMetSequencePuppiModifiedMET *
  process.fullPatMetSequenceModifiedMET *
  process.icPfMetProducer *
  process.icPuppiMetProducer
)


################################################################
# Simulation only: GenParticles, GenJets, PileupInfo
################################################################

process.icGenSequence = cms.Sequence()

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.icGenParticleProducer = producers.icGenParticleProducer.clone(
    input   = cms.InputTag("prunedGenParticles","","PAT"),
    includeMothers = cms.bool(True),
    includeDaughters = cms.bool(True),
    includeStatusFlags = cms.bool(True)
)
if isEmbed: process.icGenParticleProducer.input   = cms.InputTag("prunedGenParticles","","MERGE")


process.icGenParticleProducerFromLHEParticles = producers.icGenParticleFromLHEParticlesProducer.clone()


process.load("RecoJets.Configuration.GenJetParticles_cff")
process.genParticlesForJets.ignoreParticleIDs = cms.vuint32(
    1000022, 2000012, 2000014,
    2000016, 1000039, 5000039,
    4000012, 9900012, 9900014,
    9900016, 39, 12, 14, 16
)
process.genParticlesForJets.src = cms.InputTag("packedGenParticles")
process.load("RecoJets.JetProducers.ak4GenJets_cfi")
process.ak4GenJetsNoNuBSM  =  process.ak4GenJets.clone()
process.ak4GenJetsNoNuBSM.src=cms.InputTag("packedGenParticles")

process.selectedGenJets = cms.EDFilter("GenJetRefSelector",
    src = cms.InputTag("ak4GenJetsNoNuBSM"),
    cut = cms.string("pt > 10.0")
)

process.icGenJetProducer = producers.icGenJetProducer.clone(
    branch  = cms.string("genJetsReclustered"),
    input   = cms.InputTag("selectedGenJets"),
    inputGenParticles = cms.InputTag("prunedGenParticles","","PAT"),
    requestGenParticles = cms.bool(False),
    isSlimmed  = cms.bool(True)
)
if isEmbed: process.icGenJetProducer.inputGenParticles = cms.InputTag("prunedGenParticles","","MERGE")


process.icGenJetProducerFromSlimmed = producers.icGenJetProducer.clone(
    branch = cms.string("genJets"),
    input = cms.InputTag("slimmedGenJets"),
    inputGenParticles=cms.InputTag("genParticles"),
    requestGenParticles = cms.bool(False),
    isSlimmed = cms.bool(True)
)

if isEmbed: process.icGenJetProducerFromSlimmed.input = cms.InputTag("selectedGenJets")

process.icPileupInfoProducer = producers.icPileupInfoProducer.clone()
process.icPileupInfoProducer.input=cms.InputTag("slimmedAddPileupInfo")


if not isData:
    process.icGenSequence += (
        process.icGenParticleProducer+
        process.ak4GenJetsNoNuBSM+
        process.selectedGenJets+
        process.icGenJetProducer+
        process.icGenJetProducerFromSlimmed
    )
    if not isEmbed: process.icGenSequence += (process.icPileupInfoProducer)
    if doHT:
        process.icGenSequence += (
            process.icGenParticleProducerFromLHEParticles
        )

################################################################
# L1 Objects
################################################################

v_doBXloop = False

process.icL1EGammaProducer = cms.EDProducer('ICL1TObjectProducer<l1t::EGamma>',
    branch = cms.string("L1EGammas"),
    input = cms.InputTag("caloStage2Digis","EGamma"),
    doBXloop = cms.bool(v_doBXloop)
)

process.icL1TauProducer = cms.EDProducer('ICL1TObjectProducer<l1t::Tau>',
    branch = cms.string("L1Taus"),
    input = cms.InputTag("caloStage2Digis","Tau"),
    doBXloop = cms.bool(v_doBXloop)
)

process.icL1MuonProducer = cms.EDProducer('ICL1TObjectProducer<l1t::Muon>',
    branch = cms.string("L1Muons"),
    input = cms.InputTag("gmtStage2Digis","Muon"),
    doBXloop = cms.bool(v_doBXloop)
)


# ################################################################
# # Trigger
# ################################################################

from PhysicsTools.PatAlgos.tools.trigTools import *
process.icTriggerSequence = cms.Sequence()
process.icTriggerObjectSequence = cms.Sequence()


process.patTriggerPath = cms.Path()

switchOnTrigger(process, path = 'patTriggerPath',  outputModule = '')


process.icTriggerPathProducer = producers.icTriggerPathProducer.clone(
    branch = cms.string("triggerPaths"),
    input  = cms.InputTag("TriggerResults","","HLT"),
    inputIsStandAlone = cms.bool(True),
    inputPrescales = cms.InputTag("patTrigger")
)

if isData:
    process.icTriggerSequence += cms.Sequence(
        process.icTriggerPathProducer
    )

# mt paths

process.icIsoMu24ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu24"),
    hltPath = cms.string("HLT_IsoMu24_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icIsoMu24erObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu24er"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )


process.icIsoMu27ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu27"),
    hltPath = cms.string("HLT_IsoMu27_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icIsoMu20Tau27ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu20Tau27"),
    hltPath = cms.string("HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

# mt monitoring paths

process.icDoubleEl24ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsDoubleEl24"),
    hltPath = cms.string("HLT_DoubleEle24_eta2p1_WPTight_Gsf_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icDoubleMu20ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsDoubleMu20"),
    hltPath = cms.string("HLT_DoubleIsoMu20_eta2p1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )


# et paths

process.icEle27ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsEle27"),
    hltPath = cms.string("HLT_Ele27_WPTight_Gsf_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icEle32ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsEle32"),
    hltPath = cms.string("HLT_Ele32_WPTight_Gsf_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icEle32L1DoubleEGObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsEle32L1DoubleEG"),
    hltPath = cms.string("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icEle35ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsEle35"),
    hltPath = cms.string("HLT_Ele35_WPTight_Gsf_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icEle24Tau30ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsEle24Tau30"),
    hltPath = cms.string("HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

# tt paths

process.icDoubleTightIsoTau35ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsDoubleTightIsoTau35"),
    hltPath = cms.string("HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icDoubleMediumIsoTau40ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsDoubleMediumIsoTau40"),
    hltPath = cms.string("HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icDoubleTightIsoTau40ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsDoubleTightIsoTau40"),
    hltPath = cms.string("HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icMediumChargedIsoPFTau180HighPtRelaxedIsoProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsMediumChargedIsoPFTau180HighPtRelaxedIso"),
    hltPath = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icMediumChargedIsoPFTau180HighPtRelaxedIso1prProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsMediumChargedIsoPFTau180HighPtRelaxedIso1pr"),
    hltPath = cms.string("HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_1pr_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

# tt monitoring paths

process.icMu24TightIsoTightIDTau35ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsMu24TightIsoTightIDTau35"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icMu24MediumIsoTau35ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsMu24MediumIsoTau35"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_CrossL1"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icMu24TightIsoTau35ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsMu24TightIsoTau35"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau35_Trk1_eta2p1_Reg_CrossL1"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

# vbf paths

process.icVBFDoubleLooseChargedIsoPFTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsVBFDoubleLooseChargedIsoPFTau20"),
      hltPath = cms.string("HLT_VBF_DoubleLooseChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
      )

process.icVBFDoubleMediumChargedIsoPFTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsVBFDoubleMediumChargedIsoPFTau20"),
      hltPath = cms.string("HLT_VBF_DoubleMediumChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
      )

process.icVBFDoubleTightChargedIsoPFTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsVBFDoubleTightChargedIsoPFTau20"),
      hltPath = cms.string("HLT_VBF_DoubleTightChargedIsoPFTau20_Trk1_eta2p1_Reg_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

# vbf monitoring paths

process.icIsoMu24LooseIsoTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu24LooseIsoTau20"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_LooseChargedIsoPFTau20_SingleL1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icIsoMu24MediumIsoTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu24MediumIsoTau20"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_MediumChargedIsoPFTau20_SingleL1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

process.icIsoMu24TightIsoTau20ObjectProducer = producers.icTriggerObjectProducer.clone(
    input   = cms.InputTag("selectedPatTrigger"),
    branch = cms.string("triggerObjectsIsoMu24TightIsoTau20"),
    hltPath = cms.string("HLT_IsoMu24_eta2p1_TightChargedIsoPFTau20_SingleL1_v"),
    inputIsStandAlone = cms.bool(True),
    storeOnlyIfFired = cms.bool(False)
    )

# em paths

process.icMu23Ele12ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu23Ele12"),
      hltPath = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu23Ele12DZObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu23Ele12DZ"),
      hltPath = cms.string("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu8Ele23ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu8Ele23"),
      hltPath = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu8Ele23DZObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu8Ele23DZ"),
      hltPath = cms.string("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu12Ele23ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu12Ele23"),
      hltPath = cms.string("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu12Ele23DZObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu12Ele23DZ"),
      hltPath = cms.string("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

# em monitoring paths

process.icEle24Ele12ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsEle24Ele12"),
      hltPath = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icEle24Ele12DZObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsEle24Ele12DZ"),
      hltPath = cms.string("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

# embed sel triggers

process.icMu17Mu8ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu17Mu8"),
      hltPath = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu17Mu8DZObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu17Mu8DZ"),
      hltPath = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)

process.icMu17Mu8DZmass8ObjectProducer = producers.icTriggerObjectProducer.clone(
      input   = cms.InputTag("selectedPatTrigger"),
      branch = cms.string("triggerObjectsMu17Mu8DZmass8"),
      hltPath = cms.string("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v"),
      inputIsStandAlone = cms.bool(True),
      storeOnlyIfFired = cms.bool(False)
)


process.icTriggerObjectSequence += cms.Sequence(
    process.icIsoMu24ObjectProducer+
    process.icIsoMu24erObjectProducer+
    process.icIsoMu27ObjectProducer+
    process.icIsoMu20Tau27ObjectProducer+
    process.icDoubleEl24ObjectProducer+
    process.icDoubleMu20ObjectProducer+
    process.icEle27ObjectProducer+
    process.icEle32ObjectProducer+
    process.icEle32L1DoubleEGObjectProducer+
    process.icEle35ObjectProducer+
    process.icEle24Tau30ObjectProducer+
    process.icDoubleTightIsoTau35ObjectProducer+
    process.icDoubleMediumIsoTau40ObjectProducer+
    process.icDoubleTightIsoTau40ObjectProducer+
    process.icMediumChargedIsoPFTau180HighPtRelaxedIsoProducer+
    process.icMediumChargedIsoPFTau180HighPtRelaxedIso1prProducer+
    process.icMu24TightIsoTightIDTau35ObjectProducer+
    process.icMu24MediumIsoTau35ObjectProducer+
    process.icMu24TightIsoTau35ObjectProducer+
    process.icVBFDoubleLooseChargedIsoPFTau20ObjectProducer+
    process.icVBFDoubleMediumChargedIsoPFTau20ObjectProducer+
    process.icVBFDoubleTightChargedIsoPFTau20ObjectProducer+
    process.icIsoMu24LooseIsoTau20ObjectProducer+
    process.icIsoMu24MediumIsoTau20ObjectProducer+
    process.icIsoMu24TightIsoTau20ObjectProducer+
    process.icMu23Ele12ObjectProducer+
    process.icMu23Ele12DZObjectProducer+
    process.icMu8Ele23ObjectProducer+
    process.icMu8Ele23DZObjectProducer+
    process.icMu12Ele23ObjectProducer+
    process.icMu12Ele23DZObjectProducer+
    process.icEle24Ele12ObjectProducer+
    process.icEle24Ele12DZObjectProducer+
    process.icMu17Mu8ObjectProducer+
    process.icMu17Mu8DZObjectProducer+
    process.icMu17Mu8DZmass8ObjectProducer
    )

for name in process.icTriggerObjectSequence.moduleNames():
    mod = getattr(process, name)
    if isEmbed: mod.inputTriggerResults = cms.InputTag("TriggerResults", "","SIMembedding")


## Need to unpack filterLabels on slimmedPatTrigger then make selectedPatTrigger
process.patTriggerUnpacker = cms.EDProducer("PATTriggerObjectStandAloneUnpacker",
    patTriggerObjectsStandAlone = cms.InputTag("slimmedPatTrigger"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    unpackFilterLabels = cms.bool(True)
)

if isEmbed: process.patTriggerUnpacker.triggerResults = cms.InputTag("TriggerResults", "", "SIMembedding")

process.selectedPatTrigger = cms.EDFilter(
    'PATTriggerObjectStandAloneSelector',
    cut = cms.string('!filterLabels.empty()'),
    src = cms.InputTag('patTriggerUnpacker')
)
process.icTriggerSequence += cms.Sequence(
    process.patTriggerUnpacker +
    process.selectedPatTrigger
)


################################################################
# EventInfo
################################################################
#Load the MET filters here
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')

process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi') # outdated on miniAOD / rerun
baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])
process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist,
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
    )

if opts.LHETag: lheTag = opts.LHETag
else: lheTag = 'externalLHEProducer'

data_type = ""
if isData: data_type = "RECO"
elif isEmbed: data_type = "MERGE"
else: data_type = "PAT"

## add prefiring weights
from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

process.icEventInfoProducer = producers.icEventInfoProducer.clone(
    includePrefireWeights = cms.bool(not bool(isData) and not bool(isEmbed)),
    includeJetRho       = cms.bool(True),
    includeLHEWeights   = cms.bool(doLHEWeights),
    includeGenWeights   = cms.bool(doLHEWeights),
    includenpNLO        = cms.bool(includenpNLO),
    includeHTXS         = cms.bool(opts.includeHTXS),
    includeEmbeddingWeights = cms.bool(bool(isEmbed)),
    includeHT           = cms.bool(False),
    lheProducer         = cms.InputTag(lheTag),
    inputJetRho         = cms.InputTag("fixedGridRhoFastjetAll"),
    includeLeptonRho    = cms.bool(True),
    inputLeptonRho      = cms.InputTag("fixedGridRhoFastjetAll"),
    includeVertexCount  = cms.bool(True),
    inputVertices       = vtxLabel,
    includeCSCFilter    = cms.bool(False),
    inputCSCFilter      = cms.InputTag("BeamHaloSummary"),
    includeFiltersFromTrig = cms.bool(True),
    inputfiltersfromtrig = cms.InputTag("TriggerResults","",data_type),
    filters             = cms.PSet(
        badChargedHadronFilter  = cms.InputTag("BadChargedCandidateFilter"),
        badMuonFilter          = cms.InputTag("BadPFMuonFilter"),
        ecalBadCalibReducedMINIAODFilter = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),
    ),
    filtersfromtrig     = cms.vstring(
        "Flag_goodVertices","Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter"
        "Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter",
        "Flag_BadChargedCandidateFilter","Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter","ecalBadCalibReducedMINIAODFilter",
        )
)

process.icEventInfoSequence = cms.Sequence(
    process.ecalBadCalibReducedMINIAODFilter+
    process.BadPFMuonFilter+
    process.BadChargedCandidateFilter+
    process.prefiringweight+
    process.icEventInfoProducer
)

################################################################
# TauSpinner
################################################################

process.icTauSpinnerProducer = cms.EDProducer("ICTauSpinnerProducer",
    branch                  = cms.string("tauspinner"),
    input                   = cms.InputTag("prunedGenParticles"),
    theta                   = cms.string("0,0.25,0.5,-0.25,0.375")
)

if opts.tauSpinner:
    process.icTauSpinnerSequence = cms.Sequence(
        process.icTauSpinnerProducer
    )
else: process.icTauSpinnerSequence = cms.Sequence()

################################################################
# HTXS NNLOPS
################################################################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
    inputPruned = cms.InputTag("prunedGenParticles"),
    inputPacked = cms.InputTag("packedGenParticles"),
)
process.genParticles2HepMC = cms.EDProducer("GenParticles2HepMCConverter",
    genParticles = cms.InputTag("mergedGenParticles"),
    genEventInfo = cms.InputTag("generator"),
    signalParticlePdgIds = cms.vint32(25),
)
process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
    HepMCCollection = cms.InputTag('genParticles2HepMC','unsmeared'),
    LHERunInfo = cms.InputTag('externalLHEProducer'),
    ProductionMode = cms.string('AUTO'),
)
process.icHtxsSequence = cms.Sequence()
if opts.includeHTXS:
    process.icHtxsSequence = cms.Sequence(
        process.mergedGenParticles *
        process.genParticles2HepMC *
        process.rivetProducerHTXS
    )

################################################################
# Event
################################################################
process.icEventProducer = producers.icEventProducer.clone()


process.p = cms.Path(
    process.icSelectionSequence+
    process.pfParticleSelectionSequence+
    process.icVertexSequence+
    process.icElectronSequence+
    process.icMuonSequence+
    process.icTauSequence+
    process.icMetSequence+
    process.icPFJetSequence+
    process.icPFSequence+
    process.icPhotonSequence+
    process.icGenSequence+
    process.icL1EGammaProducer+
    process.icL1TauProducer+
    process.icL1MuonProducer+
    process.icTriggerSequence+
    process.icTriggerObjectSequence+
    process.icTauSpinnerSequence+
    process.icHtxsSequence+
    process.icEventInfoSequence+
    process.icEventProducer
)


process.schedule = cms.Schedule(process.p)

