import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")
process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.INFO.limit = 100000
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
    #'file:/tmp/aidan/SinglePhoton_700_GRR70V1-v1_101_2_tec.root'
    #'file:/tmp/aidan/0C2732B7-C4C5-E311-A8E5-00266CFFA120.root'
    'file:/tmp/aidan/Spring14miniaod__ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__MINIAODSIM__PU20bx25_POSTLS170_V5-v1__403D9883-C407-E411-80FA-003048CF6050.root'
    #'file:/tmp/aidan/ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__00000__067C0233-6ED1-E311-8C07-0025902008EC.root'
    #'/store/user/sharper/SinglePhoton/Run2012D_SinglePhoSkim_700_GRR70V1-v1_RECO/99931f3705ad96cae4b9189e0f7e85d6/SinglePhoton_700_GRR70V1-v1_100_2_TCZ.root'
])

# PFMET Type 1 (JEC) correction
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

# select global tag and pfJetMET correction automatically from datasetpath when using multicrab
noDataset = True
import sys
for arg in sys.argv:
    if arg.startswith("-CMSSW.datasetpath"):
        dataset = arg[19:]
        noDataset =  False
        break
if noDataset:
    dataset = 'none'

process.GlobalTag.globaltag = 'GR_R_70_V1::All'
print "Global Tag is ", process.GlobalTag.globaltag

process.out = cms.OutputModule("PoolOutputModule",
    ##process.FEVTSIMEventContent,
    fileName = cms.untracked.string('outfile.root')
)

process.options = cms.untracked.PSet(
    #fileMode = cms.untracked.string('NOMERGE')
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
    #input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('outfile.root')
)

process.hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

## # Primary vertex filter and no scraping events
## # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCollisionsDataAnalysis
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32 (4),
                                           maxAbsZ = cms.double (24),
                                           maxd0 = cms.double (2)
                                           )
process.primaryVertexPath = cms.Path(process.primaryVertexFilter)

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

## from RecoJets.Configuration.CaloTowersRec_cff import *
## from RecoJets.JetProducers.CaloJetParameters_cfi import *
## from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

## The next three lines are for rho computation (energy density, highly correlated to PU), see here :
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRecipesFor2011#FastJet_based_pile_up_isolation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gedGsfElectrons')
#process.muIsoSequence = setupPFMuonIso(process, 'muons')

process.load("UserCode.HEEPSkims.gsfcheckertree_cfi")
process.gsfcheckerjob.electron1EtMin = 0.
process.gsfcheckerjob.electron2EtMin = 0.
process.gsfcheckerjob.muonPtMin = 0.
process.otherStuff = cms.Sequence( process.kt6PFJets )

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.MessageLogger.suppressError = cms.untracked.vstring ('ecalLaserCorrFilter')

process.p1 = cms.Path(process.otherStuff * process.hltPhysicsDeclared * process.eeBadScFilter * process.ecalLaserCorrFilter * process.noscraping * process.primaryVertexFilter * process.pfParticleSelectionSequence * process.eleIsoSequence * process.producePFMETCorrections * process.HEEPTree)

