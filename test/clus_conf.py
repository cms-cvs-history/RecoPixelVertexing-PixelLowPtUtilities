import FWCore.ParameterSet.Config as cms

process = cms.Process("ClusterShape")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.Simulation_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load("SimGeneral.TrackingAnalysis.trackingParticles_cfi")
process.load("SimTracker.TrackAssociation.TrackAssociatorByHits_cfi")

process.load("RecoLocalTracker.Configuration.RecoLocalTracker_cff")

process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("RecoPixelVertexing.PixelLowPtUtilities.MinBiasTracking_cff")

###############################################################################
# Categories and modules
process.CategoriesAndModules = cms.PSet(
    categories = cms.untracked.vstring('MinBiasTracking'),
    debugModules = cms.untracked.vstring('*')
)

###############################################################################
# Message logger
process.MessageLogger = cms.Service("MessageLogger",
     process.CategoriesAndModules,
     cerr = cms.untracked.PSet(
         threshold = cms.untracked.string('DEBUG'),
#        threshold = cms.untracked.string('ERROR'),
         DEBUG = cms.untracked.PSet(
             limit = cms.untracked.int32(0)
         )
     ),
     destinations = cms.untracked.vstring('cerr')
)

###############################################################################
# Source
process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
       # /RelValMinBias/CMSSW_3_1_0_pre2_IDEAL_30X_v1/GEN-SIM-DIGI-RAW-HLTDEBUG
       '/store/relval/CMSSW_3_1_0_pre2/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0000/2470A74A-4103-DE11-9008-0030487A18F2.root',
       '/store/relval/CMSSW_3_1_0_pre2/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0000/AAC15191-4103-DE11-8071-0016177CA778.root',
       '/store/relval/CMSSW_3_1_0_pre2/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0000/ACD6837F-4103-DE11-BB54-0030487A18A4.root',
       '/store/relval/CMSSW_3_1_0_pre2/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0000/B6DE97A8-4203-DE11-A9E6-001D09F276CF.root',
       '/store/relval/CMSSW_3_1_0_pre2/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/IDEAL_30X_v1/0001/C207D770-9703-DE11-B4F7-001617DBD5B2.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

###############################################################################
# Cluster shape
process.clusterShape = cms.EDAnalyzer("ClusterShapeExtractor",
    trackProducer  = cms.string('allTracks'),
    hasSimHits     = cms.bool(True),
    hasRecTracks   = cms.bool(False),
    associateStrip      = cms.bool(True),
    associatePixel      = cms.bool(True),
    associateRecoTracks = cms.bool(False),
    ROUList = cms.vstring(
      'g4SimHitsTrackerHitsTIBLowTof', 'g4SimHitsTrackerHitsTIBHighTof',
      'g4SimHitsTrackerHitsTIDLowTof', 'g4SimHitsTrackerHitsTIDHighTof',
      'g4SimHitsTrackerHitsTOBLowTof', 'g4SimHitsTrackerHitsTOBHighTof',
      'g4SimHitsTrackerHitsTECLowTof', 'g4SimHitsTrackerHitsTECHighTof',
      'g4SimHitsTrackerHitsPixelBarrelLowTof',
      'g4SimHitsTrackerHitsPixelBarrelHighTof',
      'g4SimHitsTrackerHitsPixelEndcapLowTof',
      'g4SimHitsTrackerHitsPixelEndcapHighTof')
)

###############################################################################
# Paths
process.simu  = cms.Path(process.mix
                       * process.trackingParticles
                       * process.offlineBeamSpot)
process.digi  = cms.Path(process.trDigi)
process.lreco = cms.Path(process.trackerlocalreco
                       * process.clusterShape)

###############################################################################
# Global tag
process.GlobalTag.globaltag = 'IDEAL_30X::All'

###############################################################################
# Workaround
process.siPixelClusters.src = 'simSiPixelDigis'
process.siStripZeroSuppression.RawDigiProducersList = cms.VPSet(cms.PSet(
        RawDigiProducer = cms.string('simSiStripDigis'),
        RawDigiLabel = cms.string('VirginRaw')
    ))
process.siStripClusters.DigiProducersList = cms.VPSet(cms.PSet(
        DigiLabel = cms.string('ZeroSuppressed'),
        DigiProducer = cms.string('simSiStripDigis')
    ))

###############################################################################
# Schedule
process.schedule = cms.Schedule(process.simu,
                                process.digi,
                                process.lreco)
