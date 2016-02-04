import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("MssmHbb")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(-1)

##  Using MINIAOD. GlobalTag just in case jet re-clustering, L1 trigger filter  etc is needed to be done
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
#process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '74X_dataRun2_Prompt_v1')
process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = '74X_dataRun2_v5')
process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
for pset in process.GlobalTag.toGet.value():
    pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
## fix for multi-run processing
process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

output_file = 'output_new_jec.root'
## TFileService
process.TFileService = cms.Service("TFileService",
	fileName = cms.string(output_file)
)

## ============ TRIGGER FILTER ===============
## Enable below at cms.Path if needed (not recommended for Monte Carlo samples)
process.triggerSelection = cms.EDFilter( "TriggerResultsFilter",
    triggerConditions = cms.vstring(
                          "HLT_DoubleJetsC100_DoubleBTagCSV0p9_DoublePFJetsC100MaxDeta1p6_v*",
                          "HLT_DoubleJetsC100_DoubleBTagCSV0p85_DoublePFJetsC160_v*",
                          "HLT_DoubleJetsC112_DoubleBTagCSV0p9_DoublePFJetsC112MaxDeta1p6_v*",
                          "HLT_DoubleJetsC112_DoubleBTagCSV0p85_DoublePFJetsC172_v*",
                          "HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v*",
                          "HLT_QuadJet45_DoubleBTagCSV0p67_v*",
                          "HLT_QuadJet45_TripleBTagCSV0p67_v*",
                          "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v*",
                          "HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v*",
                          "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v*",
                          "HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v*",
                          
    ),
    hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
    l1tResults = cms.InputTag( "" ),
    l1tIgnoreMask = cms.bool( False ),
    l1techIgnorePrescales = cms.bool( False ),
    daqPartitions = cms.uint32( 1 ),
    throw = cms.bool( True )
)

## ============ EVENT FILTER COUNTER ===============
## Filter counter (maybe more useful for MC)
process.TotalEvents    = cms.EDProducer("EventCountProducer")
process.FilteredEvents = cms.EDProducer("EventCountProducer")

## ============ PRIMARY VERTEX FILTER ===============
process.primaryVertexFilter = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlineSlimmedPrimaryVertices"), # primary vertex collection name
   cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"), # ndof>thr=4 corresponds to sum(track_weigths) > (thr+3)/2 = 3.5 so typically 4 good tracks
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

## ============ NEW JEC ================

# add .txt files with latest JEC
# Not sure that AK4PFchs should be used. Maybe AK4PFPuppi ???
L1JEC = cms.string("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/Ntuplizer/test/Summer15_25nsV7_DATA_L1FastJet_AK4PFchs.txt");
L2JEC = cms.string("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/Ntuplizer/test/Summer15_25nsV7_DATA_L2Relative_AK4PFchs.txt");
L3JEC = cms.string("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/Ntuplizer/test/Summer15_25nsV7_DATA_L3Absolute_AK4PFchs.txt");
ResJetJEC = cms.string("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/Ntuplizer/test/Summer15_25nsV7_DATA_L2L3Residual_AK4PFchs.txt");

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
                                                                             src = cms.InputTag("slimmedJets"),
                                                                             levels = ['L1FastJet',
                                                                                       'L2Relative',
                                                                                       'L3Absolute'],
                                                                             payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.patJetsReapplyJEC = patJetsUpdated.clone(
                                                         jetSource = cms.InputTag("slimmedJets"),
                                                         jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
                                                         )

#process.p = cms.Sequence( process.patJetCorrFactorsReapplyJEC + process. patJetsReapplyJEC )

## ============  THE NTUPLIZER!!!  ===============
process.MssmHbb     = cms.EDAnalyzer("Ntuplizer",
    MonteCarlo      = cms.bool(False),
    UseFullName     = cms.bool(False),
    TotalEvents     = cms.InputTag("TotalEvents"),
    FilteredEvents  = cms.InputTag("FilteredEvents"),
    PatJets         = cms.VInputTag(
                                    cms.InputTag("slimmedJetsPuppi","","PAT"),
                                    cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked","SubJets","PAT")
                                    ), 
    PatMETs         = cms.VInputTag(
                                    cms.InputTag("slimmedMETs","","PAT"),
                                    cms.InputTag("slimmedMETsPuppi","","PAT")
                                    ), 
    PatMuons        = cms.VInputTag(
                                    cms.InputTag("slimmedMuons","","PAT")
                                    ), 
    PrimaryVertices = cms.VInputTag(
                                    cms.InputTag("offlineSlimmedPrimaryVertices","","PAT")
                                    ), 
    BTagAlgorithms = cms.vstring   (
                                    "pfCombinedInclusiveSecondaryVertexV2BJetTags",
                                    "combinedSecondaryVertexBJetTags",
                                    "pfJetBProbabilityBJetTags",
                                    "pfJetProbabilityBJetTags",
                                    "pfTrackCountingHighPurBJetTags",
                                    "pfTrackCountingHighEffBJetTags",
                                    "pfSimpleSecondaryVertexHighEffBJetTags",
                                    "pfSimpleSecondaryVertexHighPurBJetTags",
                                    "pfCombinedSecondaryVertexV2BJetTags",
                                    "pfCombinedSecondaryVertexSoftLeptonBJetTags",
                                    "pfCombinedMVABJetTags",
                                   ),
    BTagAlgorithmsAlias = cms.vstring   (
                                         "btag_csvivf",
                                         "btag_csv",
                                         "btag_jetbprob",
                                         "btag_jetprob",
                                         "btag_tchp",
                                         "btag_tche",
                                         "btag_svhe",
                                         "btag_svhp",
                                         "btag_csvv2",
                                         "btag_csvlep",
                                         "btag_csvmva",
                                        ),
    TriggerResults  = cms.VInputTag(cms.InputTag("TriggerResults","","HLT")),
    TriggerPaths    = cms.vstring  (
    ## I recommend using the version number explicitly to be able to compare 
    ## however for production one has to be careful that all versions are included.
    ## Thinking of a better solution...
    											"HLT_DoubleJetsC100_DoubleBTagCSV0p9_DoublePFJetsC100MaxDeta1p6_v",
    											"HLT_DoubleJetsC100_DoubleBTagCSV0p85_DoublePFJetsC160_v",
    											"HLT_DoubleJetsC112_DoubleBTagCSV0p9_DoublePFJetsC112MaxDeta1p6_v",
    											"HLT_DoubleJetsC112_DoubleBTagCSV0p85_DoublePFJetsC172_v",
                                    #                                   "HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v",
                                    #"HLT_QuadJet45_DoubleBTagCSV0p67_v",
                                    #"HLT_QuadJet45_TripleBTagCSV0p67_v",
                                    #"HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq200_v",
                                    #"HLT_QuadPFJet_DoubleBTagCSV_VBF_Mqq240_v",
                                    #"HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq460_v",
                                    #"HLT_QuadPFJet_SingleBTagCSV_VBF_Mqq500_v",
                                   ),
    TriggerObjectStandAlone  = cms.VInputTag(
                                             cms.InputTag("selectedPatTrigger","","PAT"),
                                             ),
    TriggerObjectLabels    = cms.vstring  (
    											"hltL1sL1DoubleJetC100",
    											"hltDoubleJetsC100",
    											"hltDoublePFJetsC100",
    											"hltDoublePFJetsC100MaxDeta1p6",
    											"hltDoublePFJetsC160",
    											"hltDoubleBTagCSV0p85",
    											"hltDoubleBTagCSV0p9",
    											"hltL1sL1DoubleJetC112",
    											"hltDoubleJetsC112",
    											"hltDoublePFJetsC112",
    											"hltDoublePFJetsC112MaxDeta1p6",
    											"hltDoublePFJetsC172",
                                   ),
#    L1ExtraJets     = cms.VInputTag(
#                                    cms.InputTag("l1extraParticles","Central","RECO"),
#                                    cms.InputTag("l1extraParticles","Forward","RECO"),
#                                    cms.InputTag("l1extraParticles","Tau","RECO")
#                                    ),
#    L1ExtraMuons    = cms.VInputTag(
#                                    cms.InputTag("l1extraParticles","","RECO")
#                                    ),
)

process.p = cms.Path(
                      process.TotalEvents *
                      process.primaryVertexFilter *
                      process.triggerSelection * 
                      process.FilteredEvents *
                      (process.patJetCorrFactorsReapplyJEC +
                      process.patJetsReapplyJEC) *
                      process.MssmHbb
                    )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [


       '/store/data/Run2015D/BTagCSV/MINIAOD/05Oct2015-v1/30000/1A46877D-766F-E511-B79A-0025905B8590.root'
] );


secFiles.extend( [
               ] )

