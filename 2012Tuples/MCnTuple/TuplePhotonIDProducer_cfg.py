import FWCore.ParameterSet.Config as cms

process = cms.Process("PhotonIDProc")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
#process.load("MagneticField.Engine.volumeBasedMagneticField_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("RecoEgamma.EgammaPhotonProducers.photonSequence_cff")
process.load("HLTrigger.HLTfilters.triggerResultsFilter_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.load("HLTrigger.HLTanalyzers.HLTAnalyser_cfi")

#process.GlobalTag.globaltag = 'MC_53_V6::All' # this works with 536 produced MC
process.GlobalTag.globaltag = 'MC_53_V15::All'
# process.GlobalTag.globaltag = 'START53_V19PR::All'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# '/store/data/Run2012C/SinglePhoton/AOD/PromptReco-v1/000/198/913/E8AA10C8-78CE-E111-BBEC-BCAEC5329720.root')	
#'file://///uscms_data/d2/bdiamond/threephoton2012/MCntuple/CMSSW_5_3_6/src/photon_53Triple.root')
# '/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_10_1_khd.root')
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_10_1_khd.root',
#'file://///uscms_data/d2/bdiamond/threephoton2012/MC_AOD/CMSSW_5_3_10/src/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU.root')
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_10_1_BNZ.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_11_1_elG.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_12_1_7b4.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_13_1_viK.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_14_2_qK8.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_15_2_a3G.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_16_2_tzW.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_17_2_rqP.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_18_1_uRn.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_19_2_1cL.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_1_2_0IE.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_20_2_zG5.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_2_1_Qlx.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_3_2_tG6.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_4_1_Re3.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_5_1_aLW.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_6_1_dgG.root',
 '/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_7_1_Ybn.root',
'/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_8_2_Rpt.root',
'/store/user/bdiamond/threephoton2012/signalMC/05SepM250/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_9_2_znK.root')
#'/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_11_1_uOt.root')
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_12_1_JEb.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_13_1_WYV.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_14_1_3gA.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_15_1_R5a.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_16_1_r8Q.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_17_1_8Ca.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_18_1_wSA.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_19_1_LXk.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_1_2_dhs.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_20_1_5Ku.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_2_1_e5I.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_3_1_gdj.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_4_1_7wz.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_5_1_jnR.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_6_1_16G.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_7_1_6Mw.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_8_1_gHD.root',
#'dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/bdiamond/threephoton2012/signalMC/06AugM200mlm/Hadronizer_MgmMatchTuneZ2star_8TeV_madgraph_cff_py_GEN_FASTSIM_HLT_PU_9_1_SXH.root')
    #fileNames = cms.untracked.vstring('dcap://cmsdcap.hep.wisc.edu:22125/pnfs/hep.wisc.edu/store/user/mbanderson/PhotonJet20-200/PhotonJet20-200-0000.root')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.HLTSel = cms.EDFilter( "TriggerResultsFilter",
                               triggerConditions = cms.vstring(
                                                  'HLT_Photon30_*'),
                               hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),
                               l1tResults = cms.InputTag( "gtDigis" ),
                               l1tIgnoreMask = cms.bool( False ),
                               l1techIgnorePrescales = cms.bool( False ),
                               daqPartitions = cms.uint32( 1 ),
                               throw = cms.bool( False )
)
# process.Out = cms.OutputModule("PoolOutputModule",
#     outputCommands = cms.untracked.vstring('drop *', 
#         'keep edmHepMCProduct_*_*_*', 
#         'keep recoBasicClusters_*_*_*', 
#         'keep recoSuperClusters_*_*_*', 
#         'keep *_PhotonIDProd_*_*', 
#         'keep *_PhotonIDProd_*_*', 
#         'keep recoPhotons_*_*_*'),
#     fileName = cms.untracked.string('Photest.root')
# )

process.photonIDAna = cms.EDAnalyzer("PhotonIDSimpleAnalyzer",
    outputFile  = cms.string('PhoIDHistsMC250full.root'),

    # Variables that must be passed before a photon candidate (a SuperCluster)
    #  gets placed into the histograms.  Basic, simple cuts.
    # Minimum Et
    minPhotonEt     = cms.double(10.0),
    # Minimum and max abs(eta)
    minPhotonAbsEta = cms.double(0.0),
    maxPhotonAbsEta = cms.double(3.0),
    # Minimum R9 = E(3x3) / E(SuperCluster)
    minPhotonR9     = cms.double(0.0),
    # Maximum HCAL / ECAL Energy
    maxPhotonHoverE = cms.double(0.2),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),#check what you need here HLT or REDIGI3..
    triggerEventTag  = cms.InputTag("hltTriggerSummaryAOD","","HLT"),#check what you need here HLT or REDIGI3..
#    trigEventTag = cms.InputTag("hltTriggerSummaryAOD","","HLT"),  
    hltlabel          = cms.string("HLT"),   #check what you need here HLT or REDIGI3..
    # Optionally produce a TTree of photons (set to False or True).
    # This slows down the analyzer, and if running
    # over 100,000+ events, this can create a large ROOT file
    createPhotonTTree  = cms.bool(False)
)

#process.p = cms.Path(process.photonSequence*process.photonIDSequence*process.photonIDAna)
process.p = cms.Path(process.HLTSel*process.photonIDAna)
# process.e = cms.EndPath(process.Out)

process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

# process.PoolSource.fileNames = [
#     '/store/relval/CMSSW_3_2_5/RelValZEE/GEN-SIM-RECO/MC_31X_V5-v1/0011/84E3AE4D-738E-DE11-A5E4-003048D374F2.root',
# '/store/relval/CMSSW_3_2_5/RelValZEE/GEN-SIM-RECO/MC_31X_V5-v1/0011/7886FE89-738E-DE11-B893-001D09F241F0.root',
# '/store/relval/CMSSW_3_2_5/RelValZEE/GEN-SIM-RECO/MC_31X_V5-v1/0011/60D71215-828E-DE11-8BFE-000423D98804.root',
# '/store/relval/CMSSW_3_2_5/RelValZEE/GEN-SIM-RECO/MC_31X_V5-v1/0011/4AA439F6-728E-DE11-AF2D-000423D98EA8.root'
# ] 
