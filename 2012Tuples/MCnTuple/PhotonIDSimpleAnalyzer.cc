// -*- C++ -*-
//
// Package:    PhotonIDSimpleAnalyzer
// Class:      PhotonIDSimpleAnalyzer
// 
/**\class PhotonIDSimpleAnalyzer PhotonIDSimpleAnalyzer.cc RecoEgamma/PhotonIdentification/test/PhotonIDSimpleAnalyzer.cc

 Description: Generate various histograms for cuts and important 
              photon ID parameters using a data sample of photons in QCD events.

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  J. Stilley, A. Askew
//  Editing Author:  M.B. Anderson
//
//         Created:  Fri May 9 11:03:51 CDT 2008
// $Id: PhotonIDSimpleAnalyzer.cc,v 1.10 2010/01/14 17:29:32 nancy Exp $
//
///////////////////////////////////////////////////////////////////////
//                    header file for this analyzer                  //
///////////////////////////////////////////////////////////////////////
#include "RecoEgamma/PhotonIdentification/plugins/PhotonIDSimpleAnalyzer.h"
#include "RecoEgamma/PhotonIdentification/plugins/TrigToolsFuncs.h"

///////////////////////////////////////////////////////////////////////
//                        CMSSW includes                             //
///////////////////////////////////////////////////////////////////////
//#include "DataFormats/EgammaCandidates/interface/PhotonIDFwd.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonID.h"
//#include "DataFormats/EgammaCandidates/interface/PhotonIDAssociation.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
// #include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

///////////////////////////////////////////////////////////////////////
//                      Root include files                           //
///////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
//                           Constructor                             //
///////////////////////////////////////////////////////////////////////
PhotonIDSimpleAnalyzer::PhotonIDSimpleAnalyzer(const edm::ParameterSet& ps)
{
  hlTriggerResults_=ps.getParameter<edm::InputTag>("HLTriggerResults");
  hltlabel_ =ps.getParameter<string>("hltlabel");
  triggerEventTag_ = ps.getParameter<edm::InputTag>("triggerEventTag");
  // Read Parameters from configuration file
  // output filename
  outputFile_   = ps.getParameter<std::string>("outputFile");
  // Read variables that must be passed to allow a 
  //  supercluster to be placed in histograms as a photon.
  minPhotonEt_     = ps.getParameter<double>("minPhotonEt");
  minPhotonAbsEta_ = ps.getParameter<double>("minPhotonAbsEta");
  maxPhotonAbsEta_ = ps.getParameter<double>("maxPhotonAbsEta");
  minPhotonR9_     = ps.getParameter<double>("minPhotonR9");
  maxPhotonHoverE_ = ps.getParameter<double>("maxPhotonHoverE");

// testing ParticleFlow Isolation for photons
 isolator.initializePhotonIsolation(kTRUE);
 isolator. setConeSize(0.3);


  // Read variable to that decidedes whether
  // a TTree of photons is created or not
  createPhotonTTree_ = ps.getParameter<bool>("createPhotonTTree");

  // open output file to store histograms
  rootFile_ = TFile::Open(outputFile_.c_str(),"RECREATE");
}

///////////////////////////////////////////////////////////////////////
//                            Destructor                             //
///////////////////////////////////////////////////////////////////////
PhotonIDSimpleAnalyzer::~PhotonIDSimpleAnalyzer()
{

// do anything here that needs to be done at desctruction time
// (e.g. close files, deallocate resources etc.)

  delete rootFile_;

}

///////////////////////////////////////////////////////////////////////
//    method called once each job just before starting event loop    //
///////////////////////////////////////////////////////////////////////
void 
PhotonIDSimpleAnalyzer::beginJob()
{
	  triggernames     = &all_triggers; 
	  triggerprescales = &all_triggerprescales; 
	  ifTriggerpassed  = &all_ifTriggerpassed;
  // go to *OUR* rootfile
  rootFile_->cd();


  PhotonIDTree = new TTree("PhotonIDTree","Photon Tree");
  PhotonIDTree->Branch("Run", &Run, "Run/I");
  PhotonIDTree->Branch("Event", &Event, "Event/I");
  PhotonIDTree->Branch("LumiSec", &LumiSec,"LumiSec/I");
  PhotonIDTree->Branch("nVert", &nVert, "nVert/I");  
  PhotonIDTree->Branch("vx", vx,"vx[nVert]/D");
  PhotonIDTree->Branch("vy", vy,"vy[nVert]/D");
  PhotonIDTree->Branch("vz", vz,"vz[nVert]/D");
  PhotonIDTree->Branch("gen_nVert", &gen_nVert, "gen_nVert/I");  
  PhotonIDTree->Branch("gen_vx", gen_vx,"gen_vx[gen_nVert]/D");
  PhotonIDTree->Branch("gen_vy", gen_vy,"gen_vy[gen_nVert]/D");
  PhotonIDTree->Branch("gen_vz", gen_vz,"gen_vz[gen_nVert]/D");
  PhotonIDTree->Branch("PU_NumInteractions", &PU_NumInteractions, "PU_NumInteractions/I");  

  
  PhotonIDTree->Branch("nPho", &nPho, "nPho/I");
  PhotonIDTree->Branch("rho", &rho,"rho/F");
//  PhotonIDTree->Branch("photoNcrys", photonNcrys,"photonNcrys[nPho]/I");
//  PhotonIDTree->Branch("photon_swissCross", photon_swissCross, "photon_swissCross[nPho]/F");
//  PhotonIDTree->Branch("photon_e2e9", photon_e2e9,"photon_e2e9[nPho]/F");
//  PhotonIDTree->Branch("photon_e6e2", photon_e6e2,"photon_e6e2[nPho]/F");
  PhotonIDTree->Branch("photonSCeta", photonSCeta, "photonSCeta[nPho]/F");
  PhotonIDTree->Branch("photonSCphi", photonSCphi, "photonSCphi[nPho]/F");
  PhotonIDTree->Branch("photonSCX", photonSCX, "photonSCX[nPho]/F");
  PhotonIDTree->Branch("photonSCY", photonSCY, "photonSCY[nPho]/F");
  PhotonIDTree->Branch("photonSCZ", photonSCZ, "photonSCZ[nPho]/F");
  PhotonIDTree->Branch("photonSCE", photonSCE, "photonSCE[nPho]/F");
  PhotonIDTree->Branch("photonSCetawidth", photonSCetawidth,"photonSCetawidth[nPho]/F");
  PhotonIDTree->Branch("photonSCphiwidth", photonSCphiwidth,"photonSCphiwidth[nPho]/F");
  PhotonIDTree->Branch("photonet", photonet, "photonet[nPho]/F");
  PhotonIDTree->Branch("photon_physeta", photon_physeta,"photon_physeta[nPho]/F");
  PhotonIDTree->Branch("photon_physphi", photon_physphi,"photon_physphi[nPho]/F");
  PhotonIDTree->Branch("photone", photone,"photone[nPho]/F");
  PhotonIDTree->Branch("photonhadTowOverEm", photonhadTowOverEm,"photonhadTowOverEm[nPho]/F"); 
  PhotonIDTree->Branch("photonsigmaIetaIeta", photonsigmaIetaIeta,"photonsigmaIetaIeta[nPho]/F");
  PhotonIDTree->Branch("photonchargedHadronIso", photonchargedHadronIso,"photonchargedHadronIso[nPho]/F");
  PhotonIDTree->Branch("photonneutralHadronIso", photonneutralHadronIso,"photonneutralHadronIso[nPho]/F");
  PhotonIDTree->Branch("photonphotonIso", photonphotonIso,"photonphotonIso[nPho]/F");
  PhotonIDTree->Branch("passelectronveto", passelectronveto,"passelectronveto[nPho]/O");
  PhotonIDTree->Branch("r9", r9,"r9[nPho]/F");
  PhotonIDTree->Branch("ntriggers",&ntriggers,"ntriggers/I");
  PhotonIDTree->Branch("triggernames","vector<std::string>",&triggernames);
  PhotonIDTree->Branch("triggerprescales","vector<int>",&triggerprescales);
  PhotonIDTree->Branch("ifTriggerpassed","vector<bool>",&ifTriggerpassed);
  

    
  // Book Histograms
  // PhotonID Histograms
  h_isoEcalRecHit_ = new TH1F("photonEcalIso",          "Ecal Rec Hit Isolation", 300, 0, 300);
  h_isoHcalRecHit_ = new TH1F("photonHcalIso",          "Hcal Rec Hit Isolation", 300, 0, 300);
  h_trk_pt_solid_  = new TH1F("photonTrackSolidIso",    "Sum of track pT in a cone of #DeltaR" , 300, 0, 300);
  h_trk_pt_hollow_ = new TH1F("photonTrackHollowIso",   "Sum of track pT in a hollow cone" ,     300, 0, 300);
  h_ntrk_solid_    = new TH1F("photonTrackCountSolid",  "Number of tracks in a cone of #DeltaR", 100, 0, 100);
  h_ntrk_hollow_   = new TH1F("photonTrackCountHollow", "Number of tracks in a hollow cone",     100, 0, 100);
  h_ebetagap_         = new TH1F("photonInEBEtagap",          "Ecal Barrel eta gap flag",  2, -0.5, 1.5);
  h_ebphigap_         = new TH1F("photonInEBEtagap",          "Ecal Barrel phi gap flag",  2, -0.5, 1.5);
  h_eeringGap_         = new TH1F("photonInEERinggap",          "Ecal Endcap ring gap flag",  2, -0.5, 1.5);
  h_eedeeGap_         = new TH1F("photonInEEDeegap",          "Ecal Endcap dee gap flag",  2, -0.5, 1.5);
  h_ebeeGap_       = new TH1F("photonInEEgap",          "Ecal Barrel/Endcap gap flag",  2, -0.5, 1.5);
  h_r9_            = new TH1F("photonR9",               "R9 = E(3x3) / E(SuperCluster)", 300, 0, 3);

  // Photon Histograms
  h_photonEt_      = new TH1F("photonEt",     "Photon E_{T}",  200,  0, 200);
  h_photonEta_     = new TH1F("photonEta",    "Photon #eta",   800, -4,   4);
  h_photonPhi_     = new TH1F("photonPhi",    "Photon #phi",   628, -1.*TMath::Pi(), TMath::Pi());
  h_hadoverem_     = new TH1F("photonHoverE", "Hadronic over EM", 200, 0, 1);

  // Photon's SuperCluster Histograms
  h_photonScEt_       = new TH1F("photonScEt",  "Photon SuperCluster E_{T}", 200,  0, 200);
  h_photonScEta_      = new TH1F("photonScEta", "Photon #eta",               800, -4,   4);
  h_photonScPhi_      = new TH1F("photonScPhi", "Photon #phi",628, -1.*TMath::Pi(), TMath::Pi());
  h_photonScEtaWidth_ = new TH1F("photonScEtaWidth","#eta-width",            100,  0,  .1);

  // Composite or Other Histograms
  h_photonInAnyGap_   = new TH1F("photonInAnyGap",     "Photon in any gap flag",  2, -0.5, 1.5);
  h_nPassingPho_      = new TH1F("photonLoosePhoton", "Total number photons (0=NotPassing, 1=Passing)", 2, -0.5, 1.5);
  h_nPassEM_          = new TH1F("photonLooseEM", "Total number photons (0=NotPassing, 1=Passing)", 2, -0.5, 1.5);
  h_nPho_             = new TH1F("photonCount",        "Number of photons passing cuts in event",  10,  0,  10);

  // Create a TTree of photons if set to 'True' in config file
  if ( createPhotonTTree_ ) {
    tree_PhotonAll_     = new TTree("TreePhotonAll", "Reconstructed Photon");
    tree_PhotonAll_->Branch("recPhoton", &recPhoton.isolationEcalRecHit, "isolationEcalRecHit/F:isolationHcalRecHit:isolationSolidTrkCone:isolationHollowTrkCone:nTrkSolidCone:nTrkHollowCone:isEBGap:isEEGap:isEBEEGap:r9:et:eta:phi:hadronicOverEm");
  }
}

///////////////////////////////////////////////////////////////////////
//                method called to for each event                    //
///////////////////////////////////////////////////////////////////////
void
PhotonIDSimpleAnalyzer::analyze(const edm::Event& evt, const edm::EventSetup& es)
{

  Run = evt.id().run();
  Event = evt.id().event();
  LumiSec = evt.id().luminosityBlock();
  using namespace std;
  using namespace edm;

// Grab Generator Info
  edm::Handle< GenEventInfoProduct > GenInfoHandle;
  evt.getByLabel( "generator", GenInfoHandle );
// Below is an example of getting GenInfo
//   double qScale = GenInfoHandle->qScale();
//   double pthat = ( GenInfoHandle->hasBinningValues() ? 
//                   (GenInfoHandle->binningValues())[0] : 0.0);
//  cout << " qScale = " << qScale << " pthat = " << pthat << endl;
  edm::Handle< HepMCProduct > EvtHandle ;
  evt.getByLabel( "generator", EvtHandle ) ;
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
//   cout << endl << " ########################## " << endl;
//   cout << " Gen event = " << Evt->event_number() << " Event ID = " << Event << endl;

//setup for conversion safe electron veto
   edm::Handle<reco::BeamSpot> bsHandle;
   evt.getByLabel("offlineBeamSpot", bsHandle);
   const reco::BeamSpot &beamspot = *bsHandle.product();

   edm::Handle<reco::ConversionCollection> hConversions;
   evt.getByLabel("allConversions", hConversions);

   edm::Handle<reco::GsfElectronCollection> hElectrons;
   evt.getByLabel("gsfElectrons", hElectrons);


  // grab PF Collection
  Handle<reco::PFCandidateCollection> pfParticlesColl;
  evt.getByLabel("particleFlow",pfParticlesColl);
    
  // grab PF jet rho
  edm::Handle<double> rH;
  evt.getByLabel("kt6PFJets","rho",rH);
  rho = *rH;
  
  // grab Vertex Collection
  edm::Handle<reco::VertexCollection> Vertices;
//  const reco::Vertex* primVtx = 0;
//  try {
//  evt.getByLabel(vtxCollectionTag_,Vertices);
  evt.getByLabel("offlinePrimaryVerticesWithBS",Vertices);
   nVert=0;
   reco::VertexCollection::const_iterator it;
   for(it = (*Vertices).begin(); it != (*Vertices).end(); it++){
//       if(!primVtx) primVtx = &*it;
    	vx[nVert] = it->x();
    	vy[nVert] = it->y();
    	vz[nVert] = it->z();
       if (nVert==0)
       {
//       std::cout << "event : " << Event << "  vertex : " << it->x() << ", " << it->y() << ", " << it->z() << std::endl;
//	   std::cout << "nVert : " << nVert << endl;

	   }
    nVert++;
    }
//std::cout << "nVert = " << nVert <<endl;
 HepMC::GenVertex* NHyperPionDecVtx = 0 ;
   gen_nVert=0;
//   currentCh_barcode = 0;
   current_barcode = 0;
 for ( HepMC::GenEvent::vertex_const_iterator
          vit=Evt->vertices_begin(); vit!=Evt->vertices_end(); vit++ )
  {

  for ( HepMC::GenVertex::particles_out_const_iterator
              pout=(*vit)->particles_out_const_begin();
            pout!=(*vit)->particles_out_const_end(); pout++ )
      {
          if ( (*pout)->pdg_id() == 54) // && (*pout)->status() == 2 ) 
          {	    
              if ( (*pout)->end_vertex() != 0 )
              {
                  NHyperPionDecVtx = (*pout)->end_vertex() ;
//					(*pout)->print() ;
                  if (current_barcode != NHyperPionDecVtx->barcode())
	              {
  	 				gen_vx[gen_nVert] = (*vit)->position().x()/10.0;
    				gen_vy[gen_nVert] = (*vit)->position().y()/10.0;
    				gen_vz[gen_nVert] = (*vit)->position().z()/10.0;
//  	 				std::cout << "gen_event : " << Evt->event_number() << "  gen_vertex : " << (*vit)->position().x()/10.0 << ", " << (*vit)->position().y()/10.0 << ", " << (*vit)->position().z()/10.0 << std::endl;
//  	 				std::cout << "gen_event : " << Evt->event_number() << "  gen_vertex : " << gen_vx[gen_nVert] << ", " << gen_vy[gen_nVert] << ", " << gen_vz[gen_nVert] << std::endl;
//					std::cout << "gen_nVert : " << gen_nVert << endl;
			  	    gen_nVert++;

   				  }
	  			  current_barcode = NHyperPionDecVtx->barcode();
                  break ;
              }
          }
//          if ( (*pout)->pdg_id() == 55) // && (*pout)->status() == 2 ) 
//           {	    
//               if ( (*pout)->end_vertex() != 0 )
//               {
//                   CHyperPionDecVtx = (*pout)->end_vertex() ;
// 					(*pout)->print() ;
//                   if (current_barcode != CHyperPionDecVtx->barcode())
// 	              {
//   	 				gen_vx[gen_nVert] = (*vit)->position().x()/10.0;
//     				gen_vy[gen_nVert] = (*vit)->position().y()/10.0;
//     				gen_vz[gen_nVert] = (*vit)->position().z()/10.0;
// //  	 				std::cout << "gen_event : " << Evt->event_number() << "  gen_vertex : " << (*vit)->position().x()/10.0 << ", " << (*vit)->position().y()/10.0 << ", " << (*vit)->position().z()/10.0 << std::endl;
// //  	 				std::cout << "gen_event : " << Evt->event_number() << "  gen_vertex : " << gen_vx[gen_nVert] << ", " << gen_vy[gen_nVert] << ", " << gen_vz[gen_nVert] << std::endl;
// //					std::cout << "gen_nVert : " << gen_nVert << endl;
// 			  	    gen_nVert++;
// 
//    				  }
// 	  			  current_barcode = NHyperPionDecVtx->barcode();
//                   break ;
//               }
//          }
      }
      if ( NHyperPionDecVtx != 0 )
      {
          break ; // break the initial loop over vertices
      }
      
//        std::cout << "gen_event : " << Evt->event_number() << "  gen_vertex : " << (*vit)->position().x() << ", " << (*vit)->position().y() << ", " << (*vit)->position().z() << std::endl;
//  	   gen_nVert++;
     }
     
  if ( NHyperPionDecVtx == 0 ) 
  {
      cout << " There is NO Neutral Hyper Pion in this event ! " << endl ;
      return ;
  }
  
  if ( NHyperPionDecVtx !=0 )
  {
//      cout << " " << endl ;
//      cout << " Neutral Hyper Pion decay found at the vertex " << NHyperPionDecVtx->barcode() <<" (barcode)" << endl ;

      vector<HepMC::GenParticle*> Children;

      for ( HepMC::GenVertex::particles_out_const_iterator NHPin = 
              NHyperPionDecVtx->particles_out_const_begin(); 
	      NHPin != NHyperPionDecVtx->particles_out_const_end(); 
	      NHPin++ ) 
      { 
         Children.push_back(*NHPin);
         
      }
//      cout << " Number of Neutral Hyper Pion (immediate) children = " << Children.size() << endl ;
      for (unsigned int ic=0; ic<Children.size(); ic++ )
      {
//          Children[ic]->print() ;   
      }
  }
// std::cout << "start "  <<endl;

// select and store stable descendants of the Hyperpion
//   
vector<HepMC::GenParticle*> StableNHPDes ;
for ( HepMC::GenVertex::particle_iterator
         des=NHyperPionDecVtx->particles_begin(HepMC::descendants);
	 des!=NHyperPionDecVtx->particles_end(HepMC::descendants); des++ )
   {
//      if ( (*des)->status() == 1 ) (*des)->print() ;
      if ( (*des)->status() == 1 ) StableNHPDes.push_back(*des) ;
   }
	


  
//   }
//   catch(...) {
//     std::cout << "VertexCollection is not available!!! " <<  std::endl;
//   }
  // grab photons
  Handle<reco::PhotonCollection> photonColl;
  evt.getByLabel("photons", "", photonColl);

  Handle<edm::ValueMap<Bool_t> > loosePhotonQual;
  evt.getByLabel("PhotonIDProd", "PhotonCutBasedIDLoose", loosePhotonQual);
  Handle<edm::ValueMap<Bool_t> > looseEMQual;
  evt.getByLabel("PhotonIDProd","PhotonCutBasedIDLooseEM",looseEMQual);
  // grab PhotonId objects  
//   Handle<reco::PhotonIDAssociationCollection> photonIDMapColl;
//   evt.getByLabel("PhotonIDProd", "PhotonAssociatedID", photonIDMapColl);
  
  // create reference to the object types we are interested in
  const reco::PhotonCollection *photons = photonColl.product();  
  const edm::ValueMap<Bool_t> *phoMap = loosePhotonQual.product();
  const edm::ValueMap<Bool_t> *lEMMap = looseEMQual.product();
  int photonCounter = 0;
  int idxpho=0;
  reco::PhotonCollection::const_iterator pho;

  nPho=0;
  for (pho = (*photons).begin(); pho!= (*photons).end() && nPho<100; pho++){      
    
    edm::Ref<reco::PhotonCollection> photonref(photonColl, idxpho);
    //reco::PhotonIDAssociationCollection::const_iterator photonIter = phoMap->find(photonref);
    //const reco::PhotonIDRef &phtn = photonIter->val;
    //const reco::PhotonRef &pho = photonIter->key;
	// std::cout << "event " << Event << " photon vertex : " << pho->vx() << ", " << pho->vy() << ", " << pho->vz()<< endl;  
    float photonEt       = pho->et();
    float superClusterEt = (pho->superCluster()->energy())/(cosh(pho->superCluster()->position().eta()));
    Bool_t LoosePhotonQu = (*phoMap)[photonref];
    h_nPassingPho_->Fill(LoosePhotonQu);
    Bool_t LooseEMQu = (*lEMMap)[photonref];
    h_nPassEM_->Fill(LooseEMQu);
    // Only store photon candidates (SuperClusters) that pass some simple cuts
    bool passCuts = (              photonEt > minPhotonEt_     ) &&
                    (      fabs(pho->eta()) > minPhotonAbsEta_ ) &&
                    (      fabs(pho->eta()) < maxPhotonAbsEta_ ) &&
                    (             pho->r9() > minPhotonR9_     ) &&
                    ( pho->hadronicOverEm() < maxPhotonHoverE_ ) ;

    passelectronveto[nPho] = !ConversionTools::hasMatchedPromptElectron(pho->superCluster(), hElectrons, hConversions, beamspot.position());

	
    if ( passCuts )
    {
      ///////////////////////////////////////////////////////
      //                fill histograms                    //
      ///////////////////////////////////////////////////////
      // PhotonID Variables
      h_isoEcalRecHit_->Fill(pho->ecalRecHitSumEtConeDR04());
      h_isoHcalRecHit_->Fill(pho->hcalTowerSumEtConeDR04());
      h_trk_pt_solid_ ->Fill(pho->trkSumPtSolidConeDR04());
      h_trk_pt_hollow_->Fill(pho->trkSumPtHollowConeDR04());
      h_ntrk_solid_->   Fill(pho->nTrkSolidConeDR04());
      h_ntrk_hollow_->  Fill(pho->nTrkHollowConeDR04());
      h_ebetagap_->        Fill(pho->isEBEtaGap());
      h_ebphigap_->        Fill(pho->isEBPhiGap());
      h_eeringGap_->        Fill(pho->isEERingGap()); 
      h_eedeeGap_->        Fill(pho->isEEDeeGap()); 
      h_ebeeGap_->      Fill(pho->isEBEEGap());
      h_r9_->           Fill(pho->r9());

      // Photon Variables
      h_photonEt_->  Fill(photonEt);
      h_photonEta_-> Fill(pho->eta());
      h_photonPhi_-> Fill(pho->phi());
      h_hadoverem_-> Fill(pho->hadronicOverEm());

      // Photon's SuperCluster Variables
      // eta is with respect to detector (not physics) vertex,
      // thus Et and eta are different from photon.
      h_photonScEt_->      Fill(superClusterEt);
      h_photonScEta_->     Fill(pho->superCluster()->position().eta());
      h_photonScPhi_->     Fill(pho->superCluster()->position().phi());
      h_photonScEtaWidth_->Fill(pho->superCluster()->etaWidth());
      
      // It passed photon cuts, mark it

      ///////////////////////////////////////////////////////
      //                fill TTree (optional)              //
      ///////////////////////////////////////////////////////
      if ( createPhotonTTree_ ) {
	recPhoton.isolationEcalRecHit    = pho->ecalRecHitSumEtConeDR04();
	recPhoton.isolationHcalRecHit    = pho->hcalTowerSumEtConeDR04();
	recPhoton.isolationSolidTrkCone  = pho->trkSumPtSolidConeDR04();
	recPhoton.isolationHollowTrkCone = pho->trkSumPtHollowConeDR04();
	recPhoton.nTrkSolidCone          = pho->nTrkSolidConeDR04();
	recPhoton.nTrkHollowCone         = pho->nTrkHollowConeDR04();
	recPhoton.isEBEtaGap                = pho->isEBEtaGap();
	recPhoton.isEBPhiGap                = pho->isEBPhiGap();
	recPhoton.isEERingGap                = pho->isEERingGap();
	recPhoton.isEEDeeGap                = pho->isEEDeeGap();
	recPhoton.isEBEEGap              = pho->isEBEEGap();
 //       recPhoton.r9                     = pho->r9();
        recPhoton.et                     = pho->et();
        recPhoton.eta                    = pho->eta();
        recPhoton.phi                    = pho->phi();
        recPhoton.hadronicOverEm         = pho->hadronicOverEm();
        recPhoton.sigmaIetaIeta          = pho->sigmaIetaIeta();

        // Fill the tree (this records all the recPhoton.* since
        // tree_PhotonAll_ was set to point at that.
        tree_PhotonAll_->Fill();
      }

      // Record whether it was near any module gap.
      // Very convoluted at the moment.
      bool inAnyGap = pho->isEBEEGap() || (pho->isEB()&&pho->isEBEtaGap()) ||(pho->isEB()&&pho->isEBPhiGap())  || (pho->isEE()&&pho->isEERingGap()) || (pho->isEE()&&pho->isEEDeeGap()) ;
      if (inAnyGap) {
        h_photonInAnyGap_->Fill(1.0);
      } else {
        h_photonInAnyGap_->Fill(0.0);
      }

      photonCounter++;

    } 
    else
    {
      // This didn't pass photon cuts, mark it
   
    }
    idxpho++;
    
    //All of this is pretty trivial, just getting the information
    //from the reco::Photon itself.
    
    photonSCeta[nPho]=pho->superCluster()->position().eta();
    photonSCphi[nPho]=pho->superCluster()->position().phi();
    photonSCX[nPho]=pho->superCluster()->position().x();
    photonSCY[nPho]=pho->superCluster()->position().y();
    photonSCZ[nPho]=pho->superCluster()->position().z();
    photonSCE[nPho]=pho->superCluster()->energy();
    photonSCetawidth[nPho]=pho->superCluster()->etaWidth();
    photonSCphiwidth[nPho]=pho->superCluster()->phiWidth();
    photonet[nPho]=pho->et();
    photon_physeta[nPho]=pho->eta();
    photon_physphi[nPho]=pho->phi();
    photone[nPho]=pho->energy();
    photonhadTowOverEm[nPho]=pho->hadTowOverEm();
    photonsigmaIetaIeta[nPho]=pho->sigmaIetaIeta();
    r9[nPho]=pho->r9();

    
    //Get the ParticleFlow photon isolation
    //fGetIsolation(const reco::Photon * photon, const reco::PFCandidateCollection* pfParticlesColl,reco::VertexRef vtx, edm::Handle< reco::VertexCollection > vertices)
	  reco::VertexRef myVtxRef(Vertices, 0);
    isolator.fGetIsolation(&*pho,pfParticlesColl.product(), myVtxRef, Vertices);
      photonchargedHadronIso[nPho] = isolator.getIsolationCharged();
	  photonneutralHadronIso[nPho] = isolator.getIsolationNeutral();
	  photonphotonIso[nPho] = isolator.getIsolationPhoton();
//     cout<<"PF  :  "<<isolator.getIsolationCharged()<<" : "<<isolator.getIsolationPhoton()<<" : "<<isolator.getIsolationNeutral()<<endl;
    

  nPho++;
  } // End Loop over photons
  h_nPho_->Fill(photonCounter);
///////////////////////////////////////////

edm::InputTag PileupSrc_("addPileupInfo");
  Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  evt.getByLabel(PileupSrc_, PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;

  // (then, for example, you can do)
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {

    //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
	PU_NumInteractions = PVI->getPU_NumInteractions();
  }

//Then when you get to your analyzer, you have to pick up on both the trigger information, and whether or not the list has changed

      Handle<TriggerResults> HLTR;
      evt.getByLabel(hlTriggerResults_,HLTR);

      Handle<trigger::TriggerEvent> triggerEventHandle;
      evt.getByLabel(triggerEventTag_,triggerEventHandle);

      
      //check
      assert(HLTR->size()==hltConfig_.size());
      

     if (HLTR.isValid()){ 
     //Clear the decks, these are data members of the Analyzer object
     all_triggerprescales.clear();
     all_ifTriggerpassed.clear();

     const edm::TriggerNames &triggerNames_ = evt.triggerNames(*HLTR);
     hlNames_=triggerNames_.triggerNames();
     
     vector<int> idx;

     //Recall:  This is the ntriggers value that we found from our loop in beginRun
     for(int i = 0; i< ntriggers;i++){
       all_triggerprescales.push_back(0);
       all_ifTriggerpassed.push_back(0);
       
       //idx stores the index of the particular trigger
       idx.push_back(triggerNames_.triggerIndex(all_triggers[i]));


     }//for(int i = 0; i< ntriggers;i++)

     Int_t hsize = Int_t(HLTR->size());
     for(int i=0;i<ntriggers;i++){
       if(idx[i] < hsize){
         all_ifTriggerpassed[i]=HLTR->accept(idx[i]);
         all_triggerprescales[i]=hltConfig_.prescaleValue( evt, es, all_triggers[i]);
       }
     }    
     }//if HLTR is Valid

  	//////////////////////        
  if (nPho>0) 
//    cout << "vx = " << vx[1] << endl;
//    cout << "r9 = " << r9[1] << endl;
    PhotonIDTree->Fill(); 
}

///////////////////////////////////////////////////////////////
void PhotonIDSimpleAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{  
  //this piece of code is run every time the run is changed, and resets all the triggers which are present if the trigger list is different
  bool changed(true); //this keeps track of whether the trigger list has changed
  if(hltConfig_.init(iRun,iSetup,hltlabel_,changed)){  // HLTlabel_ is grabbed from the cfg file
    // if init returns TRUE, initialisation has succeeded!
    if(changed){
      // The HLT config has actually changed wrt the previous Run, hence rebook your
      // histograms or do anything else dependent on the revised HLT config
      photon_triggers_in_run.clear();
      unsigned int ntriggers = hltConfig_.size();
      
      // Loop over all available triggers
      for(unsigned int t=0;t<ntriggers;++t){
	std::string hltname(hltConfig_.triggerName(t));
        string string_search ("HLT_Photon");
        
        //search the trigger name for string_search. 
        size_t found = hltname.find(string_search);
       
        if(found!=string::npos ){
          photon_triggers_in_run.push_back(hltname);
        }
      }//loop over ntriggers
      //This has to be clean for every run as Menu get changed frequently
      all_triggerprescales.clear();
      all_ifTriggerpassed.clear();
      all_triggers.clear(); 

  
      //SAVE PHOTON TRIGGER INFO
      for(int x = 0; x< (int)photon_triggers_in_run.size();x++){
        bool found = false;

        for(int i = 0; i< (int)all_triggers.size();i++){
          if(all_triggers[i]==photon_triggers_in_run[x]) found = true;
        }//loop all triggers

        if(!found)
          all_triggers.push_back(photon_triggers_in_run[x]); 
      }//loop photon triggers
           
      
    }//loop over all available triggers
  }else{
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    std::cout << " HLT config extraction failure with name " << hlTriggerResults_ << std::endl;
    // In this case, all access methods will return empty values!
  }
  ntriggers = (int)all_triggers.size();
}//ending method
//////////////////////////////////////////////////////////////////////
//                                                                  //
/////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
//    method called once each job just after ending the event loop   //
///////////////////////////////////////////////////////////////////////
void 
PhotonIDSimpleAnalyzer::endJob()
{

  // go to *OUR* root file and store histograms
  rootFile_->cd();

/* Commenting out old Histograms so we get PhotonIDTree only.
  // PhotonID Histograms
  h_isoEcalRecHit_->Write();
  h_isoHcalRecHit_->Write();
  h_trk_pt_solid_-> Write();
  h_trk_pt_hollow_->Write();
  h_ntrk_solid_->   Write();
  h_ntrk_hollow_->  Write();
  h_ebetagap_->     Write();
  h_ebphigap_->     Write();
  h_eeringGap_->     Write();
  h_eedeeGap_->     Write();
  h_ebeeGap_->   Write();
  h_r9_->        Write();

  // Photon Histograms
  h_photonEt_->  Write();
  h_photonEta_-> Write();
  h_photonPhi_-> Write();
  h_hadoverem_-> Write();

  // Photon's SuperCluster Histograms
  h_photonScEt_->      Write();
  h_photonScEta_->     Write();
  h_photonScPhi_->     Write();
  h_photonScEtaWidth_->Write();

  // Composite or Other Histograms
  h_photonInAnyGap_->Write();
  h_nPassingPho_->   Write();
  h_nPassEM_->       Write();
  h_nPho_->          Write();
*/

  // Write the new PhotonIDTree
  PhotonIDTree->Write();
  // Write the root file (really writes the TTree)
  rootFile_->Write();
  rootFile_->Close();

}

//define this as a plug-in
// DEFINE_FWK_MODULE(PhotonIDSimpleAnalyzer);
