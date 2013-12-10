#ifndef RecoEgamma_PhotonIdentification_PhotonIDSimpleAnalyzer_H
#define RecoEgamma_PhotonIdentification_PhotonIDSimpleAnalyzer_H

/**\class PhotonIDSimpleAnalyzer

 Description: Analyzer to make a load of histograms for the improvement of the PhotonID object

 Implementation:
     \\\author: J. Stilley, A. Askew May 2008
*/
//

//

// system include files
#include <memory>

// user include files

// next two are for pileup
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"


#include <string>
#include "TH1.h"
#include "TTree.h"

// PFIsolation
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "EGamma/EGammaAnalysisTools/interface/PFIsolationEstimator.h"



class TFile;

//
// class declaration
//
class PhotonIDSimpleAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PhotonIDSimpleAnalyzer( const edm::ParameterSet& );
      ~PhotonIDSimpleAnalyzer();


      virtual void analyze( const edm::Event&, const edm::EventSetup& );
      virtual void beginJob();
      virtual void endJob();
      virtual void beginRun( const edm::Run&, const edm::EventSetup& );

      HLTConfigProvider hltConfig_;
      
 private:
      //HLTConfigProvider hltConfig_;
      std::vector<std::string> photon_triggers_in_run;
      std::vector<std::string> all_triggers;
      std::vector<int> all_triggerprescales;
      std::vector<bool> all_ifTriggerpassed;
      std::vector<std::string> *triggernames;
      std::vector<int> *triggerprescales;
      std::vector<bool> *ifTriggerpassed;
      int ntriggers;


      std::vector<std::string>  hlNames_;
      edm::TriggerNames triggerNames_;  // TriggerNames class
      edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults
      edm::InputTag triggerEventTag_;
      edm::InputTag trigEventTag_;
      std::string hltlabel_;

      std::string outputFile_;   // output file
      double minPhotonEt_;       // minimum photon Et
      double minPhotonAbsEta_;   // min and
      double maxPhotonAbsEta_;   // max abs(eta)
      double minPhotonR9_;       // minimum R9 = E(3x3)/E(SuperCluster)
      double maxPhotonHoverE_;   // maximum HCAL / ECAL 
      bool   createPhotonTTree_; // Create a TTree of photon variables
      TTree *PhotonIDTree;
      Int_t Run;
      Int_t Event;
      Int_t LumiSec;
      
      Int_t nPho;
      Float_t rho;
      
      Float_t photonSCeta[100];
      Float_t photonSCphi[100];
      Float_t photonSCE[100];
      Float_t photonSCetawidth[100];
      Float_t photonSCphiwidth[100];
      Float_t photonSCX[100];
      Float_t photonSCY[100];
      Float_t photonSCZ[100];
      Float_t photonet[100];
      Float_t photon_physeta[100];
      Float_t photon_physphi[100];
      Float_t photone[100];
      Float_t photonhadTowOverEm[100];
      Float_t photonsigmaIetaIeta[100];
      PFIsolationEstimator isolator;
      Float_t photonchargedHadronIso[100];
      Float_t photonneutralHadronIso[100];
      Float_t photonphotonIso[100];
      Bool_t passelectronveto[100];
      Float_t r9[100];

      Int_t nVert;
      Int_t PU_NumInteractions;
	  Double_t vx[100];
	  Double_t vy[100];      
	  Double_t vz[100];    
	                
      Int_t current_barcode;
      Int_t gen_nVert;
	  Double_t gen_vx[100];
	  Double_t gen_vy[100];      
	  Double_t gen_vz[100];    
 
      // Will be used for creating TTree of photons.
      // These names did not have to match those from a phtn->...
      // but do match for clarity.
      struct struct_recPhoton {
        float isolationEcalRecHit;
	float isolationHcalRecHit;
	float isolationSolidTrkCone;
	float isolationHollowTrkCone;
	float nTrkSolidCone;
	float nTrkHollowCone;
        float isEBEtaGap;
        float isEBPhiGap;
	float isEERingGap;
	float isEEDeeGap;
	float isEBEEGap;
	float et;
	float eta;
	float phi;
        float hadronicOverEm;
        float sigmaIetaIeta;
      } ;
      struct_recPhoton recPhoton;

      // root file to store histograms
      TFile*  rootFile_;

      // data members for histograms to be filled

      // PhotonID Histograms
      TH1F* h_isoEcalRecHit_;
      TH1F* h_isoHcalRecHit_;
      TH1F* h_trk_pt_solid_;
      TH1F* h_trk_pt_hollow_;
      TH1F* h_ntrk_solid_;
      TH1F* h_ntrk_hollow_;
      TH1F* h_ebetagap_;
      TH1F* h_ebphigap_;
      TH1F* h_eeringGap_;
      TH1F* h_eedeeGap_;
      TH1F* h_ebeeGap_;
      TH1F* h_r9_;

      // Photon Histograms
      TH1F* h_photonEt_;
      TH1F* h_photonEta_;
      TH1F* h_photonPhi_;
      TH1F* h_hadoverem_;

      // Photon's SuperCluster Histograms
      TH1F* h_photonScEt_;       
      TH1F* h_photonScEta_;      
      TH1F* h_photonScPhi_;
      TH1F* h_photonScEtaWidth_;

      // Composite or Other Histograms
      TH1F* h_photonInAnyGap_;
      TH1F* h_nPassingPho_;
      TH1F* h_nPassEM_;
      TH1F* h_nPho_;

      // TTree
      TTree* tree_PhotonAll_;
};
#endif
