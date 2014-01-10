#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
//#include <TDCacheFile.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TString.h>

#define SIZEOF_ARRAY( a ) (sizeof( a ) / sizeof( a[ 0 ] ))
#define NCUT  18    // number of individual cuts
// To Add cuts increment NCUT and expand candcuts array accordingly. Also double check cuts index
#define NSETS  8    // number of cut sets (loose cuts, reco cuts, EV -1 cut, HE -1 cut, etc)
#define ELEMAX 100 // array length for all photons in an event

// Conditions to run code
#define NLOOSE 2
#define NRECO  1

// #define

std::map<TString, TH1F*> hName;           // Map for histograms

void CreateHistogram(const char* name,   const char* title,
                                   const char* xTitle, const char* yTitle,
                                   Int_t       nBinsX, Double_t    xLow, Double_t xUp)
{
  TH1F* h = new TH1F(name, title, nBinsX, xLow, xUp);

  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);

//  h->Sumw2(); //associate error bars

  hName[name] = h;
}

void ShiftVertex(Float_t &eta, Float_t &phi, Float_t &photonet, Float_t photone, 
					Double_t vx, Double_t vy, Double_t vz,
					Float_t SCX, Float_t SCY, Float_t SCZ){
	float_t theta = atan2(sqrt(( SCX - vx )*( SCX - vx ) + ( SCY - vy )*( SCY - vy )), ( SCZ - vz ));
    eta = -1. * log ( tan( theta/2.0 ));
	phi = atan2(( SCY - vy ), (SCX - vx )); 
	photonet = photone * sin(theta);
}

Float_t InvMass(Float_t phiA,Float_t etaA,Float_t photonetA,
				Float_t phiB,Float_t etaB,Float_t photonetB){
				
	  Float_t pxA = cos(phiA) * photonetA;
	  Float_t pyA = sin(phiA) * photonetA;
	  Float_t photonthetaA = 2. * atan(exp(-1.*etaA));
	  Float_t pzA = photonetA * cos(photonthetaA) / sin(photonthetaA);
	  
	  Float_t pxB = cos(phiB) * photonetB;
	  Float_t pyB = sin(phiB) * photonetB;
	  Float_t photonthetaB = 2. * atan(exp(-1.*etaB));
	  Float_t pzB = photonetB * cos(photonthetaB) / sin(photonthetaB);

	Float_t EleAE = sqrt(pxA*pxA + pyA*pyA + pzA*pzA);
	Float_t EleBE = sqrt(pxB*pxB + pyB*pyB + pzB*pzB);		   			 
	Float_t DotProdAB = pxA*pxB + pyA*pyB + pzA*pzB;
	Float_t MassSqAB = 2*(EleAE*EleBE - DotProdAB);
	Float_t MassAB = sqrt(MassSqAB);
	return MassAB;
}


using namespace std;
void ThreePhotonPlotsv3(void){

// Adding some switches for convenience
   Bool_t          useprimaryvertex = true;
   Bool_t          outputdata = true;  // creates a txt file of invariant mass values
   Bool_t          sampleoutput = false; // just use 50 events for sample

// Below I was seeing the effect of doubling the weight of events with the wrong vertex chosen.
// Currently hard coded with 15 wrong entries / 473 right entries
//   float		   wrongvtxwgt = 2;  // 1 - wrong / right vertices
//   float		   rightvtxwgt = 1 + (1-wrongvtxwgt)*15.0/473.0;  // explained below
// If have X entries. X = 1*wrong + 1*right. To keep x constant,
// but weight right double we have x = C*right + D*wrong. Therefore C = 1 + (1-D) wrong / right
// cout << "right weight = " << rightvtxwgt << "  wrong weight = " << wrongvtxwgt << endl;
 
  ofstream InvMassdata;
  if (useprimaryvertex == true && outputdata == true) {
     InvMassdata.open ("OUTPUTHIST_InvMassUnbinned.txt");
  } else if (useprimaryvertex == false && outputdata == true){
       InvMassdata.open ("OUTPUTHIST_InvMassUnbinnedPrim-CalcDiff.txt");
  }

  TChain *PhotonIDTree = new TChain("PhotonIDTree");
//   PhotonIDTree->Add("PhoIDHistsA.root");
//   PhotonIDTree->Add("PhoIDHistsB.root");
//   PhotonIDTree->Add("PhoIDHistsC.root");
//   PhotonIDTree->Add("PhoIDHistsD.root");
//    PhotonIDTree->Add("PhoIDHistsMC250.root");
    PhotonIDTree->Add("INPUTMCTUPLE.root");


   Int_t           Run;
   Int_t           Event;
   Int_t           LumiSec;
   Int_t           nPho;
   Int_t           PU_NumInteractions;
   Float_t		   rho;
   Float_t         photonSCeta[ELEMAX];   //[nPho]
   Float_t         photonSCphi[ELEMAX];   //[nPho]
   Float_t         photonSCX[ELEMAX];   //[nPho]
   Float_t         photonSCY[ELEMAX];   //[nPho]
   Float_t         photonSCZ[ELEMAX];   //[nPho]
   Float_t         photonSCE[ELEMAX];
   Float_t         photonSCetawidth[ELEMAX];
   Float_t         photonSCphiwidth[ELEMAX];
   Float_t         photonet[ELEMAX];   //[nPho]
   Float_t         photon_physeta[ELEMAX];   //[nPho]
   Float_t         photon_physphi[ELEMAX];   //[nPho]
   Float_t         photone[ELEMAX];   //[nPho]
   Float_t         photonhadTowOverEm[ELEMAX];
   Float_t         photonsigmaIetaIeta[ELEMAX];
   Float_t         photonchargedHadronIso[ELEMAX];
   Float_t         photonneutralHadronIso[ELEMAX];
   Float_t         photonphotonIso[ELEMAX];
   Bool_t          passelectronveto[ELEMAX];
   Float_t         r9[ELEMAX];
   
   Int_t           nVert;
   Double_t         vx[ELEMAX];
   Double_t         vy[ELEMAX];
   Double_t         vz[ELEMAX];
   Int_t           gen_nVert;
   Double_t         gen_vx[ELEMAX];
   Double_t         gen_vy[ELEMAX];
   Double_t         gen_vz[ELEMAX];

   double 			Weight; // weight from PU reweighting
   double 			w; // total of all weights
   bool           cuts[ELEMAX][NCUT];
   float 		  RhoCorrPfNeutralHadronIso[ELEMAX];
   float 		  RhoCorrPfChargedHadronIso[ELEMAX];
   float 		  RhoCorrPfPhotonIso[ELEMAX];

   
   PhotonIDTree->SetBranchAddress("Run", &Run);
   PhotonIDTree->SetBranchAddress("Event", &Event);
   PhotonIDTree->SetBranchAddress("LumiSec", &LumiSec);
   PhotonIDTree->SetBranchAddress("nPho", &nPho);
   PhotonIDTree->SetBranchAddress("rho", &rho);
   PhotonIDTree->SetBranchAddress("photonSCeta", photonSCeta);
   PhotonIDTree->SetBranchAddress("photonSCphi", photonSCphi);
   PhotonIDTree->SetBranchAddress("photonSCE", photonSCE);
   PhotonIDTree->SetBranchAddress("photonSCetawidth", photonSCetawidth);
   PhotonIDTree->SetBranchAddress("photonSCphiwidth", photonSCphiwidth);
   PhotonIDTree->SetBranchAddress("photonSCX", photonSCX);
   PhotonIDTree->SetBranchAddress("photonSCY", photonSCY);
   PhotonIDTree->SetBranchAddress("photonSCZ", photonSCZ);
   PhotonIDTree->SetBranchAddress("photonet", photonet);
   PhotonIDTree->SetBranchAddress("photon_physeta", photon_physeta);
   PhotonIDTree->SetBranchAddress("photon_physphi", photon_physphi);
   PhotonIDTree->SetBranchAddress("photone", photone);
   PhotonIDTree->SetBranchAddress("photonhadTowOverEm", photonhadTowOverEm);
   PhotonIDTree->SetBranchAddress("photonsigmaIetaIeta", photonsigmaIetaIeta);
   PhotonIDTree->SetBranchAddress("photonchargedHadronIso", photonchargedHadronIso);
   PhotonIDTree->SetBranchAddress("photonneutralHadronIso", photonneutralHadronIso);
   PhotonIDTree->SetBranchAddress("photonphotonIso", photonphotonIso);
   PhotonIDTree->SetBranchAddress("passelectronveto", passelectronveto);
   PhotonIDTree->SetBranchAddress("r9", r9);

   PhotonIDTree->SetBranchAddress("Weight", &Weight);
   PhotonIDTree->SetBranchAddress("nVert", &nVert);
   PhotonIDTree->SetBranchAddress("PU_NumInteractions", &PU_NumInteractions);
   PhotonIDTree->SetBranchAddress("vx", vx);
   PhotonIDTree->SetBranchAddress("vy", vy);
   PhotonIDTree->SetBranchAddress("vz", vz);
   PhotonIDTree->SetBranchAddress("gen_nVert", &gen_nVert);
   PhotonIDTree->SetBranchAddress("gen_vx", gen_vx);
   PhotonIDTree->SetBranchAddress("gen_vy", gen_vy);
   PhotonIDTree->SetBranchAddress("gen_vz", gen_vz);
	
    TFile *aa = new TFile("OUTPUTPLOTS.root","RECREATE");

    TH1F *nGammaCan = new TH1F("nGammaCan","Number of photons per candidate event",15,0,15);
    TH1F *pu = new TH1F("pu","pu",50,0,50);
    TH1F *puUnweighted = new TH1F("puUnweighted","puUnweighted",50,0,50);

    //Kitchen sink plots
    TH1F *LeadPhotonPT = new TH1F("LeadPhotonPT","Leading Photon p_{T}",250,0,250);
    TH1F *NextPhotonPT = new TH1F("NextPhotonPT","Next Photon p_{T}",250,0,250);
    TH1F *TrailPhotonPT = new TH1F("TrailPhotonPT","Trailing Photon p_{T}",250,0,250);


    TH1F *calcEta = new TH1F("calcEta","calculated eta",250,-1.5,1.5);
    TH1F *calcPhi = new TH1F("calcPhi","calculated Phi",250,-3.2,3.2);
    TH1F *calcEt = new TH1F("calcEt","calculated Et",250,0,600);
    TH1F *Eta = new TH1F("Eta","Physics eta",250,-1.5,1.5);
    TH1F *Phi = new TH1F("Phi","Physics Phi",250,-3.2,3.2);
    TH1F *Et = new TH1F("Et","Physics Et",250,0,600);


    TH1F *PrimaryVertex = new TH1F("PrimaryVertex","Inv. Mass from Primary Vertex",250,230,230);
    TH1F *ShiftedVertex = new TH1F("ShiftedVertex","Inv. Mass from Random Vertex",250,230,270);
//    TH1F *PrimaryVertexWeighted = new TH1F("PrimaryVertexWeighted","Weighted Inv. Mass",250,230,270);

    TH1F *Diffvx = new TH1F("Diffvx","vx - gen_vx",250,-10,-10);
    TH1F *Diffvy = new TH1F("Diffvy","vy - gen_vy",250,-10,-10);
    TH1F *Diffvz = new TH1F("Diffvz","vz - gen_vz",250,-10,-10);

    TH1F *DiffEta = new TH1F("DiffEta","Physics eta - calc eta",250,-3,-3);
    TH1F *DiffPhi = new TH1F("DiffPhi","Physics phi - calc phi",250,3.2,3.2);
    TH1F *DiffEt = new TH1F("DiffEt","Physics et - calc et",250,-1,-1);
    TH1F *DiffEtzoom = new TH1F("DiffEtzoom","Physics et - calc et",250,-1,-1);
    TH1F *DiffEtaOverPhys = new TH1F("DiffEtaOverPhys","Physics eta - calc eta",250,-5,-5);
    TH1F *DiffPhiOverPhys = new TH1F("DiffPhiOverPhys","Physics phi - calc phi",250,5,5);
    TH1F *DiffEtOverPhys = new TH1F("DiffEtOverPhys","Physics et - calc et",250,-1,-1);



	vector<TString> SetName;
    SetName.push_back("loose");
    SetName.push_back("reco");
    SetName.push_back("EV_cutexcluded");
    SetName.push_back("HE_cutexcluded");
    SetName.push_back("NHIso_cutexcluded");
    SetName.push_back("CHIso_cutexcluded");
    SetName.push_back("PHIso_cutexcluded");
    SetName.push_back("SigIeta_cutexcluded");

    
//     TH1F *Hr9 = new TH1F("r9","r9",200,0,6);
//     TH1F *Hrho = new TH1F("rho","rho",200,-50,50);

for ( size_t i=0; i<SetName.size(); i++ ){
    CreateHistogram("Hpasselectronveto_" + SetName[i], "passelectronveto_" + SetName[i] + " (1=no tracker seed present)", "","", 3,0,3); 
    CreateHistogram("HphotonhadTowOverEm_" + SetName[i], "HphotonhadTowOverEm_" + SetName[i] + " (cut at 0.05)", "","",150,0,.06);
    CreateHistogram("HphotonneutralHadronIso_" + SetName[i], "photonneutralHadronIso_" + SetName[i] + " (cut at 3.5 + 0.04*pho_Pt)", "","",100,0,10); 
    CreateHistogram("HphotonchargedHadronIso_" + SetName[i], "photonchargedHadronIso_" + SetName[i] + " (cut at 2.6)", "","",100,0,10);
    CreateHistogram("HphotonphotonIso_" + SetName[i], "photonphotonIso_" + SetName[i] + " (cut at 1.3 + 0.005*pho_Pt)", "","",100,0,10);
    CreateHistogram("HphotonsigmaIetaIeta_" + SetName[i], "HphotonsigmaIetaIeta_" + SetName[i] + " (cut at 0.012)", "","",100,0,0.015);
    CreateHistogram("InvMass_" + SetName[i], "Invariant mass_" + SetName[i], "","",600,0,600);
    CreateHistogram("InvMass_" + SetName[i]+"unweighted", "Invariant mass_" + SetName[i]+"unweighted", "","",600,0,600);
	}

    vector <Int_t> RUN;
    vector <Int_t> EVENT;
    
    int matchvert = 0;
    Float_t nentries;

	if  (sampleoutput){
    	nentries = 50;
    }	
	else { 
    	nentries = PhotonIDTree->GetEntries();
    }
    
    cout << "nentries: " << nentries << endl;
    //Loop over Events
    for (int entry=0; entry<nentries; ++entry){
	
    // loop over cuts, set to false on both arrays
    
      PhotonIDTree->GetEntry(entry);
      if (entry%10000==0) cout << "entry: " << entry << " / " << nentries << endl;
      
	  w = Weight;
      
//		cout << "nVert: " << nVert << endl;
//		cout << "PU_NumInteractions: " << PU_NumInteractions << endl;
		pu->Fill(PU_NumInteractions,w);
		puUnweighted->Fill(PU_NumInteractions);

    // logic 1=pass, 0=fail, -1=I don't care

// 	 [NCUT]-> pt,electronveto,singleTowerHeCut,looseNHIso,looseCHIso,loosePHIso,looseSigmaIEtaIEta,tightNHIso,tightCHIso,tightPHIso,tightSigmaIEtaIEta,r9,fabs(SCeta), unused, unused}
//    therefore:            {{pt,ev,HE,---loose---,---tight---,r9,EB,-----unused----} 
 int candcuts[NSETS][NCUT]= {{ 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [0] loose photon
   							 { 1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,}, // [1] recophoton
   							 { 1,-1, 1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [2] -1 EV cut
   							 { 1, 1,-1, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [3] -1 HE cut
   							 { 1, 1, 1,-1, 1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [4] -1 NHIso cut
   							 { 1, 1, 1, 1,-1, 1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [5] -1 CHIso cut
   							 { 1, 1, 1, 1, 1,-1, 1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1,}, // [6] -1 PHIso cut
   							 { 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1, 1, 1,-1,-1,-1,-1,-1} };// [7] -1 SigIeta cut


// iset is loose cuts, reco cuts, 1st -1 cut, 2nd -1 cut, etc
// initialize candpass array to all true
// candpass is array of photons in event stating pass or fail for each set of cuts.
  Int_t ncandpass[NSETS]={0};
  bool candpass[NSETS][ELEMAX];
    

// initialize candpass array for each event
  for (int iset=0;iset<NSETS;iset++) {
    for (int iele=0;iele<ELEMAX;iele++){
      candpass[iset][iele]= true;
    }
  }

// initialize cuts array for each event
  for (int icut=0;icut<NCUT;icut++) {
    for (int iele=0;iele<ELEMAX;iele++){
      cuts[iele][icut]= 0;
    }
  }
  
      Int_t RecoEle[100];
      Int_t nRecoEle=0;


      //Loop Over Photons
      for (int pho=0;pho<nPho;++pho){      
      	   	
	   		//begin defining rho correction areas
		float EA_CH = 0; 
		float EA_NE = 0;
		float EA_PH = 0;
	    
		if(fabs(photon_physeta[pho]) < 1.0){
		  EA_CH = 0.012;
		  EA_NE = 0.030;
		  EA_PH = 0.148;
		}
		
		if(fabs(photon_physeta[pho]) > 1.0	
		   && fabs(photon_physeta[pho]) < 1.479
		   ){
		  EA_CH = 0.010;
		  EA_NE = 0.057;
		  EA_PH = 0.130;
		}

		RhoCorrPfNeutralHadronIso[pho] = ((photonneutralHadronIso[pho]) - (rho)*EA_NE);
   	    RhoCorrPfChargedHadronIso[pho] = ((photonchargedHadronIso[pho]) - (rho)*EA_CH);
 	    RhoCorrPfPhotonIso[pho] = ((photonphotonIso[pho]) - (rho)*EA_PH);
// 	
		//end rho correction areas	
		//Fill rho corrected histograms
// 		H_RCPfNHIso->Fill(RhoCorrPfNeutralHadronIso,w);
// 		H_RCPfCHIso->Fill(RhoCorrPfChargedHadronIso,w);
// 		H_RCPfPHIso->Fill(RhoCorrPfPhotonIso,w);				
		

        cuts[pho][0] = (photonet[pho] > 20.);
        cuts[pho][1] = passelectronveto[pho];
        cuts[pho][2] = (photonhadTowOverEm[pho] < 0.05);
        cuts[pho][3] = (RhoCorrPfNeutralHadronIso[pho] < (3.5 + 0.04 * photonet[pho]));
        cuts[pho][4] = (RhoCorrPfChargedHadronIso[pho] < 2.6);
        cuts[pho][5] = (RhoCorrPfPhotonIso[pho] < (1.3 + 0.005 * photonet[pho]));
        cuts[pho][6] = (photonsigmaIetaIeta[pho] < 0.012);
        cuts[pho][7] = (RhoCorrPfChargedHadronIso[pho] < 0.7);
        cuts[pho][8] = (RhoCorrPfPhotonIso[pho] < (0.5 + 0.005 * photonet[pho]));
        cuts[pho][9] = (RhoCorrPfNeutralHadronIso[pho] < (0.4 + 0.04 * photonet[pho]));
        cuts[pho][10] = (photonsigmaIetaIeta[pho] < 0.011);
        cuts[pho][11] = (r9[pho]>0.);
        cuts[pho][12] = (fabs(photonSCeta[pho])<1.44);
      
//Loop over sets
for (int iset=0;iset<NSETS;iset++) {
        
		// do candcuts and fill results into candpass
		for (int icut=0;(icut<NCUT)&&(candpass[iset][pho]);icut++) {
		  if (candcuts[iset][icut]==1) {
			if (cuts[pho][icut]!=1) candpass[iset][pho]=false;
		  } else if (candcuts[iset][icut]==0) {
			if (cuts[pho][icut]!=0) candpass[iset][pho]=false;
		  }
		} // close icut loop (we ignore a cut if it's set to -1)
		
		// count photons
		if (candpass[iset][pho]) {
		  ncandpass[iset]++;
  		    // now fill histos for all photons passing cuts (regardless of photons/event)
		  if (iset == 0) {
			// fill loose photons
		  } else if (iset == 1) {
		  	RecoEle[nRecoEle] = pho;
		    nRecoEle++;
			// fill reco photons
// 			Hpasselectronveto->Fill(passelectronveto[pho],w);
// 			HphotonhadTowOverEm->Fill(photonhadTowOverEm[pho],w);
// 			HphotonsigmaIetaIeta->Fill(photonsigmaIetaIeta[pho],w);
 		  } // etc...
		} // endif candpass cases

    } //end loop over all photons
   } //end iset loop


//Loop over sets
for (int iset=0;iset<NSETS;iset++) {
// now we have candpass array filled for each event, look only at three photon events
//    if ((ncandpass[iset]>2)||((ncandpass[iset]==2)&&((ncandpass[1]-ncandpass[iset])>0))) {
     if (ncandpass[iset]>2) {
     // good three photon event
	if (iset==0){ //record Run and Event since we have 3 possible photons (only once per cut set)
      nGammaCan->Fill(nRecoEle,w);
//  	Hrho->Fill(rho,w);
      RUN.push_back(Run); 
	  EVENT.push_back(Event);
//	  cout << "THREE GOOD: event# " << entry << "  nRecoEle " << nRecoEle << "  ncandpass[1] " << ncandpass[1] << "  nloose " << ncandpass[0] << endl;
	}
    // find the three photons
//    int threephoton[3];
    int threephoton[3] = {-1,-1,-1};
    float threephotonpt[4] = {-1, -1, -1, -1};
	// find highest pt passing photon
	for (int i=0;i<nRecoEle;i++) {
		int iph = RecoEle[i];
		if (candpass[iset][iph]) {
		  if (photonet[iph]>threephotonpt[0]) {
			threephotonpt[0] = photonet[iph];
			threephoton[0] = iph; 
//			if (threephoton[0]==-1){cout << "Eek! I didn't find the highest passing three photon!" << endl ;} // warning if not filled
		  }
		}
	} // close loop over photons
	// find 2nd highest pt passing photon
	for (int i=0;i<nRecoEle;i++) {
		int iph = RecoEle[i];
		if ((candpass[iset][iph])&&(iph!=threephoton[0])) {
		  if (photonet[iph]>threephotonpt[1]) {
			threephotonpt[1] = photonet[iph];
			threephoton[1] = iph;
		  }
		}
	} // close loop over photons
	// find 3rd highest pt "candidate" photon OR highest pt reco that is not passing
	for (int i=0;i<nRecoEle;i++) {
		int iph = RecoEle[i];
		if ((iph!=threephoton[0])&&(iph!=threephoton[1])) {
		  if (ncandpass[iset]>2) {
			if (candpass[iset][iph]) { 
			  if (photonet[iph]>threephotonpt[2]) {
				threephotonpt[2] = photonet[iph];
				threephoton[2] = iph;
//		        if (iset==0){cout << "3 case: event# " << entry << "   index " << iph << "   pt " << photonet[iph] << "   ev " << passelectronveto[iph] << "   sigmaIetaIeta " << photonsigmaIetaIeta[iph] << endl;} 
			  }
			}
		  }
		  else if ((ncandpass[1]-2)>0) {
			 if (candpass[1][iph]) { 
			  if (photonet[iph]>threephotonpt[2]) {
				threephotonpt[2] = photonet[iph];
				threephoton[2] = iph;
//		        if (iset==0){cout << "2+1 case: event# " << entry << "   index " << iph << "   pt " << photonet[iph] << "   ev " << passelectronveto[iph] << "   sigmaIetaIeta " << photonsigmaIetaIeta[iph] << endl;} 
			  }
			 }
		  }
		} // endif not first two highest pt photons
	} // close loop over RecoPhotons
		
	// begin plot filling for "candidate" events
    if (iset == 0) {
		// fill loose photon event info
		LeadPhotonPT->Fill(threephotonpt[0],w);
		NextPhotonPT->Fill(threephotonpt[1],w);
		TrailPhotonPT->Fill(threephotonpt[2],w); 
  	} else if (iset == 2) {
		// fill -1 EV photons
  	} // etc...
  	
  	float_t etaA = photon_physeta[threephoton[0]];
	float_t phiA = photon_physphi[threephoton[0]];
	float_t etA  = photonet[threephoton[0]];  	

  	float_t etaB = photon_physeta[threephoton[1]];
	float_t phiB = photon_physphi[threephoton[1]];
	float_t etB  = photonet[threephoton[1]];  	

  	float_t etaC = photon_physeta[threephoton[2]];
	float_t phiC = photon_physphi[threephoton[2]];
	float_t etC  = photonet[threephoton[2]];  	
	
	float_t MassAB = InvMass(phiA,etaA,etA,phiB,etaB,etB);
	float_t MassAC = InvMass(phiA,etaA,etA,phiC,etaC,etC);
	float_t MassBC = InvMass(phiB,etaB,etB,phiC,etaC,etC);
	  
	hName["InvMass_" + SetName[iset]]->Fill(MassAB,w);
	hName["InvMass_" + SetName[iset]]->Fill(MassAC,w);
	hName["InvMass_" + SetName[iset]]->Fill(MassBC,w);

	hName["InvMass_" + SetName[iset] + "unweighted"]->Fill(MassAB);
	hName["InvMass_" + SetName[iset] + "unweighted"]->Fill(MassAC);
	hName["InvMass_" + SetName[iset] + "unweighted"]->Fill(MassBC);
	
	bool vertmatch = true;

	if ( iset==0 ) {
		Diffvx->Fill(vx[0]-gen_vx[0],w);
		Diffvy->Fill(vy[0]-gen_vy[0],w);
		Diffvz->Fill(vz[0]-gen_vz[0],w);
	
// Below was just for seeing effects of doubling weight of wrong vertex selection.
// 			for (int it=1; it<nVert;it++) 
// 			  { if (fabs(vz[0]-gen_vz[0])>fabs(vz[it]-gen_vz[0])) {
// 				  vertmatch = false;
// 				  }
// 			  }
// 			if (vertmatch == true)
// 			  {
// 			  PrimaryVertexWeighted->Fill(MassAB,rightvtxwgt);
// 			  PrimaryVertexWeighted->Fill(MassAC,rightvtxwgt);
// 			  PrimaryVertexWeighted->Fill(MassBC,rightvtxwgt);
// 			  }
// 			else 
// 			  {
// 			  PrimaryVertexWeighted->Fill(MassAB,wrongvtxwgt);
// 			  PrimaryVertexWeighted->Fill(MassAC,wrongvtxwgt);
// 			  PrimaryVertexWeighted->Fill(MassBC,wrongvtxwgt);
// 			  }	
	  PrimaryVertex->Fill(MassAB,w);
	  PrimaryVertex->Fill(MassAC,w);
	  PrimaryVertex->Fill(MassBC,w);
	  
	if (outputdata == true){
		  InvMassdata << MassAB << " " << w << "\n";
		  InvMassdata << MassAC << " " << w << "\n";
		  InvMassdata << MassBC << " " << w << "\n";		
	}

	Eta->Fill(etaA,w);
	Eta->Fill(etaB,w);
	Eta->Fill(etaC,w);
	Phi->Fill(phiA,w);
	Phi->Fill(phiB,w);
	Phi->Fill(phiC,w);
	Et->Fill(etA,w);
	Et->Fill(etB,w);
	Et->Fill(etC,w);	
	}  		  

//   	  if ( iset==0 ) {
//   	  cout << "eta: " << etaA << "   phi: " << phiA << "   photonet: " << etA << "\n"; 
//   	  }
  	if (useprimaryvertex == true) {
//	if ( iset==0 && r9[threephoton[0]]<0.94 && r9[threephoton[1]]<0.94 && r9[threephoton[2]]<0.94) {
	if ( iset==0) {
	Diffvx->Fill(vx[0]-gen_vx[0],w);
	Diffvy->Fill(vy[0]-gen_vy[0],w);
	Diffvz->Fill(vz[0]-gen_vz[0],w);
	vertmatch = true;
	
		for (int it=1; it<nVert;it++) {
			if (fabs(vz[0]-gen_vz[0])>fabs(vz[it]-gen_vz[0])) {
			vertmatch = false;
//			cout << "Event " << entry  << " with "<< nVert << " vertices " << 
//			"   gen_vz[0]= "<< gen_vz[0] <<
//			"   vz[0]= " << vz[0] << 
//			"   vz["<< it <<"]= " << vz[it] << endl <<
//			"fabs(vz[0]-gen_vz[0]) = " << fabs(vz[0]-gen_vz[0]) <<
//			"   fabs(vz[it]-gen_vz[0]) = " << fabs(vz[it]-gen_vz[0]) <<
//			 endl << endl;
			}
	 	}
		if (vertmatch == true)
		 {
		matchvert++;
		}	
	}
	}

  	if (useprimaryvertex == false) {
//	if ( iset==0 && r9[threephoton[0]]<0.94 && r9[threephoton[1]]<0.94 && r9[threephoton[2]]<0.94) {
	if ( iset==0) {
//	for (int it=0; it<nVert;it++) { 
//	for (int randvert=0;randvert<nVert;randvert++) { 
   	srand(time(NULL));
 	  int randvert = (rand() % nVert);
//	  int randvert = 0;
//  	  cout << "randvert " << 
      ShiftVertex(etaA, phiA, etA, photone[threephoton[0]], 
					vx[randvert], vy[randvert], vz[randvert],
					photonSCX[threephoton[0]], photonSCY[threephoton[0]], photonSCZ[threephoton[0]]);
      ShiftVertex(etaB, phiB, etB, photone[threephoton[1]], 
					vx[randvert], vy[randvert], vz[randvert],
					photonSCX[threephoton[1]], photonSCY[threephoton[1]], photonSCZ[threephoton[1]]);
      ShiftVertex(etaC, phiC, etC, photone[threephoton[2]], 
					vx[randvert], vy[randvert], vz[randvert],
					photonSCX[threephoton[2]], photonSCY[threephoton[2]], photonSCZ[threephoton[2]]);
	 

// Next block of code is just to check out the vertex info. Remove before full run.
// 	  if ( iset==0 ) {
// 	  cout << "vx: " << vx[0] << "  vy: " << vy[0] << "  vz: " << vz[0] << "\n"; 
// 	  cout << "eta: " << etaA << "   phi: " << phiA << "   photonet: " << etA << "\n";  
// 	  }
	  
	  
	  
	  float_t MassABshift = InvMass(phiA,etaA,etA,phiB,etaB,etB);
	  float_t MassACshift = InvMass(phiA,etaA,etA,phiC,etaC,etC);
	  float_t MassBCshift = InvMass(phiB,etaB,etB,phiC,etaC,etC);
	  
	  cout << "MassAB   =  " << MassAB << endl;
if (outputdata == true){
	  InvMassdata << MassAB << "\n";
	  InvMassdata << MassAC << "\n";
	  InvMassdata << MassBC << "\n";				  
}
		ShiftedVertex->Fill(MassABshift,w);
		ShiftedVertex->Fill(MassACshift,w);
		ShiftedVertex->Fill(MassBCshift,w);
	
	
	
	calcEta->Fill(etaA,w);
	calcEta->Fill(etaB,w);
	calcEta->Fill(etaC,w);
	calcPhi->Fill(phiA,w);
	calcPhi->Fill(phiB,w);
	calcPhi->Fill(phiC,w);
	calcEt->Fill(etA,w);
	calcEt->Fill(etB,w);
	calcEt->Fill(etC,w);
		
double diffphiA = photon_physphi[threephoton[0]]-phiA;
double diffphiB = photon_physphi[threephoton[1]]-phiB;
double diffphiC = photon_physphi[threephoton[2]]-phiC;
	
	DiffEta->Fill(photon_physeta[threephoton[0]]-etaA,w);
	DiffEta->Fill(photon_physeta[threephoton[1]]-etaB,w);
	DiffEta->Fill(photon_physeta[threephoton[2]]-etaC,w);
	DiffPhi->Fill(min(fabs(diffphiA), (2*3.14159 - fabs(diffphiA))),w);
	DiffPhi->Fill(min(fabs(diffphiB), (2*3.14159 - fabs(diffphiB))),w);
	DiffPhi->Fill(min(fabs(diffphiC), (2*3.14159 - fabs(diffphiC))),w);
	DiffEt->Fill(photonet[threephoton[0]]-etA,w);
	DiffEt->Fill(photonet[threephoton[1]]-etB,w);
	DiffEt->Fill(photonet[threephoton[2]]-etC,w);
	DiffEtzoom->Fill(photonet[threephoton[0]]-etA,w);
	DiffEtzoom->Fill(photonet[threephoton[1]]-etB,w);
	DiffEtzoom->Fill(photonet[threephoton[2]]-etC,w);
// 

//	if (((photonet[threephoton[2]]-etC)/photonet[threephoton[2]])<-0.8){
// 	if ((min(fabs(diffphiC), (2*3.14159 - fabs(diffphiC))))>1.5){
// 		cout << "I have DiffPhi >1.5 I'm event: " << entry << endl;
//  	  cout << "vx: " << vx[0] << "  vy: " << vy[0] << "  vz: " << vz[0] << "\n"; 
//  	  cout << "SCx: " << photonSCX[threephoton[2]] << "  SCy: " << photonSCY[threephoton[2]] << "  SCz: " << photonSCZ[threephoton[2]] << "\n"; 
//  	  cout << "Phys eta: " << photon_physeta[threephoton[2]] << " Phys phi: " << photon_physphi[threephoton[2]] << "  Phys photonet: " << photonet[threephoton[2]] << "\n";  
//  	  cout << "eta: " << etaC << "   phi: " << phiC << "   photonet: " << etC << "\n";  
// 	  cout << "r9: " << r9[threephoton[2]] << "  diffphi:  " << (min(fabs(diffphiC), (2*3.14159 - fabs(diffphiC)))) << "\n"<<  "\n";  
// 	  
// 		}

	DiffEtaOverPhys->Fill((photon_physeta[threephoton[0]]-etaA)/photon_physeta[threephoton[0]],w);
	DiffEtaOverPhys->Fill((photon_physeta[threephoton[1]]-etaB)/photon_physeta[threephoton[1]],w);
	DiffEtaOverPhys->Fill((photon_physeta[threephoton[2]]-etaC)/photon_physeta[threephoton[2]],w);
	DiffPhiOverPhys->Fill((photon_physphi[threephoton[0]]-phiA)/photon_physphi[threephoton[0]],w);
	DiffPhiOverPhys->Fill((photon_physphi[threephoton[1]]-phiB)/photon_physphi[threephoton[1]],w);
	DiffPhiOverPhys->Fill((photon_physphi[threephoton[2]]-phiC)/photon_physphi[threephoton[2]],w);
	DiffEtOverPhys->Fill((photonet[threephoton[0]]-etA)/photonet[threephoton[0]],w);
	DiffEtOverPhys->Fill((photonet[threephoton[1]]-etB)/photonet[threephoton[1]],w);
	DiffEtOverPhys->Fill((photonet[threephoton[2]]-etC)/photonet[threephoton[2]],w);
	 }

//	  }
  	}
//  	  if ( iset==0 ) {
//  	  cout << "------        " << "\n\n\n";
//   	  }
  	//  plot filling for "candidate" photons in the event (pt ordered)
	for (int i=0;i<3;i++) {
 			hName["Hpasselectronveto_" + SetName[iset]]->Fill(passelectronveto[threephoton[i]],w);
 			hName["HphotonhadTowOverEm_" + SetName[iset]]->Fill(photonhadTowOverEm[threephoton[i]],w);
  			hName["HphotonneutralHadronIso_" + SetName[iset]]->Fill(RhoCorrPfNeutralHadronIso[threephoton[i]],w);
 			hName["HphotonchargedHadronIso_" + SetName[iset]]->Fill(RhoCorrPfChargedHadronIso[threephoton[i]],w);
 			hName["HphotonphotonIso_" + SetName[iset]]->Fill(RhoCorrPfPhotonIso[threephoton[i]],w);
  			hName["HphotonsigmaIetaIeta_" + SetName[iset]]->Fill(photonsigmaIetaIeta[threephoton[i]],w);
    } // end loop over 3 candidates

   }// endif three passing or two + reco 
      
   } //end iset loop

  }//entries
    
    cout << "Total events with >2 loose or 2+1 photons = "<< EVENT.size() << endl;
	cout << "Total events matching gen vertices = " << matchvert << endl;
	cout << "Total events selecting other vertices = " <<  EVENT.size() - matchvert << endl;
	cout << "Percent of time selecting correct vertex = " << 100*matchvert/EVENT.size() << " %" << endl;	
	aa->Write();
    aa->Close();
if (outputdata == true){InvMassdata.close();}
}	
