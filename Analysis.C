#include <TLorentzVector.h>
#include <vector>
#include <iostream>
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "../TMVA_Regression.h"
#include <math.h>
#include <algorithm> 
#include <string>
#include "../plugin/puweicalc.h"
#include "../untuplizerv8.h" 
#include "../plugin/Jets.h"
#include "../plugin/ElectronSelections_v2.h"
#include "../plugin/MuonSelections_v2.h"
#include "../plugin/PhotonSelections_v.h"
#include "../plugin/RoccoR.cc" 
#include "../plugin/ScalingFactors.h"
#include "external/BTagCalibrationStandalone.cpp"
#include "external/BTagCalibrationStandalone.h"

struct HZgCand {
	Int_t bjet[2];
	vector<float> reg_fac;
	Int_t pho;
	Int_t ev;
};

void mcMu(TreeReader &data, std::vector<int> &accepted){
  accepted.clear();
  Int_t    nMC = data.GetInt("nMC");
  Float_t* mcEta = data.GetPtrFloat("mcEta");
  Int_t*    mcPID = data.GetPtrInt("mcPID");
  Int_t*    mcMomPID = data.GetPtrInt("mcMomPID");
  
  for (int i = 0; i < nMC; i++){
	if (fabs(mcPID[i]) != 13) continue;
	if (fabs(mcEta[i]) > 2.4) continue;
	accepted.push_back(i);
    }
}

void xAna(const char* inpath, TString outpath = "minitree.root", Int_t channel=0, Float_t xs=889.3, Float_t lumi=35.9, Int_t aMCatNLO=0, Int_t DoRegression=0, 
bool mucorr = false){

	TreeReader data(inpath);

	TFile* fo = TFile::Open(outpath.Data(),"RECREATE");
	if (!fo || fo->IsZombie())  FATAL("TFile::Open() failed");
	fo->cd();

	TFile *fe_trg1 = new TFile("external/SFs_Leg1_Ele23_HZZSelection_Tag35.root");
	TFile *fe_trg2 = new TFile("external/SFs_NoRefit.root");
	TFile *fe_dZ = new TFile("external/eleDz_Ele23_Ele12_vsAbsEta_NonDZ_DZ_Moriond17_SFs.root");
	TFile *fe_LowGSF = new TFile("external/eleReco_HZZ_Moriond17_SFs.root");
	TFile *fe_GSF = new TFile("external/eleReco_Moriond17_SFs.root");
	TFile *fe_LowID  = new TFile("external/eleMVA_HZZ_Moriond17_SFs.root");
	TFile *fg_ID  = new TFile("external/phoMVA_WP90_Moriond17_SFs.root");
	TFile *fe_ID  = new TFile("external/eleMVA_WP90_Moriond17_SFs.root");
	TFile *fm_trg1_0to09 = new TFile("external/sf_Mu8Leg_Eta0to09.root");
	TFile *fm_trg1_09to12 = new TFile("external/sf_Mu8Leg_Eta09to12.root");
	TFile *fm_trg1_12to21  = new TFile("external/sf_Mu8Leg_Eta12to21.root");
	TFile *fm_trg1_21to24 = new TFile("external/sf_Mu8Leg_Eta21to24.root");
	TFile *fm_trg2_0to09 = new TFile("external/sf_Mu17Leg_Eta0to09.root");
	TFile *fm_trg2_09to12 = new TFile("external/sf_Mu17Leg_Eta09to12.root");
	TFile *fm_trg2_12to21  = new TFile("external/sf_Mu17Leg_Eta12to21.root");
	TFile *fm_trg2_21to24 = new TFile("external/sf_Mu17Leg_Eta21to24.root");
	TFile *fm_HZZ = new TFile("external/muSelectionAndRecoSF_HZZ_Moriond17.root");
	
	TH2F  *he_GSF = (TH2F*) fe_GSF->Get("EGamma_SF2D");
	TH2F  *he_LowGSF = (TH2F*) fe_LowGSF->Get("EGamma_SF2D");
	TH2F  *he_LowID  = (TH2F*) fe_LowID->Get("EGamma_SF2D");
	TH2F  *he_ID  = (TH2F*) fe_ID->Get("EGamma_SF2D");
	TH2F  *he_trg1 = (TH2F*) fe_trg1->Get("EGamma_SF2D");
	TH2F  *he_trg2 = (TH2F*) fe_trg2->Get("EGamma_SF2D");
	TH2F  *he_dZ = (TH2F*) fe_dZ->Get("hEta1_Eta2_SF");
	TH2F  *hg_ID  = (TH2F*) fg_ID->Get("EGamma_SF2D");
	TH1F  *hm_trg1_0to09  = (TH1F*) fm_trg1_0to09->Get("scale_factor");
	TH1F  *hm_trg2_0to09  = (TH1F*) fm_trg2_0to09->Get("scale_factor");
	TH1F  *hm_trg1_09to12  = (TH1F*) fm_trg1_09to12->Get("scale_factor");
	TH1F  *hm_trg2_09to12  = (TH1F*) fm_trg2_09to12->Get("scale_factor");
	TH1F  *hm_trg1_12to21  = (TH1F*) fm_trg1_12to21->Get("scale_factor");
	TH1F  *hm_trg2_12to21  = (TH1F*) fm_trg2_12to21->Get("scale_factor");
	TH1F  *hm_trg1_21to24  = (TH1F*) fm_trg1_21to24->Get("scale_factor");
	TH1F  *hm_trg2_21to24  = (TH1F*) fm_trg2_21to24->Get("scale_factor");
	TH2F  *hm_HZZSF = (TH2F*) fm_HZZ->Get("FINAL");
	TH2F  *hm_HZZSFErr = (TH2F*) fm_HZZ->Get("ERROR");

	TFile *fg_R9 = new TFile("external/R9rewei_Moriond17_AfterPreApr_v1.root");
	TGraph *gg_R9EE = (TGraph*) fg_R9->Get("transffull5x5R9EE");
	TGraph *gg_R9EB = (TGraph*) fg_R9->Get("transffull5x5R9EB"); 
	
	
	//B-tag Scale Factors
	BTagCalibration calib("csvv2", "external/CSVv2_Moriond17_B_H.csv");
	
	std::string measType = "comb";
	
	//HEAVY FLAVOR : b
	BTagCalibrationReader reader_loose_b(BTagEntry::OP_LOOSE, "central"); 
	BTagCalibrationReader reader_loose_b_up(BTagEntry::OP_LOOSE, "up"); 
	BTagCalibrationReader reader_loose_b_down(BTagEntry::OP_LOOSE, "down"); 

	BTagCalibrationReader reader_medium_b(BTagEntry::OP_MEDIUM, "central"); 
	BTagCalibrationReader reader_medium_b_up(BTagEntry::OP_MEDIUM, "up"); 
	BTagCalibrationReader reader_medium_b_down(BTagEntry::OP_MEDIUM, "down"); 

	BTagCalibrationReader reader_tight_b(BTagEntry::OP_TIGHT, "central"); 
	BTagCalibrationReader reader_tight_b_up(BTagEntry::OP_TIGHT, "up"); 
	BTagCalibrationReader reader_tight_b_down(BTagEntry::OP_TIGHT, "down"); 
	
	//C 
	BTagCalibrationReader reader_loose_c(BTagEntry::OP_LOOSE, "central"); 
	BTagCalibrationReader reader_loose_c_up(BTagEntry::OP_LOOSE, "up"); 
	BTagCalibrationReader reader_loose_c_down(BTagEntry::OP_LOOSE, "down"); 

	BTagCalibrationReader reader_medium_c(BTagEntry::OP_MEDIUM, "central"); 
	BTagCalibrationReader reader_medium_c_up(BTagEntry::OP_MEDIUM, "up"); 
	BTagCalibrationReader reader_medium_c_down(BTagEntry::OP_MEDIUM, "down"); 

	BTagCalibrationReader reader_tight_c(BTagEntry::OP_TIGHT, "central"); 
	BTagCalibrationReader reader_tight_c_up(BTagEntry::OP_TIGHT, "up"); 
	BTagCalibrationReader reader_tight_c_down(BTagEntry::OP_TIGHT, "down"); 
	
	//LIGHT FLAVOR
	BTagCalibrationReader reader_loose_l(BTagEntry::OP_LOOSE, "central"); 
	BTagCalibrationReader reader_loose_l_up(BTagEntry::OP_LOOSE, "up"); 
	BTagCalibrationReader reader_loose_l_down(BTagEntry::OP_LOOSE, "down"); 

	BTagCalibrationReader reader_medium_l(BTagEntry::OP_MEDIUM, "central"); 
	BTagCalibrationReader reader_medium_l_up(BTagEntry::OP_MEDIUM, "up"); 
	BTagCalibrationReader reader_medium_l_down(BTagEntry::OP_MEDIUM, "down"); 

	BTagCalibrationReader reader_tight_l(BTagEntry::OP_TIGHT, "central"); 
	BTagCalibrationReader reader_tight_l_up(BTagEntry::OP_TIGHT, "up"); 
	BTagCalibrationReader reader_tight_l_down(BTagEntry::OP_TIGHT, "down"); 
	
	//READER - CENTRAL
	reader_loose_b.load(calib, BTagEntry::FLAV_B, measType);
	reader_loose_c.load(calib, BTagEntry::FLAV_C, measType);
	reader_loose_l.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_medium_b.load(calib, BTagEntry::FLAV_B, measType);
	reader_medium_c.load(calib, BTagEntry::FLAV_C, measType);
	reader_medium_l.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_tight_b.load(calib, BTagEntry::FLAV_B, measType);
	reader_tight_c.load(calib, BTagEntry::FLAV_C, measType);
	reader_tight_l.load(calib, BTagEntry::FLAV_UDSG, measType);

	//READER - UP
	reader_loose_b_up.load(calib, BTagEntry::FLAV_B, measType);
	reader_loose_c_up.load(calib, BTagEntry::FLAV_C, measType);
	reader_loose_l_up.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_medium_b_up.load(calib, BTagEntry::FLAV_B, measType);
	reader_medium_c_up.load(calib, BTagEntry::FLAV_C, measType);
	reader_medium_l_up.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_tight_b_up.load(calib, BTagEntry::FLAV_B, measType);
	reader_tight_c_up.load(calib, BTagEntry::FLAV_C, measType);
	reader_tight_l_up.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	//READER - DOWN
	reader_loose_b_down.load(calib, BTagEntry::FLAV_B, measType);
	reader_loose_c_down.load(calib, BTagEntry::FLAV_C, measType);
	reader_loose_l_down.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_medium_b_down.load(calib, BTagEntry::FLAV_B, measType);
	reader_medium_c_down.load(calib, BTagEntry::FLAV_C, measType);
	reader_medium_l_down.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	reader_tight_b_down.load(calib, BTagEntry::FLAV_B, measType);
	reader_tight_c_down.load(calib, BTagEntry::FLAV_C, measType);
	reader_tight_l_down.load(calib, BTagEntry::FLAV_UDSG, measType);
	
	//EVENT COUNTERS
	Int_t TotalEvent = 0;
	Int_t PV = 0;
	Int_t Trigger = 0;
	Int_t GloSel = 0;
	Int_t Ele = 0;
	Int_t Mu = 0;
	Int_t Pho = 0;
	Int_t mDilepton = 0;
	Int_t Dijet = 0;
	Int_t BTaggedJets = 0;
	Int_t mHiggs = 0;
	Int_t mBJetHP = 0;
	Int_t mBJetMP = 0;
	Int_t mBJetLP = 0;
	
	//YIELD COUNTERS
	long double  filter = 0, trigger = 0, globsel = 0, electron = 0, muon = 0, dilepton = 0, dijet = 0, btagjets = 0, photon = 0;
	long double higgsM = 0, bjetLP = 0, bjetMP = 0, bjetHP = 0;

	RoccoR  rc("../plugin/rcdata.2016.v3");
	
	//VARIABLES
	Int_t run, Lumi, nVtx, nMC;
	Int_t category, totalEvents;
	Long64_t event; 
	
	Float_t puWei = 1.; 
	Float_t puWeiUp = 1.;
	Float_t puWeiDown = 1.;
	Float_t mcWei = 1.;
	Float_t genWei = 1.;
	Float_t TotalWei = 1.;
	
	//LEPTONS
	Float_t llepPt, tlepPt, llepEta, tlepEta, llepPhi, tlepPhi, llepSCEta, tlepSCEta, llepSCPhi, tlepSCPhi, llepR9, tlepR9, llepPhoMVA, tlepPhoMVA;
	Float_t Mll, mcMll;
	Int_t llepCharge, tlepCharge;
	
	//JETS
	Float_t jet1Pt, jet2Pt, jet1Eta, jet2Eta, jet1Phi, jet2Phi, jet1En, jet2En;
	Int_t jet1BTagType, jet2BTagType;
	Float_t jet1BTagDiscr, jet2BTagDiscr;
	Float_t Mjj;
	
	//PHOTON
	Float_t PhoEt, PhoEta, PhoPhi, PhoSCEta, PhoR9, PhoMVA, PhoSSMVA, PhoMass, PhoPFChIso;

	//THREE-BODY MASS
	Float_t mH;
	
	//DELTA-Rs
	Float_t dRjj, dEtajj, dRlj, dRljMin, dRljMax, dRjg, dRjgMax, dRjgMin, dRlg, dRlgMax, dRlgMin;

	//TLORENTZVECTORS
	TLorentzVector Zll, Lepton[2], bJet[2], Photon, Higgs;
	
	//SCALING FACTORS
	//Nominal 
	Double_t TotalSF, EleTrgSF, MuTrgSF, EtaSF, EleVetoSF, PhotonSF, CSVSF, EleIDSF, HZZMuIDSF, GSF;
	//Up
	Double_t EleTrgSFUp, MuTrgSFUp, EtaSFUp, EleVetoSFUp, PhotonSFUp, CSVSFUp, EleIDSFUp, HZZMuIDSFUp, GSFUp;
	//Down
	Double_t EleTrgSFDown, MuTrgSFDown, EtaSFDown, EleVetoSFDown, PhotonSFDown, CSVSFDown, EleIDSFDown, HZZMuIDSFDown, GSFDown;	
	//SF Uncertainties
	Double_t UnEleTrgSF, UnMuTrgSF, UnEtaSF, UnEleVetoSF, UnPhotonSF, UnEleIDSF, UnHZZMuIDSF, UnGSF;
	
	
	TTree *tree = new TTree("outTree", "minitree for ZH->ll+Zg->ll+jjg");
	
	tree->Branch("totalEvents", &totalEvents, "totalEvents/I");
	tree->Branch("run", &run, "run/I");
	tree->Branch("Lumi", &Lumi, "Lumi/I");
	tree->Branch("event", &event, "event/L");
	tree->Branch("nMC", &nMC, "nMC/I");
	tree->Branch("nVtx", &nVtx, "nVtx/I");
	tree->Branch("category", &category, "category/I");
	tree->Branch("puWei", &puWei, "puWei/F");	
	tree->Branch("puWeiUp", &puWeiUp, "puWeiUp/F");	
	tree->Branch("puWeiDown", &puWeiDown, "puWeiDown/F");	
	tree->Branch("mcWei", &mcWei, "mcWei/F");
	tree->Branch("genWei", &genWei, "genWei/F");
	tree->Branch("TotalWei", &TotalWei, "TotalWei/F");	
	
	//LEPTONS
	tree->Branch("llepPt", &llepPt, "llepPt/F");
	tree->Branch("tlepPt", &tlepPt, "tlepPt/F");
	tree->Branch("llepEta", &llepEta, "llepEta/F");	
	tree->Branch("tlepEta", &tlepEta, "tlepEta/F");
	tree->Branch("llepPhi", &llepPhi, "llepPhi/F");
	tree->Branch("tlepPhi", &tlepPhi, "tlepPhi/F");
	tree->Branch("llepSCEta", &llepSCEta, "llepSCEta/F");	
	tree->Branch("tlepSCEta", &tlepSCEta, "tlepSCEta/F");
	tree->Branch("llepSCPhi", &llepSCPhi, "llepSCPhi/F");
	tree->Branch("tlepSCPhi", &tlepSCPhi, "tlepSCPhi/F");
	tree->Branch("llepR9", &llepR9, "llepR9/F");	
	tree->Branch("tlepR9", &tlepR9, "tlepR9/F");
	tree->Branch("llepPhoMVA", &llepPhoMVA, "llepPhoMVA/F");	
	tree->Branch("tlepPhoMVA", &tlepPhoMVA, "tlepPhoMVA/F");
	tree->Branch("Mll", &Mll, "Mll/F");	
	tree->Branch("mcMll", &mcMll, "mcMll/F");
	tree->Branch("llepCharge", &llepCharge, "llepCharge/I");
	tree->Branch("tlepCharge", &tlepCharge, "tlepCharge/I");

	//JETS
	tree->Branch("jet1Pt", &jet1Pt, "jet1Pt/F");
	tree->Branch("jet2Pt", &jet2Pt, "jet2Pt/F");
	tree->Branch("jet1Eta", &jet1Eta, "jet1Eta/F");
	tree->Branch("jet2Eta", &jet2Eta, "jet2Eta/F");
	tree->Branch("jet1Phi", &jet1Phi, "jet1Phi/F");
	tree->Branch("jet2Phi", &jet2Phi, "jet2Phi/F");
	tree->Branch("jet1En", &jet1En, "jet1En/F");
	tree->Branch("jet2En", &jet2En, "jet2En/F");
	tree->Branch("jet1BTagType", &jet1BTagType, "jet1BTagType/I");
	tree->Branch("jet2BTagType", &jet2BTagType, "jet2BTagType/I");
	tree->Branch("jet1BTagDiscr", &jet1BTagDiscr, "jet1BTagDiscr/F");
	tree->Branch("jet2BTagDiscr", &jet2BTagDiscr, "jet2BTagDiscr/F");
	tree->Branch("Mjj", &Mjj, "Mjj/F");
	
	//PHOTON
	tree->Branch("PhoEt", &PhoEt, "PhoEt/F");
	tree->Branch("PhoEta", &PhoEta, "PhoEta/F");
	tree->Branch("PhoSCEta", &PhoSCEta, "PhoSCEta/F");
	tree->Branch("PhoPhi", &PhoPhi, "PhoPhi/F");
	tree->Branch("PhoMass", &PhoMass, "PhoMass/F");
	tree->Branch("PhoR9", &PhoR9, "PhoR9/F");
	tree->Branch("PhoMVA", &PhoMVA, "PhoMVA/F");
	tree->Branch("PhoSSMVA", &PhoSSMVA, "PhoSSMVA/F");
	tree->Branch("PhoPFChIso", &PhoPFChIso, "PhoPFChIso/F");
	
	//THREE-BODY MASS
	tree->Branch("mH", &mH, "mH/F");
	
	//DELTA-Rs
	tree->Branch("dRjj", &dRjj, "dRjj/F");
	tree->Branch("dEtajj", &dEtajj, "dEtajj/F");
	tree->Branch("dRlj", &dRlj, "dRlj/F");
	tree->Branch("dRlg", &dRlg, "dRlg/F");
	tree->Branch("dRjg", &dRjg, "dRjg/F");
	tree->Branch("dRljMin", &dRljMin, "dRljMin/F");
	tree->Branch("dRljMax", &dRljMax, "dRljMax/F");
	tree->Branch("dRjgMin", &dRjgMin, "dRjgMin/F");
	tree->Branch("dRjgMax", &dRjgMax, "dRjgMax/F");
	tree->Branch("dRlgMin", &dRlgMin, "dRlgMin/F");
	tree->Branch("dRlgMax", &dRlgMax, "dRlgMax/F");

	//SCALING FACTORS
	//Nominal
	tree->Branch("TotalSF", &TotalSF, "TotalSF/D");
	tree->Branch("EleTrgSF", &EleTrgSF, "EleTrgSF/D");
	tree->Branch("MuTrgSF", &MuTrgSF, "MuTrgSF/D");
	tree->Branch("EtaSF", &EtaSF, "EtaSF/D");
	tree->Branch("EleVetoSF", &EleVetoSF, "EleVetoSF/D");
	tree->Branch("PhotonSF", &PhotonSF, "PhotonSF/D");	
	tree->Branch("CSVSF", &CSVSF, "CSVSF/D");	
	tree->Branch("EleIDSF", &EleIDSF, "EleIDSF/D");	
	tree->Branch("HZZMuIDSF", &HZZMuIDSF, "HZZMuIDSF/D");	
	tree->Branch("GSF", &GSF, "GSF/D");
	//Up
	tree->Branch("EleTrgSFUp", &EleTrgSFUp, "EleTrgSFUp/D");	
	tree->Branch("MuTrgSFUp", &MuTrgSFUp, "MuTrgSFUp/D");	
	tree->Branch("EtaSFUp", &EtaSFUp, "EtaSFUp/D");	
	tree->Branch("EleVetoSFUp", &EleVetoSFUp, "EleVetoSFUp/D");	
	tree->Branch("PhotonSFUp", &PhotonSFUp, "PhotonSFUp/D");	
	tree->Branch("CSVSFUp", &CSVSFUp, "CSVSFUp/D");	
	tree->Branch("EleIDSFUp", &EleIDSFUp, "EleIDSFUp/D");	
	tree->Branch("HZZMuIDSFUp", &HZZMuIDSFUp, "HZZMuIDSFUp/D");		
	tree->Branch("GSFUp", &GSFUp, "GSFUp/D");	
	//Down
	tree->Branch("EleTrgSFDown", &EleTrgSFDown, "EleTrgSFDown/D");	
	tree->Branch("MuTrgSFDown", &MuTrgSFDown, "MuTrgSFDown/D");	
	tree->Branch("EtaSFDown", &EtaSFDown, "EtaSFDown/D");	
	tree->Branch("EleVetoSFDown", &EleVetoSFDown, "EleVetoSFDown/D");	
	tree->Branch("PhotonSFDown", &PhotonSFDown, "PhotonSFDown/D");	
	tree->Branch("CSVSFDown", &CSVSFDown, "CSVSFDown/D");	
	tree->Branch("EleIDSFDown", &EleIDSFDown, "EleIDSFDown/D");	
	tree->Branch("HZZMuIDSFDown", &HZZMuIDSFDown, "HZZMuIDSFDown/D");	
	tree->Branch("GSFDown", &GSFDown, "GSFDown/D");
	//Uncertainties
	tree->Branch("UnEleTrgSF", &UnEleTrgSF, "UnEleTrgSF/D");
	tree->Branch("UnMuTrgSF", &UnMuTrgSF, "UnMuTrgSF/D");
	tree->Branch("UnEtaSF", &UnEtaSF, "UnEtaSF/D");
	tree->Branch("UnEleVetoSF", &UnEleVetoSF, "UnEleVetoSF/D");
	tree->Branch("UnPhotonSF", &UnPhotonSF, "UnPhotonSF/D");
	tree->Branch("UnEleIDSF", &UnEleIDSF, "UnEleIDSF/D");	
	tree->Branch("UnHZZMuIDSF", &UnHZZMuIDSF, "UnHZZMuIDSF/D");
	tree->Branch("UnGSF", &UnGSF, "UnGSF/D");	
	
	//PILEUP REWEIGHTING
	PUWeightCalculator puCalc;
	PUWeightCalculator puCalcUp;
	PUWeightCalculator puCalcDown;
	
	if (data.HasMC()){
      puCalc.Init("80X_puwei/36p0_invfb/summer16/PU_histo_13TeV_GoldenJSON_69200nb.root");
      puCalcUp.Init("80X_puwei/36p0_invfb/summer16/PU_histo_13TeV_GoldenJSON_71300nb.root");
      puCalcDown.Init("80X_puwei/36p0_invfb/summer16/PU_histo_13TeV_GoldenJSON_69000nb.root");
    }

	if (aMCatNLO == 1){
		Int_t totalEvents = 0;
		for (Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){
		 data.GetEntry(ev);
		 if (data.HasMC()){
			  float genWeight = data.GetFloat("genWeight");
			  if (genWeight > 0) totalEvents++; else totalEvents--;
			}
        }
      mcWei = (totalEvents != 0) ? (xs*lumi/totalEvents) : 1.;
	  cout << totalEvents << endl;
	}
	else mcWei = xs*lumi/data.GetEntriesFast();
	cout << xs << " " << lumi << " " << endl;
	
	for(Long64_t ev = 0; ev < data.GetEntriesFast(); ev++){
		if (ev % 5000000 == 0){fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());}
		
	HZgCand HZgCand = {};
	HZgCand.ev = ev;
	
	totalEvents = data.GetEntriesFast();
	run   = data.GetInt("run");
	Lumi = data.GetInt("lumis");
	event = data.GetLong64("event");
	nVtx  = data.GetInt("nVtx");

	TotalEvent = data.GetEntriesFast();
	data.GetEntry(ev);

	//VERTEX FILTER AND TRIGGER
	Bool_t isPVGood = data.GetBool("isPVGood");
	if (isPVGood == false) continue;
	PV++;
	filter 		+= mcWei;

	Long64_t HLTEleMuX = (Long64_t) data.GetLong64("HLTEleMuX");
	ULong64_t   HLTEleMuXIsPrescaled  = (ULong64_t)data.GetLong64("HLTEleMuXIsPrescaled");

	if(channel==0){
		if( (HLTEleMuX >>  40 & 1) == 0 && (HLTEleMuX >>  5 & 1) == 0) continue;
	}
	if (channel == 1){
		if ((HLTEleMuX >> 14 & 1) == 0 && (HLTEleMuX >> 15 & 1) == 0 && ((HLTEleMuX >> 41) & 1) == 0 && ((HLTEleMuX >> 42) & 1) == 0) continue;
	}
	 Trigger++;
	 trigger	+= mcWei;
	
	Int_t  nMC      		= 0;
	Int_t* mcPID    		= NULL;
	Int_t* mcMomPID 		= NULL;
	Int_t* mcGMomPID 		= NULL;
	Int_t* jetPartonID 		= NULL;
	float* mcMass   		= NULL;
	float* mcPt     		= NULL;
	float* mcEta    		= NULL;
	float* mcPhi    		= NULL;
	float* mcE	  			= NULL;
	UShort_t* mcStatusFlag 	= NULL;
	Float_t* jetP4Smear		= NULL;
	Float_t* jetP4SmearUp	= NULL;
	Float_t* jetP4SmearDo	= NULL;

	if (data.HasMC()){
		nMC          	= data.GetInt("nMC");
		mcPID        	= data.GetPtrInt("mcPID");
		mcMomPID     	= data.GetPtrInt("mcMomPID");
		mcGMomPID    	= data.GetPtrInt("mcGMomPID");
		mcMass       	= data.GetPtrFloat("mcMass");
		mcPt         	= data.GetPtrFloat("mcPt");
		mcEta        	= data.GetPtrFloat("mcEta");
		mcPhi        	= data.GetPtrFloat("mcPhi");
		mcE		   		= data.GetPtrFloat("mcE");
		mcStatusFlag 	= (UShort_t*) data.GetPtrShort("mcStatusFlag");
		jetPartonID 	= data.GetPtrInt("jetPartonID"); 
		jetP4Smear 		= data.GetPtrFloat("jetP4Smear");
		jetP4SmearUp	= data.GetPtrFloat("jetP4SmearUp");
		jetP4SmearDo 	= data.GetPtrFloat("jetP4SmearDo");
	}
	
	Int_t  nPho     = data.GetInt("nPho");	
	//float* phoEt    = data.GetPtrFloat("phoEt");
	float* phoEt      = data.GetPtrFloat("phoCalibEt"); 
	float* phoEta     = data.GetPtrFloat("phoEta");
	float* phoPhi     = data.GetPtrFloat("phoPhi");
	float* phoSCEta   = data.GetPtrFloat("phoSCEta");
	float* phoSCPhi   = data.GetPtrFloat("phoSCPhi");
	float* phoR9      = data.GetPtrFloat("phoR9Full5x5");
	float* phoIDMVA   = data.GetPtrFloat("phoIDMVA");
	float* phoPFChIso = data.GetPtrFloat("phoPFChIso");
	
	Int_t  nJet     = data.GetInt("nJet");
	Float_t* jetPt_ = data.GetPtrFloat("jetPt");
	Float_t* jetEta = data.GetPtrFloat("jetEta");
	Float_t* jetPhi = data.GetPtrFloat("jetPhi");
	Float_t* jetEn_ = data.GetPtrFloat("jetEn");
	Float_t* jetBtag = data.GetPtrFloat("jetCSV2BJetTags");

	float* elePt     = data.GetPtrFloat("eleCalibPt");
	float* eleEta    = data.GetPtrFloat("eleEta");
	float* elePhi    = data.GetPtrFloat("elePhi");
	float* eleSCEta  = data.GetPtrFloat("eleSCEta");//use for categorization
	float* eleSCPhi  = data.GetPtrFloat("eleSCPhi");//use for categorization
	float* eleR9     = data.GetPtrFloat("eleR9");//use for categorization
	Int_t* eleCharge = data.GetPtrInt("eleCharge");
	UShort_t* eleIDbit = (UShort_t*) data.GetPtrShort("eleIDbit");

	Int_t  nMu         = data.GetInt("nMu");
	float* muPt_	   = data.GetPtrFloat("muPt"); 
	float* muEta       = data.GetPtrFloat("muEta");
	float* muPhi       = data.GetPtrFloat("muPhi");
	Int_t* muCharge    = data.GetPtrInt("muCharge");//for roch corr
	Int_t* muType      = data.GetPtrInt("muType");//for roch corr
	Int_t* muTrkLayers = data.GetPtrInt("muTrkLayers");//for roch corr
	
	//GLOBAL SELECTIONS
	if (nPho < 1 || nJet < 2) continue;
	GloSel++;
	globsel		+= mcWei;
	
	//LEPTONS
	vector<int> eleID;
	vector<int> ZeeID;
	vector<int> muID;
	vector<int> Zmm;
	vector<float> mucorrPt; mucorrPt.clear();
	vector<float> muPt; muPt.clear();
	Int_t nLep = 0;
	Int_t nmcLep = 0;

	//ROCHESTER CORRECTION
	Float_t corrPt;
	for (Int_t i = 0; i < nMu ; i++){
		Float_t rand1, rand2;
		TRandom g;
		rand1 = g.Uniform(0,1);
		rand2 = g.Uniform(0,1);

		double SF = 0;
		vector<int> imcMu;
		int countmc = 0;
		if (data.HasMC()) mcMu(data, imcMu);
		TLorentzVector rochMu;
		rochMu.SetPtEtaPhiM(muPt_[i], muEta[i], muPhi[i], 0.1057);
		TLorentzVector mcMu;
	
		if (data.HasMC()){
			//mc objects
			Int_t    nMC = data.GetInt("nMC");
			Float_t* mcPt = data.GetPtrFloat("mcPt");
			Float_t* mcEta = data.GetPtrFloat("mcEta");
			Float_t* mcPhi = data.GetPtrFloat("mcPhi");
			
			for (unsigned int mc = 0; mc < imcMu.size(); mc++){
				if (countmc > 0) break;
				TLorentzVector mcMutemp;
				mcMutemp.SetPtEtaPhiM(mcPt[imcMu[mc]], mcEta[imcMu[mc]], mcPhi[imcMu[mc]], 0.1057);
				Float_t dRmc = mcMutemp.DeltaR(rochMu);
				if (dRmc > 0.1) continue;
				mcMu = mcMutemp;
				countmc++; 	
			}
				if (countmc > 0) SF = rc.kScaleFromGenMC(muCharge[i], muPt_[i], muEta[i], muPhi[i], muTrkLayers[i], mcMu.Pt(), rand1, 0, 0);
				if (countmc == 0) SF = rc.kScaleAndSmearMC(muCharge[i], muPt_[i], muEta[i], muPhi[i], muTrkLayers[i], rand1, rand2, 0, 0);
		}
	
		else SF = rc.kScaleDT(muCharge[i], muPt_[i], muEta[i], muPhi[i], 0, 0); 
				muPt_[i] = muPt_[i]*SF;
				mucorrPt.push_back(muPt_[i]);
	}
		muPt.assign(mucorrPt.begin(), mucorrPt.end()); 	  

	vector<TLorentzVector> multz; multz.clear();
	ElectronID16NonTrgMVA(data, 98, ZeeID, multz); //ELECTRON MVA ID 
	HZZmuID16forZg(data, muID, Zmm, muPt, channel, multz); //HZZ ID FOR MUON CHANNEL

	//ELECTRON LOOP
	if(channel == 0){
		if (ZeeID.size() < 2) continue;
		if (elePt[ZeeID[0]] < 25.) continue;
		if (elePt[ZeeID[1]] < 15.) continue;

		llepSCEta 		= eleSCEta[ZeeID[0]];
		llepR9    		= eleR9[ZeeID[0]];
		llepCharge    	= eleCharge[ZeeID[0]];
		llepEta 		= Lepton[0].Eta();
		llepPt 			= Lepton[0].Pt();
		llepPhi			= Lepton[0].Phi();

		tlepSCEta		= eleSCEta[ZeeID[1]];
		tlepR9    		= eleR9[ZeeID[1]];
		tlepCharge    	= eleCharge[ZeeID[1]];
		tlepEta 		= Lepton[1].Eta();
		tlepPt 			= Lepton[1].Pt();
		tlepPhi 		= Lepton[1].Phi();

		Lepton[0].SetPtEtaPhiM(elePt[ZeeID[0]], eleEta[ZeeID[0]], elePhi[ZeeID[0]], 0.511*0.001);
		Lepton[1].SetPtEtaPhiM(elePt[ZeeID[1]], eleEta[ZeeID[1]], elePhi[ZeeID[1]], 0.511*0.001);

		Ele++;
		electron	+= mcWei;
	}
	
	//MUON LOOP
	else if (channel == 1){ 
		if (Zmm.size() < 2) continue;
		if (muPt[Zmm[0]] < 20.) continue;
		if (muPt[Zmm[1]] < 10.) continue;

		llepPt  	= Lepton[0].Pt();
		llepEta  	= Lepton[0].Eta();
		llepPhi 	= Lepton[0].Phi();
		llepCharge 	= muCharge[Zmm[0]];
		
		tlepPt  	= Lepton[1].Pt();
		tlepEta  	= Lepton[1].Eta();
		tlepPhi 	= Lepton[1].Phi();
		tlepCharge 	= muCharge[Zmm[1]];
		
		Lepton[0].SetPtEtaPhiM(muPt[Zmm[0]], muEta[Zmm[0]], muPhi[Zmm[0]], 105.7*0.001);
		Lepton[1].SetPtEtaPhiM(muPt[Zmm[1]], muEta[Zmm[1]], muPhi[Zmm[1]], 105.7*0.001);
	  
		Mu++;
		muon		+= mcWei;
    }
	
	//DILEPTON INVARIANT MASS
	Zll	= Lepton[0] + Lepton[1];
	if (Zll.M() < 50.) continue;
	Mll = Zll.M();
	mDilepton++;
	dilepton	+= mcWei;

	//JET SELECTIONS
	vector<int> jetList; jetList.clear(); 
	vector<int> jet1List; jet1List.clear();
	vector<int> jet2List; jet2List.clear();
	vector<int> acc_jetpair; acc_jetpair.clear();
	
	vector<float> jetPt, jetEn; jetPt.clear(); jetEn.clear();
	
	//SMEARING CORRECTIONS
	for (int i = 0; i < nJet; i++){
		float jetP4Smear_ = 1.;
		if(data.HasMC()) jetP4Smear_ = jetP4SmearDo[i];
		jetPt.push_back(jetPt_[i]*jetP4Smear_);
		jetEn.push_back(jetEn_[i]*jetP4Smear_);
	}

	for (int n = 0; n < nJet; n++){
		if (!JetPFLooseId(data, n)) continue;
		if (fabs(jetEta[n]) > 2.4) continue;
		if (deltaR(jetEta[n], jetPhi[n], Lepton[0].Eta(), Lepton[0].Phi()) < 0.4) continue; 
		if (deltaR(jetEta[n], jetPhi[n], Lepton[1].Eta(), Lepton[1].Phi()) < 0.4) continue;
		jetList.push_back(n); 
	}
	if (jetList.size() < 2) continue;

	for (size_t i=0;i<jetList.size();i++){
		for(size_t j=i+1;j<jetList.size();j++){
			jet1List.push_back(jetList[i]);
			jet2List.push_back(jetList[j]);
		}
	}
	for (size_t i=0;i<jet1List.size();i++){
		if(jetPt[jet1List[i]]<jetPt[jet2List[i]]){
				int tmp; tmp = jet1List[i];
				jet1List[i] = jet2List[i];
				jet2List[i] = tmp;
			}
		}
		
	vector<float> jjDR;	jjDR.clear(); 	
	vector<float> bjRegfactor[2]; 
	bjRegfactor[0].clear(); 
	bjRegfactor[1].clear();
	for (size_t n=0;n<jet1List.size();n++) jjDR.push_back(deltaR(jetEta[jet1List[n]],jetPhi[jet1List[n]],jetEta[jet2List[n]],jetPhi[jet2List[n]]));
	for (size_t n=0;n<jet1List.size();n++){
		if		(DoRegression==0) {bjRegfactor[0].push_back(1);} //NoReg
		else if (DoRegression==1) {bjRegfactor[0].push_back(TMVA_15var_jetGenJet_nu_leading_3_1(data,jet1List[n],jjDR[n]));} //Official 15var separate
		else if (DoRegression==2) {bjRegfactor[0].push_back(TMVA_15plus3_jetGenJet_nu_leading_2_27(data,jet1List[n],jjDR[n]));} //15+3 (separate)
		else if (DoRegression==3) {bjRegfactor[0].push_back(TMVA_15var_jetGenJet_nu_2_27(data,jet1List[n]));} //Offical 15var 	
	}
	for (size_t n=0;n<jet2List.size();n++){
		if		(DoRegression==0) {bjRegfactor[1].push_back(1);} //NoReg
		else if (DoRegression==1) {bjRegfactor[1].push_back(TMVA_15var_jetGenJet_nu_trailing_3_1(data,jet2List[n],jjDR[n]));} //Official 15var separate
		else if (DoRegression==2) {bjRegfactor[1].push_back(TMVA_15plus3_jetGenJet_nu_trailing_2_27(data,jet2List[n],jjDR[n]));} //15+3 (separate)
		else if (DoRegression==3) {bjRegfactor[1].push_back(TMVA_15var_jetGenJet_nu_2_27(data,jet2List[n]));} //Offical 15var 
	}
	for (size_t n=0; n<jet1List.size(); n++){	
		if (jetPt[jet1List[n]]*bjRegfactor[0][n] < 25.) continue; 
		if (jetPt[jet2List[n]]*bjRegfactor[1][n] < 25.) continue;
		
		TLorentzVector jet[2];
		jet[0].SetPtEtaPhiE(jetPt[jet1List[n]]*bjRegfactor[0][n], jetEta[jet1List[n]], jetPhi[jet1List[n]], jetEn[jet1List[n]]*bjRegfactor[0][n]);
		jet[1].SetPtEtaPhiE(jetPt[jet2List[n]]*bjRegfactor[1][n], jetEta[jet2List[n]], jetPhi[jet2List[n]], jetEn[jet2List[n]]*bjRegfactor[1][n]);
		TLorentzVector jjM = jet[0] + jet[1];
		if (jjM.M() < 50. || jjM.M() > 150.) continue;  
		acc_jetpair.push_back(n);
	}
	if (acc_jetpair.size() < 1) continue;
	Dijet++;
	dijet		+= mcWei;
	
	//Highest sum of b-tagging scores
	float jetSumbTag= 0;
	int sel_jetPair = -1;
	for(size_t nPair=0; nPair<acc_jetpair.size(); nPair++){
		float jetSumbTag_tmp = jetBtag[jet1List[acc_jetpair[nPair]]] + jetBtag[jet2List[acc_jetpair[nPair]]];
		if (jetSumbTag < jetSumbTag_tmp){
			jetSumbTag = jetSumbTag_tmp;
			sel_jetPair = acc_jetpair[nPair];
		}
	}
	if (sel_jetPair==-1) continue; 
	
	vector<float> bjRegfactor_; bjRegfactor_.clear();
	bjRegfactor_.push_back(bjRegfactor[0][sel_jetPair]);
	bjRegfactor_.push_back(bjRegfactor[1][sel_jetPair]);
	HZgCand.reg_fac = bjRegfactor_;
	HZgCand.bjet[0] = jet1List[sel_jetPair];
	HZgCand.bjet[1] = jet2List[sel_jetPair];

	TLorentzVector bJet[2];
	bJet[0].SetPtEtaPhiE(jetPt[HZgCand.bjet[0]]*HZgCand.reg_fac[0], jetEta[HZgCand.bjet[0]], jetPhi[HZgCand.bjet[0]], jetEn[HZgCand.bjet[0]]*HZgCand.reg_fac[0]);
	bJet[1].SetPtEtaPhiE(jetPt[HZgCand.bjet[1]]*HZgCand.reg_fac[1], jetEta[HZgCand.bjet[1]], jetPhi[HZgCand.bjet[1]], jetEn[HZgCand.bjet[1]]*HZgCand.reg_fac[1]);
	TLorentzVector diJet = bJet[0] + bJet[1];
	Mjj = diJet.M();
	
	//B-tag selections;
	// 0: Untagged
	// 1: Loose
	// 2: Medium
	// 3: Tight
	int btagCat[2];
	btagCat[0] = JetBtagCSVv2_Cat(data, jet1List[sel_jetPair]);
	btagCat[1] = JetBtagCSVv2_Cat(data, jet2List[sel_jetPair]);
	
	if(btagCat[0]==0 || btagCat[1]==0) continue; //Reject untagged jets!
	BTaggedJets++;
	btagjets	+= mcWei;

	//PHOTON SELECTIONS
	vector<int> phoID;
	PhotonIDMVA2016(data, phoID);
	if (phoID.size() == 0) continue;
	if(phoEt[phoID[0]] < 15.)continue;
	Pho++;
	photon		+= mcWei;

	int phoIndex = -999;
	Int_t nSelPho  = 0;
	bool selpho = false;
	bool haspho = false;
	Int_t isPromptPhoton = 0;
	
	TLorentzVector tmppho;
	for (size_t j = 0; j < phoID.size(); ++j){
		
		if(haspho == true) break;
		if (phoEt[phoID[j]] < 15.) continue;
		
		tmppho.SetPtEtaPhiM(phoEt[phoID[j]], phoEta[phoID[j]], phoPhi[phoID[j]], 0.);
		
		if (tmppho.DeltaR(Lepton[0]) < 0.4 || tmppho.DeltaR(Lepton[1]) < 0.4) continue;
		if (tmppho.DeltaR(bJet[0]) < 0.4 || tmppho.DeltaR(bJet[1]) < 0.4) continue;
		
		Higgs = tmppho + diJet;
		
		if (Higgs.M() > 180.|| Higgs.M() < 100.) continue;
		if (tmppho.Pt()/Higgs.M() < 15./110.) continue;
		
		Photon = tmppho;
		haspho = true;
		phoIndex =  phoID[j];
	} 
	if (phoIndex == -999) continue;
	mHiggs++;
	
	PhoEt      		= Photon.Pt();
	PhoEta     		= Photon.Eta();
	PhoPhi     		= Photon.Phi();
	PhoSCEta   		= phoSCEta[phoIndex];
	PhoR9      		= phoR9[phoIndex];
    PhoMVA     		= phoIDMVA[phoIndex];
	PhoPFChIso 		= phoPFChIso[phoIndex];
	PhoMass	   		= Photon.M();
	mH		   		= Higgs.M();
	dRjj       		= bJet[0].DeltaR(bJet[1]);
	dEtajj     		= fabs(bJet[1].Eta() - bJet[0].Eta());
	dRlj	   		= bJet[0].DeltaR(Lepton[0]);
	dRljMin    		= TMath::Min(bJet[0].DeltaR(Lepton[0]), bJet[1].DeltaR(Lepton[1]));
	dRljMin    		= TMath::Max(bJet[0].DeltaR(Lepton[0]), bJet[1].DeltaR(Lepton[1]));
	dRjg	   		= Photon.DeltaR(bJet[0]);
	dRjgMin    		= TMath::Min(Photon.DeltaR(bJet[0]), Photon.DeltaR(bJet[1]));
	dRjgMax    		= TMath::Max(Photon.DeltaR(bJet[0]), Photon.DeltaR(bJet[1]));
	dRlg	   		= Photon.DeltaR(Lepton[0]);
	dRlgMin	   		= TMath::Min(Photon.DeltaR(Lepton[0]), Photon.DeltaR(Lepton[1]));
	dRlgMax	   		= TMath::Max(Photon.DeltaR(Lepton[0]), Photon.DeltaR(Lepton[1]));
	jet1Pt 	   		= jetPt[HZgCand.bjet[0]]*HZgCand.reg_fac[0];
	jet2Pt 	   		= jetPt[HZgCand.bjet[1]]*HZgCand.reg_fac[1];
	jet1Eta    		= jetEta[HZgCand.bjet[0]];
	jet2Eta    		= jetEta[HZgCand.bjet[1]];
	jet1Phi    		= jetPhi[HZgCand.bjet[0]];
	jet2Phi    		= jetPhi[HZgCand.bjet[1]]; 
	jet1En     		= jetEn[HZgCand.bjet[0]]*HZgCand.reg_fac[0];
	jet2En     		= jetEn[HZgCand.bjet[1]]*HZgCand.reg_fac[1];
	jet1BTagType 	= btagCat[0];
	jet2BTagType 	= btagCat[1];
	jet1BTagDiscr 	= jetBtag[jet1List[sel_jetPair]];
	jet2BTagDiscr 	= jetBtag[jet2List[sel_jetPair]];
	
	if (data.HasMC()){
		  float* puTrue = data.GetPtrFloat("puTrue");
		  puWei     = (float) puCalc.GetWeight(run, puTrue[1]); // in-time PU
		  puWeiUp   = (float) puCalcUp.GetWeight(run, puTrue[1]); // in-time PU
		  puWeiDown = (float) puCalcDown.GetWeight(run, puTrue[1]); // in-time PU
		  float generatorWeight = data.GetFloat("genWeight");
		  if(aMCatNLO == 1) genWei = (generatorWeight > 0) ? 1. : -1.; else genWei = 1.;
    }
	
	//SCALING FACTORS	
	//B-tagging SFs
	double SF_loose_b = 1.; 		double SF_medium_b = 1.; 		double SF_tight_b = 1.;
	double SF_loose_c = 1.; 		double SF_medium_c = 1.; 		double SF_tight_c = 1.;
	double SF_loose_l = 1.; 		double SF_medium_l = 1.; 		double SF_tight_l = 1.;	
	double SF_loose_b_up = 1.; 		double SF_medium_b_up = 1.; 	double SF_tight_b_up = 1.;
	double SF_loose_c_up = 1.; 		double SF_medium_c_up = 1.; 	double SF_tight_c_up = 1.;
	double SF_loose_l_up = 1.; 		double SF_medium_l_up = 1.; 	double SF_tight_l_up = 1.;
	double SF_loose_b_down = 1.; 	double SF_medium_b_down = 1.; 	double SF_tight_b_down = 1.;
	double SF_loose_c_down = 1.; 	double SF_medium_c_down = 1.; 	double SF_tight_c_down = 1.;
	double SF_loose_l_down = 1.; 	double SF_medium_l_down = 1.; 	double SF_tight_l_down = 1.;
	
	double TotalCSVWeight = 1.; 	double TotalCSVWeight_up = 1.; double TotalCSVWeight_down = 1.; 
	
	for(int iJet=0; iJet<2; iJet++){
		double csv =  jetBtag[HZgCand.bjet[iJet]];
		double JetPt = jetPt[HZgCand.bjet[iJet]]*HZgCand.reg_fac[iJet];
		double JetEta = jetEta[HZgCand.bjet[iJet]];
		
		if(csv < 0.) csv = 0.; // comb as a measurement type has bounds of [0,1] only as a discriminating cut!
		if(csv > 1.) csv = 1.;
	
		if (fabs(jetPartonID[iJet]) == 5 ){ //Heavy-Flavour (B-jets)		  
		
			double iCSV2_loose_b 		= reader_loose_b.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_medium_b 		= reader_medium_b.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_tight_b 		= reader_tight_b.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
				
			double iCSV2_loose_b_up 	= reader_loose_b_up.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_medium_b_up 	= reader_medium_b_up.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_tight_b_up 	= reader_tight_b_up.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);	
				
			double iCSV2_loose_b_down 	= reader_loose_b_down.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_medium_b_down 	= reader_medium_b_down.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);
			double iCSV2_tight_b_down 	= reader_tight_b_down.eval(BTagEntry::FLAV_B, JetEta, JetPt, csv);	

		if(iCSV2_loose_b!=0 || iCSV2_medium_b!=0 || iCSV2_tight_b!=0 || iCSV2_loose_b_up!=0 || iCSV2_medium_b_up!=0 || iCSV2_tight_b_up!=0 || 
		iCSV2_loose_b_down!=0 || iCSV2_medium_b_down!=0 || iCSV2_tight_b_down!=0){
				
					SF_loose_b 			*= iCSV2_loose_b;
					SF_medium_b 		*= iCSV2_medium_b;
					SF_tight_b 			*= iCSV2_tight_b;
					
					SF_loose_b_up 		*= iCSV2_loose_b_up;
					SF_medium_b_up 		*= iCSV2_medium_b_up;
					SF_tight_b_up 		*= iCSV2_tight_b_up;
					
					SF_loose_b_down 	*= iCSV2_loose_b_down;
					SF_medium_b_down 	*= iCSV2_medium_b_down;
					SF_tight_b_down 	*= iCSV2_tight_b_down;
				}
		}
		else if(fabs(jetPartonID[iJet]) == 4 ){  //Charm jets
	
			double iCSV2_loose_c 		= reader_loose_c.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_medium_c 		= reader_medium_c.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_tight_c 		= reader_tight_c.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			
			double iCSV2_loose_c_up 	= reader_loose_c_up.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_medium_c_up 	= reader_medium_c_up.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_tight_c_up 	= reader_tight_c_up.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			
			double iCSV2_loose_c_down 	= reader_loose_c_down.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_medium_c_down 	= reader_medium_c_down.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
			double iCSV2_tight_c_down 	= reader_tight_c_down.eval(BTagEntry::FLAV_C, JetEta, JetPt, csv);
		
			if(iCSV2_loose_c!=0 || iCSV2_medium_c!=0 || iCSV2_tight_c!=0 || iCSV2_loose_c_up!=0 || iCSV2_medium_c_up!=0 || iCSV2_tight_c_up!=0 || 
			iCSV2_loose_c_down!=0 || iCSV2_medium_c_down!=0 || iCSV2_tight_c_down!=0){
					
					SF_loose_c 			*= iCSV2_loose_c;
					SF_medium_c 		*= iCSV2_medium_c;
					SF_tight_c 			*= iCSV2_tight_c;
					
					SF_loose_c_up 		*= iCSV2_loose_c_up;
					SF_medium_c_up 		*= iCSV2_medium_c_up;
					SF_tight_c_up 		*= iCSV2_tight_c_up;
					
					SF_loose_c_down 	*= iCSV2_loose_c_down;
					SF_medium_c_down 	*= iCSV2_medium_c_down;
					SF_tight_c_down 	*= iCSV2_tight_c_down;
				}
			}
		else { //Light-Flavour jets
	
			double iCSV2_loose_l 		= reader_loose_l.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_medium_l 		= reader_medium_l.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_tight_l 		= reader_tight_l.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);

			double iCSV2_loose_l_up 	= reader_loose_l_up.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_medium_l_up 	= reader_medium_l_up.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_tight_l_up 	= reader_tight_l_up.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);

			double iCSV2_loose_l_down 	= reader_loose_l_down.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_medium_l_down 	= reader_medium_l_down.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);
			double iCSV2_tight_l_down 	= reader_tight_l_down.eval(BTagEntry::FLAV_UDSG, JetEta, JetPt, csv);

			if(iCSV2_loose_l!=0 || iCSV2_medium_l!=0 || iCSV2_tight_l!=0 || iCSV2_loose_l_up!=0 || iCSV2_medium_l_up!=0 || iCSV2_tight_l_up!=0 ||
			iCSV2_loose_l_down!=0 || iCSV2_medium_l_down!=0 || iCSV2_tight_l_down!=0){
					
					SF_loose_l 			*= iCSV2_loose_l;
					SF_medium_l 		*= iCSV2_medium_l;
					SF_tight_l 			*= iCSV2_tight_l;
					
					SF_loose_l_up 		*= iCSV2_loose_l_up;
					SF_medium_l_up 		*= iCSV2_medium_l_up;
					SF_tight_l_up 		*= iCSV2_tight_l_up;
					
					SF_loose_l_down 	*= iCSV2_loose_l_down;
					SF_medium_l_down 	*= iCSV2_medium_l_down;
					SF_tight_l_down 	*= iCSV2_tight_l_down;
				}
			}
		TotalCSVWeight = (SF_loose_b*SF_loose_c*SF_loose_l)*(SF_medium_b*SF_medium_c*SF_medium_l)*(SF_tight_b*SF_tight_c*SF_tight_l);
		TotalCSVWeight_up = (SF_loose_b_up*SF_loose_c_up*SF_loose_l_up)*(SF_medium_b_up*SF_medium_c_up*SF_medium_l_up)*(SF_tight_b_up*SF_tight_c_up*SF_tight_l_up);
		TotalCSVWeight_down = (SF_loose_b_down*SF_loose_c_down*SF_loose_l_down)*(SF_medium_b_down*SF_medium_c_down*SF_medium_l_down)*(SF_tight_b_down*SF_tight_c_down*SF_tight_l_down);
	}
	
		CSVSF		= TotalCSVWeight;
		CSVSFUp 	= TotalCSVWeight_up;
		CSVSFDown 	= TotalCSVWeight_down;
		
	//Photon SFs
	PhotonMVAIDSFs(hg_ID, PhoEt, PhoSCEta, PhotonSF, UnPhotonSF, PhotonSFUp, PhotonSFDown);
	
	//Electron Veto SFs
	EleVetoSFs(PhoSCEta, EleVetoSF, UnEleVetoSF, EleVetoSFUp, EleVetoSFDown);
	
	//Electron Trigger SFs
	if(channel==0){ 
		EleTriggerSFs(he_trg1, he_trg2, llepSCEta, tlepSCEta, llepPt, tlepPt, EleTrgSF, EleTrgSFUp, EleTrgSFDown, UnEleTrgSF);
		EleID(he_ID, llepSCEta, tlepSCEta, llepPt, tlepPt, EleIDSF, UnEleIDSF, EleIDSFUp, EleIDSFDown);
		EtaSFs(he_dZ, EtaSF, EtaSFUp, EtaSFDown, UnEtaSF);
		G_SF(he_GSF, llepSCEta, tlepSCEta, GSF, UnGSF, GSFUp, GSFDown);
		TotalSF = PhotonSF*EleVetoSF*EtaSF*EleTrgSF*EleIDSF*GSF*CSVSF;
	}
	if(channel==1){
		MuonTriggerSFs(hm_trg1_0to09, hm_trg1_09to12, hm_trg1_12to21, hm_trg1_21to24, hm_trg2_0to09, hm_trg2_09to12, hm_trg2_12to21, hm_trg2_21to24,
		llepPt, tlepPt, llepEta, tlepEta, MuTrgSF, UnMuTrgSF, MuTrgSFUp, MuTrgSFDown);
		HZZMuID(hm_HZZSF, hm_HZZSFErr, llepEta, tlepEta, llepPt, tlepPt, HZZMuIDSF, UnHZZMuIDSF, HZZMuIDSFUp, HZZMuIDSFDown);
		TotalSF = PhotonSF*EleVetoSF*MuTrgSF*HZZMuIDSF*CSVSF; 
	}
	
	if (data.HasMC()) TotalWei = TotalSF*puWei*genWei*mcWei; else TotalWei = 1.;
	
	higgsM += TotalWei;
	
	//CATEGORIZATIONS
	//These categorizations were designed in such a way there is no overlapping of categories
	//LOW PURITY CATEGORY: both jets are loose tagged 
	if(btagCat[0]==1 && btagCat[1]==1) category = 0; 

	//MEDIUM PURITY CATEGORY: atleast 1 bjet is medium tagged, others are tagged as loose||medium
	if (((btagCat[0]!=3 || btagCat[1]!=3) && (btagCat[0]==2 || btagCat[1]==2))) category = 1; 
	
	//HIGH PURITY CATEGORY: atleast 1 bjet is tight tagged, others are tagged as loose||medium||tight
	if ((btagCat[0]==3 || btagCat[1]==3)) category = 2; 

	
	if(category == 0){
		mBJetLP++;
		bjetLP		+= TotalWei;	
	}
	if(category == 1){
		mBJetMP++;
		bjetMP		+= TotalWei;
	}
	if(category == 2){
		mBJetHP++;
		bjetHP		+= TotalWei;
	}
	
	tree->Fill();
	
	phoID.clear(); 
	eleID.clear();
	muID.clear();
	ZeeID.clear();
	Zmm.clear();
	
	}
	
	cout << "    UNWEIGHTED EVENTS    " << endl;
	cout << "Total no. of events: " << totalEvents << endl;
	cout << "Scraping Filter: " << PV << endl;
	cout << "Trigger: " << Trigger << endl;
	cout << "nPho > 0 || nJet > 1: " << GloSel << endl;
	cout << "Electron: " << Ele << endl;
	cout << "Muon: " << Mu << endl;
	cout << "Dilepton: " << mDilepton << endl;
	cout << "Dijet: " << Dijet << endl;
	cout << "B-tagged jets: " << BTaggedJets << endl;
	cout << "Photon MVA ID: " << Pho << endl;
	cout << "Higgs Candidates: " << mHiggs << endl;
	cout << "Low Purity: " << mBJetLP << endl;
	cout << "Medium Purity: " << mBJetMP << endl;
	cout << "High Purity: " << mBJetHP << endl;
	cout << "-----------------------------" << endl;
	cout << "    WEIGHTED EVENTS    " << endl;
	cout << "Total Yield: " << xs*lumi<< endl;
	cout << "Scraping Filter: " << filter << endl;
	cout << "Trigger: " << trigger << endl;
	cout << "nPho > 0 || nJet > 1: " << globsel << endl;
	cout << "Electron: " << electron << endl;
	cout << "Muon: " << muon << endl;
	cout << "Dilepton: " << dilepton << endl;
	cout << "Dijet: " << dijet << endl;
	cout << "B-tagged jets: " << btagjets << endl;
	cout << "Photon MVA ID: " << photon << endl;
	cout << "Higgs Candidates: " << higgsM << endl;
	cout << "Low Purity: " << bjetLP << endl;
	cout << "Medium Purity: " << bjetMP << endl;
	cout << "High Purity: " << bjetHP << endl;

	fo->cd();
	tree->Write("", TObject::kOverwrite);
	
	delete tree;
	delete fo; 
	
	delete fg_R9;
	delete fe_trg1;
	delete fe_trg2;
	delete fe_dZ;
	delete fe_LowGSF;
	delete fe_LowID;
	delete fg_ID;
	delete fm_trg1_0to09;
	delete fm_trg1_09to12; 
	delete fm_trg1_12to21;
	delete fm_trg1_21to24;
	delete fm_trg2_0to09;
	delete fm_trg2_09to12;
	delete fm_trg2_12to21;
	delete fm_trg2_21to24;
	delete fm_HZZ;
}