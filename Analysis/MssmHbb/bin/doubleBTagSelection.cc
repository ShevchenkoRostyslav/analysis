#include <string>
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "stdlib.h"

#include "Analysis/MssmHbb/interface/MssmHbb.h"
#include "Analysis/MssmHbb/interface/json.h"
#include "Analysis/MssmHbb/interface/BTagCalibrationStandalone.h"
#include "Analysis/MssmHbb/interface/Weights.h"

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;

double weight(TH2F *, const float &, const float &, const float &, const float &);
double weight2(TH2F *, TH2F *, const float &, const float &, const float &, const float &);
double ptRewfunction(const float &);
double ptRew2D(TH2F *, const float &, const float &);
double dEtaRew(const float &);
double bTagRew(TH2F *, const float&, const float &);

// =============================================================================================
int main(int argc, char * argv[])
{

	   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

	   // Input files list
	   //std::string inputList = "rootFileListBTagCSV.txt";
	   std::string inputList = "/nfs/dust/cms/user/shevchen/samples/miniaod/76X/MSSMHbb/SUSYGluGluToBBHToBB_M-300_TuneCUETP8M1_13TeV-pythia8.txt";
	   // Initialisation of MssmHbb class
	   MssmHbb analysis(inputList);
	   // Process selected JSON file
	   if(!analysis.isMC()) analysis.processJsonFile("Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt");

	   // Add std::vector<std::string> of the Trigger Objects that you would like to apply.
	   // Also Trigger Results name will be stored, according to the trigger objects names
	   // By default - LowM MssmHbb trigger objects is used
	   std::vector<std::string> triggerObjects = {"hltL1sL1DoubleJetC100","hltDoubleJetsC100","hltDoubleBTagCSV0p85","hltDoublePFJetsC160"};
	   if(!analysis.isMC()) analysis.addTriggerObjects();

	   //Select which sets of cuts should be applied:for Low mass selection = true; for High M = false.
	   if(analysis.isMC()) analysis.setLowMSelection(true);
	   //Create Object to calculate the weights:
	   Weights weightCalc(analysis.getLowMSelection());
	   std::map<std::string, double> weight;

	   //.................................Input corrections.................................
	   TFile *fBTagEff = 0, *fPtEff = 0;
	   TH1F *btagEff0p9 = 0, *btagEff0p9_1p4 = 0, *btagEff1p4_2p5 = 0;
	   TH2F *btag2DEff0p85 = 0, *hPtEff = 0;

	   if(analysis.getLowMSelection()){
		   //Online BTag Trigger Efficiency produced by Ye Chen
		   fBTagEff = new TFile("input_corrections/RelOnlineBTagCSV0p9Eff_PtEta.root","read");
		   btagEff0p9 = (TH1F*) fBTagEff ->Get("heh4");	 // eta <0.9
		   btagEff0p9_1p4 = (TH1F*) fBTagEff ->Get("heh3"); // 1.4 > eta >0.9
		   btagEff1p4_2p5 = (TH1F*) fBTagEff ->Get("heh2"); // 2.5 > eta > 1.4

		   //Online Pt trigger efficiency:
		   fPtEff = new TFile("input_corrections/TwoDPtLowMassEfficiency.root","read");
		   hPtEff = (TH2F*) fPtEff ->Get("TwoDEff_Num"); // 2D

	   }
	   else {

		   // For high mass trigger only 2D efficiency were provided
		   fBTagEff = new TFile("input_corrections/TwoDBTagCSV0p85_2D_PtEta.root");
		   btag2DEff0p85 = (TH2F*) fBTagEff->Get("h2ehn");

		   fPtEff = new TFile("input_corrections/TwoDPtHighMassEfficiency.root","read");
		   hPtEff = (TH2F*) fPtEff ->Get("TwoDEff_Num");

	   }

	   // Add PileUp reweighting
	   std::map<std::string, TFile*> fPileUpData;
	   std::map<std::string, TH1F* > hPileUpData;
	   fPileUpData["central"] = new TFile("input_corrections/PileUp_2015Dec_central.root","read");
	   hPileUpData["central"] = (TH1F*) fPileUpData["central"]->Get("pileup");
	   fPileUpData["down"] = new TFile("input_corrections/PileUp_2015Dec_up.root","read");
	   hPileUpData["down"] = (TH1F*) fPileUpData["down"]->Get("pileup");
	   fPileUpData["up"] = new TFile("input_corrections/PileUp_2015Dec_down.root","read");
	   hPileUpData["up"] = (TH1F*) fPileUpData["up"]->Get("pileup");
	   TFile *fPileUpMC = new TFile("input_corrections/MC_Fall15_PU25_V1.root","read");
	   TH1F *hPileUpMC = (TH1F*) fPileUpMC->Get("pileup");


	   //Add BTagCalibration calculators needed for Offline BTag SF:
	   BTagCalibration calib("csvv2", "input_corrections/SFbLib.csv");
	   BTagCalibrationReader reader(&calib,               // calibration instance
	                                BTagEntry::OP_TIGHT,  // operating point
	                                "mujets",               // measurement type
	                                "central");           // systematics type
	   BTagCalibrationReader reader_up(&calib, BTagEntry::OP_TIGHT, "mujets", "up");  // sys up
	   BTagCalibrationReader reader_down(&calib, BTagEntry::OP_TIGHT, "mujets", "down");  // sys down

   int nCand = 0;
   std::vector <int> leadJetsNumber;

   std::map < std::string, TH1F *> h1;
   h1["h_jet_N"]         = new TH1F("h_jet_N" , "", 20, 0., 20.);

	h1["h_jet1_Pt"]       = new TH1F("h_jet1_Pt" , "", 100, 0., 1000.);
	h1["h_jet1_Eta"]      = new TH1F("h_jet1_Eta", "", 30, -5, 5.);
	h1["h_jet1_Phi"]      = new TH1F("h_jet1_Phi", "", 100, -3.2, 3.2);
	h1["h_jet1_BTag"]     = new TH1F("h_jet1_BTag", "", 100, 0., 1.);
	h1["h_jet1_BTagZoom"] = new TH1F("h_jet1_BTagZoom", "", 100, 0.94, 1.);

   	h1["h_jet2_Pt"]       = new TH1F("h_jet2_Pt" , "", 100, 0., 1000.);
   	h1["h_jet2_Eta"]      = new TH1F("h_jet2_Eta", "", 30, -5, 5.);
  	h1["h_jet2_Phi"]      = new TH1F("h_jet2_Phi", "", 100, -3.2, 3.2);
  	h1["h_jet2_BTag"]     = new TH1F("h_jet2_BTag", "", 100, 0., 1.);
  	h1["h_jet2_BTagZoom"] = new TH1F("h_jet2_BTagZoom", "", 100, 0.94, 1.);

	h1["h_obj_M"]  		  = new TH1F("h_obj_M12","",100,0.,1000.);
	h1["h_obj_Pt"]		  = new TH1F("h_obj_Pt","",100,0.,1000.);
	h1["h_obj_Eta"]		  = new TH1F("h_obj_Eta","",100,-5.,5.);
	h1["h_obj_Phi"]       = new TH1F("h_obj_Phi","",100,-3.2,3.2);

	h1["h_djet_Pt12"]     = new TH1F("h_djet_Pt12","",100,0.,200);
	h1["h_djet_Eta12"]    = new TH1F("h_djet_Eta12","",100,-5.,5.);
	h1["h_djet_Phi12"]    = new TH1F("h_djet_Phi12","",100,-3.2,3.2);
	h1["h_djet_PtTOT"]	  = new TH1F("h_djet_PtTOT","",100,0,1);

	h1["h_jet_dEta"]	  = new TH1F("h_jet_dEta","",30,0,3);

	//Plots to study reweighting:
	h1["h_jet1_PtNoRew"]      = new TH1F("h_jet1_PtNoRew", "", 100, 0., 1000.);
	h1["h_jet1_Ptrew1"]       = new TH1F("h_jet1_Ptrew1" , "", 100, 0., 1000.);
	h1["h_jet1_Ptrew2"]       = new TH1F("h_jet1_Ptrew2" , "", 100, 0., 1000.);
	h1["h_jet1_PtplusDEta"]   = new TH1F("h_jet1_PtplusDEta" , "", 100, 0., 1000.);
	h1["h_jet1_All"]          = new TH1F("h_jet1_All" , "", 100, 0., 1000.);
	h1["h_jet1_All2"]         = new TH1F("h_jet1_All2" , "", 100, 0., 1000.);

	h1["h_jet1_qq_Pt"]		  = new TH1F("h_jet1_qq_Pt", "", 100, 0., 1000.);
	h1["h_jet1_bb_Pt"]		  = new TH1F("h_jet1_bb_Pt", "", 100, 0., 1000.);
	h1["h_jet1_cc_Pt"]		  = new TH1F("h_jet1_cc_Pt", "", 100, 0., 1000.);
	h1["h_jet1_qb_Pt"]		  = new TH1F("h_jet1_qb_Pt", "", 100, 0., 1000.);
	h1["h_jet1_qc_Pt"]		  = new TH1F("h_jet1_qc_Pt", "", 100, 0., 1000.);
	h1["h_jet1_cb_Pt"]		  = new TH1F("h_jet1_cb_Pt", "", 100, 0., 1000.);

   // Analysis of events
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {

      analysis.event(i);

      //Select only good JSon runs
      if (!analysis.isMC() && !analysis.selectJson()) continue;

      //Trigger Selection
      if( !analysis.isMC() && !analysis.getTriggerResult()) continue;

      //Match offline Jets to online Objects
      if (!analysis.isMC()) analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());

      //Define Jet Collection
      auto offlineJets = analysis.collection<Jet>("Jets");

      if (offlineJets -> size() < 2) continue;


      // Selection of double b-tag sample:

      leadJetsNumber.clear();
      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
      {
      	Jet jet = offlineJets -> at(iJet);
		if(jet.idLoose() ) leadJetsNumber.push_back(iJet);

      }
      if( leadJetsNumber.size() < 2) continue;
       // Match Offline jet objects to trigger objects
      Jet fLeadJet = offlineJets -> at(leadJetsNumber[0]);
      Jet sLeadJet = offlineJets -> at(leadJetsNumber[1]);
      //std::cout<<"Flavours: "<<fLeadJet.flavour("Hadron")<<" "<<std::endl;
      if( fLeadJet.pt() < 100 || sLeadJet.pt() < 100 ) continue;
      if( std::abs(fLeadJet.eta()) > 2.2 || std::abs(sLeadJet.eta()) > 2.2 ) continue;
      if( std::abs(fLeadJet.eta() - sLeadJet.eta()) > 1.6 ) continue;
      if( fLeadJet.deltaR(sLeadJet) <= 1) continue;
      if( fLeadJet.btag() < 0.935 || sLeadJet.btag() < 0.935 ) continue;

      h1["h_jet1_PtNoRew"] -> Fill(fLeadJet.pt());
      if(!analysis.isMC() && !analysis.OnlineSelection(fLeadJet,sLeadJet) ) continue;

     if(analysis.isMC()) {

   	  //.........................Calculate weights...............
   	  // general weights which is independent of selection type OR already includes
   	  // selection depends criterias inside.
   	  weight["FactorPt"] = weightCalc.FactorizationPtWeight(fLeadJet.pt(), sLeadJet.pt());
   	  weight["dEta"]     = weightCalc.dEtaWeight(abs(fLeadJet.eta() - sLeadJet.eta()));
   	  weight["2DPt"]     = weightCalc.TwoDPtWeight(hPtEff,fLeadJet.pt(),sLeadJet.pt());
   	  weight["Lumi"]     = weightCalc.LumiWeight(2318.278306,analysis.luminosity());

   	  weight["PileUpCentral"] = weightCalc.PileUpWeight(hPileUpData["central"],hPileUpMC,analysis.nTruePileup());
   	  weight["PileUpDown"]    = weightCalc.PileUpWeight(hPileUpData["down"],hPileUpMC,analysis.nTruePileup());
   	  weight["PileUpUp"]      = weightCalc.PileUpWeight(hPileUpData["up"],hPileUpMC,analysis.nTruePileup());

   	  //Selection depending weights
   	  if(analysis.getLowMSelection()){
       	  weight["BTag"] = weightCalc.BTagWeight(btagEff0p9,btagEff0p9_1p4,btagEff1p4_2p5,fLeadJet.pt(),fLeadJet.eta()) *
       			  	  	   weightCalc.BTagWeight(btagEff0p9,btagEff0p9_1p4,btagEff1p4_2p5,sLeadJet.pt(),sLeadJet.eta());
   	  }
   	  else {
       	  weight["BTag"] = weightCalc.BTagWeight(btag2DEff0p85, fLeadJet.pt(),fLeadJet.eta())*
       			  	  	   weightCalc.BTagWeight(btag2DEff0p85, sLeadJet.pt(), sLeadJet.eta());
   	  }

   	  h1["h_jet1_Ptrew1"] -> Fill(fLeadJet.pt(),weight["FactorPt"]);
   	  h1["h_jet1_Ptrew2"] -> Fill(fLeadJet.pt(),weight["2DPt"]);
   	  h1["h_jet1_PtplusDEta"] -> Fill(fLeadJet.pt(),weight["2DPt"]*weight["dEta"]);
     	 h1["h_jet1_All"] -> Fill(fLeadJet.pt(),weight["2DPt"]*weight["dEta"]*weight["BTag"]);
     }

     double reweight;
     if(!analysis.isMC()) reweight = 1.;
     else reweight = weight["2DPt"]*weight["dEta"]*weight["BTag"]*weight["PileUpCentral"]*weight["Lumi"];
     //else reweight = weight2(PtEff,btagEff,fLeadJet.pt(),fLeadJet.eta(),sLeadJet.pt(),sLeadJet.eta());

     h1["h_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     h1["h_jet1_Eta"] -> Fill(fLeadJet.eta(),reweight);
     h1["h_jet1_Phi"] -> Fill(fLeadJet.phi(),reweight);
     h1["h_jet1_BTag"] -> Fill(fLeadJet.btag(),reweight);
     h1["h_jet1_BTagZoom"] -> Fill(fLeadJet.btag(),reweight);

     h1["h_jet2_Pt"] -> Fill(sLeadJet.pt(),reweight);
     h1["h_jet2_Eta"] -> Fill(sLeadJet.eta(),reweight);
     h1["h_jet2_Phi"] -> Fill(sLeadJet.phi(),reweight);
     h1["h_jet2_BTag"] -> Fill(sLeadJet.btag(),reweight);
     h1["h_jet2_BTagZoom"] -> Fill(sLeadJet.btag(),reweight);

     h1["h_jet_dEta"] -> Fill(std::abs(fLeadJet.eta() - sLeadJet.eta()),reweight);

     TLorentzVector obj;
     obj = fLeadJet.p4() + sLeadJet.p4();

     h1["h_obj_Pt"] -> Fill(obj.Pt(),reweight);
     h1["h_obj_Eta"] -> Fill(obj.Eta(),reweight);
     h1["h_obj_Phi"] -> Fill(obj.Phi(),reweight);
     h1["h_obj_M"] -> Fill(obj.M(),reweight);

     h1["h_djet_Pt12"] -> Fill(fLeadJet.pt() - sLeadJet.pt(),reweight);
     h1["h_djet_Eta12"] -> Fill( fLeadJet.eta() - sLeadJet.eta(),reweight);
     h1["h_djet_Phi12"] -> Fill( fLeadJet.phi() - sLeadJet.phi(),reweight);
     h1["h_djet_PtTOT"] -> Fill( (fLeadJet.pt() - sLeadJet.pt() )/ (fLeadJet.pt() + sLeadJet.pt()) ,reweight);

     if(fLeadJet.flavour("Hadron") == 0 && sLeadJet.flavour("Hadron") == 0 )
     {
     	h1["h_jet1_qq_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }
     else if (fLeadJet.flavour("Hadron") == 4 && sLeadJet.flavour("Hadron") == 4)
     {
     	h1["h_jet1_cc_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }
     else if (fLeadJet.flavour("Hadron") == 5 && sLeadJet.flavour("Hadron") == 5)
     {
     	h1["h_jet1_bb_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }
     else if ( (fLeadJet.flavour("Hadron") == 0 || sLeadJet.flavour("Hadron") == 0) && (sLeadJet.flavour("Hadron") == 5 || fLeadJet.flavour("Hadron") == 5 ) )
     {
     	h1["h_jet1_qb_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }
     else if ( (fLeadJet.flavour("Hadron") == 4 || sLeadJet.flavour("Hadron") == 4) && ( sLeadJet.flavour("Hadron") == 0 || fLeadJet.flavour("Hadron") == 0 ) )
     {
     	h1["h_jet1_qc_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }
     else if ( (fLeadJet.flavour("Hadron") == 5 || sLeadJet.flavour("Hadron") == 5) && (sLeadJet.flavour("Hadron") == 4 || fLeadJet.flavour("Hadron") == 4 ) )
     {
     	h1["h_jet1_cb_Pt"] -> Fill(fLeadJet.pt(),reweight);
     }

      ++nCand;

   }
   std::cout<<"Number of Candidates: "<<nCand<<std::endl;
   std::string fileName = "output.root";
   TFile * outFile = new TFile(fileName.c_str(),"RECREATE");

   //double lumiData2015DPrRec3 = 553.149951;
//   double lumiData2015DPrRec4 = 1001.179266;
   for ( auto& ih1 : h1 )
   {
//   	if(analysis.isMC() && fileName.find("QCD") != std::string::npos ) ih1.second->Scale(lumiData2015DPrRec4/analysis.luminosity());
//   	else if (analysis.isMC() && fileName.find("TT") != std::string::npos) 	{ih1.second->Scale(lumiData2015DPrRec4/analysis.luminosity(850.));}
   	ih1.second -> Write();
   }
   outFile -> Close();

   std::cout<<" Output File: "<<fileName<<" was stored!"<<std::endl;

   return 0;
//
}

double weight(TH2F * BTagH, const float &pt1, const float &eta1, const float & pt2, const float & eta2)
{
  double weightPt1 = ptRewfunction(pt1);
	double weightPt2 = ptRewfunction(pt2);
	double dEta = std::abs(eta1 - eta2);
	double weightdEta = dEtaRew(dEta);
	double weightBTag1 = bTagRew(BTagH,pt1,eta1);
	double weightBTag2 = bTagRew(BTagH,pt2,eta2);
	double result = weightPt1 * weightPt2 * weightdEta * weightBTag1 * weightBTag2;
	return result;
}

double weight2(TH2F * PtH, TH2F * BTagH, const float &pt1, const float &eta1, const float & pt2, const float & eta2)
{
	double weightPt1 = ptRew2D(PtH,pt1,eta1);
	double weightPt2 = ptRew2D(PtH,pt2,eta2);
	double dEta = std::abs(eta1 - eta2);
	double weightdEta = dEtaRew(dEta);
	double weightBTag1 = bTagRew(BTagH,pt1,eta1);
	double weightBTag2 = bTagRew(BTagH,pt2,eta2);
	double result = weightPt1 * weightPt2 * weightdEta * weightBTag1 * weightBTag2;
	//double result = weightBTag1 * weightBTag2;
	return result;
}

double ptRewfunction(const float &pt)
{
	return 1./(1.+exp( -8.87968e-02 * (pt - 1.03879e+02) )) * 1./ (1.+exp( -1.88620e-01 * (pt - 1.02891e+02) ));
}
double ptRew2D(TH2F *PtH, const float & pt,const float &eta)
{
	if(pt <500) return PtH -> Interpolate(pt,eta);
	else return 1;
}
double dEtaRew(const float & dEta)
{
//	return 1./(1.+exp( 1.15733e+02 * (dEta - 1.59862e+00) )) * 1./ (1.+exp( -2.06004e-02 * (dEta + 2.22684e+02) ));
	return 1./(1.+exp( 1.41132e+02 * (dEta - 1.60010e+00) )) * 1./ (1.+exp( -2.06352e-02 * (dEta + 2.27738e+02) ));
}
double bTagRew(TH2F *histo, const float &pt, const float &eta)
{
	if(pt <1000) return histo -> Interpolate(pt,eta);
	else return 1;
}
