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

int main(int argc, char * argv[])
{
   	TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   	// Input files list
	std::string inputList = "/afs/desy.de/user/c/chayanit/private/MSSMHBB_AnaRosty/CMSSW_7_6_3/src/Analysis/MssmHbb/bin/BTagCSVData_Run2015CandD_25ns-16Dec2015.txt";
   
	// Initialisation of MssmHbb class
	MssmHbb analysis(inputList);
	// Process selected JSON file
	//if(!analysis.isMC()) analysis.processJsonFile("Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt");
	if(!analysis.isMC()) analysis.processJsonFile("Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt");

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
	TFile *fBTagEff = 0, *fPtEff = 0, *fHtWeight = 0;
	TH1F *btagEff0p9 = 0, *btagEff0p9_1p4 = 0, *btagEff1p4_2p5 = 0, *HtRatio = 0;
	TH2F *btag2DEff0p85 = 0, *hPtEff = 0;

	if(analysis.getLowMSelection()){
	  //Online BTag Trigger Efficiency produced by Ye Chen
	  fBTagEff = new TFile("input_corrections/RelOnlineBTagCSV0p9Eff_PtEta.root","read");
	  btagEff0p9 = (TH1F*) fBTagEff ->Get("heh4"); // eta <0.9
	  btagEff0p9_1p4 = (TH1F*) fBTagEff ->Get("heh3"); // 1.4 > eta >0.9
	  btagEff1p4_2p5 = (TH1F*) fBTagEff ->Get("heh2"); // 2.5 > eta > 1.4

	  //Online Pt trigger efficiency:
	  fPtEff = new TFile("input_corrections/TwoDPtLowMassEfficiency.root","read");
	  hPtEff = (TH2F*) fPtEff ->Get("TwoDEff_Num"); // 2D

	  // Add Ht reweighting:
	  fHtWeight = new TFile("input_corrections/HtRatio.root","read");
	  HtRatio = (TH1F*) fHtWeight->Get("hRatio");
	}
	else {

	  // For high mass trigger only 2D efficiency were provided
	  fBTagEff = new TFile("input_corrections/TwoDBTagCSV0p85_2D_PtEta.root");
	  btag2DEff0p85 = (TH2F*) fBTagEff->Get("h2ehn");
	  fPtEff = new TFile("input_corrections/TwoDPtHighMassEfficiency.root","read");
	  hPtEff = (TH2F*) fPtEff ->Get("TwoDEff_Num");

	  fHtWeight = new TFile("input_corrections/HtRatio.root","read"); // Useless for the timebinning.
	  HtRatio = (TH1F*) fHtWeight->Get("hRatio");

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

	//For the third jet with lower WP
	BTagCalibrationReader reader_j3(&calib, BTagEntry::OP_LOOSE, "mujets", "central");  // sys up
	BTagCalibrationReader reader_j3_up(&calib, BTagEntry::OP_LOOSE, "mujets", "up");  // sys up
	BTagCalibrationReader reader_j3_down(&calib, BTagEntry::OP_LOOSE, "mujets", "down");  // sys down

	//Setup output file name
	//name can me specified explicitly with method: createOutputFile(fileName);
	std::string fileName = "/nfs/dust/cms/user/chayanit/TripleBTagReverseSelection_76X";
//	std::string fileName = "./TripleBTagReverseSelection_76X";
	analysis.SetupStandardOutputFile(fileName);
//   analysis.createOutputFile("./TripleBTagReverseSelection_76X_cmsdas16.root");

	//Setup Branches
	analysis.setBranches();

	int counter = 0;
	bool goodLeadingJets = true;
	// Analysis of events

	std::cout<<"Number of Entries: "<<analysis.size()<<std::endl;
	for ( int i = 0 ; i < analysis.size() ; ++i )
	  {
	    analysis.event(i);
	    analysis.setTotalNumberOfEvents(analysis.size());

	    //Select only good JSon runs
	    if (!analysis.isMC() && !analysis.selectJson()) continue;
	    
	    //Trigger Selection
	    if( !analysis.isMC() && !analysis.getTriggerResult()) continue;

	    //Used to clean values of the variables for the new event
	    analysis.cleanVariables();

      	    //mHat correction to signal MC:
      	    if(analysis.isSignalMc() && !analysis.mHatCut(0.7)) continue;

	    //Define Jet Collection
	    auto offlineJets = analysis.collection<Jet>("Jets");

	    if (offlineJets -> size() < 3) continue;

	    //Match offline Jets to online Objects
	    if (!analysis.isMC()) analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());

	    //Set Total Number of Jets
	    analysis.setNjets(offlineJets->size());
	    //Define Vertex collection
	    auto offlinePrimaryVertices = analysis.collection<Vertex>("Vertices");

	    counter = 0;
	    goodLeadingJets = false;
	    Jet LeadJet[20];
	    // Selection of double b-tag sample:
	    for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
		Jet jet = offlineJets -> at(iJet);

		//Calculate Ht
		analysis.addHt(jet.pt());

		//Only jets that pass Loose identification will be considered
		if(!jet.idLoose()) break;
		counter++;

		//Only first 5 jets is taken into account
		//To reduce size of the output file
		if(counter > 5) break;
		LeadJet[counter - 1] = jet;

		analysis.setJetCounter(counter-1);
		analysis.setJetVariables(jet);

		//Selection cuts for first three leading jets
		if(counter == 1 || counter == 2 || counter == 3 ){

		  if(counter == 1 || counter == 2){
		    if(jet.pt() < analysis.Pt1Cut()) break;
		    if(abs(jet.eta()) > 2.2) break;
		    if(jet.btag() < analysis.BTag1Cut()) break;

		    if(counter == 2){
		      if(abs(LeadJet[0].eta() - LeadJet[1].eta()) > analysis.dEtaCut()) break;

		      // Online selection is applied only at the first 2 leading jets
		      if(!analysis.isMC() && !analysis.OnlineSelection(LeadJet[0],LeadJet[1])) break;
		    }
		    if(analysis.isMC()) analysis.calculateBTagSF(reader,reader_up,reader_down);
		  }
		  //Selection cuts for the third leading jet
		  if(counter == 3){
		    if(jet.pt() < analysis.Pt3Cut() ) break;
		    if(abs(jet.eta()) > 2.2 ) break;
		    if( jet.btag() > analysis.BTag3Cut() ) break; // Reverse Loose WP

		    // Delta R cuts
		    if(LeadJet[0].deltaR(LeadJet[1]) <= 1) break;
		    if(LeadJet[0].deltaR(LeadJet[2]) <= 1) break;
		    if(LeadJet[1].deltaR(LeadJet[2]) <= 1) break;

		    if(analysis.isMC()) analysis.calculateBTagSF(reader_j3,reader_j3_up,reader_j3_down);

		    goodLeadingJets = true;
		  }
		}
	      }
	    if(!goodLeadingJets) continue;
	
	    //Vertices study
            for(int iVtx = 0; iVtx < offlinePrimaryVertices->size(); ++iVtx){
	      Vertex vertex = offlinePrimaryVertices->at(iVtx);
	      analysis.addVertexInfo(vertex);
            }

            // Method that calculates dEta, dPhi between frist two leading jets
            // Also all parameters of di-jet object is calculated
            analysis.calculateJetVariables();

            if(analysis.isMC()){

	      //.........................Calculate weights...............
	      // general weights which is independent of selection type OR already includes
	      // selection depends criterias inside.
	      weight["FactorPt"] = weightCalc.FactorizationPtWeight(LeadJet[0].pt(), LeadJet[1].pt());
	      weight["dEta"]     = weightCalc.dEtaWeight(abs(LeadJet[0].eta() - LeadJet[1].eta()));
	      weight["2DPt"]     = weightCalc.TwoDPtWeight(hPtEff,LeadJet[0].pt(),LeadJet[1].pt());
	      weight["Lumi"]     = weightCalc.LumiWeight(2318.278306,analysis.luminosity());
	      weight["Ht"]       = weightCalc.HtWeight(HtRatio,analysis.getHt());
	      //TODO: PileUp reweighting is wrong (truePileUp distribution should be used)!!!
	      weight["PileUpCentral"] = weightCalc.PileUpWeight(hPileUpData["central"],hPileUpMC,analysis.nTruePileup());
	      weight["PileUpDown"]    = weightCalc.PileUpWeight(hPileUpData["down"],hPileUpMC,analysis.nTruePileup());
	      weight["PileUpUp"]      = weightCalc.PileUpWeight(hPileUpData["up"],hPileUpMC,analysis.nTruePileup());

	      //Selection depending weights
	      if(analysis.getLowMSelection()){
		weight["BTag"] = weightCalc.BTagWeight(btagEff0p9,btagEff0p9_1p4,btagEff1p4_2p5,LeadJet[0].pt(),LeadJet[0].eta()) *
		  weightCalc.BTagWeight(btagEff0p9,btagEff0p9_1p4,btagEff1p4_2p5,LeadJet[1].pt(),LeadJet[1].eta());
	      }
	      else {
		weight["BTag"] = weightCalc.BTagWeight(btag2DEff0p85, LeadJet[0].pt(),LeadJet[0].eta())*
		  weightCalc.BTagWeight(btag2DEff0p85, LeadJet[1].pt(), LeadJet[1].eta());
	      }

	      //
	      analysis.setPtFactorizationWeight(weight["FactorPt"]);
	      analysis.setdEtaWeight(weight["dEta"]);
	      analysis.setBTagWeight(weight["BTag"]);
	      analysis.setPt2DWeight(weight["2DPt"]);
	      analysis.setLumiWeight(weight["Lumi"]);
	      analysis.setHtWeight(weight["Ht"]);
	      analysis.setPileUpWeight(weight["PileUpCentral"],weight["PileUpUp"],weight["PileUpDown"]);

	      //Calculation of the flavour composition, based on HadronFlavours
	      analysis.calculateFlavourComposition();
            }

	    analysis.fillTree();

	  }
	std::cout<<"Number of Candidates: "<<analysis.getNumberOfCandidates()<<std::endl;
	analysis.writeTree();

	return 0;
	//
}


