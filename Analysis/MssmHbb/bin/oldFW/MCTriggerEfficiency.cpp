/*
 * MCTriggerEfficiency.cpp
 *
 *  Created on: Apr 1, 2016
 *      Author: shevchen
 */

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

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;

int main(int argc, char **argv) {

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
	   analysis.addTriggerObjects();

	   //Select which sets of cuts should be applied:for Low mass selection = true; for High M = false.
	   analysis.setLowMSelection(true);

	   std::string fileName = "/nfs/dust/cms/user/shevchen/output/DoubleBTagSelection_76X_EfficiencyCalculation";
	   analysis.SetupStandardOutputFile(fileName);

	   std::cout<<"Number of Entries: "<<analysis.size()<<std::endl;

	   std::map<std::string, TH1D* > h;
	   h["Eff_noMatching_pt1"] = new TH1D("Eff_noMatching_pt1","Trigger Efficiency without matching vs Pt1",100,0.,500);
	   h["Eff_noMatching_pt2"] = new TH1D("Eff_noMatching_pt2","Trigger Efficiency without matching vs Pt2",100,0.,500);
	   h["Eff_noMatching_BTag1"] = new TH1D("Eff_noMatching_BTag1","Trigger Efficiency without matching vs BTag1",100,0.,1.);
	   h["Eff_noMatching_BTag2"] = new TH1D("Eff_noMatching_BTag2","Trigger Efficiency without matching vs BTag2",100,0.,1.);
	   h["Eff_noMatching_dEta"] = new TH1D("Eff_noMatching_dEta","Trigger Efficiency without matching vs Pt1",100,0.,3.);

	   h["Num_noMatching_pt1"] = new TH1D("Num_noMatching_pt1","Trigger Efficiency without matching vs Pt1",100,0.,500);
	   h["Num_noMatching_pt2"] = new TH1D("Num_noMatching_pt2","Trigger Efficiency without matching vs Pt2",100,0.,500);
	   h["Num_noMatching_BTag1"] = new TH1D("Num_noMatching_BTag1","Trigger Efficiency without matching vs BTag1",100,0.,1.);
	   h["Num_noMatching_BTag2"] = new TH1D("Num_noMatching_BTag2","Trigger Efficiency without matching vs BTag2",100,0.,1.);
	   h["Num_noMatching_dEta"] = new TH1D("Num_noMatching_dEta","Trigger Efficiency without matching vs Pt1",100,0.,3.);

	   h["Den_noMatching_pt1"] = new TH1D("Den_noMatching_pt1","Trigger Efficiency without matching vs Pt1",100,0.,500);
	   h["Den_noMatching_pt2"] = new TH1D("Den_noMatching_pt2","Trigger Efficiency without matching vs Pt2",100,0.,500);
	   h["Den_noMatching_BTag1"] = new TH1D("Den_noMatching_BTag1","Trigger Efficiency without matching vs BTag1",100,0.,1.);
	   h["Den_noMatching_BTag2"] = new TH1D("Den_noMatching_BTag2","Trigger Efficiency without matching vs BTag2",100,0.,1.);
	   h["Den_noMatching_dEta"] = new TH1D("Den_noMatching_dEta","Trigger Efficiency without matching vs Pt1",100,0.,3.);

	   h["Eff_Matching_pt1"] = new TH1D("Eff_Matching_pt1","Trigger Efficiency with matching vs Pt1",100,0.,500);
	   h["Eff_Matching_pt2"] = new TH1D("Eff_Matching_pt2","Trigger Efficiency with matching vs Pt2",100,0.,500);
	   h["Eff_Matching_BTag1"] = new TH1D("Eff_Matching_BTag1","Trigger Efficiency with matching vs BTag1",100,0.,1.);
	   h["Eff_Matching_BTag2"] = new TH1D("Eff_Matching_BTag2","Trigger Efficiency with matching vs BTag2",100,0.,1.);
	   h["Eff_Matching_dEta"] = new TH1D("Eff_Matching_dEta","Trigger Efficiency with matching vs Pt1",100,0.,3.);

	   h["Num_Matching_pt1"] = new TH1D("Num_Matching_pt1","Trigger Efficiency with matching vs Pt1",100,0.,500);
	   h["Num_Matching_pt2"] = new TH1D("Num_Matching_pt2","Trigger Efficiency with matching vs Pt2",100,0.,500);
	   h["Num_Matching_BTag1"] = new TH1D("Num_Matching_BTag1","Trigger Efficiency with matching vs BTag1",100,0.,1.);
	   h["Num_Matching_BTag2"] = new TH1D("Num_Matching_BTag2","Trigger Efficiency with matching vs BTag2",100,0.,1.);
	   h["Num_Matching_dEta"] = new TH1D("Num_Matching_dEta","Trigger Efficiency with matching vs Pt1",100,0.,3.);

	   h["Den_Matching_pt1"] = new TH1D("Den_Matching_pt1","Trigger Efficiency with matching vs Pt1",100,0.,500);
	   h["Den_Matching_pt2"] = new TH1D("Den_Matching_pt2","Trigger Efficiency with matching vs Pt2",100,0.,500);
	   h["Den_Matching_BTag1"] = new TH1D("Den_Matching_BTag1","Trigger Efficiency with matching vs BTag1",100,0.,1.);
	   h["Den_Matching_BTag2"] = new TH1D("Den_Matching_BTag2","Trigger Efficiency with matching vs BTag2",100,0.,1.);
	   h["Den_Matching_dEta"] = new TH1D("Den_Matching_dEta","Trigger Efficiency with matching vs Pt1",100,0.,3.);



	   int counter = 0;

	   for ( int i = 0 ; i < analysis.size() ; ++i )
	   {

	      analysis.event(i);

	      //Define Jet Collection
	      auto offlineJets = analysis.collection<Jet>("Jets");

	      if (offlineJets -> size() < 2) continue;

	      //Match offline Jets to online Objects
	      analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());
	      Jet LeadJet[20];
	      counter = 0;

	      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
	      	Jet jet = offlineJets -> at(iJet);

	      	//Only jets that pass Loose identification will be considered
			if(!jet.idLoose()) continue;
			counter++;
			if(counter > 3) continue;

			LeadJet[counter - 1] = jet;

	      }

	      if(counter < 2) continue;

	      //No matching
	      //As a function of pt1
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_noMatching_pt1"]->Fill(LeadJet[0].pt());
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_noMatching_pt1"]->Fill(LeadJet[0].pt());

	      //As a function of pt2
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_noMatching_pt2"]->Fill(LeadJet[1].pt());
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_noMatching_pt2"]->Fill(LeadJet[1].pt());

	      //As a function of btag1
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_noMatching_BTag1"]->Fill(LeadJet[0].btag());
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_noMatching_BTag1"]->Fill(LeadJet[0].btag());

	      //As a function of btag2
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_noMatching_BTag2"]->Fill(LeadJet[1].btag());
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_noMatching_BTag2"]->Fill(LeadJet[1].btag());

	      //As a function of dEta
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935) h["Den_noMatching_dEta"]->Fill(std::abs(LeadJet[0].eta() - LeadJet[1].eta()));
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && analysis.getTriggerResult()) h["Num_noMatching_dEta"]->Fill(std::abs(LeadJet[0].eta() - LeadJet[1].eta()));


	      //With matching
	      //As a function of pt1
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_Matching_pt1"]->Fill(LeadJet[0].pt());
	      if(analysis.OnlineSelection(LeadJet[0],LeadJet[1]) && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_Matching_pt1"]->Fill(LeadJet[0].pt());

	      //As a function of pt2
	      if(LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_Matching_pt2"]->Fill(LeadJet[1].pt());
	      if(analysis.OnlineSelection(LeadJet[0],LeadJet[1]) && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_Matching_pt2"]->Fill(LeadJet[1].pt());

	      //As a function of btag1
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_Matching_BTag1"]->Fill(LeadJet[0].btag());
	      if(analysis.OnlineSelection(LeadJet[0],LeadJet[1]) && LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_Matching_BTag1"]->Fill(LeadJet[0].btag());

	      //As a function of btag2
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6) h["Den_Matching_BTag2"]->Fill(LeadJet[1].btag());
	      if(analysis.OnlineSelection(LeadJet[0],LeadJet[1]) && LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && std::abs(LeadJet[0].eta() - LeadJet[1].eta()) < 1.6 && analysis.getTriggerResult()) h["Num_Matching_BTag2"]->Fill(LeadJet[1].btag());

	      //As a function of dEta
	      if(LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935) h["Den_Matching_dEta"]->Fill(std::abs(LeadJet[0].eta() - LeadJet[1].eta()));
	      if(analysis.OnlineSelection(LeadJet[0],LeadJet[1]) && LeadJet[0].pt() > 100 && LeadJet[1].pt() > 100 && LeadJet[0].btag() > 0.935 && LeadJet[1].btag() > 0.935 && analysis.getTriggerResult()) h["Num_Matching_dEta"]->Fill(std::abs(LeadJet[0].eta() - LeadJet[1].eta()));


	   }

	   h["Eff_noMatching_pt1"]->Divide(h["Num_noMatching_pt1"],h["Den_noMatching_pt1"],1,1,"b");
	   h["Eff_noMatching_pt2"]->Divide(h["Num_noMatching_pt2"],h["Den_noMatching_pt2"],1,1,"b");
	   h["Eff_noMatching_BTag1"]->Divide(h["Num_noMatching_BTag1"],h["Den_noMatching_BTag1"],1,1,"b");
	   h["Eff_noMatching_BTag2"]->Divide(h["Num_noMatching_BTag2"],h["Den_noMatching_BTag2"],1,1,"b");
	   h["Eff_noMatching_dEta"]->Divide(h["Num_noMatching_dEta"],h["Den_noMatching_dEta"],1,1,"b");

	   h["Eff_Matching_pt1"]->Divide(h["Num_Matching_pt1"],h["Den_Matching_pt1"],1,1,"b");
	   h["Eff_Matching_pt2"]->Divide(h["Num_Matching_pt2"],h["Den_Matching_pt2"],1,1,"b");
	   h["Eff_Matching_BTag1"]->Divide(h["Num_Matching_BTag1"],h["Den_Matching_BTag1"],1,1,"b");
	   h["Eff_Matching_BTag2"]->Divide(h["Num_Matching_BTag2"],h["Den_Matching_BTag2"],1,1,"b");
	   h["Eff_Matching_dEta"]->Divide(h["Num_Matching_dEta"],h["Den_Matching_dEta"],1,1,"b");

	   for( const auto & ih : h){
		   ih.second->Write();
	   }

	   analysis.writeTree();


	return 0;
}



