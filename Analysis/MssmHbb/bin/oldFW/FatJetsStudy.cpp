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
	   std::string inputList = "/nfs/dust/cms/user/shevchen/samples/miniaod/76X/MSSMHbb/SUSYGluGluToBBHToBB_M-700_TuneCUETP8M1_13TeV-pythia8.txt";
	   // Initialisation of MssmHbb class
	   MssmHbb analysis(inputList);
	   // Process selected JSON file
	   if(!analysis.isMC()) analysis.processJsonFile("Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt");

	   // Add std::vector<std::string> of the Trigger Objects that you would like to apply.
	   // Also Trigger Results name will be stored, according to the trigger objects names
	   // By default - LowM MssmHbb trigger objects is used
	   std::vector<std::string> triggerObjects = {"hltL1sL1DoubleJetC100","hltDoubleJetsC100","hltDoubleBTagCSV0p85","hltDoublePFJetsC160"};
	   analysis.addTriggerObjects();

	   //Select which sets of cuts should be applied:for Low mass selection = true; for High M = false.z
	   analysis.setLowMSelection(true);

	   std::string fileName = "/nfs/dust/cms/user/shevchen/output/DoubleBTagSelection_76X_FatJets";
	   analysis.SetupStandardOutputFile(fileName);

	   std::cout<<"Number of Entries: "<<analysis.size()<<std::endl;

	   std::map<std::string, TH1D* > h;
	   h["Njets"] = new TH1D("Njets","",20,0,20);

	   h["LeadJet_pt1"] 	= new TH1D("LeadJet_pt1","",100,0,1500);
	   h["LeadJet_pt2"] 	= new TH1D("LeadJet_pt2","",100,0,1500);
	   h["LeadJet_eta1"] 	= new TH1D("LeadJet_eta1","",100,-2.2,2.2);
	   h["LeadJet_eta2"] 	= new TH1D("LeadJet_eta2","",100,-2.2,2.2);
	   h["LeadJet_btag1"] 	= new TH1D("LeadJet_btag1","",100,0,1);
	   h["LeadJet_btag2"] 	= new TH1D("LeadJet_btag2","",100,0,1);
	   h["LeadJet_dR12"] 	= new TH1D("LeadJet_dR12","",100,0,6);

	   h["FatJet_pt1"] 		= new TH1D("FatJet_pt1","",100,0,1500);
	   h["FatJet_pt2"] 		= new TH1D("FatJet_pt2","",100,0,1500);
	   h["FatJet_eta1"] 	= new TH1D("FatJet_eta1","",100,-2.2,2.2);
	   h["FatJet_eta2"] 	= new TH1D("FatJet_eta2","",100,-2.2,2.2);
	   h["FatJet_dR12"] 	= new TH1D("FatJet_dR12","",100,0,6);

	   h["Obj12_pt"] 		= new TH1D("Obj12_pt","",100,0,1500);
	   h["Obj12_eta"] 		= new TH1D("Obj12_eta","",100,-3.,3.);
	   h["Obj12_M"] 		= new TH1D("Obj12_M","",100,0,1500);

	   h["FatObj12_pt"] 	= new TH1D("FatObj12_pt","",100,0,1500);
	   h["FatObj12_eta"] 	= new TH1D("FatObj12_eta","",100,-3.,3.);
	   h["FatObj12_M"] 		= new TH1D("FatObj12_M","",100,0,1500);

	   h["SmallJet_N"] 		= new TH1D("SmallJet_N","",20,0,20);
	   h["SmallJet_pt"]		= new TH1D("SmallJet_pt","",100,0,1500);
	   h["SmallJet_eta"]	= new TH1D("SmallJet_eta","",100,-5.,5.);
	   h["SmallJet_btag"]	= new TH1D("SmallJet_btag","",100,0,1);
	   h["SmallJet_dR"]		= new TH1D("SmallJet_dR","",100,0,6);
	   h["SmallJet_dR_zoom"]= new TH1D("SmallJet_dR_zoom","",100,0,2);
	   h["SmallJet_dR_mid"] = new TH1D("SmallJet_dR_mid","",100,-1.1,1.1);

	   //gen jets:
	   analysis.addTree<Jet>("GenJets","MssmHbb/Events/slimmedGenJets");

	   //gen particles
	   analysis.addTree<GenParticle>("GenParticles","MssmHbb/Events/prunedGenParticles");

	   h["GenParticle_M"]	= new TH1D("GenParticle_M","",100,0,1500);
	   h["GenParticle_pt"]	= new TH1D("GenParticle_pt","",100,0,1500);
	
	   h["GenParticles_Sum_M"]   = new TH1D("GenParticles_Sum_M","",100,0,1500);
       h["GenParticles_Sum_pt"]  = new TH1D("GenParticles_Sum_pt","",100,0,1500);

	   h["GenJet_M"]		= new TH1D("GenJet_M","",100,0.,1500.);
	   h["GenJet_pt"]		= new TH1D("GenJet_pt","",100,0.,1500);




	   int counter = 0;
	   int n_small_jets = 0;
	   bool goodLeadingJets = false;


	   for ( int i = 0 ; i < analysis.size() ; ++i )
	   {

	      analysis.event(i);

	      //Define Jet Collection
	      auto offlineJets = analysis.collection<Jet>("Jets");

	      //Define Gen Jet collection
	      auto genJets = analysis.collection<Jet>("GenJets"	);

	      //Define Gen Particles collection
	      auto genParticles = analysis.collection<GenParticle>("GenParticles");

	      if (offlineJets -> size() < 2) continue;

	      //Match offline Jets to online Objects
	      analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());
	      Jet LeadJet[20];
	      TLorentzVector FatJet1, FatJet2, FatObj12, Obj12;
	      n_small_jets = 0;
	      counter = 0;
	      goodLeadingJets = false;

	      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
	      	Jet jet = offlineJets -> at(iJet);

	      	//Only jets that pass Loose identification will be considered
			if(!jet.idLoose()) continue;
			counter++;

			LeadJet[counter - 1] = jet;

			if(counter == 1 || counter == 2 ){
				if(jet.pt() < 100 ) break;
				if(abs(jet.eta()) > 2.2) break;
				if(jet.btag() < analysis.BTag1Cut()) break;

				if(counter == 2){
					if(LeadJet[0].deltaR(LeadJet[1]) <= 1.1) break;
					if(abs(LeadJet[0].eta() - LeadJet[1].eta()) > analysis.dEtaCut()) break;
					FatJet1 = LeadJet[0].p4();
					FatJet2 = LeadJet[1].p4();
					goodLeadingJets = true;
				}
			}
			else{
				Jet closest_lead_jet;
				if (jet.deltaR(LeadJet[0]) < jet.deltaR(LeadJet[1])) closest_lead_jet = LeadJet[0];
				else closest_lead_jet = LeadJet[1];

				if(jet.pt() > 10){

					if (jet.deltaR(LeadJet[0]) < 1.1 || jet.deltaR(LeadJet[1]) < 1.1){
						h["SmallJet_pt"] ->Fill(jet.pt());
						h["SmallJet_eta"]->Fill(jet.eta());
						h["SmallJet_btag"]->Fill(jet.btag());
						h["SmallJet_dR"] ->Fill(jet.deltaR(closest_lead_jet));
						h["SmallJet_dR_zoom"] ->Fill(jet.deltaR(closest_lead_jet));
						n_small_jets ++;
					}

					if (jet.deltaR(LeadJet[0]) < 1.1 && jet.deltaR(LeadJet[1]) > 1.1) {
						FatJet1 += jet.p4();
					}
					else if (jet.deltaR(LeadJet[1]) < 1.1 && jet.deltaR(LeadJet[0]) > 1.1){
						FatJet2 += jet.p4();
					}
					else if (jet.deltaR(LeadJet[0])	< 1.1 && jet.deltaR(LeadJet[1]) < 1.1){
						h["SmallJet_dR_mid"]->Fill(jet.deltaR(LeadJet[0]) - jet.deltaR(LeadJet[1]));
						if(jet.deltaR(LeadJet[0]) < jet.deltaR(LeadJet[1])) FatJet1 += jet.p4();
						else FatJet2 += jet.p4();
					}
				}

			}

	      }

	      if(!goodLeadingJets) continue;
	      if(counter < 2) continue;

	      //Fat Jets
	      FatObj12 = FatJet1 + FatJet2;
	      h["Njets"]->Fill(offlineJets->size());
	      h["FatJet_pt1"]->Fill(FatJet1.Pt());
	      h["FatJet_pt2"]->Fill(FatJet2.Pt());
	      h["FatJet_eta1"]->Fill(FatJet1.Eta());
	      h["FatJet_eta2"]->Fill(FatJet2.Eta());
	      h["FatJet_dR12"]->Fill(FatJet1.DeltaR(FatJet2));

	      h["FatObj12_M"]->Fill(FatObj12.M());
	      h["FatObj12_pt"]->Fill(FatObj12.Pt());
	      h["FatObj12_eta"]->Fill(FatObj12.Eta());

	      Obj12 = LeadJet[0].p4() + LeadJet[1].p4();
		  h["LeadJet_pt1"] -> Fill(LeadJet[0].pt());
		  h["LeadJet_pt2"] -> Fill(LeadJet[1].pt());
		  h["LeadJet_eta1"] -> Fill(LeadJet[0].eta());
		  h["LeadJet_eta2"] -> Fill(LeadJet[1].eta());
		  h["LeadJet_btag1"] -> Fill(LeadJet[0].btag());
		  h["LeadJet_btag2"] -> Fill(LeadJet[1].btag());
		  h["LeadJet_dR12"]->Fill(LeadJet[0].deltaR(LeadJet[1]));

	      h["Obj12_M"]->Fill(Obj12.M());
	      h["Obj12_pt"]->Fill(Obj12.Pt());
	      h["Obj12_eta"]->Fill(Obj12.Eta());

	      h["SmallJet_N"]->Fill(n_small_jets);

	      TLorentzVector genJetH, genParticleH, genJJH;
	      int higgsCounter = 0;
	      genJetH = genJets->at(0).p4() + genJets->at(1).p4();

		  h["GenJet_M"]  ->Fill(genJetH.M());
		  h["GenJet_pt"] ->Fill(genJetH.Pt());
	      for(int iGen=0;iGen<genParticles->size();++iGen){
			GenParticle part = genParticles->at(iGen);	
			if(part.pdgId() == 36) {
				higgsCounter++;
				genParticleH = part.p4();
			}
			if(part.higgsDaughter() && part.pdgId() != 36 ){
				genJJH += part.p4();	
			}
	      }
		
	      h["GenParticle_M"] -> Fill(genParticleH.M());
	      h["GenParticle_pt"] -> Fill(genParticleH.Pt());
		
	      h["GenParticles_Sum_M"] -> Fill(genJJH.M());
          h["GenParticles_Sum_pt"] -> Fill(genJJH.Pt());
	   }

	   for( const auto & ih : h){
		   ih.second->Write();
	   }

	   analysis.writeTree();


	return 0;
}



