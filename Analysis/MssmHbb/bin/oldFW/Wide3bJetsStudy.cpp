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
	   std::string inputList = "/nfs/dust/cms/user/shevchen/samples/miniaod/76X/MSSMHbb/SUSYGluGluToBBHToBB_M-1100_TuneCUETP8M1_13TeV-pythia8.txt";
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

	   std::string fileName = "/nfs/dust/cms/user/shevchen/output/TripleBTagSelection_76X_FatJets";
	   analysis.SetupStandardOutputFile(fileName);

	   std::cout<<"Number of Entries: "<<analysis.size()<<std::endl;

	   std::map<std::string, TH1D* > h;
	   h["Njets"] = new TH1D("Njets","",20,0,20);

	   h["LeadJet_pt1"] 	= new TH1D("LeadJet_pt1","",100,0,1500);
	   h["LeadJet_pt2"] 	= new TH1D("LeadJet_pt2","",100,0,1500);
	   h["LeadJet_pt3"] 	= new TH1D("LeadJet_pt3","",100,0,1500);
	   h["LeadJet_eta1"] 	= new TH1D("LeadJet_eta1","",100,-2.2,2.2);
	   h["LeadJet_eta2"] 	= new TH1D("LeadJet_eta2","",100,-2.2,2.2);
	   h["LeadJet_eta3"] 	= new TH1D("LeadJet_eta3","",100,-2.2,2.2);
	   h["LeadJet_btag1"] 	= new TH1D("LeadJet_btag1","",100,0,1);
	   h["LeadJet_btag2"] 	= new TH1D("LeadJet_btag2","",100,0,1);
	   h["LeadJet_btag3"] 	= new TH1D("LeadJet_btag3","",100,0,1);
	   h["LeadJet_dR12"] 	= new TH1D("LeadJet_dR12","",100,0,6);
	   h["LeadJet_dR13"] 	= new TH1D("LeadJet_dR13","",100,0,6);
	   h["LeadJet_dR23"] 	= new TH1D("LeadJet_dR23","",100,0,6);

	   h["NumberOf3dJets"] = new TH1D("NumberOf3dJets","",20,0,20);
	   h["NewLeadJet_pt3"] 		= new TH1D("NewLeadJet_pt3","",100,0,1500);
	   h["NewLeadJet_eta3"] 	= new TH1D("NewLeadJet_eta3","",100,-2.2,2.2);
	   h["NewLeadJet_btag3"] 	= new TH1D("NewLeadJet_btag3","",100,0,1);

	   h["FatJet_pt1"] 		= new TH1D("FatJet_pt1","",100,0,1500);
	   h["FatJet_pt2"] 		= new TH1D("FatJet_pt2","",100,0,1500);
	   h["FatJet_pt3"] 		= new TH1D("FatJet_pt3","",100,0,1500);
	   h["FatJet_eta1"] 	= new TH1D("FatJet_eta1","",100,-2.2,2.2);
	   h["FatJet_eta2"] 	= new TH1D("FatJet_eta2","",100,-2.2,2.2);
	   h["FatJet_eta3"] 	= new TH1D("FatJet_eta3","",100,-2.2,2.2);
	   h["FatJet_btag3"] 	= new TH1D("FatJet_btag3","",100,0,1);
	   h["FatJet_dR12"] 	= new TH1D("FatJet_dR12","",100,0,6);
	   h["FatJet_dR13"] 	= new TH1D("FatJet_dR13","",100,0,6);
	   h["FatJet_dR23"] 	= new TH1D("FatJet_dR23","",100,0,6);

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

	   h["GenHMass_"]

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

	      if (offlineJets -> size() < 3) continue;

	      //Match offline Jets to online Objects
	      analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());
	      Jet LeadJet[20];
	      std::vector<Jet> NewLeadJet;
	      TLorentzVector FatObj12, Obj12;
	      std::vector<int> third_jet_counter;
	      n_small_jets = 0;
	      counter = 0;
	      goodLeadingJets = false;
	      double btag3Max = 0.8;

	      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
	      	Jet jet = offlineJets -> at(iJet);
	      	if(!jet.idLoose()) continue;
	      	++counter;

	      	LeadJet[counter - 1] = jet;

	      	//Check first two jets:
	      	if(counter == 1 || counter == 2){
	      		if(jet.pt() < analysis.Pt1Cut() ) break;
	      		if(abs(jet.eta()) > 2.2) break;
	      		if(jet.btag() < analysis.BTag1Cut()) break;
	      		NewLeadJet.push_back(jet);
	      		if(counter == 2){
	      			if(abs(LeadJet[0].eta() - LeadJet[1].eta()) > analysis.dEtaCut()) break;
	      			if(LeadJet[0].deltaR(LeadJet[1]) <= 1.1) break;
	      		}
	      	}

	      	//All jets except first 2
	      	if(counter > 2){
	      		//Select third Main jet
	      		if(jet.pt() > 30 && std::abs(jet.eta()) < 2.2 && jet.btag() > analysis.BTag3Cut() && LeadJet[0].deltaR(jet) > 1.1 && LeadJet[1].deltaR(jet) > 1.1 ){
	      			NewLeadJet.push_back(jet);
	      			btag3Max = NewLeadJet.at(2).btag();
	      			third_jet_counter.push_back(counter);
	      			goodLeadingJets = true;
	      		}
	      	}

	      }

	      if(!goodLeadingJets) continue;
	      if(counter < 3) continue;

	      counter = 0;
	      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
		      	Jet jet = offlineJets -> at(iJet);
		      	if(!jet.idLoose()) continue;
		      	++counter;

		      	//All jets except first 2 and 3d Main jet
		      	if(counter > 2 && counter != third_jet_counter.at(0)){
		      		//Construct Wide Jets
	      			Jet closest_lead_jet;
	      			//Add jets to first Main Jet
	      			if (jet.deltaR(NewLeadJet.at(0)) < 1.1 && jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(NewLeadJet.at(1)) && jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(NewLeadJet.at(2))){
	      				NewLeadJet.at(0).p4( NewLeadJet.at(0).p4() + jet.p4() );
	      				closest_lead_jet = NewLeadJet.at(0);
	      			}
	      			//Add jets to second Main Jet
	      			if (jet.deltaR(NewLeadJet.at(1)) < 1.1 && jet.deltaR(NewLeadJet.at(1)) < jet.deltaR(NewLeadJet.at(0)) && jet.deltaR(NewLeadJet.at(1)) < jet.deltaR(NewLeadJet.at(2))){
	      				NewLeadJet.at(1).p4( NewLeadJet.at(1).p4() + jet.p4() );
	      				closest_lead_jet = NewLeadJet.at(1);
	      			}
	      			//Add jets to third Main Jet
	      			if (jet.deltaR(NewLeadJet.at(2)) < 1.1 && jet.deltaR(NewLeadJet.at(2)) < jet.deltaR(NewLeadJet.at(0)) && jet.deltaR(NewLeadJet.at(2)) < jet.deltaR(NewLeadJet.at(1))){
	      				NewLeadJet.at(2).p4( NewLeadJet.at(2).p4() + jet.p4() );
	      				closest_lead_jet = NewLeadJet.at(2);
	      				//Choose btag for third jet
						if(jet.btag() > btag3Max){
							NewLeadJet.at(2).btag(jet.btag());
							btag3Max = jet.btag();
						}
	      			}
	      			//For histograms:
	      			if (jet.deltaR(NewLeadJet.at(0)) < 1.1 || jet.deltaR(NewLeadJet.at(1)) < 1.1 || jet.deltaR(NewLeadJet.at(2)) < 1.1){
		      			h["SmallJet_pt"] ->Fill(jet.pt());
						h["SmallJet_eta"]->Fill(jet.eta());
						h["SmallJet_btag"]->Fill(jet.btag());
						h["SmallJet_dR"] ->Fill(jet.deltaR(closest_lead_jet));
	      			}
	      	}
      }


	      //Fat Jets
	      FatObj12 = NewLeadJet.at(0).p4() + NewLeadJet.at(1).p4();
	      h["Njets"]->Fill(offlineJets->size());
	      h["FatJet_pt1"]->Fill(NewLeadJet.at(0).pt());
	      h["FatJet_pt2"]->Fill(NewLeadJet.at(1).pt());
	      h["FatJet_pt3"]->Fill(NewLeadJet.at(2).pt());
	      h["FatJet_eta1"]->Fill(NewLeadJet.at(0).eta());
	      h["FatJet_eta2"]->Fill(NewLeadJet.at(1).eta());
	      h["FatJet_eta3"]->Fill(NewLeadJet.at(2).eta());
		  h["FatJet_btag3"] -> Fill(NewLeadJet.at(2).btag());
	      h["FatJet_dR12"]->Fill(NewLeadJet.at(0).deltaR(NewLeadJet.at(1)));
	      h["FatJet_dR13"]->Fill(NewLeadJet.at(0).deltaR(NewLeadJet.at(2)));
	      h["FatJet_dR23"]->Fill(NewLeadJet.at(1).deltaR(NewLeadJet.at(2)));

		  h["NumberOf3dJets"] -> Fill(NewLeadJet.size()-2);
		  h["NewLeadJet_pt3"] -> Fill(NewLeadJet.at(2).pt());
		  h["NewLeadJet_eta3"] 	->Fill(NewLeadJet.at(2).eta());
		  h["NewLeadJet_btag3"] -> Fill(NewLeadJet.at(2).btag());
	      h["FatObj12_M"]->Fill(FatObj12.M());
	      h["FatObj12_pt"]->Fill(FatObj12.Pt());
	      h["FatObj12_eta"]->Fill(FatObj12.Eta());

	      Obj12 = LeadJet[0].p4() + LeadJet[1].p4();
		  h["LeadJet_pt1"] -> Fill(LeadJet[0].pt());
		  h["LeadJet_pt2"] -> Fill(LeadJet[1].pt());
		  h["LeadJet_pt3"] -> Fill(LeadJet[2].pt());
		  h["LeadJet_eta1"] -> Fill(LeadJet[0].eta());
		  h["LeadJet_eta2"] -> Fill(LeadJet[1].eta());
		  h["LeadJet_eta3"] -> Fill(LeadJet[2].eta());
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


/*
 	      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
	      {
	      	Jet jet = offlineJets -> at(iJet);

	      	//Only jets that pass Loose identification will be considered
			if(!jet.idLoose()) continue;
			counter++;

			if(counter == 1 || counter == 2 || counter == 3 ){
				if(counter != 3){
					if(jet.pt() < 100 ) break;
					if(abs(jet.eta()) > 2.2) break;
					if(jet.btag() < analysis.BTag1Cut()) break;
					LeadJet[counter - 1] = jet;
				}
				else {
					LeadJet[counter - 1] = jet;
//					if(jet.btag() < analysis.BTag3Cut()) break;
//					if(abs(jet.eta()) > 2.2 ) break;
//					if(jet.pt() < 30) break;
//					if(LeadJet[0].deltaR(LeadJet[2]) <= 1.1) break;
//					if(LeadJet[1].deltaR(LeadJet[2]) <= 1.1) break;
					if(abs(LeadJet[0].eta() - LeadJet[1].eta()) > analysis.dEtaCut()) break;
					if(LeadJet[0].deltaR(LeadJet[1]) <= 1.1) break;

					if(jet.btag() > analysis.BTag3Cut() && jet.pt() > 30 && std::abs(jet.eta()) < 2.2 && LeadJet[0].deltaR(LeadJet[1]) > 1.1 && LeadJet[0].deltaR(jet) > 1.1 && LeadJet[1].deltaR(jet) > 1.1){
						LeadJet[2] = jet;
					}
					FatJet1 = LeadJet[0].p4();
					FatJet2 = LeadJet[1].p4();
//					goodLeadingJets = true;
				}
			}
			if(counter >= 3 && jet.pt() > 30 && std::abs(jet.eta()) < 2.2 && jet.btag() > analysis.BTag3Cut()){
				if(jet.deltaR(LeadJet[1]) > 1.1 && jet.deltaR(LeadJet[0]) > 1.1){
					NewLeadJet.push_back(jet);
					FatJet3 = NewLeadJet.at(0).p4();
					btag3 = NewLeadJet.at(0).btag();
					goodLeadingJets = true;
				}
			}



			if (counter > 3){
				Jet closest_lead_jet;
				if (jet.deltaR(LeadJet[0]) < jet.deltaR(LeadJet[1]) && jet.deltaR(LeadJet[0]) < jet.deltaR(LeadJet[2])) closest_lead_jet = LeadJet[0];
				else if (jet.deltaR(LeadJet[1]) < jet.deltaR(LeadJet[0]) && jet.deltaR(LeadJet[1]) < jet.deltaR(LeadJet[2])) closest_lead_jet = LeadJet[1];
				else if (jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(LeadJet[0]) && jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(LeadJet[1])) closest_lead_jet = NewLeadJet.at(0);

				if(jet.pt() > 30){

					if (jet.deltaR(LeadJet[0]) < 1.1 || jet.deltaR(LeadJet[1]) < 1.1 || jet.deltaR(NewLeadJet.at(0)) < 1.1){
						h["SmallJet_pt"] ->Fill(jet.pt());
						h["SmallJet_eta"]->Fill(jet.eta());
						h["SmallJet_btag"]->Fill(jet.btag());
						h["SmallJet_dR"] ->Fill(jet.deltaR(closest_lead_jet));
						h["SmallJet_dR_zoom"] ->Fill(jet.deltaR(closest_lead_jet));
						n_small_jets ++;
					}

					if (jet.deltaR(LeadJet[0]) < 1.1 && jet.deltaR(LeadJet[0]) < jet.deltaR(LeadJet[1]) && jet.deltaR(LeadJet[0]) < jet.deltaR(NewLeadJet.at(0))) {
						FatJet1 += jet.p4();
					}
					else if (jet.deltaR(LeadJet[1]) < 1.1 && jet.deltaR(LeadJet[1]) < jet.deltaR(LeadJet[0]) && jet.deltaR(LeadJet[1]) < jet.deltaR(NewLeadJet.at(0))){
						FatJet2 += jet.p4();
					}
					else if (jet.deltaR(NewLeadJet.at(0)) < 1.1 && jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(LeadJet[0]) && jet.deltaR(NewLeadJet.at(0)) < jet.deltaR(LeadJet[1])){
						FatJet3 += jet.p4();
						if(jet.btag() > btag3){
							NewLeadJet.at(0).btag(jet.btag());
							btag3 = jet.btag();
						}
				}

			}
			}

	      }
	     */

