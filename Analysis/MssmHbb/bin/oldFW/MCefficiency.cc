#include <string>
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>
#include <stdlib.h>

#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"
#include "TRandom1.h"

#include "Analysis/MssmHbb/interface/TriggerEfficiency.h"
#include "Analysis/MssmHbb/interface/json.h"

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;


// =============================================================================================
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // Input files list MonteCarloStudies
   std::string inputList = "rootFileList.txt";//"/nfs/dust/cms/user/shevchen/samples/miniaod/JetHT/JetHTFileList.txt";
   TriggerEfficiency analysis(inputList,"MssmHbb/Events/EventInfo");
   
   if(analysis.isMC())
   {
     // cross sections
     analysis.crossSections("MssmHbb/Metadata/CrossSections");
     analysis.listCrossSections();

     //generator filter
     analysis.generatorFilter("MssmHbb/Metadata/GeneratorFilter");
     analysis.listGeneratorFilter();
   }

   analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsPuppi");
   if(analysis.isMC()) analysis.addTree<Jet> ("genJets","MssmHbb/Events/slimmedGenJets");

   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   analysis.addTriggerObjects();  
   
   //std::string triggerLogicName = "HLT_PFJet60_v2";
   
   int nCand = 0;
   
   std::map < std::string, TH1F *> h1;
   h1["h_jet_N"]           		= new TH1F("h_jet_N" , "", 20, 0., 20.);
	
   h1["h_Only1jet_jet1_Pt"]     = new TH1F("h_Only1jet_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_thirdJet_Pt"]			= new TH1F("h_thirdJet_Pt","",100,0.,1000.);
   h1["h_fourthJet_Pt"]			= new TH1F("h_fourthJet_Pt","",100,0.,1000.);\
   h1["dPhi"]					= new TH1F("dPhi","",100,-2.*TMath::Pi(),2.*TMath::Pi());

   //Histograms for Efficiency
   h1["hpf40_eff"]				= new TH1F("hpf40_eff","",200,0.,500.);
   h1["hpf40_effNumerator"]  	= new TH1F("hpf40_effNumerator","",200,0.,500.);
   h1["hpf40_effDenumerator"]	= new TH1F("hpf40_effDenumerator","",200,0.,500.);
   
   h1["hpf60_eff"]				= new TH1F("hpf60_eff","",200,0.,500.);
   h1["hpf60_effNumerator"]  	= new TH1F("hpf60_effNumerator","",200,0.,500.);
   h1["hpf60_effDenumerator"]	= new TH1F("hpf60_effDenumerator","",200,0.,500.);
   
   h1["hpf80_eff"]				= new TH1F("hpf80_eff","",200,0.,500.);
   h1["hpf80_effNumerator"]  	= new TH1F("hpf80_effNumerator","",200,0.,500.);
   h1["hpf80_effDenumerator"]	= new TH1F("hpf80_effDenumerator","",200,0.,500.);
   
   h1["hpf100_effMCmodel_Num"]		= new TH1F("hpf100_effMCmodel_Num","",200,0.,500.);
   h1["hpf100_effMCmodel_Denum"]	= new TH1F("hpf100_effMCmodel_Denum","",200,0.,500.);

   h1["hpf40for60_eff"]				= new TH1F("hpf40for60_eff","",200,0.,500.);
   h1["hpf40for60_effNumerator"]  	= new TH1F("hpf40for60_effNumerator","",200,0.,500.);
   h1["hpf40for60_effDenumerator"]	= new TH1F("hpf40for60_effDenumerator","",200,0.,500.);
   
   h1["hpf60for80_eff"]				= new TH1F("hpf60for80_eff","",200,0.,500.);
   h1["hpf60for80_effNumerator"]  	= new TH1F("hpf60for80_effNumerator","",200,0.,500.);
   h1["hpf60for80_effDenumerator"]	= new TH1F("hpf60for80_effDenumerator","",200,0.,500.);
   
   h1["T&P_probe40"]				= new TH1F("T&P_probe40","",200,0.,500.);
   h1["T&P_tag40"]					= new TH1F("T&P_tag40","",200,0.,500.);
   h1["T&Ppf40_eff"]				= new TH1F("T&Ppf40_eff","",200,0.,500.);
   
   h1["T&P_probe60"]				= new TH1F("T&P_probe60","",200,0.,500.);
   h1["T&P_tag60"]					= new TH1F("T&P_tag60","",200,0.,500.);
   h1["T&Ppf60_eff"]				= new TH1F("T&Ppf60_eff","",200,0.,500.);
   
   h1["T&P_probe80"]				= new TH1F("T&P_probe80","",200,0.,500.);
   h1["T&P_tag80"]					= new TH1F("T&P_tag80","",200,0.,500.);
   h1["T&Ppf80_eff"]				= new TH1F("T&Ppf80_eff","",200,0.,500.);
   
   h1["T&P_probe100"]				= new TH1F("T&P_probe100","",200,0.,500.);
   h1["T&P_tag100"]					= new TH1F("T&P_tag100","",200,0.,500.);
   h1["T&Ppf100_eff"]				= new TH1F("T&Ppf100_eff","",200,0.,500.);
   h1["N_pf80_l1"] 				    = new TH1F("N_pf80_l1","",10,0,10);
   h1["N_pf80_l2"]                  = new TH1F("N_pf80_l2","",10,0,10);
   h1["N_pf80_l3"]                  = new TH1F("N_pf80_l3","",10,0,10);

   if(!analysis.isMC()) analysis.processJsonFile("goodJson.txt");
   std::cout<<"Total Entries: "<<analysis.size()<<std::endl;
   
   const Candidate  *l2p5MPF40Cand, *l3MPF40Cand, *l1MPF60Cand, *l2p5MPF60Cand, *l3MPF60Cand, *l1MPF80Cand, *l2p5MPF80Cand, *l3MPF80Cand, 
   				    *J2l2p5MPF40Cand, *J2l3MPF40Cand, *J2l1MPF60Cand, *J2l2p5MPF60Cand, *J2l3MPF60Cand, *J2l1MPF80Cand, *J2l2p5MPF80Cand, *J2l3MPF80Cand;
   
	TRandom1 rand;
	int rand_val = 0;
	bool PtVeto = false, goodLeadingJets = false;
	int counter;
   // Analysis of events
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
   		
      analysis.event(i);
      if (!analysis.isMC() && !analysis.selectJson()) continue;
      //if( !analysis.triggerResult("HLT_PFJet40_v") && !analysis.triggerResult("HLT_PFJet60_v") && !analysis.triggerResult("HLT_PFJet80_v")) continue;
      
      // Selection of double jet sample:
      auto jets = analysis.collection<Jet>("Jets");
      analysis.match<Jet,TriggerObject>("Jets",analysis.triggerObjectNames());

      counter = 0;
      goodLeadingJets = true;
      PtVeto = true;
      Jet LeadJet[2];
      if(jets->size() < 2) continue;
      for( int iJet = 0; iJet < jets -> size(); ++iJet)
      {
      	Jet jet = jets -> at(iJet);
		if(!jet.idLoose()) continue;
		counter++;
		if(counter == 1 || counter == 2){
			LeadJet[counter - 1] = jet;
			if(jet.pt() < 90) goodLeadingJets = false;
			if(abs(jet.eta()) > 2) goodLeadingJets = false;
			if(!goodLeadingJets) break;
		}
		else{
			if(jet.pt() > 40) PtVeto = false;
			if(!PtVeto) break;
		}
	  }
      if(counter < 2) continue;
      if(!PtVeto) continue;
      if(abs(LeadJet[0].phi() - LeadJet[1].phi()) < 2.7 || abs(LeadJet[0].phi() - LeadJet[1].phi()) > TMath::Pi() ) goodLeadingJets = false;
      if((LeadJet[0].pt() - LeadJet[1].pt()) / (LeadJet[0].pt() + LeadJet[1].pt()) > 0.2 ) goodLeadingJets = false;
      if (!goodLeadingJets) continue;
      h1["dPhi"] -> Fill(abs(LeadJet[0].phi() - LeadJet[1].phi()));
      // Size ob Trigger Objects
      h1["N_pf80_l1"] -> Fill(analysis.collection<TriggerObject>("hltL1sL1SingleJet52")->size());
      h1["N_pf80_l2"] -> Fill(analysis.collection<TriggerObject>("hltSingleCaloJet50")->size());
      h1["N_pf80_l3"] -> Fill(analysis.collection<TriggerObject>("hltSinglePFJet80")->size());

      //Trigger matching:
      //l1MPF40Cand = LeadJet[0].matched("hltL1sL1ZeroBias");
      l2p5MPF40Cand = LeadJet[0].matched("hltSingleCaloJet10");
      l3MPF40Cand = LeadJet[0].matched("hltSinglePFJet40");
      
      l1MPF60Cand = LeadJet[0].matched("hltL1sL1SingleJet36");
      l2p5MPF60Cand = LeadJet[0].matched("hltSingleCaloJet40");
      l3MPF60Cand = LeadJet[0].matched("hltSinglePFJet60");
      
      l1MPF80Cand = LeadJet[0].matched("hltL1sL1SingleJet52");
      l2p5MPF80Cand = LeadJet[0].matched("hltSingleCaloJet50");
      l3MPF80Cand = LeadJet[0].matched("hltSinglePFJet80");
	  //Prepare Pt efficiency calculation:
	  
	  if(analysis.triggerResult("HLT_PFJet40_v") && l2p5MPF40Cand && l3MPF40Cand && l2p5MPF40Cand->pt() > 100. && l3MPF40Cand->pt() > 100. ) h1["hpf40_effNumerator"] -> Fill(LeadJet[0].pt());
	  if(analysis.triggerResult("HLT_PFJet40_v") ) h1["hpf40_effDenumerator"] -> Fill(LeadJet[0].pt());
	  
	  if(analysis.triggerResult("HLT_PFJet60_v") && l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand && l1MPF60Cand->pt() > 100. && l2p5MPF60Cand->pt() > 100. && l3MPF60Cand->pt() > 100. && LeadJet[0].pt() > 80) h1["hpf60_effNumerator"] -> Fill(LeadJet[0].pt());
	  if(analysis.triggerResult("HLT_PFJet60_v") && LeadJet[0].pt() > 80) h1["hpf60_effDenumerator"] -> Fill(LeadJet[0].pt());
	  
	  if(analysis.triggerResult("HLT_PFJet80_v") &&  l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && l1MPF80Cand->pt() > 100. && l2p5MPF80Cand->pt() > 100. && l3MPF80Cand->pt() > 100. && LeadJet[0].pt() > 80) h1["hpf80_effNumerator"] -> Fill(LeadJet[0].pt());
	  if(analysis.triggerResult("HLT_PFJet80_v") && LeadJet[0].pt() > 80) h1["hpf80_effDenumerator"] -> Fill(LeadJet[0].pt());
	  
	  if(analysis.triggerResult("HLT_PFJet80_v") &&  l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && l1MPF80Cand->pt() > 100. && l2p5MPF80Cand->pt() > 100. && l3MPF80Cand->pt() > 100. && LeadJet[0].pt() > 80) h1["hpf100_effMCmodel_Num"] -> Fill(LeadJet[0].pt());
	  if(LeadJet[0].pt() > 80) h1["hpf100_effMCmodel_Denum"] -> Fill(LeadJet[0].pt());

	  if(analysis.triggerResult("HLT_PFJet40_v")  && l2p5MPF40Cand && l3MPF40Cand && l2p5MPF40Cand->pt() > 40. && l3MPF40Cand->pt() > 60. ) h1["hpf40for60_effNumerator"] -> Fill(LeadJet[0].pt());
	  if(analysis.triggerResult("HLT_PFJet40_v") ) h1["hpf40for60_effDenumerator"] -> Fill(LeadJet[0].pt());
	  
	  if(analysis.triggerResult("HLT_PFJet60_v") && l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand && l1MPF60Cand->pt() > 52. && l2p5MPF60Cand->pt() > 50. && l3MPF60Cand->pt() > 80. ) h1["hpf60for80_effNumerator"] -> Fill(LeadJet[0].pt());
	  if(analysis.triggerResult("HLT_PFJet60_v") ) h1["hpf60for80_effDenumerator"] -> Fill(LeadJet[0].pt());
	 
	  if(l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand) h1["h_Only1jet_jet1_Pt"]  -> Fill(LeadJet[0].pt());
      
      rand_val = round( rand.Gaus(0.5,0.5));
      Jet tagJet, probeJet;
      if (rand_val < 1)
      {
      	tagJet = LeadJet[0];
      	probeJet = LeadJet[1];
      }
      else
      {
      	tagJet = LeadJet[1];
      	probeJet = LeadJet[0];
      }
      
      //l1MPF40Cand = tagJet.matched("hltL1sL1ZeroBias");	
      l2p5MPF40Cand = tagJet.matched("hltSingleCaloJet10");
      l3MPF40Cand = tagJet.matched("hltSinglePFJet40");
      
      l1MPF60Cand = tagJet.matched("hltL1sL1SingleJet36");	
      l2p5MPF60Cand = tagJet.matched("hltSingleCaloJet40");
      l3MPF60Cand = tagJet.matched("hltSinglePFJet60");
      
      l1MPF80Cand = tagJet.matched("hltL1sL1SingleJet52");
      l2p5MPF80Cand = tagJet.matched("hltSingleCaloJet50");
      l3MPF80Cand = tagJet.matched("hltSinglePFJet80");
      
      //Trigger matching:
      //J2l1MPF40Cand = probeJet.matched("hltL1sL1ZeroBias");	
      J2l2p5MPF40Cand = probeJet.matched("hltSingleCaloJet10");
      J2l3MPF40Cand = probeJet.matched("hltSinglePFJet40");
      
      J2l1MPF60Cand = probeJet.matched("hltL1sL1SingleJet36");	
      J2l2p5MPF60Cand = probeJet.matched("hltSingleCaloJet40");
      J2l3MPF60Cand = probeJet.matched("hltSinglePFJet60");
      
      J2l1MPF80Cand = probeJet.matched("hltL1sL1SingleJet52");
      J2l2p5MPF80Cand = probeJet.matched("hltSingleCaloJet50");
      J2l3MPF80Cand = probeJet.matched("hltSinglePFJet80");
      
      if( analysis.triggerResult("HLT_PFJet60_v") && l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand && J2l1MPF60Cand && J2l2p5MPF60Cand && J2l3MPF60Cand) h1["T&P_probe60"] -> Fill(probeJet.pt()); //probe
      if( analysis.triggerResult("HLT_PFJet60_v") && l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand) h1["T&P_tag60"] -> Fill(probeJet.pt()); // tag
      
      if( analysis.triggerResult("HLT_PFJet40_v") && l2p5MPF40Cand && l3MPF40Cand && J2l2p5MPF40Cand && J2l3MPF40Cand) h1["T&P_probe40"] -> Fill(probeJet.pt()); //probe
      if( analysis.triggerResult("HLT_PFJet40_v") && l2p5MPF40Cand && l3MPF40Cand) h1["T&P_tag40"] -> Fill(probeJet.pt()); // tag
      
      if( analysis.triggerResult("HLT_PFJet80_v") && l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && J2l1MPF80Cand && J2l2p5MPF80Cand && J2l3MPF80Cand) h1["T&P_probe80"] -> Fill(probeJet.pt()); //probe
      if( analysis.triggerResult("HLT_PFJet80_v") && l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand) h1["T&P_tag80"] -> Fill(probeJet.pt()); // tag
    
      if( analysis.triggerResult("HLT_PFJet80_v") && l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && l1MPF80Cand->pt() > 100 && l2p5MPF80Cand -> pt() > 100. && l3MPF80Cand->pt() > 100.)
      {
      	if(J2l1MPF80Cand && J2l2p5MPF80Cand && J2l3MPF80Cand && J2l1MPF80Cand->pt() > 100. && J2l2p5MPF80Cand->pt() > 100. && J2l3MPF80Cand -> pt() > 100.) h1["T&P_probe100"] -> Fill(probeJet.pt()); //probe
      	h1["T&P_tag100"] -> Fill(probeJet.pt()); // tag
      }

      if(counter >= 3) h1["h_thirdJet_Pt"]->Fill(jets->at(2).pt());
      if(counter >= 4) h1["h_fourthJet_Pt"] -> Fill(jets->at(3).pt());

	  
           ++nCand;

   }
   h1["hpf40_eff"] -> Divide(h1["hpf40_effNumerator"],h1["hpf40_effDenumerator"],1.,1.,"b");
   h1["hpf60_eff"] -> Divide(h1["hpf60_effNumerator"],h1["hpf60_effDenumerator"],1.,1.,"b");
   h1["hpf80_eff"] -> Divide(h1["hpf80_effNumerator"],h1["hpf80_effDenumerator"],1.,1.,"b");
   
   h1["hpf40for60_eff"] -> Divide(h1["hpf40for60_effNumerator"],h1["hpf40for60_effDenumerator"],1.,1.,"b");
   h1["hpf60for80_eff"] -> Divide(h1["hpf60for80_effNumerator"],h1["hpf60for80_effDenumerator"],1.,1.,"b");
   
   h1["T&Ppf40_eff"] -> Divide(h1["T&P_probe40"],h1["T&P_tag40"],1.,1.,"b");
   h1["T&Ppf60_eff"] -> Divide(h1["T&P_probe60"],h1["T&P_tag60"],1.,1.,"b");
   h1["T&Ppf80_eff"] -> Divide(h1["T&P_probe80"],h1["T&P_tag80"],1.,1.,"b");
   h1["T&Ppf100_eff"] -> Divide(h1["T&P_probe100"],h1["T&P_tag100"],1.,1.,"b");
   
   std::cout<<"Number of Candidates: "<<nCand<<std::endl;
   std::string fileName = "output.root";
   fileName = analysis.tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
   if(analysis.isMC()) fileName = "outputFiles/ResultsDoubleSelection" + boost::filesystem::basename(fileName) + ".root";
   else fileName = "outputFiles/ResultsDoubleSelection" + boost::filesystem::basename(fileName) + ".root";
   //fileName = "output.root";
   TFile * outFile = new TFile(fileName.c_str(),"RECREATE");
   outFile->mkdir("HLT_PFJet60");
   outFile->mkdir("HLT_PFJet80");
   double lumiData2015DPrRec4 = 1001.179266; 
   for ( auto& ih1 : h1 )
   {
   	if(analysis.isMC()) ih1.second->Scale(lumiData2015DPrRec4/analysis.luminosity());
   	if(ih1.first.find("hpf60") != std::string::npos)
   	outFile -> cd("HLT_PFJet60");
   	if(ih1.first.find("hpf80") != std::string::npos)
   	outFile -> cd("HLT_PFJet80");
   	ih1.second -> Write();
   	outFile -> cd("");
   }
  
   outFile -> Close();
   
   
   std::cout<<"File: "<<fileName<<" was created"<<std::endl;
   
   return 0;
//
}
