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

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;
double weight(const float &pt1, const float &pt2, const float &dEta, const float &btag1, const float &btag2);

// =============================================================================================
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   bool jetht = 0;
   //For btag efficiency:
   TFile * fileBTagEff = new TFile("~/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/MssmHbb/bin/outputFiles/bTagEfficiency2D_pteta.root");
   TH2F *btagEff = (TH2F*) fileBTagEff ->Get("h2ehn");

   // Input files list MonteCarloStudies
   std::string inputList = "rootFileListBTagCSV.txt";
   MssmHbb analysis(inputList);
   
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
	std::vector<std::string> jetTriggerObjects;
   if(!jetht) {
   
   analysis.addTree<TriggerObject> ("hltL1sL1DoubleJetC100","MssmHbb/Events/selectedPatTrigger/hltL1sL1DoubleJetC100");
   analysis.addTree<TriggerObject> ("hltDoubleJetsC100","MssmHbb/Events/selectedPatTrigger/hltDoubleJetsC100");
   analysis.addTree<TriggerObject> ("hltDoublePFJetsC100","MssmHbb/Events/selectedPatTrigger/hltDoublePFJetsC100");
   analysis.addTree<TriggerObject> ("hltDoubleBTagCSV0p9","MssmHbb/Events/selectedPatTrigger/hltDoubleBTagCSV0p9");
   analysis.addTree<TriggerObject> ("hltDoublePFJetsC100MaxDeta1p6","MssmHbb/Events/selectedPatTrigger/hltDoublePFJetsC100MaxDeta1p6");
   analysis.triggerResults("MssmHbb/Events/TriggerResults");
   std::string triggerLogicName = "HLT_DoubleJetsC100_DoubleBTagCSV0p9_DoublePFJetsC100MaxDeta1p6_v";
   
   jetTriggerObjects.push_back("hltL1sL1DoubleJetC100");
   jetTriggerObjects.push_back("hltDoubleJetsC100");
   jetTriggerObjects.push_back("hltDoubleBTagCSV0p9");
   jetTriggerObjects.push_back("hltDoublePFJetsC100");
   jetTriggerObjects.push_back("hltDoublePFJetsC100MaxDeta1p6");
   }
   else {
   analysis.addTree<TriggerObject> ("hltL1sL1SingleJet52","MssmHbb/Events/selectedPatTrigger/hltL1sL1SingleJet52");
   analysis.addTree<TriggerObject> ("hltSinglePFJet80","MssmHbb/Events/selectedPatTrigger/hltSinglePFJet80");
   analysis.addTree<TriggerObject> ("hltSingleCaloJet50","MssmHbb/Events/selectedPatTrigger/hltSingleCaloJet50");
   
   jetTriggerObjects.push_back("hltL1sL1SingleJet52");
   jetTriggerObjects.push_back("hltSinglePFJet80");
   jetTriggerObjects.push_back("hltSingleCaloJet50");
   }

   
   if(!analysis.isMC()) analysis.processJsonFile("goodJson.txt");
   
   int nCand = 0;
   std::vector <int> leadJetsNumber;
   
   std::map < std::string, TH1F *> h1;
   h1["h_jet_N"]         = new TH1F("h_jet_N" , "", 20, 0., 20.);
	
	h1["h_jet1_Pt"]       = new TH1F("h_jet1_Pt" , "", 100, 0., 1000.);
	h1["h_jet1_Eta"]      = new TH1F("h_jet1_Eta", "", 100, -5, 5.);
	h1["h_jet1_Phi"]      = new TH1F("h_jet1_Phi", "", 100, -3.2, 3.2);
	h1["h_jet1_BTag"]     = new TH1F("h_jet1_BTag", "", 100, 0., 1.);
	
   	h1["h_jet2_Pt"]       = new TH1F("h_jet2_Pt" , "", 100, 0., 1000.);
   	h1["h_jet2_Eta"]      = new TH1F("h_jet2_Eta", "", 100, -5, 5.);
  	h1["h_jet2_Phi"]      = new TH1F("h_jet2_Phi", "", 100, -3.2, 3.2);
  	h1["h_jet2_BTag"]     = new TH1F("h_jet2_BTag", "", 100, 0., 1.);
  	
	h1["h_obj_M"]  		  = new TH1F("h_obj_M12","",100,0.,1000.);	
	
	h1["h_jet_dEta"]	  = new TH1F("h_jet_dEta","",100,0,3);
	                
   
   // Analysis of events
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
   	
      analysis.event(i);
      
      if (!analysis.isMC() && !analysis.selectJson()) continue;
      auto offlineJets = analysis.collection<Jet>("Jets");
      if (offlineJets -> size() < 2) continue;
      if(!analysis.isMC()) analysis.match<Jet,TriggerObject>("Jets",jetTriggerObjects);
      
      
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
      
      if( fLeadJet.pt() < 80 || sLeadJet.pt() < 80 ) continue;
      if( std::abs(fLeadJet.eta()) > 2.2 || std::abs(sLeadJet.eta()) > 2.2 ) continue;
      if( fLeadJet.deltaR(sLeadJet) <= 1) continue;
      //if( fLeadJet.btag() < 0.941 || sLeadJet.btag() < 0.941 ) continue;
      
      unsigned int trigger_selection = 0;
      for( const auto &it : jetTriggerObjects)
      {
      	if(fLeadJet.matched(it) && sLeadJet.matched(it )) trigger_selection++;
      }
      if(trigger_selection != jetTriggerObjects.size()) continue;
      if(!jetht){
      bool status = true;
      	if(fLeadJet.matched("hltDoublePFJetsC100MaxDeta1p6") && sLeadJet.matched("hltDoublePFJetsC100MaxDeta1p6"))
		{
			const Candidate *onlineJet1 = fLeadJet.matched("hltDoublePFJetsC100MaxDeta1p6");
			const Candidate *onlineJet2 = sLeadJet.matched("hltDoublePFJetsC100MaxDeta1p6");
			if( std::abs(onlineJet1->eta() - onlineJet2-> eta()) > 1.6) status = false; //Check Pt and dEta condition	else return false;
		}
		if(!status) continue;
      } 
      
     double reweight;
     if(!analysis.isMC()) reweight = 1.;
     else reweight = weight(fLeadJet.pt(),sLeadJet.pt(),std::abs(fLeadJet.eta() - sLeadJet.eta()), btagEff -> Interpolate(fLeadJet.pt(),fLeadJet.eta()), btagEff -> Interpolate(sLeadJet.pt(),sLeadJet.eta()));
      
     h1["h_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     h1["h_jet1_Eta"] -> Fill(fLeadJet.eta(),reweight);
     h1["h_jet1_Phi"] -> Fill(fLeadJet.phi(),reweight);
     h1["h_jet1_BTag"] -> Fill(fLeadJet.btag(),reweight);
     
     h1["h_jet2_Pt"] -> Fill(sLeadJet.pt(),reweight);
     h1["h_jet2_Eta"] -> Fill(sLeadJet.eta(),reweight);
     h1["h_jet2_Phi"] -> Fill(sLeadJet.phi(),reweight);
     h1["h_jet2_BTag"] -> Fill(sLeadJet.btag(),reweight);
     
     h1["h_jet_dEta"] -> Fill(std::abs(fLeadJet.eta() - sLeadJet.eta()),reweight);
     
     TLorentzVector obj;
     obj = fLeadJet.p4() + sLeadJet.p4();
     

     h1["h_obj_M"] -> Fill(obj.M(),reweight);
	 
	     
      ++nCand;

   }
   std::cout<<"Number of Candidates: "<<nCand<<std::endl;
   std::string fileName = "output.root";
   fileName = analysis.tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
   fileName = "outputFiles/SelectionDataComp" + boost::filesystem::basename(fileName) + ".root";
   TFile * outFile = new TFile(fileName.c_str(),"RECREATE");
   
   //double lumiData2015DPrRec3 = 553.149951;
   double lumiData2015DPrRec4 = 1001.179266; 
   for ( auto& ih1 : h1 )
   {
   	if(analysis.isMC()) ih1.second->Scale(lumiData2015DPrRec4/analysis.luminosity());	
   	ih1.second -> Write();
   }
   outFile -> Close();
   
   std::cout<<" Output File: "<<fileName<<" was stored!"<<std::endl; 
   
   return 0;
//
}

double weight(const float &pt1, const float &pt2, const float &dEta, const float &btag1, const float &btag2)
{
	double weightPt2 = 1.;
	double weightPt1 = 1.;
	//weightPt1 = 1./(1.+exp( -8.76241e-02 * (pt1 - 1.05187e+02) )) * 1./ (1.+exp( -1.88955e-01 * (pt1 - 1.02451e+02) ));
	//weightPt2 = 1./(1.+exp( -8.76241e-02 * (pt2 - 1.05187e+02) )) * 1./ (1.+exp( -1.88955e-01 * (pt2 - 1.02451e+02) ));
	double weightdEta = 1./(1.+exp( 1.15733e+02 * (dEta - 1.59862e+00) )) * 1./ (1.+exp( -2.06004e-02 * (dEta + 2.22684e+02) ));
	double result;
	if( pt1 >= 1000.) result = weightPt1 * weightPt2 * weightdEta;
	else result = weightPt1 * weightPt2 * weightdEta * btag1 * btag2;
	return result;
}
