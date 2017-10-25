#include <string>
#include <iostream>
#include <vector>
#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TEfficiency.h"

#include "Analysis/MssmHbb/interface/MssmHbb.h"
#include "Analysis/MssmHbb/interface/json.h"

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;

double weight(const float &,const float &);
double ptRew2D(TH2F *, const float &, const float &);
double bTagRew(TH2F *, const float&, const float &);

// =============================================================================================
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // Input files list MonteCarloStudies
   std::string inputList = "rootFileList.txt";
   MssmHbb analysis(inputList,"MssmHbb/Events/EventInfo");
   //For Pt,eta efficiency:
   TFile * filePtEff = new TFile("~/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/MssmHbb/bin/outputFiles/ResultsJetHT2015DPromptReco4.root");
   TH2F *PtEff = (TH2F*) filePtEff ->Get("h2_eff80");
   //For btag efficiency:
   TFile * fileBTagEff = new TFile("~/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/MssmHbb/bin/outputFiles/bTagEfficiency2D_pteta_zomin.root");
   TH2F *btagEff = (TH2F*) fileBTagEff ->Get("h2ehn");
   
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
   // HLT_PFJet60 object
   analysis.addTree<TriggerObject> ("hltL1sL1SingleJet36","MssmHbb/Events/selectedPatTrigger/hltL1sL1SingleJet36");
   analysis.addTree<TriggerObject> ("hltSinglePFJet60","MssmHbb/Events/selectedPatTrigger/hltSinglePFJet60");
   analysis.addTree<TriggerObject> ("hltSingleCaloJet40","MssmHbb/Events/selectedPatTrigger/hltSingleCaloJet40");
   analysis.triggerResults("MssmHbb/Events/TriggerResults");

   // HLT_PFJet80 objects
   analysis.addTree<TriggerObject> ("hltL1sL1SingleJet52","MssmHbb/Events/selectedPatTrigger/hltL1sL1SingleJet52");
   analysis.addTree<TriggerObject> ("hltSinglePFJet80","MssmHbb/Events/selectedPatTrigger/hltSinglePFJet80");
   analysis.addTree<TriggerObject> ("hltSingleCaloJet50","MssmHbb/Events/selectedPatTrigger/hltSingleCaloJet50");
   
   std::vector<std::string> jetTriggerObjects;
   
   jetTriggerObjects.push_back("hltL1sL1SingleJet36");
   jetTriggerObjects.push_back("hltSinglePFJet60");
   jetTriggerObjects.push_back("hltSingleCaloJet40");
   
   jetTriggerObjects.push_back("hltL1sL1SingleJet52");
   jetTriggerObjects.push_back("hltSinglePFJet80");
   jetTriggerObjects.push_back("hltSingleCaloJet50");
   
   
   //std::string triggerLogicName = "HLT_PFJet60_v2";
   
   int nCand = 0;
   std::vector <int> leadJetsNumber;
   
   
   std::map < std::string, TH1F *> h1;
   h1["h_jet_N"]           		= new TH1F("h_jet_N" , "", 20, 0., 20.);
	
   h1["h_Only1jet_jet1_Pt"]     = new TH1F("h_Only1jet_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_Only1jet_jet1_Eta"]    = new TH1F("h_Only1jet_jet1_Eta", "", 100, -5, 5.);
   h1["h_Only1jet_jet1_Phi"]    = new TH1F("h_Only1jet_jet1_Phi", "", 100, -3.2, 3.2);
   
   h1["h_Only1jet_jet1_PtBTag"] = new TH1F("h_Only1jet_jet1_PtBTag","",50,0.,1000.);
   h1["h_Only1jet_jet1_PtAll"] = new TH1F("h_Only1jet_jet1_PtAll","",50,0.,1000.);
   
   
   h1["h_TwoJets_jet1_Pt"]      = new TH1F("h_TwoJets_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_TwoJets_jet1_Eta"]     = new TH1F("h_TwoJets_jet1_Eta", "", 100, -5, 5.);
   h1["h_TwoJets_jet1_Phi"]     = new TH1F("h_TwoJets_jet1_Phi", "", 100, -3.2, 3.2);
   h1["h_TwoJets_jet2_Pt"]      = new TH1F("h_TwoJets_jet2_Pt" , "", 100, 0., 1000.);
   h1["h_TwoJets_jet2_Eta"]     = new TH1F("h_TwoJets_jet2_Eta", "", 100, -5, 5.);
   h1["h_TwoJets_jet2_Phi"]     = new TH1F("h_TwoJets_jet2_Phi", "", 100, -3.2, 3.2);
   
   h1["h_TwoBtagJets_jet1_Pt"]      = new TH1F("h_TwoBtagJets_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_TwoBtagJets_jet1_Eta"]     = new TH1F("h_TwoBtagJets_jet1_Eta", "", 100, -5, 5.);
   h1["h_TwoBtagJets_jet1_Phi"]     = new TH1F("h_TwoBtagJets_jet1_Phi", "", 100, -3.2, 3.2);
   h1["h_TwoBtagJets_jet2_Pt"]      = new TH1F("h_TwoBtagJets_jet2_Pt" , "", 100, 0., 1000.);
   h1["h_TwoBtagJets_jet2_Eta"]     = new TH1F("h_TwoBtagJets_jet2_Eta", "", 100, -5, 5.);
   h1["h_TwoBtagJets_jet2_Phi"]     = new TH1F("h_TwoBtagJets_jet2_Phi", "", 100, -3.2, 3.2);
   
   h1["h_jet1_BTag"]       		 = new TH1F("h_jet1_BTag", "", 100, 0., 1.);
   
   //Pt spectra study
   h1["h_PtStudy_ptCut_jet1_Pt"] = new TH1F("h_PtStudy_ptCut_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_PtStudy_etaCut_jet1_Pt"]= new TH1F("h_PtStudy_etaCut_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_PtStudy_dRCut_jet1_Pt"] = new TH1F("h_PtStudy_dRCut_jet1_Pt" , "", 100, 0., 1000.);
   h1["h_PtStudy_AllCut_jet1_Pt"]= new TH1F("h_PtStudy_AllCut_jet1_Pt" , "", 100, 0., 1000.);

   //Histograms for PF60
   h1["hpf60_l1j1_Pt"]     		= new TH1F("hpf60_l1j1_Pt","",100,0.,500.);
   h1["hpf60_l1j1_Eta"]    		= new TH1F("hpf60_l1j1_Eta","",100,-5.,5.);
   h1["hpf60_l1j1_Phi"]    		= new TH1F("hpf60_l1j1_Phi", "",100,-3.2,3.2);

   h1["hpf60_l2p5j1_Pt"]   		= new TH1F("hpf60_l2p5j1_Pt","",100,0.,500.);
   h1["hpf60_l2p5j1_Eta"]  		= new TH1F("hpf60_l2p5j1_Eta","",100,-5.,5.);
   h1["hpf60_l2p5j1_Phi"]  		= new TH1F("hpf60_l2p5j1_Phi", "",100,-3.2,3.2);

   h1["hpf60_l3j1_Pt"]       	= new TH1F("hpf60_l3j1_Pt","",100,0.,500.);
   h1["hpf60_l3j1_Eta"]      	= new TH1F("hpf60_l3j1_Eta","",100,-5.,5.);
   h1["hpf60_l3j1_Phi"]      	= new TH1F("hpf60_l3j1_Phi", "",100,-3.2,3.2);
   
   //Histograms for PF80
   h1["hpf80_l1j1_Pt"]     		= new TH1F("hpf80_l1j1_Pt","",100,0.,500.);
   h1["hpf80_l1j1_Eta"]    		= new TH1F("hpf80_l1j1_Eta","",100,-5.,5.);
   h1["hpf80_l1j1_Phi"]    		= new TH1F("hpf80_l1j1_Phi", "",100,-3.2,3.2);

   h1["hpf80_l2p5j1_Pt"]   		= new TH1F("hpf80_l2p5j1_Pt","",100,0.,500.);
   h1["hpf80_l2p5j1_Eta"]  		= new TH1F("hpf80_l2p5j1_Eta","",100,-5.,5.);
   h1["hpf80_l2p5j1_Phi"]  		= new TH1F("hpf80_l2p5j1_Phi", "",100,-3.2,3.2);

   h1["hpf80_l3j1_Pt"]       	= new TH1F("hpf80_l3j1_Pt","",100,0.,500.);
   h1["hpf80_l3j1_Eta"]      	= new TH1F("hpf80_l3j1_Eta","",100,-5.,5.);
   h1["hpf80_l3j1_Phi"]      	= new TH1F("hpf80_l3j1_Phi", "",100,-3.2,3.2);
   
   //Histograms for Efficiency
   h1["hpf60_eff"]				= new TH1F("hpf60_eff","",100,0.,500.);
   h1["hpf60_effNumerator"]  	= new TH1F("hpf60_effNumerator","",100,0.,500.);
   h1["hpf60_effDenumerator"]	= new TH1F("hpf60_effDenumerator","",100,0.,500.);
   
   h1["hpf80_eff"]				= new TH1F("hpf80_eff","",100,0.,500.);
   h1["hpf80_effNumerator"]  	= new TH1F("hpf80_effNumerator","",100,0.,500.);
   h1["hpf80_effDenumerator"]	= new TH1F("hpf80_effDenumerator","",100,0.,500.);
   
   //Difference betwwen 2015D and C
   
   h1["hpf60_effL1"]				= new TH1F("hpf60_effL1","",100,0.,500.);
   h1["hpf60_effDenumeratorL1"]  	= new TH1F("hpf60_effDenumeratorL1","",100,0.,500.);
   
   h1["hpf60_effL2p5"]				= new TH1F("hpf60_effL2p5","",100,0.,500.);
   h1["hpf60_effDenumeratorL2p5"]  	= new TH1F("hpf60_effDenumeratorL2p5","",100,0.,500.);
   
   h1["hpf60_effL3"]				= new TH1F("hpf60_effL3","",100,0.,500.);
   h1["hpf60_effDenumeratorL3"]  	= new TH1F("hpf60_effDenumeratorL3","",100,0.,500.);
   
   
   // 160Pt treshold lvl
   h1["hpf60_eff160"]				= new TH1F("hpf60_eff160","",100,0.,500.);
   h1["hpf60_effNumerator160"]  	= new TH1F("hpf60_effNumerator160","",100,0.,500.);
   
   h1["hpf80_eff160"]				= new TH1F("hpf80_eff160","",100,0.,500.);
   h1["hpf80_effNumerator160"]  	= new TH1F("hpf80_effNumerator160","",100,0.,500.);
   
   // dEta efficiency
   h1["hpf80_l3Dij1_Pt"]			= new TH1F("hpf80_l3Dij1_Pt","",100,0.,500.);
   h1["hpf80_l3Dij1_Eta"]      		= new TH1F("hpf80_l3Dij1_Eta","",100,-5.,5.);
   h1["hpf80_l3Dij1_Phi"]      		= new TH1F("hpf80_l3Dij1_Phi", "",100,-3.2,3.2);
   
   h1["hpf80_l3Dij2_Pt"]			= new TH1F("hpf80_l3Dij2_Pt","",100,0.,500.);
   h1["hpf80_l3Dij2_Eta"]      		= new TH1F("hpf80_l3Dij2_Eta","",100,-5.,5.);
   h1["hpf80_l3Dij2_Phi"]      		= new TH1F("hpf80_l3Dij2_Phi", "",100,-3.2,3.2);
   
   h1["hpf80_l3Dij_dEta"]			= new TH1F("hpf80_l3Dij_dEta","",300,0.,3.);
   h1["hpf80_dEta_effNumerator"] 	= new TH1F("hpf80_dEta_effNumerator","",300,0.,3.);
   h1["hpf80_dEta_effDenumerator"] 	= new TH1F("hpf80_dEta_effDenumerator","",300,0.,3.);
   h1["hpf80_dEta_eff"] 			= new TH1F("hpf80_dEta_eff","",300,0.,3.);
   
   // Slices
   h1["h_slice_80_90"]				= new TH1F("h_slice_80_90","",100,0.,200.);
   h1["h_slice_90_100"]				= new TH1F("h_slice_90_100","",100,0.,200.);
   h1["h_slice_100_120"]			= new TH1F("h_slice_100_120","",100,0.,200.);
   h1["h_slice_170_180"]			= new TH1F("h_slice_170_180","",100,50.,250.);
   h1["h_slice_180_190"]			= new TH1F("h_slice_180_190","",100,50.,250.);
   h1["h_slice_280_290"]			= new TH1F("h_slice_280_290","",100,200.,350.);
   h1["h_slice_270_280"]			= new TH1F("h_slice_270_280","",100,200.,350.);
   
   h1["h_slice_80_90_dif"]				= new TH1F("h_slice_80_90_dif","",100,-70.,70.);
   h1["h_slice_90_100_dif"]				= new TH1F("h_slice_90_100_dif","",100,-70.,70.);
   h1["h_slice_100_120_dif"]			= new TH1F("h_slice_100_120_dif","",100,-70.,70.);
   h1["h_slice_170_180_dif"]			= new TH1F("h_slice_170_180_dif","",60,-70.,70.);
   h1["h_slice_180_190_dif"]			= new TH1F("h_slice_180_190_dif","",40,-70.,70.);
   h1["h_slice_280_290_dif"]			= new TH1F("h_slice_280_290_dif","",20,-70.,70.);
   h1["h_slice_270_280_dif"]			= new TH1F("h_slice_270_280_dif","",20,-70.,70.);
   
   TH2F *h2_offL3 = new TH2F("h2_offL3","",100,0.,500.,100,0.,500);
   
   if(!analysis.isMC()) analysis.processJsonFile("goodJson.txt");
	   	   
   // Analysis of events
   //std::cout<<" Entries: "<<analysis.luminosity()<<" "<<analysis.size()<<" "<<analysis.crossSection()<<" "<<analysis.size()/analysis.crossSection()<<std::endl; 
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
   		
      analysis.event(i);

      if (!analysis.isMC() && !analysis.selectJson()) continue;
      if (!analysis.isMC() && !analysis.triggerResult("HLT_PFJet60_v") && !analysis.triggerResult("HLT_PFJet80_v")) continue;
      
      // Selection of double b-tag sample:
      auto jets = analysis.collection<Jet>("Jets");
      if(jets -> size() < 1 ) continue;
      if(!analysis.isMC()) analysis.match<Jet,TriggerObject>("Jets",jetTriggerObjects);
      
      leadJetsNumber.clear();
      for( int iJet = 0; iJet < jets -> size(); ++iJet)
      {
      	Jet jet = jets -> at(iJet);
		if(jet.idLoose() ) leadJetsNumber.push_back(iJet);
		}
      if( leadJetsNumber.size() < 1) continue;
      
      Jet fLeadJet = jets -> at(leadJetsNumber[0]);
      if(abs(fLeadJet.eta()) > 2.5) continue;
         
      //Trigger matching:
      //Two Leading jets matchs online objects
      const Candidate * l1MPF60Cand = fLeadJet.matched("hltL1sL1SingleJet36");
      const Candidate * l2p5MPF60Cand = fLeadJet.matched("hltSingleCaloJet40");	  
      const Candidate * l3MPF60Cand = fLeadJet.matched("hltSinglePFJet60");
      
      //PFJet60 efficiency:
      if(analysis.triggerResult("HLT_PFJet60_v") && l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand && l1MPF60Cand->pt() > 100. && l2p5MPF60Cand->pt() > 100. && l3MPF60Cand->pt() > 100. ) h1["hpf60_effNumerator"] -> Fill(fLeadJet.pt());
	  if(analysis.triggerResult("HLT_PFJet60_v") ) h1["hpf60_effDenumerator"] -> Fill(fLeadJet.pt());
	  if(l1MPF60Cand && l2p5MPF60Cand && l3MPF60Cand && l1MPF60Cand->pt() > 160. && l2p5MPF60Cand->pt() > 160. && l3MPF60Cand->pt() > 160. ) h1["hpf60_effNumerator160"] -> Fill(fLeadJet.pt());
      
      const Candidate * l1MPF80Cand = fLeadJet.matched("hltL1sL1SingleJet52");	  
      const Candidate * l2p5MPF80Cand = fLeadJet.matched("hltSingleCaloJet50");
      const Candidate * l3MPF80Cand = fLeadJet.matched("hltSinglePFJet80");
	  
	  //PFJet80 efficiency:
	  if(analysis.triggerResult("HLT_PFJet80_v") &&  l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && l1MPF80Cand->pt() > 100. && l2p5MPF80Cand->pt() > 100. && l3MPF80Cand->pt() > 100. ) h1["hpf80_effNumerator"] -> Fill(fLeadJet.pt());
	  if(l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && l1MPF80Cand->pt() > 160. && l2p5MPF80Cand->pt() > 160. && l3MPF80Cand->pt() > 160. ) h1["hpf80_effNumerator160"] -> Fill(fLeadJet.pt());
	  if(analysis.triggerResult("HLT_PFJet80_v") ) h1["hpf80_effDenumerator"] -> Fill(fLeadJet.pt());
	  
	  double reweight = 1., reweight1 = 1.;
     if(analysis.isMC())
     {
     	//reweight1 = 1./(1.+exp( -8.76241e-02 * (fLeadJet.pt() - 1.05187e+02) )) * 1./ (1.+exp( -1.88955e-01 * (fLeadJet.pt() - 1.02451e+02) ));
     	reweight1 = ptRew2D(PtEff,fLeadJet.pt(),fLeadJet.eta());
     }
	  

	  if(!analysis.isMC() && l1MPF80Cand && l2p5MPF80Cand && l3MPF80Cand && abs(fLeadJet.eta()) < 2.2 && l1MPF80Cand->pt() > 100. && l2p5MPF80Cand->pt() > 100. && l3MPF80Cand->pt() > 100.)
	  {
	  	h1["h_Only1jet_jet1_Pt"] -> Fill(fLeadJet.pt());
      	h1["h_Only1jet_jet1_Eta"] -> Fill(fLeadJet.eta());
      	h1["h_Only1jet_jet1_Phi"] -> Fill(fLeadJet.phi());
      	if(fLeadJet.btag() > 0.941) h1["h_Only1jet_jet1_PtBTag"] -> Fill(fLeadJet.pt());
      	if(fLeadJet.btag() > 0.941 && fLeadJet.pt() > 100. && abs(fLeadJet.eta()) < 2.2) h1["h_Only1jet_jet1_PtAll"] -> Fill(fLeadJet.pt());
      }
      if(analysis.isMC() && abs(fLeadJet.eta()) < 2.2)
      {
	  	h1["h_Only1jet_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight1);
      	h1["h_Only1jet_jet1_Eta"] -> Fill(fLeadJet.eta(),reweight1);
      	h1["h_Only1jet_jet1_Phi"] -> Fill(fLeadJet.phi(),reweight1);
      	if(fLeadJet.btag() > 0.941) h1["h_Only1jet_jet1_PtBTag"] -> Fill(fLeadJet.pt(), reweight1 * bTagRew(btagEff,fLeadJet.pt(),fLeadJet.eta()));
      	if(fLeadJet.btag() > 0.941 && fLeadJet.pt() > 100. && abs(fLeadJet.eta()) < 2.2) h1["h_Only1jet_jet1_PtAll"] -> Fill(fLeadJet.pt(),reweight1 * bTagRew(btagEff,fLeadJet.pt(),fLeadJet.eta()));
      }
      
      if(analysis.triggerResult("HLT_PFJet80_v") && l3MPF80Cand && l1MPF80Cand && l2p5MPF80Cand)
     {
     	if(l3MPF80Cand->pt() > 80 && l3MPF80Cand->pt() < 90)  { h1["h_slice_80_90"] -> Fill(fLeadJet.pt()); 	h1["h_slice_80_90_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 90 && l3MPF80Cand->pt() < 100) { h1["h_slice_90_100"] -> Fill(fLeadJet.pt()); 	h1["h_slice_90_100_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 100 && l3MPF80Cand->pt() < 120){ h1["h_slice_100_120"] -> Fill(fLeadJet.pt());	h1["h_slice_100_120_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 170 && l3MPF80Cand->pt() < 180){ h1["h_slice_170_180"] -> Fill(fLeadJet.pt());	h1["h_slice_170_180_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 180 && l3MPF80Cand->pt() < 190){ h1["h_slice_180_190"] -> Fill(fLeadJet.pt());	h1["h_slice_180_190_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 280 && l3MPF80Cand->pt() < 290){ h1["h_slice_280_290"] -> Fill(fLeadJet.pt());	h1["h_slice_280_290_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     	if(l3MPF80Cand->pt() > 270 && l3MPF80Cand->pt() < 280){ h1["h_slice_270_280"] -> Fill(fLeadJet.pt());	h1["h_slice_270_280_dif"] -> Fill(fLeadJet.pt() - l3MPF80Cand->pt()); }
     }

	  // Eta Efifciency
	  if(leadJetsNumber.size() < 2) continue;
      Jet sLeadJet = jets -> at(leadJetsNumber[1]);
	  if(abs(fLeadJet.eta()) > 2.5 || abs(sLeadJet.eta()) > 2.5) continue;
      //if(analysis.isMC()) reweight = weight(fLeadJet.pt(),sLeadJet.pt());
	  
	  
	  const Candidate * l3DiMPF80sCand = sLeadJet.matched("hltSinglePFJet80");	  
	  const Candidate * l1MPF80sCand = sLeadJet.matched("hltL1sL1SingleJet52");
      const Candidate * l2p5MPF80sCand = sLeadJet.matched("hltSingleCaloJet50");

	  if(l3MPF80Cand && l3DiMPF80sCand) h1["hpf80_l3Dij_dEta"] -> Fill (std::abs(l3MPF80Cand->eta() - l3DiMPF80sCand->eta()) );
	  
	  if(analysis.triggerResult("HLT_PFJet80_v") && l3MPF80Cand && l3DiMPF80sCand && l3DiMPF80sCand->pt() > 100 && l3MPF80Cand->pt() > 100 && std::abs(l3MPF80Cand->eta() - l3DiMPF80sCand->eta()) < 1.6 ) h1["hpf80_dEta_effNumerator"] -> Fill( std::abs(fLeadJet.eta() - sLeadJet.eta() ) );
	  if(analysis.triggerResult("HLT_PFJet80_v") && l3MPF80Cand && l3DiMPF80sCand && l3DiMPF80sCand->pt() > 100 && l3MPF80Cand->pt() > 100) h1["hpf80_dEta_effDenumerator"] -> Fill( std::abs(fLeadJet.eta() - sLeadJet.eta() ) ) ;
		
     if(l3MPF80Cand) h2_offL3->Fill(fLeadJet.pt(),l3MPF80Cand->pt());
     
     if(!analysis.isMC())
     {
     	if(!sLeadJet.matched("hltL1sL1SingleJet52") || !sLeadJet.matched("hltSinglePFJet80") || !sLeadJet.matched("hltSingleCaloJet50")) continue;
     	if(!fLeadJet.matched("hltL1sL1SingleJet52") || !fLeadJet.matched("hltSinglePFJet80") || !fLeadJet.matched("hltSingleCaloJet50")) continue;
     	if( fLeadJet.matched("hltL1sL1SingleJet52") == sLeadJet.matched("hltL1sL1SingleJet52") || fLeadJet.matched("hltSinglePFJet80") == sLeadJet.matched("hltSinglePFJet80") || fLeadJet.matched("hltSingleCaloJet50") == sLeadJet.matched("hltSingleCaloJet50") ) continue;
     	if(l1MPF80Cand->pt() < 100. || l2p5MPF80Cand->pt() < 100. || l3MPF80Cand->pt() < 100. || l1MPF80sCand->pt() < 100. || l2p5MPF80sCand->pt() < 100. || l3DiMPF80sCand->pt() < 100.) continue;
     }
     
     if(analysis.isMC()) reweight = ptRew2D(PtEff,fLeadJet.pt(),fLeadJet.eta()) * ptRew2D(PtEff,sLeadJet.pt(),sLeadJet.eta());
     
     h1["h_TwoJets_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     h1["h_TwoJets_jet1_Eta"] -> Fill(fLeadJet.eta(),reweight);
     h1["h_TwoJets_jet1_Phi"] -> Fill(fLeadJet.phi(),reweight);
     
     h1["h_TwoJets_jet2_Pt"] -> Fill(sLeadJet.pt(),reweight);
     h1["h_TwoJets_jet2_Eta"] -> Fill(sLeadJet.eta(),reweight);
     h1["h_TwoJets_jet2_Phi"] -> Fill(sLeadJet.phi(),reweight);
     
     if(fLeadJet.pt() > 100. && sLeadJet.pt() > 100) h1["h_PtStudy_ptCut_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     if(fLeadJet.eta() < 2.2 && sLeadJet.eta() < 2.2) h1["h_PtStudy_etaCut_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     if(fLeadJet.deltaR(sLeadJet) > 1 ) h1["h_PtStudy_dRCut_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     if( fLeadJet.pt() > 100. && sLeadJet.pt() > 100 && abs(fLeadJet.eta()) < 2.2 && abs(sLeadJet.eta()) < 2.2 && fLeadJet.deltaR(sLeadJet) > 1 && abs(fLeadJet.eta() - sLeadJet.eta()) < 1.6 ) h1["h_PtStudy_AllCut_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
     
     if(fLeadJet.btag() >= 0.9 && sLeadJet.btag() >= 0.9) 
     {
     	h1["h_TwoBtagJets_jet1_Pt"] -> Fill(fLeadJet.pt(),reweight);
    	h1["h_TwoBtagJets_jet1_Eta"] -> Fill(fLeadJet.eta(),reweight);
    	h1["h_TwoBtagJets_jet1_Phi"] -> Fill(fLeadJet.phi(),reweight);
     
     	h1["h_TwoBtagJets_jet2_Pt"] -> Fill(sLeadJet.pt(),reweight);
     	h1["h_TwoBtagJets_jet2_Eta"] -> Fill(sLeadJet.eta(),reweight);
     	h1["h_TwoBtagJets_jet2_Phi"] -> Fill(sLeadJet.phi(),reweight);
     }
        
           ++nCand;

   }
   h1["hpf60_eff"] -> Divide(h1["hpf60_effNumerator"],h1["hpf60_effDenumerator"],1.,1.,"b");
   h1["hpf80_eff"] -> Divide(h1["hpf80_effNumerator"],h1["hpf80_effDenumerator"],1.,1.,"b");
   h1["hpf60_effL1"] -> Divide(h1["hpf60_effNumerator"],h1["hpf60_effDenumeratorL1"],1.,1.,"b");
   h1["hpf60_effL2p5"] -> Divide(h1["hpf60_effNumerator"],h1["hpf60_effDenumeratorL2p5"],1.,1.,"b");
   h1["hpf60_effL3"] -> Divide(h1["hpf60_effNumerator"],h1["hpf60_effDenumeratorL3"],1.,1.,"b");
   h1["hpf60_eff160"] -> Divide(h1["hpf60_effNumerator160"],h1["hpf60_effDenumerator"],1.,1.,"b");
   h1["hpf80_eff160"] -> Divide(h1["hpf80_effNumerator160"],h1["hpf80_effDenumerator"],1.,1.,"b");
   h1["hpf80_dEta_eff"] -> Divide(h1["hpf80_dEta_effNumerator"], h1["hpf80_dEta_effDenumerator"], 1., 1., "b");
   
   std::cout<<"Number of Candidates: "<<nCand<<std::endl;
   std::string fileName = "output.root";
   fileName = analysis.tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
   fileName = "outputFiles/ResultsMC" + boost::filesystem::basename(fileName) + ".root";
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
   h2_offL3->Write();
  
   outFile -> Close();
   
   
   std::cout<<"File: "<<fileName<<" was created"<<std::endl;
   
   return 0;
//
}

double weight(const float &pt1, const float &pt2)
{
	double weightPt2 = 1.;
	double weightPt1 = 1.;
	weightPt1 = 1./(1.+exp( -8.76241e-02 * (pt1 - 1.05187e+02) )) * 1./ (1.+exp( -1.88955e-01 * (pt1 - 1.02451e+02) ));
	weightPt2 = 1./(1.+exp( -8.76241e-02 * (pt2 - 1.05187e+02) )) * 1./ (1.+exp( -1.88955e-01 * (pt2 - 1.02451e+02) ));
	
	return weightPt1 * weightPt2;
}
double ptRew2D(TH2F *PtH, const float & pt,const float &eta)
{
	if(pt <500) return PtH -> Interpolate(pt,eta);
	else return 1;
}
double bTagRew(TH2F *histo, const float &pt, const float &eta)
{
	if(pt <1000) return histo -> Interpolate(pt,eta);
	else return 1;
}
