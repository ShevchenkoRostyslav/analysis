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

double ptRewTnP(const float &, const bool &);
double ptRewComb(const float &, const bool &);
double ptRewRef(const float &, const bool &);

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;


// =============================================================================================
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   // Input files list MonteCarloStudies
   std::string inputList = "/nfs/dust/cms/user/shevchen/samples/miniaod/JetHT/JetHTFileListRun2015D05Oct.txt";
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
   analysis.setBranches();
   
   int nCand = 0;
   
   //std::map < std::string, TH1F *> h1;

   if(!analysis.isMC()) analysis.processJsonFile("goodJson.txt");
   std::cout<<"Total Entries: "<<analysis.size()<<std::endl;

	TRandom1 rand;
	TLorentzVector obj;
	int rand_val = 0;
	bool PtVetoB = false, goodLeadingJets = false;
	int counter;
   // Analysis of events
	double lumiData2015DPrRec4 = 1001.179266;
	std::string fileName = "output.root";
	fileName = analysis.tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
	if(analysis.isMC()) fileName = "/nfs/dust/cms/user/shevchen/output/ResultsTreesMC" + boost::filesystem::basename(fileName) + ".root";
	else fileName = "/nfs/dust/cms/user/shevchen/output/ResultsTrees" + boost::filesystem::basename(fileName) + ".root";
	TFile * outFile = new TFile(fileName.c_str(),"RECREATE");

	int num = 0;
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
   		
      analysis.event(i);
      if (!analysis.isMC() && !analysis.selectJson()) continue;
      if( !analysis.isMC() && !analysis.triggerResult("HLT_PFJet60_v") && !analysis.triggerResult("HLT_PFJet80_v")) continue;

      // Selection of double jet sample:
      auto jets = analysis.collection<Jet>("Jets");
      analysis.match<Jet,TriggerObject>("Jets",analysis.triggerObjectNames());

      counter = 0;
      goodLeadingJets = true;
      PtVetoB = true;

      Jet LeadJet[20];
      analysis.zeroingBranches();
      analysis.Njets80 = 0;
      if(jets->size() < 2) continue;


      for( int iJet = 0; iJet < jets -> size(); ++iJet)
      {
      	Jet jet = jets -> at(iJet);
      	if(jet.pt() >80){
      		analysis.Njets80 ++;
      	}
		if(!jet.idLoose()) continue;
		counter++;
		LeadJet[counter - 1] = jet;
		if(counter == 1 || counter == 2){

			if(jet.pt() < 90) goodLeadingJets = false;
			if(abs(jet.eta()) > 2.5) goodLeadingJets = false;
			if(!goodLeadingJets) break;
		}
		else{
			if(jet.pt() > 40) {
				analysis.ptVeto = jet.pt();
				LeadJet[counter - 1] = jet;
				//PtVetoB = false;
			}
			if(!PtVetoB) break;
		}
	  }
      if(counter < 2) continue;
      if(!PtVetoB) continue;
      analysis.dPhiFS = acos(cos(LeadJet[0].phi() - LeadJet[1].phi()));
      //if(analysis.dPhiFS < 2.7 ) goodLeadingJets = false;
      //if((LeadJet[0].pt() - LeadJet[1].pt()) / (LeadJet[0].pt() + LeadJet[1].pt()) > 0.2 ) goodLeadingJets = false;
      if (!goodLeadingJets) continue;

      //Trigger matching:
      analysis.NL1Object = analysis.collection<TriggerObject>("hltL1sL1SingleJet52")->size();
      analysis.NL2Object = analysis.collection<TriggerObject>("hltSingleCaloJet50")->size();
      analysis.NL3Object = analysis.collection<TriggerObject>("hltSinglePFJet80")->size();

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

      if(counter > 3 ) num = 4;
      if(counter == 3) num = 3;
      if(counter == 2) num = 2;

      for(int i = 0; i<num;i++){
    	  analysis.LeadPt[i] = LeadJet[i].pt();
    	  analysis.LeadEta[i] = LeadJet[i].eta();
    	  analysis.LeadPhi[i] = LeadJet[i].phi();
    	  analysis.LeadBTag[i] = LeadJet[i].btag();
    	  if(i < 2){
    	      analysis.LeadMatch60[i] = analysis.matchToPF60(LeadJet[i]);
    	      analysis.LeadMatch80[i] = analysis.matchToPF80(LeadJet[i]);
    	      analysis.LeadMatch100[i] = analysis.matchToPF100(LeadJet[i]);

    	      //Efficiency calculation
    		  if(analysis.isMC()) {
    		  	if(analysis.LeadMatch100[i] ) analysis.mcModel_Num[i] = LeadJet[i].pt();
    		  	analysis.mcModel_Denum[i] = LeadJet[i].pt();
    		  }

    	      if(analysis.LeadMatch100[i]) analysis.hpf100_Num[i] = LeadJet[i].pt();
    	      if(analysis.triggerResult("HLT_PFJet80_v")) analysis.hpf100_Denum[i] = LeadJet[i].pt();

    	      //De-randomization:
    	      if(i == 0) {
    	    	  tagJet = LeadJet[i+1];
    	    	  probeJet = LeadJet[i];
    	      } else if(i == 1){
    	    	  tagJet = LeadJet[i-1];
    	    	  probeJet = LeadJet[i];
    	      }

    	      //TnP
    	      if( analysis.matchToPF100(tagJet))
    	      {
    	      	if(analysis.matchToPF100(probeJet)) analysis.TnP100_probe[i] = probeJet.pt(); //probe
    	      	analysis.TnP100_tag[i] = probeJet.pt(); // tag
    	      }
    	  }
      }

      analysis.Njets = jets->size();

      analysis.dEtaFS = abs(analysis.LeadEta[0] - analysis.LeadEta[1]);

	  analysis.weightPtComb = ptRewComb(analysis.LeadPt[0],analysis.isMC()) *ptRewComb(analysis.LeadPt[1],analysis.isMC());
	  analysis.weightPtRef = ptRewRef(analysis.LeadPt[0],analysis.isMC()) *ptRewRef(analysis.LeadPt[1],analysis.isMC());
	  analysis.weightPtTnP = ptRewTnP(analysis.LeadPt[0],analysis.isMC()) *ptRewTnP(analysis.LeadPt[1],analysis.isMC());
      if( analysis.isMC()) {
    	  analysis.weightLumi = lumiData2015DPrRec4/analysis.luminosity();
      }

      obj = LeadJet[0].p4() + LeadJet[1].p4();
      analysis.M12 = obj.M();

      analysis.fillTree();
           ++nCand;

   }

   std::cout<<"Number of Candidates: "<<nCand<<std::endl;

   analysis.writeTree();
   outFile -> Close();
   
   
   std::cout<<"File: "<<fileName<<" was created"<<std::endl;
   
   return 0;
//
}

double ptRewTnP(const float &pt, const bool &mc)
{
	return 1./(1.+exp( -8.05781e-02 * (pt - 9.33941e+01) )) * 1./ (1.+exp( -1.46596e-01 * (pt - 1.03105e+02) ));
}

double ptRewComb(const float &pt, const bool &mc)
{
	return 1./(1.+exp( -1.06208e-01 * (pt - 1.06026e+02) )) * 1./ (1.+exp( -2.87063e-01 * (pt - 1.00032e+02) ));
}

double ptRewRef(const float &pt, const bool &mc)
{
	if(mc) return 1./(1.+exp( -8.84154e-02 * (pt - 1.00909e+02) )) * 1./ (1.+exp( -1.82032e-01 * (pt - 9.98953e+01) ));
	else return 1./(1.+exp( -9.67376e-02 * (pt - 1.02491e+02) )) * 1./ (1.+exp( -1.85832e-01 * (pt - 1.02067e+02) ));
}
