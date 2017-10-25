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

//Delcare functions
double factorizationPtWeight1D(const double &);
double dEtaWeight(const double &);
double BTagWeight(TH2F*region1,TH2F* region2,TH2F* region3, const double &pt, const double &eta);
double twoDPtWeight(TH2F *, const double &, const double &);

using namespace std;
using namespace analysis::selection;
using namespace analysis::tools;


// =============================================================================================
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms

   //Online Pt trigger efficiency:
   TFile * filePtEff = new TFile("~/cms/cmssw-analysis/CMSSW_7_5_2/src/Analysis/MssmHbb/bin/input_corrections/TwoDPtHighMassEfficiency.root");
   TH2F *hPtEff = (TH2F*) filePtEff ->Get("TwoDEff_Num"); // 2D

   // Input files list MssmHbb
   std::string inputList = "/nfs/dust/cms/user/shevchen/samples/miniaod/76X/JetHT/Run2015C_25ns-16Dec2015-v1.txt";
   TriggerEfficiency analysis(inputList); //Default second argument is "MssmHbb/Events/EventInfo"
   
   // Process selected JSON file
   if(!analysis.isMC()) analysis.processJsonFile("Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt");

   // Add std::vector<std::string> of the Trigger Objects that you would like to apply.
   // Also Trigger Results name will be stored, according to the trigger objects names
   // By default - LowM MssmHbb trigger objects is used
   std::vector<std::string> trigObj = {"hltL1sL1ZeroBias","hltSingleCaloJet10","hltSinglePFJet40",
		   	   	   	   	   	   	   	   "hltL1sL1SingleJet36","hltSingleCaloJet40","hltSinglePFJet60",
									   "hltL1sL1SingleJet52","hltSingleCaloJet50","hltSinglePFJet80"};
   if(!analysis.isMC()) analysis.addTriggerObjects(trigObj);

   //Setup output file name
   //name can me specified explicitly with method: createOutputFile(fileName);
   std::string fileName = "/nfs/dust/cms/user/shevchen/output/EfficiencyStudy_76X";
   analysis.SetupStandardOutputFile(fileName);

   //Setup Branches. Shoudl be used AFTER declaration of trigger objects
   analysis.setBranches();

   int counter = 0;
   bool goodLeadingJets = true;
   // Analysis of events

   std::cout<<"Number of Entries: "<<analysis.size()<<std::endl;
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {

      analysis.event(i);

      //Select only good JSon runs
      if (!analysis.isMC() && !analysis.selectJson()) continue;

      //Trigger Selection
      if( !analysis.isMC() && !analysis.triggerResult("HLT_PFJet60_v") && !analysis.triggerResult("HLT_PFJet80_v")) continue;

      //Used to clean values of the variables for the new event
      analysis.cleanVariables();

      //Define Jet Collection
      auto offlineJets = analysis.collection<Jet>("Jets");

//      if (offlineJets -> size() < 2) continue;

      //Match offline Jets to online Objects
      if (!analysis.isMC()) analysis.match<Jet,TriggerObject>("Jets",analysis.getTriggerObjectNames());

      //Set Total Number of Jets
      analysis.setNjets(offlineJets->size());

      counter = 0;
      goodLeadingJets = false;
      Jet LeadJet[20];
      // Selection of double jet sample:
      for( int iJet = 0; iJet < offlineJets -> size(); ++iJet)
      {
         	Jet jet = offlineJets -> at(iJet);
          	//Calculate HT
           	analysis.addHt(jet.pt());

            //Only jets that pass Loose identification will be considered
      		if(!jet.idLoose()) continue;
      		counter++;

      		//Only first 5 jets is taken into account
      		//To reduce size of the output file
      		if(counter > 5) break;
      		LeadJet[counter - 1] = jet;

      		analysis.setJetCounter(counter-1);
      		analysis.setJetVariables(jet);
      		analysis.matchToPF60(jet);
      		analysis.matchToPF80(jet);
      		analysis.matchToPF100(jet);
      		analysis.matchToPF160(jet);
      		analysis.matchToPF100L3(jet);

      		//Selection cuts for first leading jets
      		if(counter == 1)
      		{
      			if(jet.pt() < 90) break;
      			if(abs(jet.eta()) > 2.2) break;
      			goodLeadingJets = true;
      		}
      }
      if(!goodLeadingJets) continue;

      // Method that calculates dEta, dPhi between frist t2o leading jets
      // Also all parameters of di-jet object is calculated
      analysis.calculateJetVariables();

      //Flag to switch from single to double Jet topology
      bool double_topology = false;
      if(LeadJet[1].pt() >= 90 && abs(LeadJet[1].eta()) <= 2.2) double_topology = true;
      analysis.setTopology(double_topology);
      if (double_topology) analysis.matchToPF100dEta1p6(LeadJet[0],LeadJet[1]);

      //Set trigger objects variables
      analysis.setTriggerObjectVars();

      if(analysis.isMC()){
          	analysis.setLumiWeight(2182.680439,analysis.luminosity());
           	analysis.setPtFactorizationWeight(factorizationPtWeight1D(LeadJet[0].pt()) * factorizationPtWeight1D(LeadJet[1].pt()));
           	if(double_topology) {
           		analysis.setPt2DWeight(twoDPtWeight(hPtEff,LeadJet[0].pt(),LeadJet[1].pt()));
           		analysis.setdEtaWeight(abs(LeadJet[0].eta() - LeadJet[1].eta()));
           	}
      }

      analysis.fillTree();


   }
   std::cout<<"Number of Candidates: "<<analysis.getNumberOfCandidates()<<std::endl;
   analysis.writeTree();

}

double factorizationPtWeight1D(const double &pt)
{
	return 1./(1.+exp( -1.60938e-01 * (pt - 1.64206e+02) )) * 1./ (1.+exp( -1.53883e-01 * (pt + 3.50936e+01) ));
}
double dEtaWeight(const double & dEta)
{
	return 1./(1.+exp( 1.41132e+02 * (dEta - 1.60010e+00) )) * 1./ (1.+exp( -2.06352e-02 * (dEta + 2.27738e+02) ));
}
double BTagWeight(TH2F*region1,TH2F* region2,TH2F* region3, const double &pt, const double &eta){
	if(pt <1000) {
		if(eta < 0.9) return region1 -> Interpolate(pt);
		else if (eta > 0.9 &&eta < 1.4) return region2 -> Interpolate(pt);
		else return region3 -> Interpolate(pt);

	}
	else {
		if(eta < 0.9) return region1 -> Interpolate(995.);
		else if (eta > 0.9 &&eta < 1.4) return region2 -> Interpolate(995.);
		else return region3 -> Interpolate(995.);
	}
}

double twoDPtWeight(TH2F *histo, const double &pt1, const double &pt2){
	if (pt1 >= 500) {
		if (pt2 >= 500){
			return histo->Interpolate(489,489);
		}
		else return histo->Interpolate(pt2,489);
	}
	else return histo -> Interpolate(pt2,pt1);
}
