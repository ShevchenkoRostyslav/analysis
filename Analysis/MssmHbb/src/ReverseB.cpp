/*
 * ReverseB.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/ReverseB.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

ReverseB::ReverseB(const std::string & inputFilelist, const bool & lowM, const bool & test) :
					selectionDoubleB(inputFilelist,lowM,test)
{
	nJets_ = 3;
	if(lowM){
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp014DoubleWithMatching","hltDoublePFJetsC100","hltDoublePFJetsC100MaxDeta1p6"};
        pt3_ = 40.;
        btag3_ = 0.8;
	}
	else {
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp026DoubleWithMatching","hltDoublePFJetsC160"};
        pt3_ = 40.;
        btag3_ = 0.46;
	}
}

ReverseB::~ReverseB() {
	// TODO Auto-generated destructor stub
}

const bool ReverseB::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){

	//To avoid code duplication frind class were used
	if(!selectionDoubleB::leadingJetSelection(offlineJets)) return false;

	if(TEST) std::cout<<"I'm in ReverseB::leadingJetSelection"<<std::endl;

	//Selection of good Leading Jets:
	//Only jets that pass Loose identification will be considered
	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);
	Jet jet3 = offlineJets->at(2);

	//Trigger Selection
	if(!isMC() && !this->triggerResult(triggerLogicName_)) return false;

	//Online selection:
	if(!isMC() && !this->OnlineSelection(jet1,jet2)) return false;

	//Pt requirements
	if (jet3.pt() < pt3_) return false;

	//Eta requirements
	if (std::abs(jet3.eta()) > eta3_) return false;

	//BTag requirements
	if (jet3.btag() > btag3_) return false;

	//deltaR requirements
	if (jet2.deltaR(jet3) <= dR_) return false;
	if (jet1.deltaR(jet3) <= dR_) return false;

	return true;
}

void ReverseB::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){

	if(TEST) std::cout<<"I'm in ReverseB::fillHistograms"<<std::endl;
	//TODO: re-implement with selectionDoubleB friendness

	selectionDoubleB::fillHistograms(offlineJets,weight);

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);
	Jet jet3 = offlineJets->at(2);
	TLorentzVector obj12;
	obj12 = jet1.p4() + jet2.p4();

	(histo_.getHisto())["jet_pt3"]->Fill(jet3.pt(),weight);
	(histo_.getHisto())["jet_eta3"]->Fill(jet3.eta(),weight);
	(histo_.getHisto())["jet_phi3"]->Fill(jet3.phi(),weight);
	(histo_.getHisto())["jet_btag_csv3"]->Fill(jet3.btag(),weight);
	(histo_.getHisto())["jet_btag_cmva3"]->Fill(jet3.btag("btag_csvmva"),weight);

	if(TEST) std::cout<<"I'm out of ReverseB::fillHistograms"<<std::endl;
}

const double ReverseB::assignWeight(){
	double weight = 1;
	if(isMC()) {
		weight = weight_["dEta"] * weight_["Lumi"] * weight_["2DPt"] * weight_["BTag"] * weight_["PU_central"] * weight_["SFb_central"];
	}
	return weight;
}
