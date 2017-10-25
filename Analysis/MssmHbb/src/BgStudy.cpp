/*
 * BgStudy.cpp
 *
 *  Created on: Apr 23, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/BgStudy.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

BgStudy::BgStudy(const std::string & inputFilelist,const double & dataLumi, const bool & lowM, const bool & test) :
		MssmHbbSignal(inputFilelist,dataLumi,lowM,test)
{
	nJets_ = 3;
	if(lowM){
		btag1_ = -100;
		btag2_ = -100;
        btag3_ = -100;
        btagOP3_ = 1;
	}
	else {
        btag3_ = 0.46;
	}
	baseOutputName_ = "BgStudy";
}

BgStudy::~BgStudy() {
}

//Member methods

const bool BgStudy::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){

	//To avoid code duplication frind class were used
	if(!selectionDoubleB::leadingJetSelection(offlineJets)) return false;

//	if(TEST) std::cout<<"I'm in BgStudy::leadingJetSelection"<<std::endl;

	//Selection of good Leading Jets:
	//Only jets that pass Loose identification will be considered
	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);
	Jet jet3 = offlineJets->at(2);

	//Trigger Selection
//	if(!this->triggerResult(triggerLogicName_)) return false;

	//Online selection:
//	if(!this->OnlineSelection(jet1,jet2)) return false;

	//Pt requirements
	if (jet3.pt() < pt3_) return false;

	//Eta requirements
	if (std::abs(jet3.eta()) > eta3_) return false;

	//BTag requirements
	if (jet3.btag() < btag3_) return false;

	//deltaR requirements
	if (jet2.deltaR(jet3) <= dR_) return false;
	if (jet1.deltaR(jet3) <= dR_) return false;

	if(jet1.flavour() != 5) return false;
	if(jet2.flavour() != 5) return false;
	if(jet3.flavour() != 5) return false;

	return true;
}

void BgStudy::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){

//	if(TEST) std::cout<<"I'm in BgStudy::fillHistograms"<<std::endl;

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);
	Jet jet3 = offlineJets->at(2);
	TLorentzVector obj12;
	obj12 = jet1.p4() + jet2.p4();

	MssmHbbSignal::fillHistograms(offlineJets,weight);

	(histo_.getHisto())["data_Mbb_20"]->Fill(obj12.M(),weight);
	(histo_.getHisto())["data_Mbb_30"]->Fill(obj12.M(),weight);
	(histo_.getHisto())["data_Mbb_50"]->Fill(obj12.M(),weight);
}

const double BgStudy::assignWeight(){
double weight = 1;
	if(isMC()) {
		weight = weight_["Lumi"] * weight_["PU_central"] * weight_["SFb_central"] * weight_["SFl_central"];
	}
	return weight;
}



