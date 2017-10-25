/*
 * doubleB.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/selectionDoubleB.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

selectionDoubleB::selectionDoubleB(const std::string & inputFilelist, const double & dataLumi, const bool & lowM, const bool & test) :
JetAnalysisBase(inputFilelist,dataLumi,lowM,test) {
	nJets_ = 2;
        if(lowM){
        	triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
        	triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp014DoubleWithMatching","hltDoublePFJetsC100","hltDoublePFJetsC100MaxDeta1p6"};
        	pt1_ = 100.; pt2_ = 100.; pt3_ = 40;
        	eta1_ = 2.2; eta2_ = 2.2; eta3_ = 2.2;
        	btag1_ = 0.8; btag2_ = 0.8; btag3_ = 0.8;
        	btagOP1_ = 1; btagOP2_ = 1; btagOP3_ = 1;
        	dR_ = 1; dEta_ = 1.55; mHat_ = 0.7;
        }
        else {
        	triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v";
        	triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp026DoubleWithMatching","hltDoublePFJetsC160"};
        	pt1_ = 160.; pt2_ = 160.; pt3_ = 40;
        	eta1_ = 2.2; eta2_ = 2.2; eta3_ = 2.2;
        	btag1_ = 0.8; btag2_ = 0.8; btag3_ = 0.8;
        	btagOP1_ = 1; btagOP2_ = 1; btagOP3_ = 1;
        	dR_ = 1; dEta_ = 100.; mHat_ = 0.7;
        }
	baseOutputName_ = "DoubleB";
}

selectionDoubleB::~selectionDoubleB() {
	// TODO Auto-generated destructor stub
}

const bool selectionDoubleB::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){
//	if(TEST) std::cout<<"I'm in selectionDoubleB::leadingJetSelection"<<std::endl;
	//Selection of good Leading Jets:
	//Only jets that pass Loose identification will be considered

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);

	//Pt requirements
	if(!cuts_.check("pt1",jet1.pt() >= pt1_)) return false;
	if(!cuts_.check("pt2",jet2.pt() >= pt2_)) return false;

	//Eta requirements
	if(!cuts_.check("eta1",std::abs(jet1.eta()) <= eta1_)) return false;
	if(!cuts_.check("eta2",std::abs(jet2.eta()) <= eta2_)) return false;

	//dEta in case of LowM selection
	if(lowM_){
		if(!cuts_.check("deta12",std::abs(jet1.eta() - jet2.eta()) <= dEta_)) return false;
	}

	//deltaR requirements
	if(!cuts_.check("dR12",jet1.deltaR(jet2) > dR_)) return false;

	//BTag requirements
	if(!cuts_.check("btag1",jet1.btag() >= btag1_)) return false;
	if(!cuts_.check("btag2",jet2.btag() >= btag2_)) return false;



	//Trigger Selection
	if(!cuts_.check("TriggerBit",this->triggerResult(triggerLogicName_))) return false;

	//Online selection:
	if(!cuts_.check("TriggerMatching",this->OnlineSelection(jet1,jet2))) return false;



	return true;
}

void selectionDoubleB::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){
	if(TEST) std::cout<<"I'm in selectionDoubleB::fillHistograms"<<std::endl;

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);

	int njets = 0;
	for (int iJet = 0; iJet < offlineJets -> size(); ++iJet){
		if(offlineJets->at(iJet).pt() > 20) ++njets;
	}
	(histo_.getHisto())["NumberOfJets"]->Fill(njets);
//	(histo_.getHisto())["NumberOfJets"]->Fill(offlineJets->size());

	(histo_.getHisto())["jet_pt1"]->Fill(jet1.pt(),weight);
	(histo_.getHisto())["jet_pt2"]->Fill(jet2.pt(),weight);

	(histo_.getHisto())["jerResolution1"]->Fill(jet1.JerResolution(),weight);
	(histo_.getHisto())["jerResolution2"]->Fill(jet2.JerResolution(),weight);

	(histo_.getHisto())["jet12_assym"]->Fill((jet1.pt()-jet2.pt())/(jet1.pt()+jet2.pt()),weight);

	(histo_.getHisto())["jet_eta1"]->Fill(jet1.eta(),weight);
	(histo_.getHisto())["jet_eta2"]->Fill(jet2.eta(),weight);

	(histo_.getHisto())["jet_phi1"]->Fill(jet1.phi(),weight);
	(histo_.getHisto())["jet_phi2"]->Fill(jet2.phi(),weight);

	(histo_.getHisto())["jet_btag_csv1"]->Fill(jet1.btag(),weight);
	(histo_.getHisto())["jet_btag_csv2"]->Fill(jet2.btag(),weight);

	(histo_.getHisto())["jet_btag_cmva1"]->Fill(jet1.btag("btag_csvmva"),weight);
	(histo_.getHisto())["jet_btag_cmva2"]->Fill(jet2.btag("btag_csvmva"),weight);

	(histo_.getHisto())["jet_deta12"]->Fill(jet1.eta() - jet2.eta(),weight);

	(histo_.getHisto())["jet_dR12"]->Fill(jet1.deltaR(jet2),weight);

	if(isMC() && genJets_->size()>1){
		Jet genjet1 = genJets_->at(0);
		Jet genjet2 = genJets_->at(1);
		(histo_.getHisto())["Genjet_pt1"]->Fill(genjet1.pt(),weight);
		(histo_.getHisto())["Genjet_pt2"]->Fill(genjet2.pt(),weight);

		(histo_.getHisto())["Genjet_eta1"]->Fill(genjet1.eta(),weight);
		(histo_.getHisto())["Genjet_eta2"]->Fill(genjet2.eta(),weight);

		(histo_.getHisto())["Genjet_phi1"]->Fill(genjet1.phi(),weight);
        (histo_.getHisto())["Genjet_phi2"]->Fill(genjet2.phi(),weight);

		(histo_.getHisto())["jet_flavor1"]->Fill(jet1.flavour(),weight);
        (histo_.getHisto())["jet_flavor2"]->Fill(jet2.flavour(),weight);

        (histo_.getHisto2D())["Genjet_vs_off_pt1"]->Fill(genjet1.pt(),jet1.pt(),weight);
	}

	TLorentzVector obj12;
	obj12 = jet1.p4() + jet2.p4();
	(histo_.getHisto())["diJet_pt"]->Fill(obj12.Pt(),weight);
	(histo_.getHisto())["diJet_eta"]->Fill(obj12.Eta(),weight);
	(histo_.getHisto())["diJet_phi"]->Fill(obj12.Phi(),weight);
	(histo_.getHisto())["diJet_m"]->Fill(obj12.M(),weight);
	(histo_.getHisto())["diJet_m_5GeV"]->Fill(obj12.M(),weight);
	(histo_.getHisto())["diJet_m_20GeV"]->Fill(obj12.M(),weight);

	//Flavour composition
	if(isMC()){
		if(jet1.flavour() == 5 && jet2.flavour() == 5) (histo_.getHisto())["diJet_m_bb"]->Fill(obj12.M(),weight);
		else if ((jet1.flavour() == 5 && jet2.flavour() == 4) || (jet2.flavour() == 5 && jet1.flavour() == 4)) (histo_.getHisto())["diJet_m_bc"]->Fill(obj12.M(),weight);
		else if ((jet1.flavour() == 5 && jet2.flavour() == 0) || (jet2.flavour() == 5 && jet1.flavour() == 0)) (histo_.getHisto())["diJet_m_bl"]->Fill(obj12.M(),weight);
		else if (jet1.flavour() == 0 && jet2.flavour() == 0) (histo_.getHisto())["diJet_m_ll"]->Fill(obj12.M(),weight);
		else if (jet1.flavour() == 4 && jet2.flavour() == 4) (histo_.getHisto())["diJet_m_cc"]->Fill(obj12.M(),weight);
		else if ((jet1.flavour() == 4 && jet2.flavour() == 0) || (jet2.flavour() == 4 && jet1.flavour() == 0)) (histo_.getHisto())["diJet_m_cl"]->Fill(obj12.M(),weight);
		else std::logic_error("Undefined flavours at: selectionDoubleB::fillHistograms");
	}

}

const double selectionDoubleB::assignWeight(){
	double weight = 1;
	if(isMC()) {
		weight = weight_["PtEff_central"] * weight_["PU_central"] * weight_["SFb_central"] * weight_["SFl_central"];
	}
	return weight;
}
