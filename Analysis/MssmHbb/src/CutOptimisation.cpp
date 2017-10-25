/*
 * doubleB.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/CutOptimisation.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

CutOptimisation::CutOptimisation(const std::string & inputFilelist, const double & dataLumi, const bool & lowM,const int & njets, const bool & test) :
		selectionDoubleB(inputFilelist,dataLumi,lowM,test) {
	baseOutputName_ = "CutOptimisation";
	nJets_ = 3;
	JESshift_ = 0;
	JERshift_ = 0;
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
}

CutOptimisation::~CutOptimisation() {
	// TODO Auto-generated destructor stub
}


const bool CutOptimisation::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){

//	if(TEST) std::cout<<"I'm in MssmHbbSignal::leadingJetSelection"<<std::endl;

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
	if(!cuts_.check("deta12",std::abs(jet1.eta() - jet2.eta()) <= dEta_)) return false;

	//deltaR requirements
	if(!cuts_.check("dR12",jet1.deltaR(jet2) > dR_)) return false;

	//BTag requirements
//	if(!cuts_.check("btag1",jet1.btag() >= btag1_)) return false;
//	if(!cuts_.check("btag2",jet2.btag() >= btag2_)) return false;



	//Trigger Selection
	if(!cuts_.check("TriggerBit",this->triggerResult(triggerLogicName_))) return false;

	//Online selection:
	if(!cuts_.check("TriggerMatching",this->OnlineSelection(jet1,jet2))) return false;

	if(nJets_ == 3) {
		Jet jet3 = offlineJets->at(2);
		//Pt requirements
		if(!cuts_.check("pt3",jet3.pt() >= pt3_)) return false;

		//Eta requirements
		if(!cuts_.check("eta3",std::abs(jet3.eta()) <= eta3_)) return false;

		//deltaR requirements
		if(!cuts_.check("dR23",jet2.deltaR(jet3) > dR_)) return false;
		if(!cuts_.check("dR13",jet1.deltaR(jet3) > dR_)) return false;

		//BTag requirements
//		if(!cuts_.check("btag3",jet3.btag() >= btag3_)) return false;
	}

	return true;
}

void CutOptimisation::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){
	if(TEST) std::cout<<"I'm in CutOptimisation::fillHistograms"<<std::endl;

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);
	Jet jet3 = offlineJets->at(2);
	TLorentzVector obj12;
	obj12 = jet1.p4() + jet2.p4();

	(histo_.getHisto())["diJet_m"]->Fill(obj12.M(),weight);
	bool standard_cuts = jet1.btag() > btag1_ && jet2.btag() > btag2_ && jet3.btag() > btag3_;

	//Cut optimisation
	std::vector<std::string> pt_cuts = {"40","50","60","70","80","90","100","120","140","160","180","200"};
	std::vector<std::string> d_csv	 = {"0.46","0.8","0.935"};
	std::vector<std::string> d_cmva	 = {"-0.715","0.185","0.875"};
	for(const auto & s : pt_cuts) {
		if(jet1.pt() >= std::stod(s) && standard_cuts ) (histo_.getHisto())[("jet1_m12_pt>"+s)]->Fill(obj12.M(),weight);
		if(jet2.pt() >= std::stod(s) && standard_cuts ) (histo_.getHisto())[("jet2_m12_pt>"+s)]->Fill(obj12.M(),weight);
		if(jet3.pt() >= std::stod(s) && standard_cuts ) (histo_.getHisto())[("jet3_m12_pt>"+s)]->Fill(obj12.M(),weight);
	}
	for(const auto & s1 : d_csv) {
		for(const auto & s2 : d_csv) {
			for(const auto & s3 : d_csv) {
				if(jet1.btag() > std::stod(s1) && jet2.btag() > std::stod(s2) && jet3.btag() > std::stod(s3)) (histo_.getHisto())[("jet1_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
				if(jet2.btag() > std::stod(s1) && jet1.btag() > std::stod(s2) && jet3.btag() > std::stod(s3)) (histo_.getHisto())[("jet2_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
				if(jet3.btag() > std::stod(s1) && jet1.btag() > std::stod(s2) && jet2.btag() > std::stod(s3)) (histo_.getHisto())[("jet3_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
			}
		}
	}
	for(const auto & s1 : d_cmva) {
		for(const auto & s2 : d_cmva) {
			for(const auto & s3 : d_cmva) {
		if(jet1.btag("btag_csvmva") > std::stod(s1) && jet2.btag("btag_csvmva") > std::stod(s2) && jet3.btag("btag_csvmva") > std::stod(s3)) (histo_.getHisto())[("jet1_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
		if(jet2.btag("btag_csvmva") > std::stod(s1) && jet1.btag("btag_csvmva") > std::stod(s2) && jet3.btag("btag_csvmva") > std::stod(s3)) (histo_.getHisto())[("jet2_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
		if(jet3.btag("btag_csvmva") > std::stod(s1) && jet1.btag("btag_csvmva") > std::stod(s2) && jet2.btag("btag_csvmva") > std::stod(s3)) (histo_.getHisto())[("jet3_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)]->Fill(obj12.M(),weight);
			}
		}
	}



}

void CutOptimisation::runAnalysis(const std::string &json, const std::string &output, const int &size){
	if(TEST) std::cout<<"I'm in CutOptimisation::runAnalysis"<<std::endl;

	this->SetupStandardOutputFile(output);
	this->setupAnalysis(json);
	std::cout<<"Total number of events: "<<this->size()<<std::endl;
	this->makeHistograms(size);
	this->applySelection();

	this->writeHistograms();


}

const double CutOptimisation::assignWeight(){
	double weight = 1;
	std::string file_name = outputFile_->GetName();
	if(isMC()) {
		weight = weight_["PtEff_central"] * weight_["PU_central"];// * weight_["SFb_central"] * weight_["SFl_central"];
//		if(file_name.find("madgraphMLM") != std::string::npos){
//			weight = weight * weight_["Ht"];
//		}
	}
	return weight;
}

void CutOptimisation::writeHistograms(){
	std::string full_name;
	if(TEST) std::cout<<"I'm in CutOptimisation::writeHistograms"<<std::endl;
	outputFile_->cd();

	for(const auto & h: (histo_.getHisto())){
		if(h.second->GetEntries() == 0) continue; 			//skip empty histograms
		if(isMC()) h.second->Scale(weight_["Lumi"]);
		h.second->Write();
	}
	outputFile_->Close();
}
