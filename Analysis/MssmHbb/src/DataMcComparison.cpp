/*
 * doubleB.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/DataMcComparison.h"
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem ;

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

DataMcComparison::DataMcComparison(const std::string & inputFilelist, const double & dataLumi, const bool & lowM,const int & njets, const bool & test) :
		selectionDoubleB(inputFilelist,dataLumi,lowM,test) {
	baseOutputName_ = "DataMC";
	if(njets == 3)
	baseOutputName_ += "_3b";
	nJets_ = njets;
	JESshift_ = 0;
	JERshift_ = 0;
    if(lowM){
    	triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
    	triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp014DoubleWithMatching","hltDoublePFJetsC100","hltDoublePFJetsC100MaxDeta1p6"};
    	pt1_ = 100.; pt2_ = 100.; pt3_ = 40;
    	eta1_ = 2.2; eta2_ = 2.2; eta3_ = 2.2;
    	btag1_ = 0.8484; btag2_ = 0.8484; btag3_ = 0.8484;
    	btagOP1_ = 1; btagOP2_ = 1; btagOP3_ = 1;
    	dR_ = 1; dEta_ = 1.55; mHat_ = 0.7;
    }
    else {
    	triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v";
    	triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp026DoubleWithMatching","hltDoublePFJetsC160"};
    	pt1_ = 160.; pt2_ = 160.; pt3_ = 40;
    	eta1_ = 2.2; eta2_ = 2.2; eta3_ = 2.2;
    	btag1_ = 0.8484; btag2_ = 0.8484; btag3_ = 0.8484;
    	btagOP1_ = 1; btagOP2_ = 1; btagOP3_ = 1;
    	dR_ = 1; dEta_ = 100.; mHat_ = 0.7;
    }
}

DataMcComparison::~DataMcComparison() {
	// TODO Auto-generated destructor stub
}


const bool DataMcComparison::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){

	if(TEST) std::cout<<"I'm in MssmHbbSignal::leadingJetSelection"<<std::endl;

	//Selection of good Leading Jets:
	//Only jets that pass Loose identification will be considered
	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);

	//To avoid code duplication frind class were used
	if(!selectionDoubleB::leadingJetSelection(offlineJets)) return false;

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
		if(!cuts_.check("btag3",jet3.btag() >= btag3_)) return false;


	}

	return true;
}

void DataMcComparison::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){
	if(TEST) std::cout<<"I'm in DataMcComparison::fillHistograms"<<std::endl;

	jet1_ = offlineJets->at(0);
	jet2_ = offlineJets->at(1);
	if(nJets_ == 3) jet3_ = offlineJets->at(2);

	diJetObject_ = jet1_.p4() + jet2_.p4();
	jet_size__ = offlineJets->size();

	double syst_weight = weight;
	std::string syst;
	std::string variation;
	
	if(JERshift_ == 0 && JESshift_ == 0){
		(histo_.getHisto())["NumberOfJets"]->Fill(offlineJets->size());
		selectionDoubleB::fillHistograms(offlineJets,weight);
		//Central part:
		histoToFill("","",weight);
		if(isMC()){
			//Pt Efficiency
			syst = "_eff_pTonl_";
			variation = "up";
			syst_weight = weight/weight_["eff_pTonl_central"] * weight_["eff_pTonl_up"];
			histoToFill(syst,variation,syst_weight);

			variation = "down";
			syst_weight = weight/weight_["eff_pTonl_central"] * weight_["eff_pTonl_down"];
			histoToFill(syst,variation,syst_weight);

			//PU
			syst = "_PU_";
			variation = "up";
			syst_weight = weight/weight_["PU_central"] * weight_["PU_up"];
			histoToFill(syst,variation,syst_weight);

			variation = "down";
			syst_weight = weight/weight_["PU_central"] * weight_["PU_down"];
			histoToFill(syst,variation,syst_weight);

			//SFb
			syst = "_eff_b_";
			variation = "up";
			syst_weight = weight/weight_["eff_b_central"] * weight_["eff_b_up"];
			histoToFill(syst,variation,syst_weight);

			variation = "down";
			syst_weight = weight/weight_["eff_b_central"] * weight_["eff_b_down"];
			histoToFill(syst,variation,syst_weight);

			//SFl
			syst = "_eff_l_";
			variation = "up";
			syst_weight = weight/weight_["eff_l_central"] * weight_["eff_l_up"];
			histoToFill(syst,variation,syst_weight);

			variation = "down";
			syst_weight = weight/weight_["eff_l_central"] * weight_["eff_l_down"];
			histoToFill(syst,variation,syst_weight);

			//Online BTag
			syst = "_eff_bonl_";
			variation = "up";
			syst_weight = weight/weight_["eff_bonl_central"] * weight_["eff_bonl_up"];
			histoToFill(syst,variation,syst_weight);

			variation = "down";
			syst_weight = weight/weight_["eff_bonl_central"] * weight_["eff_bonl_down"];
			histoToFill(syst,variation,syst_weight);
		}
	}
	else if (JESshift_ < 0 && JERshift_ == 0){
		syst = "_scale_j_";
		variation = "down";
		histoToFill(syst,variation,weight);
	}
	else if (JESshift_ > 0 && JERshift_ == 0){
		syst = "_scale_j_";
		variation = "up";
		histoToFill(syst,variation,weight);
	}
	else if (JERshift_ > 0 && JESshift_ == 0){
		syst = "_res_j_";
		variation = "up";
		histoToFill(syst,variation,weight);
	}
	else if (JERshift_ < 0 && JESshift_ == 0){
		syst = "_res_j_";
		variation = "down";
		histoToFill(syst,variation,weight);
	}
	else {
		std::cerr<<"ERROR in DataMcComparison::fillHistograms"<<std::endl;
		exit(-1);
	}

}

void DataMcComparison::runAnalysis(const std::string &json, const std::string &output, const int &size){
	if(TEST) std::cout<<"I'm in DataMcComparison::runAnalysis"<<std::endl;

	this->SetupStandardOutputFile(output);
	this->setupAnalysis(json);
	std::cout<<"Total number of events: "<<this->size()<<std::endl;
	this->makeHistograms(size);
/*	if(isMC()){
		for(int i = 0; i < 3 ; ++i){
		JERshift_ = 0;
		if(i==0) JESshift_ = 0;
		else if (i == 1){
			JESshift_ = 1;
		}
		else {
			JESshift_ = - 1;
		}
		this->applySelection();
		}

		for(int i = 1; i < 3; ++i){
			JESshift_ = 0;
			if (i == 1){
				JERshift_ = 2;
			}
			else {
				JERshift_ = -2;
			}
			this->applySelection();
		}
	}
	else */
		this->applySelection();

	this->writeHistograms();


}

std::shared_ptr<tools::Collection<tools::Jet> > DataMcComparison::modifyJetCollection(tools::Jet & jet,
						std::shared_ptr<tools::Collection<tools::Jet> > & initialJets
						){
	Jet shiftedJet = jet;
	TLorentzVector p4Shift = jet.p4() * (1 + 2.* JESshift_ * jet.jecUncert());
	jet.p4(p4Shift);
	initialJets->add(jet);
	return initialJets;
}

const double DataMcComparison::assignWeight(){
	double weight = 1;
	std::string file_name = outputFile_->GetName();
	if(isMC()) {
//		weight = weight_["PtEff_central"] * weight_["PU_central"] * weight_["SFb_central"] * weight_["SFl_central"];
//		weight *= weight_["BTagEff_central"];
		if(file_name.find("madgraphMLM") != std::string::npos){
//			weight = weight * weight_["Ht"];
		}
	}
	return weight;
}

void DataMcComparison::histoToFill(const std::string & syst, const std::string & variation, const double & weight){

	(histo_.getHisto())["jet_b_pt1"+syst+variation]->Fill(jet1_.pt(),weight);
	(histo_.getHisto())["jet_b_pt2"+syst+variation]->Fill(jet2_.pt(),weight);

	(histo_.getHisto())["jet12_b_assym"+syst+variation]->Fill((jet1_.pt()-jet2_.pt())/(jet1_.pt()+jet2_.pt()),weight);

	(histo_.getHisto())["jet_b_eta1"+syst+variation]->Fill(jet1_.eta(),weight);
	(histo_.getHisto())["jet_b_eta2"+syst+variation]->Fill(jet2_.eta(),weight);

	(histo_.getHisto())["jet_b_phi1"+syst+variation]->Fill(jet1_.phi(),weight);
	(histo_.getHisto())["jet_b_phi2"+syst+variation]->Fill(jet2_.phi(),weight);

	(histo_.getHisto())["jet_b_btag_csv1"+syst+variation]->Fill(jet1_.btag(),weight);
	(histo_.getHisto())["jet_b_btag_csv2"+syst+variation]->Fill(jet2_.btag(),weight);

	(histo_.getHisto())["jet_b_btag_cmva1"+syst+variation]->Fill(jet1_.btag("btag_csvmva"),weight);
	(histo_.getHisto())["jet_b_btag_cmva2"+syst+variation]->Fill(jet2_.btag("btag_csvmva"),weight);

	(histo_.getHisto())["jet_b_deta12"+syst+variation]->Fill(jet1_.eta() - jet2_.eta(),weight);
	(histo_.getHisto())["jet_b_dphi12"+syst+variation]->Fill(  acos(cos(jet1_.phi() - jet2_.phi())) ,weight);

	(histo_.getHisto())["jet_b_dR12"+syst+variation]->Fill(jet1_.deltaR(jet2_),weight);

	(histo_.getHisto())["diJet_b_pt"+syst+variation]->Fill(diJetObject_.Pt(),weight);
	(histo_.getHisto())["diJet_b_eta"+syst+variation]->Fill(diJetObject_.Eta(),weight);
	(histo_.getHisto())["diJet_b_phi"+syst+variation]->Fill(diJetObject_.Phi(),weight);
	(histo_.getHisto())["diJet_b_m"+syst+variation]->Fill(diJetObject_.M(),weight);

	(histo_.getHisto())["jet_b_Ht"+syst+variation]->Fill(this->Ht(),weight);
	(histo_.getHisto())["NumberOfJets_b"+syst+variation]->Fill(jet_size__,weight);
	if(nJets_ == 3){
		(histo_.getHisto())["jet_b_pt3"+syst+variation]->Fill(jet3_.pt(),weight);
		(histo_.getHisto())["jet_b_eta3"+syst+variation]->Fill(jet3_.eta(),weight);
		(histo_.getHisto())["jet_b_phi3"+syst+variation]->Fill(jet3_.phi(),weight);
		(histo_.getHisto())["jet_b_btag_csv3"+syst+variation]->Fill(jet3_.btag(),weight);
		(histo_.getHisto())["jet_b_btag_cmva3"+syst+variation]->Fill(jet3_.btag("btag_csvmva"),weight);
		(histo_.getHisto())["jet_b_dR23"+syst+variation]->Fill(jet2_.deltaR(jet3_),weight);
		(histo_.getHisto())["jet_b_dR13"+syst+variation]->Fill(jet1_.deltaR(jet3_),weight);
		(histo_.getHisto())["jet_b_deta13"+syst+variation]->Fill(jet1_.eta()-jet3_.eta(),weight);
		(histo_.getHisto())["jet_b_deta23"+syst+variation]->Fill(jet2_.eta()-jet3_.eta(),weight);
		(histo_.getHisto())["jet_b_dphi13"+syst+variation]->Fill( acos(cos(jet1_.phi() - jet3_.phi())),weight);
		(histo_.getHisto())["jet_b_dphi23"+syst+variation]->Fill( acos(cos(jet2_.phi() - jet3_.phi())),weight);

	}

}

void DataMcComparison::writeHistograms(){
	std::string full_name;
	if(TEST) std::cout<<"I'm in DataMcComparison::writeHistograms"<<std::endl;

	outputFile_->mkdir("jet_b_pt1"," Pt1 distribution");
	outputFile_->mkdir("jet_b_pt2"," Pt2 distribution");
	outputFile_->mkdir("jet_b_pt3"," Pt3 distribution");
	outputFile_->mkdir("jet_b_Ht"," Ht distribution");
	outputFile_->mkdir("jet_b_eta1"," Eta1 distribution");
	outputFile_->mkdir("jet_b_eta2"," Eta2 distribution");
	outputFile_->mkdir("jet_b_eta3"," Eta3 distribution");
	outputFile_->mkdir("jet_b_phi1"," Phi1 distribution");
	outputFile_->mkdir("jet_b_phi2"," Phi2 distribution");
	outputFile_->mkdir("jet_b_phi3"," Phi3 distribution");
	outputFile_->mkdir("jet12_b_assym"," Pt assymetry distribution");
	outputFile_->mkdir("jet_b_btag_csv1"," BTagCSVv2 1 distribution");
	outputFile_->mkdir("jet_b_btag_csv2"," BTagCSVv2 2 distribution");
	outputFile_->mkdir("jet_b_btag_csv3"," BTagCSVv2 3 distribution");
	outputFile_->mkdir("jet_b_btag_cmva1"," BTagCMVA1 distribution");
	outputFile_->mkdir("jet_b_btag_cmva2"," BTagCMVA2 distribution");
	outputFile_->mkdir("jet_b_btag_cmva3"," BTagCMVA3 distribution");
	outputFile_->mkdir("jet_b_deta12", "dEta 12");
	outputFile_->mkdir("jet_b_dphi12", "dPhi 12");
	outputFile_->mkdir("jet_b_dR12"," dR 12");
	outputFile_->mkdir("jet_b_dR23"," dR 23");
	outputFile_->mkdir("jet_b_dR13"," dR 13");
	outputFile_->mkdir("diJet_b_pt"," Pt of di-Jet object");
	outputFile_->mkdir("diJet_b_eta"," Eta of di-Jet object");
	outputFile_->mkdir("diJet_b_phi"," Phi of di-Jet object");
	outputFile_->mkdir("diJet_b_m"," M of di-Jet object");
	outputFile_->mkdir("NumberOfJets_b"," Number of jets");
	outputFile_->mkdir("general", " General distributions");


	for(const auto & h: (histo_.getHisto())){
		if(h.second->GetEntries() == 0) continue; 			//skip empty histograms
		if(isMC()) h.second->Scale(weight_["Lumi"]);
		if(h.first.find("jet_b_pt1") != std::string::npos){
			outputFile_->cd("jet_b_pt1");
		}
		else if (h.first.find("jet_b_pt2") != std::string::npos){
			outputFile_->cd("jet_b_pt2");
		}
		else if (h.first.find("jet_b_pt3") != std::string::npos){
			outputFile_->cd("jet_b_pt3");
		}
		else if (h.first.find("jet_b_Ht") != std::string::npos){
			outputFile_->cd("jet_b_Ht");
		}
		else if (h.first.find("jet_b_eta1") != std::string::npos){
			outputFile_->cd("jet_b_eta1");
		}
		else if (h.first.find("jet_b_eta2") != std::string::npos){
			outputFile_->cd("jet_b_eta2");
		}
		else if (h.first.find("jet_b_eta3") != std::string::npos){
			outputFile_->cd("jet_b_eta3");
		}
		else if (h.first.find("jet_b_phi1") != std::string::npos){
			outputFile_->cd("jet_b_phi1");
		}
		else if (h.first.find("jet_b_phi2") != std::string::npos){
			outputFile_->cd("jet_b_phi2");
		}
		else if (h.first.find("jet_b_phi3") != std::string::npos){
			outputFile_->cd("jet_b_phi3");
		}
		else if (h.first.find("jet12_b_assym") != std::string::npos){
			outputFile_->cd("jet12_b_assym");
		}
		else if (h.first.find("jet_b_btag_csv1") != std::string::npos){
			outputFile_->cd("jet_b_btag_csv1");
		}
		else if (h.first.find("jet_b_btag_csv2") != std::string::npos){
			outputFile_->cd("jet_b_btag_csv2");
		}
		else if (h.first.find("jet_b_btag_csv3") != std::string::npos){
			outputFile_->cd("jet_b_btag_csv3");
		}
		else if (h.first.find("jet_b_btag_cmva1") != std::string::npos){
			outputFile_->cd("jet_b_btag_cmva1");
		}
		else if (h.first.find("jet_b_btag_cmva2") != std::string::npos){
			outputFile_->cd("jet_b_btag_cmva2");
		}
		else if (h.first.find("jet_b_btag_cmva3") != std::string::npos){
			outputFile_->cd("jet_b_btag_cmva3");
		}
		else if (h.first.find("jet_b_deta12") != std::string::npos){
			outputFile_->cd("jet_b_deta12");
		}
		else if (h.first.find("jet_b_dphi12") != std::string::npos){
			outputFile_->cd("jet_b_dphi12");
		}
		else if (h.first.find("jet_b_dR12") != std::string::npos){
			outputFile_->cd("jet_b_dR12");
		}
		else if (h.first.find("jet_b_dR13") != std::string::npos){
			outputFile_->cd("jet_b_dR13");
		}
		else if (h.first.find("jet_b_dR23") != std::string::npos){
			outputFile_->cd("jet_b_dR23");
		}
		else if (h.first.find("diJet_b_pt") != std::string::npos){
			outputFile_->cd("diJet_b_pt");
		}
		else if (h.first.find("diJet_b_eta") != std::string::npos){
			outputFile_->cd("diJet_b_eta");
		}
		else if (h.first.find("diJet_b_phi") != std::string::npos){
			outputFile_->cd("diJet_b_phi");
		}
		else if (h.first.find("diJet_b_m") != std::string::npos){
			outputFile_->cd("diJet_b_m");
		}
		else if (h.first.find("NumberOfJets_b") != std::string::npos){
			outputFile_->cd("NumberOfJets_b");
		}
		else{
			outputFile_->cd("general");
		}
		h.second->Write();
		outputFile_->cd("");
	}
	outputFile_->Close();
}

void DataMcComparison::SetupStandardOutputFile(const std::string & outputFileName){

		//get the full file name and path from the Tree
		std::string fullName = fs::basename(inputFilelist_);
		std::string outputName = baseOutputName_;
		if(outputFileName != "") outputName+="_"+outputFileName;
		if(lowM_) outputName += "_lowM_";
		else {outputName += "_highM_";}


		outputName += fullName;// + ".root";
		this->createOutputFile(outputName);
//	else{
//		this->createOutputFile(outputFileName);
//	}

}
