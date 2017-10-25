/*
 * doubleB.cpp
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/X750Search.h"

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

X750Search::X750Search(const std::string & inputFilelist, const double & dataLumi, const bool & lowM, const bool & test) :
TriggerEfficiency(inputFilelist,dataLumi,lowM,test) {
	nJets_ = 2;
	if(lowM_){
		pt1_ = 100;
		pt2_ = 100;
	}
	else{
		pt1_ = 160.;
		pt2_ = 160.;
	}
	eta1_ = 2.2;
	eta2_ = 2.2;
	dR_ = 1.;
	triggerLogicName_ = "HLT_PFJet60_v";
	triggerObjectName_ = {"hltL1sL1SingleJet36","hltSingleCaloJet40","hltSinglePFJet60","hltL1sL1SingleJet52","hltSingleCaloJet50","hltSinglePFJet80"};
	baseOutputName_ = "X750Eff";
}

X750Search::~X750Search() {
	// TODO Auto-generated destructor stub
}

void X750Search::applySelection(){
	//	if(TEST) std::cout<<"I'm in applySelection"<<std::endl;
		bool goodLeadingJets = false;

		double totWeight = 0;
		TLorentzVector diJetObject;
		double Ht;
		//Event loop:
		auto nevents = 0;
		if(TEST) nevents = 0.4*this->size();
		else nevents = this->size();

		std::cout<<"Events to process: "<<nevents<<std::endl;

		for(auto i = 0; i < nevents; ++ i){
			this->event(i);

			//Select only good JSon runs for Data
			if(!isMC() && !this->selectJson() ) continue;

		    //Define Jet Collection
			auto offlineJets = this->collection<Jet>("Jets");
			auto newLeadJets = std::make_shared<Collection<Jet> >();

			//Define MC specific collections:
			if(isMC()){
				genJets_ = (std::shared_ptr<Collection<Jet> >) this->collection<Jet>("GenJets");
			}

			if(offlineJets->size() < nJets_ ) continue;

			//Match offline Jets to online Objects

			this->match<Jet,TriggerObject>("Jets",this->getTriggerObjectNames());

		    //Define Vertex collection
		    auto offlinePrimaryVertices = this->collection<Vertex>("Vertices");

		    goodLeadingJets = false;
		    Ht = 0.;
		    double btagmax1 = -15;
		    int index_max1 = 0;
		    double btagmax2 = -15;
//		    int index_max2 = 0;
		    int counter = 0;
			Jet LeadJet[10];
		    //Jet Selection
		    for( int iJet = 0; iJet < offlineJets -> size(); ++iJet){
		    	Jet jet = offlineJets->at(iJet);

		    	//ht calculation
		    	if(jet.idLoose()) Ht +=jet.pt();

		    		//Selection of good Leading Jets:
		    		if(!jet.idLoose()) break;
		    		if(jet.pt() < 100) continue;
				    if(jet.btag() < 0.1) continue;
		    		if(std::abs(jet.eta()) > 2.4) continue;

		    		if(jet.btag() > btagmax1) {
		    			LeadJet[0] = jet;
		    			btagmax1 = jet.btag();
		    			index_max1 = iJet;
		    		}
		    		else if (jet.btag() >= btagmax2 && iJet != index_max1){
		    			LeadJet[1] = jet;
		    			btagmax2 = jet.btag();
//		    			index_max2 = iJet;
		    		}

		    		++counter;

		    		if(counter >= nJets_) goodLeadingJets = true;
		    }
		    if(!goodLeadingJets) continue;

		    //Leading Jet Selection:
		    if(LeadJet[0].pt() < pt1_) continue;
		    if(LeadJet[1].pt() < pt2_) continue;

		    if(LeadJet[0].btag() < 0.1) continue;
		    if(LeadJet[1].btag() < 0.1) continue;

			if (std::abs(LeadJet[0].eta()) > eta1_) continue;
			if (std::abs(LeadJet[1].eta()) > eta2_) continue;

			if (LeadJet[0].deltaR(LeadJet[1]) <= dR_) continue;

			newLeadJets = this->modifyJetCollection(LeadJet[0],newLeadJets);
			newLeadJets = this->modifyJetCollection(LeadJet[1],newLeadJets);

//			std::cout<<"Lead: Njets = "<<counter<<" Index: 1 = "<<index_max1<<" 2 = "<<index_max2<<" btag = "<<LeadJet[0].btag()<<" "<<LeadJet[1].btag()<<" pt = "<<LeadJet[0].pt()<<" "<<LeadJet[1].pt()<<std::endl;

		    //Should return true for trigger efficiency study
		    //Weights:
		    if(isMC()){
		    	//.........................Calculate weights...............
		    	//TODO: Data luminosity!!!
		    	if(!signalMC_) {
		    		weight_["Lumi"] 			=pWeight_->LumiWeight(dataLumi_,this->luminosity());
		    	}
//		    	else weight_["Lumi"] 			=pWeight_->LumiWeight(dataLumi_,this->numberFilteredGenEvents()/xsection_[returnMassPoint()]);
		    	weight_["Ht"]       		= pWeight_->HtWeight(hCorrections1D_["hHtWeight"],Ht);

		    }
		    totWeight = this->assignWeight();
		    this->fillHistograms(newLeadJets,totWeight);
			(histo_.getHisto())["NumberOfSelectedEvents"] -> Fill(offlineJets->size(),totWeight);

		}
			double nGenMHat = 0;
//		(histo_.getHisto())["TotalNumberOfGenEvents"] -> Fill(this->numberGenEvents());
		(histo_.getHisto())["NumberOfFilteredEvents"] -> Fill(this->size());
		if(signalMC_) {
			(histo_.getHisto())["xsection"] -> Fill(xsection_[returnMassPoint()]);
//			nGenMHat = this->numberFilteredGenEvents() * weight_["Lumi"];
			(histo_.getHisto())["NumberOfGenEvents_afterMHat"] -> Fill(nGenMHat);
		}
		else (histo_.getHisto())["xsection"] -> Fill(this->crossSection());

	//	this->writeHistograms();
	//	outputFile_->Close();
}

