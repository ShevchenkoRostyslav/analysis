/*
 * JetAnalysisBase.cpp
 *
 *  Created on: 17 апр. 2016 г.
 *      Author: rostyslav
 */
#include <stdlib.h>
#include <chrono>
#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"


using namespace analysis;
using namespace analysis::tools;
using namespace analysis::mssmhbb;
//using namespace std::chrono;

double signal_shape(const double & M12_gen){
	double result = 1;
	double par[6];
	par[0] = 7.08179e-01; par[1] = 7.74883e+00; par[2] = 1.10010e+03;
	par[3] = 4.48655e-02; par[4] = 8.82140e+02; par[5] = 2.12431e+02;
	double M  = M12_gen;
	double M2 = M*M;
	double G = par[1];
	double G2 = G*G;
	double E = par[2];
	double E2= E*E;
	double rel_BW = par[0] * 2 / TMath::Pi() * G2*E2 / ( (M2-E2)*(M2-E2)+ M2*M2*(G2/E2) );
	double gaus = par[3] * TMath::Gaus(M,par[4],par[5]);
	if(M >= 1000 && M <= 1200) result = rel_BW / gaus;
	else result = 0;
	return result;
}

JetAnalysisBase::JetAnalysisBase(const std::string & inputFilelist, const double & dataLumi, const bool & lowM, const bool & test) :
								 Analysis(inputFilelist,"MssmHbb/Events/EventInfo"),
								 dataLumi_(dataLumi),
								 lowM_(lowM),
								 triggerLogicName_(""),
								 nJets_(1),
								 TEST(test){
	if(this->isMC()){
		//Define MC type
		if(findStrings(inputFilelist,"susy")) signalMC_ = true;
		else signalMC_ = false;
		if(signalMC_) this->setupXSections();

		//Add specific to MC trees
		this->addTree<GenParticle>("GenParticles","MssmHbb/Events/prunedGenParticles");
		if(isMC()) this->addTree<Jet>("GenJets","MssmHbb/Events/slimmedGenJets");

		// Add MC information:
		this->crossSections("MssmHbb/Metadata/CrossSections");
		this->generatorFilter("MssmHbb/Metadata/GeneratorFilter");
		if(signalMC_) this->eventFilter("MssmHbb/Metadata/EventFilter");//,"MssmHbb/mHatFilter/EventFilter");
		else this->eventFilter("MssmHbb/Metadata/EventFilter");

		// Show MC information
		this->listCrossSections();
		this->listGeneratorFilter();

		if(isMC()){
			std::cout<<"***********************************"<<std::endl;
			if(signalMC_) std::cout<<"This is signal SUSY MC"<<std::endl;
			else std::cout<<" This is BG MC"<<std::endl;
			std::cout<<"***********************************"<<std::endl;
		}
		else {
			std::cout<<"***********************************"<<std::endl;
			std::cout<<"This is Real Data"<<std::endl;
			std::cout<<"***********************************"<<std::endl;
		}
	}
	JESshift_ = 0;
	JERshift_ = 0;
	//Muon tree
	this->addTree<Muon> ("Muons","MssmHbb/Events/slimmedMuons");

	//Trigger trees
	this->triggerResults("MssmHbb/Events/TriggerResults");

	// Tree for Jets
	this->addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");
	// Tree for Vertices
	this->addTree<Vertex> ("Vertices","MssmHbb/Events/offlineSlimmedPrimaryVertices");

	baseOutputName_ = "JetAnalysisBase";

	std::string selection_type;
	if(lowM){
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp014DoubleWithMatching","hltDoublePFJetsC100","hltDoublePFJetsC100MaxDeta1p6"};
        pt1_ = 100.; pt2_ = 100.;
        btag1_ = 0.8484; btag2_ = 0.8484;
        btagOP1_ = 1; btagOP2_ = 1;
        btag3_ = 0.8484;
        btagOP3_ = 1;	//Mid OP
        selection_type = "Low Mass";
	}
	else {
        btag3_ = 0.46;
        btagOP3_ = 0; // Loose OP
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltBTagCaloCSVp026DoubleWithMatching","hltDoublePFJetsC160"};
        pt1_ = 160.; pt2_ = 160.;
        btag1_ = 0.8484; btag2_ = 0.8484;
        btagOP1_ = 1; btagOP2_ = 1;
        selection_type = "High Mass";
	}
//	cuts_ = CutFlow(baseOutputName_,selection_type);
	mHat_ = 0.7;
}

JetAnalysisBase::JetAnalysisBase(const std::string & inputFilelist, const bool & lowM) :
								 Analysis(inputFilelist,"MssmHbb/Events/EventInfo"),
								 dataLumi_(1),
								 lowM_(lowM),
								 triggerLogicName_(""),
								 nJets_(1),
								 TEST(false){
	if(this->isMC()){
		//Define MC type
		if(findStrings(inputFilelist,"susy")) signalMC_ = true;
		else signalMC_ = false;
		if(signalMC_) this->setupXSections();

		//Add specific to MC trees
		this->addTree<GenParticle>("GenParticles","MssmHbb/Events/prunedGenParticles");
		if(isMC()) this->addTree<Jet>("GenJets","MssmHbb/Events/slimmedGenJets");

		// Add MC information:
		this->crossSections("MssmHbb/Metadata/CrossSections");
		this->generatorFilter("MssmHbb/Metadata/GeneratorFilter");
		if(signalMC_) this->eventFilter("MssmHbb/Metadata/EventFilter");//,"MssmHbb/mHatFilter/EventFilter");
		else this->eventFilter("MssmHbb/Metadata/EventFilter");

		// Show MC information
		this->listCrossSections();
		this->listGeneratorFilter();

		if(isMC()){
			std::cout<<"***********************************"<<std::endl;
			if(signalMC_) std::cout<<"This is signal SUSY MC"<<std::endl;
			else std::cout<<" This is BG MC"<<std::endl;
			std::cout<<"***********************************"<<std::endl;
		}
		else {
			std::cout<<"***********************************"<<std::endl;
			std::cout<<"This is Real Data"<<std::endl;
			std::cout<<"***********************************"<<std::endl;
		}
	}
	JESshift_ = 0;
	JERshift_ = 0;
	//Muon tree
	this->addTree<Muon> ("Muons","MssmHbb/Events/slimmedMuons");

	//Trigger trees
	this->triggerResults("MssmHbb/Events/TriggerResults");

	// Tree for Jets
	this->addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");
	// Tree for Vertices
	this->addTree<Vertex> ("Vertices","MssmHbb/Events/offlineSlimmedPrimaryVertices");

	triggerLogicName_ = "";
	triggerObjectName_ = {};
	pt1_ = -100; pt2_ = -100;
	btag1_ = -100; btag2_ = -100;
	btagOP1_ = -100; btagOP2_ = -100;
	btag3_ = -100;
	btagOP3_ = -100;	//Mid OP

	baseOutputName_ = "JetAnalysisBase";

	std::string selection_type;
	if(lowM_){
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV0p9_DoublePFJetsC100MaxDeta1p6_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltDoublePFJetsC100","hltDoubleBTagCSV0p9","hltDoublePFJetsC100MaxDeta1p6"};
        selection_type = "Low Mass";
	}
	else {
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV0p85_DoublePFJetsC160_v";
		triggerObjectName_ = {"hltL1sDoubleJetC100","hltDoubleJetsC100","hltDoubleBTagCSV0p85","hltDoublePFJetsC160"};
        selection_type = "High Mass";
	}

	mHat_ = 0.7;
}

JetAnalysisBase::~JetAnalysisBase() {
	//if(TEST) std::cout<<"I'm at JetAnalysisBase::~JetAnalysisBase"<<std::endl;
	// TODO Auto-generated destructor stub
}

void JetAnalysisBase::setupAnalysis(const std::string & json){
	if(TEST) std::cout<<"I'm in setupAnalysis"<<std::endl;
	if(!isMC()) this->processJsonFile(json);
	if(isMC())  this->loadCorrections();
	this->addTriggerObjects(triggerObjectName_);
}

void JetAnalysisBase::applySelection(){

	if(TEST) std::cout<<"I'm in applySelection"<<std::endl;

	//different counters declaration
	int TotalNumberOfGenEvents = 0;						// Total number of generated events
	int NumberOfGenEvents_afterMHat = 0;				// Total number of generated events after Mhat cut
	double NumberOfGenEvents_afterMHat_rewPU = 0;		// Total number of generated events after Mhat cut and PU reweighting
	int NumberOfEventsAfterSelection = 0;				// Total number of the events that pass selection
	double NumberOfEventsAfterSelection_weighted = 0;	// Weighted number of the events that pass selection

	bool goodLeadingJets = false;
	double totWeight = 0;
	TLorentzVector diJetObject;
	double Ht = 0;
	std::vector<double> lumis = {19717, //  B-F
								 15956}; // G-H
	//Event loop:

//	Instance of BTagScaleFactor lib:
//	BTagScaleFactor BTagLib;

//	if(isMC()) BTagScaleFactor BTagLib_("csvv2",defaultInputDir_ + "bin/input_corrections/SFbLib.csv");

	int n_true_b = 0;

	std::string selection_type;
	if(lowM_)selection_type = "Low Mass";
	else selection_type = "High Mass";
	cuts_ = CutFlow(baseOutputName_,selection_type);

	std::array<Jet,5> LeadJet;
	std::vector<std::pair<int,Jet> > TrueBJets;
	TrueBJets.reserve(2);
	std::array<std::vector<BTagScaleFactor::ScaleFactor> , 3 > BTagSFs;

	if(TEST) TotalNumberOfGenEvents = 100;
	else TotalNumberOfGenEvents = this->size();
	std::cout<<"Events to process: "<<TotalNumberOfGenEvents<<std::endl;
//	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	for(auto i = 0; i < TotalNumberOfGenEvents; ++ i){
		this->event(i);

		cuts_.check("Total");

		if(isSignalMC() && !cuts_.check("mHat",(mHat() >= 0.7 * returnMassPoint()))) continue;
		++NumberOfGenEvents_afterMHat;

		//PU reweighting for MC
		if(isMC()) {
			weight_["PU_central"] 	= pWeight_->PileUpWeight(*hCorrections1D_["hPileUpData_central"],*hCorrections1D_["hPileUpMC"],this->nTruePileup());
			NumberOfGenEvents_afterMHat_rewPU += weight_["PU_central"];
//			NumberOfGenEvents_afterMHat_rewPU = NumberOfGenEvents_afterMHat;
			(histo_.getHisto())["nTruePileup"]->Fill(this->nTruePileup());
		}

		//Select only good JSon runs for Data
		if(!isMC()) if(!cuts_.check("json",this->selectJson() )) continue;

	    //Define Jet Collection
		auto offlineJets = this->collection<Jet>("Jets");
		auto shiftedJets = std::make_shared<Collection<Jet> >();

//		std::array<Jet,5> LeadJet;
//			std::vector<std::pair<int,Jet> > TrueBJets;
//			std::array<BTagScaleFactor::ScaleFactor,3> BTagSFs;

		//Clean jet vector and reserve the memory:
		TrueBJets.clear();
		TrueBJets.reserve(2);
//		Clean BTagSFs
//		BTagSFs.fill(BTagScaleFactor::ScaleFactor());
//		Clean LeadJet
		LeadJet.fill(Jet());

		//Define MC specific collections:
		if(isMC()){
			genJets_ = (std::shared_ptr<Collection<Jet> >) this->collection<Jet>("GenJets");
			//Apply JER smearing if this is MC:
		    offlineJets->matchTo(*genJets_,3,0.2);
		    offlineJets->smearTo(*genJets_,JERshift_);
		}

		if(!cuts_.check("nJets",offlineJets->size() >= nJets_)) continue;
		//Match offline Jets to online Objects

		this->match<Jet,TriggerObject>("Jets",this->getTriggerObjectNames());

	    //Define Vertex collection
	    auto offlinePrimaryVertices = this->collection<Vertex>("Vertices");

	    goodLeadingJets = false;

	    Ht = 0.;
	    //Jet Selection
	    for( int iJet = 0; iJet < offlineJets -> size(); ++iJet){
	    	Jet jet = offlineJets->at(iJet);
	    	shiftedJets = this->modifyJetCollection(jet,shiftedJets);
	    	jet = shiftedJets->at(iJet);

	    	if(isMC() && jet.flavour() == 5) {
	    		TrueBJets.push_back( std::make_pair(iJet,jet) );
	    		++n_true_b;
	    	}

	    	//ht calculation
	    	if(jet.idLoose() && jet.pt() > 20) Ht +=jet.pt();

	    	if(iJet < nJets_){
	    		//Selection of good Leading Jets:
	    		if(!jet.idLoose()) break;
//	    		if(!this->leadingJetSelection(shiftedJets)) break;

	    	    //Clear SFb
//	    		BTagSFs[iJet].BTagScaleFactor::ScaleFactor::clear();

	    		LeadJet[iJet] = jet;

	    		//Calculate SFb
	    		if(isMC()){
	    			if(iJet < 2){
	    				BTagSFs[iJet] = BTagLib_.calculateBTagSF(jet,btagOP1_);
	    			}
	    			if(iJet == 2){
	    				BTagSFs[iJet] = BTagLib_.calculateBTagSF(jet,btagOP3_);
	    			}
	    		}

	    		if(iJet == nJets_ - 1) goodLeadingJets = true;
	    	}
	    }
	    // STUDY OF THE FLAVOUR COMPOSITION
	    int j = 0;
	    if(isMC()){
	    	for(const auto & jpair : TrueBJets){
	    		std::string name = "bJet" + std::to_string(j) + "_pt";
	    		(histo_.getHisto())[name]->Fill(jpair.second.pt());
	    		name = "bJet" + std::to_string(j) + "_eta";
	    		(histo_.getHisto())[name]->Fill(jpair.second.eta());
	    		name = "bJet" + std::to_string(j) + "_btag";
	    		(histo_.getHisto())[name]->Fill(jpair.second.btag());
	    		name = "bJet" + std::to_string(j) + "_order";
	    		(histo_.getHisto())[name]->Fill(jpair.first);
	    		++j;
	    	}
	    }


	    (histo_.getHisto())["nTrue_b"]->Fill(n_true_b);
	    n_true_b = 0;
	    if(!cuts_.check("LoosID",goodLeadingJets)) continue;
//	    if(!goodLeadingJets) continue;
	    //Should return true for trigger efficiency study
	    if(!cuts_.check("LeadingJetSelection",leadingJetSelection(shiftedJets))) continue;
	    this->Ht(Ht);
	    //Weights:

	    if(isMC()){
	    	//.........................Calculate weights...............
        	// general weights which is independent of selection type OR already includes
        	// selection depends criterias inside.
//          weight_["FactorPt"] 		= pWeight_->FactorizationPtWeight(LeadJet[0].pt(), LeadJet[1].pt());

	    	double f1 		= pWeight_->triggerCorrectionFunction(LeadJet[0].pt(), LeadJet[0].eta());
	    	double eff1		= pWeight_->PtTriggerEfficiency(LeadJet[0].pt(), LeadJet[0].eta());
	    	double f2		= pWeight_->triggerCorrectionFunction(LeadJet[1].pt(), LeadJet[1].eta());
	    	double eff2		= pWeight_->PtTriggerEfficiency(LeadJet[1].pt(), LeadJet[1].eta());

	    	double df1	 	= f1 * 0.1 * (1. - eff1);
	    	double df2	 	= f2 * 0.1 * (1. - eff2);

	    	double F 		= f1*f2;
 	    	double dF		= F * std::sqrt( pow(df1/f1,2) + pow(df2/f2,2) );

	    	//For systematics
//	    	weight_["PtTriggerEfficiency"]= pWeight_->PtTriggerEfficiency(LeadJet[0].pt(), LeadJet[0].eta()) * pWeight_->PtTriggerEfficiency(LeadJet[1].pt(), LeadJet[1].eta());
	    	weight_["dEta"]     		= pWeight_->dEtaWeight(std::abs(LeadJet[0].eta() - LeadJet[1].eta()));
	    	//TODO: Data luminosity!!!

//	    	std::string file_name = outputFile_->GetName();
//	    	if(file_name.find("madgraphMLM") != std::string::npos && ( file_name.find("bEnriched") != std::string::npos || file_name.find("TuneCUETP8M1") != std::string::npos )){
	    	weight_["Ht"]  = pWeight_->HtWeight(hCorrections1D_["hHtWeight"],Ht);
//	    	}
//	    	else weight_["Ht"] = 1;

	    	//2SIGMA variation!!!!!!
	    	weight_["PU_down"]    	= pWeight_->PileUpWeight(*hCorrections1D_["hPileUpData_down"],*hCorrections1D_["hPileUpMC"],this->nTruePileup());
	    	weight_["PU_up"]      	= pWeight_->PileUpWeight(*hCorrections1D_["hPileUpData_up"],*hCorrections1D_["hPileUpMC"],this->nTruePileup());
//	    	std::cout<<"wt: "<<weight_["PU_central"]<<" "<<weight_["PU_down"]<<" "<<weight_["PU_up"]<<std::endl;
//	    	std::cout<<"wt: "<<weight_["PU_central"]<<" "<<weight_["PU_down"]/weight_["PU_central"]<<" "<<weight_["PU_up"]/weight_["PU_central"]<<"\n"<<std::endl;

	    	weight_["PtEff_central"]	= F;
	    	weight_["PtEff_up"]       	= F + 2.*dF;
	    	weight_["PtEff_down"]		= F - 2.*dF;

	    	std::pair<double,double> OnlineSFbs_j1, OnlineSFbs_j2;
	    	OnlineSFbs_j1 = pWeight_->BTagTriggerEfficiency(LeadJet[0].pt(),dataLumi_);
	    	OnlineSFbs_j2 = pWeight_->BTagTriggerEfficiency(LeadJet[1].pt(),dataLumi_);

	    	weight_["BTagEff_central"] = OnlineSFbs_j1.first * OnlineSFbs_j2.first;
	    	weight_["BTagEff_up"]      = (OnlineSFbs_j1.first + 2. * OnlineSFbs_j1.second) * (OnlineSFbs_j2.first + 2. * OnlineSFbs_j2.second);
	    	weight_["BTagEff_down"]    = (OnlineSFbs_j1.first - 2. * OnlineSFbs_j1.second) * (OnlineSFbs_j2.first - 2. * OnlineSFbs_j2.second);

	    	weight_["SFl_central"] = 1;
	    	weight_["SFl_up"] = 1;
	    	weight_["SFl_down"] = 1;

	    	weight_["SFb_central"] = 1;
	    	weight_["SFb_up"] = 1;
	    	weight_["SFb_down"] = 1;

	    	weight_["dR"] = 1;
	    	if(baseOutputName_ == "bbx") {
	    		weight_["dR"] = pWeight_->dRWeight(LeadJet[0].deltaR(LeadJet[1]));
	    		weight_["M12"] = pWeight_->M12Weight((LeadJet[0].p4() + LeadJet[1].p4()).M());
	    		weight_["Ht_bbx"] = pWeight_->Ht_bbxWeight(*hCorrections1D_["hHtWeight_bbx"],Ht);
	    	}

	    	combineBTagSFs(BTagSFs,lumis,dataLumi_);

//	    	double Lumi_GtoH = 16382;
//	    	double Lumi_BtoF = 20233;
//	    	double Lumi_tot = Lumi_GtoH + Lumi_BtoF;
	    	//..........................SF weight...............
//	    	for(const auto & sf : BTagSFs){
//	    		for(int i = 0; i<nJets_;++i){
//	    			double SF_GtoH = 0, SF_BtoF = 0, eSF_GtoH = 0, eSF_BtoF = 0, SF = 1;
//	    			for(const auto& sf : BTagSFs[i]){
//	    				if(sf.flavour != 0){ //only for b/c
//	    					SF *= sf.central * Lumi_BtoF / ;
//	    				}
//	    				else {
//
//	    				}
//	    			}
//
//
//	    			if(sf[i].flavour != 0){ //only for b/c
//	    				weight_["SFb_central"] 	*= sf[i].central;
//	    				weight_["SFb_up"]			*= (2.*(sf[i].up - sf[i].central) + sf[i].central);
//	    				weight_["SFb_down"]		*= (2.*(sf[i].down - sf[i].central) + sf[i].central);
//	    			}
//	    			else {
//	    				weight_["SFl_central"] 	*= sf[i].central;
//	    				weight_["SFl_up"]			*= (2.*(sf[i].up - sf[i].central) + sf[i].central);
//	    				weight_["SFl_down"]		*= (2.*(sf[i].down - sf[i].central) + sf[i].central);
//	    			}
//	    		}
//	    	}
//	    	combineBTagSFs(BTagSFs);
        	  //Selection depending weights
//        	  if(lowM_){
//            	  weight_["BTag"] = pWeight_->BTagWeight(hCorrections1D_["hRelBTagEff0p9"],hCorrections1D_["hRelBTagEff0p9_1p4"],hCorrections1D_["hRelBTagEff1p4_2p5"],LeadJet[0].pt(),LeadJet[0].eta()) *
//            			  	  	   pWeight_->BTagWeight(hCorrections1D_["hRelBTagEff0p9"],hCorrections1D_["hRelBTagEff0p9_1p4"],hCorrections1D_["hRelBTagEff1p4_2p5"],LeadJet[1].pt(),LeadJet[1].eta());
//        	  }
//        	  else {
//            	  weight_["BTag"] = pWeight_->BTagWeight(hCorrections2D_["hRelBTagEff2D"], LeadJet[0].pt(),LeadJet[0].eta())*
//            			  	  	   pWeight_->BTagWeight(hCorrections2D_["hRelBTagEff2D"], LeadJet[1].pt(), LeadJet[1].eta());
//        	  }

				if(signalMC_) {
					// PDF uncertainty
					pdfUncertainties(this->genScale(), this->pdf());
				}

	    }
//	    weight_["M12"] = pWeight_->M12Weight(*hCorrections1D_["hM12Weight_bbx"],(LeadJet[0].p4() + LeadJet[1].p4()).M());

	    ++NumberOfEventsAfterSelection;

	    totWeight = this->assignWeight();
	    this->fillHistograms(shiftedJets,totWeight);
		NumberOfEventsAfterSelection_weighted += totWeight;
	}

	double xsection = 0;
	std::string file_name = outputFile_->GetName();
	if(isMC()){
		if(!signalMC_) {
			xsection = this->crossSection();
			if(file_name.find("bEnriched") != std::string::npos || file_name.find("BGenFilter") != std::string::npos || file_name.find("TTJets") != std::string::npos)
				xsection = this->crossSection("myCrossSection");
		}
		else {
			xsection = xsection_[returnMassPoint()];
		}
		weight_["Lumi"] 			=pWeight_->LumiWeight(dataLumi_,NumberOfGenEvents_afterMHat_rewPU / xsection);
//		if(file_name.find("bEnriched") != std::string::npos || file_name.find("BGenFilter") != std::string::npos || file_name.find("QCD") != std::string::npos){
//			weight_["Lumi"] 			=pWeight_->LumiWeight(dataLumi_,1. / xsection); // To be able to Scale histograms after!!!
//		}
	}
//	std::cout<<"WTF: "<<weight_["Lumi"]<<std::endl;
	(histo_.getHisto())["xsection"]->Fill(xsection);
	(histo_.getHisto())["lumi_weight"]->Fill(weight_["Lumi"]);
	(histo_.getHisto())["TotalNumberOfGenEvents"]->Fill(TotalNumberOfGenEvents) ;
	(histo_.getHisto())["NumberOfGenEvents_afterMHat"]->Fill(NumberOfGenEvents_afterMHat);
	(histo_.getHisto())["NumberOfGenEvents_afterMHat_rewPU"]->Fill(NumberOfGenEvents_afterMHat_rewPU);

	cuts_.Print();

	if(isMC()){
		std::cout<<"***********************************"<<std::endl;
		std::cout<<"# generated events: "<<TotalNumberOfGenEvents<<std::endl;
		if(isSignalMC()) std::cout<<"# generated events after MHat cut: "<<NumberOfGenEvents_afterMHat<<std::endl;
		std::cout<<"# generated events after MHat cut and PU rew : "<<NumberOfGenEvents_afterMHat_rewPU<<std::endl;
		std::cout<<"# selected events: "<<NumberOfEventsAfterSelection<<std::endl;
		std::cout<<"# selected weighted evetns: "<<NumberOfEventsAfterSelection_weighted<<std::endl;
		std::cout<<"Cross section: "<<xsection<<" Effective Lumi: "<<NumberOfGenEvents_afterMHat_rewPU / xsection<<std::endl;
		std::cout<<"***********************************"<<std::endl;
	}
	else{
		std::cout<<"***********************************"<<std::endl;
		std::cout<<"Total number of Events: "<<TotalNumberOfGenEvents<<std::endl;
		std::cout<<"Number of selected Events: "<<NumberOfEventsAfterSelection<<std::endl;
		std::cout<<"***********************************"<<std::endl;
	}

//	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
//	std::cout<<"EVENT LOOP tooks "<<std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()<<" ms."<<std::endl;

//	this->writeHistograms();
//	outputFile_->Close();
}
void JetAnalysisBase::addTriggerObjects(const std::vector<std::string> &triggerObjectName, const std::string & path)
{
	if(TEST) std::cout<<"I'm in JetAnalysisBase::addTriggerObjects"<<std::endl;

	if(triggerObjectName.size() == 0)
	{
		std::cerr<<"Error: Empty vector of triggerObjectNames were sepcified. Interupt!"<<std::endl;
		exit(1);
	}

	triggerObjectName_ = triggerObjectName;
	//Add trees with different Trigger Objects and specify Trigger Logic name
	for(const auto & triggerObject : triggerObjectName)
	{
		this->addTree<TriggerObject>(triggerObject,(path + triggerObject).c_str());
	}
}

void JetAnalysisBase::loadCorrections(){
	if(TEST) std::cout<<"I'm in JetAnalysisBase::loadCorrections"<<std::endl;
	TH1::AddDirectory(0);
//	if(lowM_){
//		// Add Ht reweighting:
//		std::string file_name = outputFile_->GetName();
////		if(file_name.find("BGenFilter") != std::string::npos ){
////			fCorrections_["fHtWeight"] = std::make_unique<TFile>("input_corrections/DataMC_HT_Correction_lowM_BGenFilter.root","read");
////			hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
////		}
//		std::string file_ht_name;
//		if(file_name.find("bEnriched") != std::string::npos || file_name.find("BGenFilter") != std::string::npos){
////			fCorrections_["fHtWeight"] = std::make_unique<TFile>("input_corrections/bbx_HT_Correction_lowM_QCD_b_Ht.root","read");
////			Use the same reweighting function for all topologies
//			if(file_name.find("bbx") == std::string::npos && nJets_ == 2) file_ht_name = "bin/input_corrections/DataMC_HT_Correction_lowM_QCD_not_rew.root";
//			else if (file_name.find("bbx") == std::string::npos && nJets_ == 3) file_ht_name = "bin/input_corrections/DataMC_HT_Correction_lowM_QCD_not_rew.root";
//			if (file_name.find("bbx") != std::string::npos) file_ht_name = "bin/input_corrections/DataMC_HT_Correction_lowM_QCD_not_rew.root";
//			fCorrections_["fHtWeight"] = std::make_unique<TFile>((defaultInputDir_ + file_ht_name).c_str(),"read");
//			hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
//			std::cout<<"HT corrections for bEnriched + BGenFilter samples were loaded"<<std::endl;
//		}
//		else if (file_name.find("TuneCUETP8M1") != std::string::npos) {
//			fCorrections_["fHtWeight"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/DataMC_HT_Correction_lowM_TuneCUETP8M1.root").c_str(),"read");
//			hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
//		}
//
//	}
//	else {
//		// Add Ht reweighting:
//		std::string file_name = outputFile_->GetName();
//		if(file_name.find("BGenFilter") != std::string::npos ){
//			fCorrections_["fHtWeight"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/DataMC_HT_Correction_highM_BGenFilter.root").c_str(),"read");
//			hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
//		}
//		else if (file_name.find("TuneCUETP8M1") != std::string::npos) {
//			fCorrections_["fHtWeight"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/DataMC_HT_Correction_highM_TuneCUETP8M1.root").c_str(),"read");
//			hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
//		}
//	}
	fCorrections_["fPileUpData_central"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/Data_Pileup_2016_ReReco.root").c_str(),"read");
	CheckZombie(*fCorrections_["fPileUpData_central"].get());
	CheckZombieObjectInTFile(*fCorrections_["fPileUpData_central"].get(),"pileup");
	fCorrections_["fPileUpData_up"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/Data_Pileup_2016_ReReco_up2Sigma.root").c_str(),"read");
	CheckZombie(*fCorrections_["fPileUpData_up"].get());
	CheckZombieObjectInTFile(*fCorrections_["fPileUpData_up"].get(),"pileup");
	fCorrections_["fPileUpData_down"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/Data_Pileup_2016_ReReco_down2Sigma.root").c_str(),"read");
	CheckZombie(*fCorrections_["fPileUpData_down"].get());
	CheckZombieObjectInTFile(*fCorrections_["fPileUpData_down"].get(),"pileup");
	hCorrections1D_["hPileUpData_central"] = (TH1D*) fCorrections_["fPileUpData_central"]->Get("pileup");
	hCorrections1D_["hPileUpData_up"] = (TH1D*) fCorrections_["fPileUpData_up"]->Get("pileup");
	hCorrections1D_["hPileUpData_down"] = (TH1D*) fCorrections_["fPileUpData_down"]->Get("pileup");
	fCorrections_["fPileUpMC"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/MC_Autumn16_PU25ns_V1.root").c_str(),"read");
	CheckZombie(*fCorrections_["fPileUpMC"].get());
	CheckZombieObjectInTFile(*fCorrections_["fPileUpMC"].get(),"pileup");
	hCorrections1D_["hPileUpMC"] = (TH1D*) fCorrections_["fPileUpMC"] -> Get("pileup");
	//Add Ht weights
	fCorrections_["fHtWeight"] = std::make_unique<TFile>((defaultInputDir_ + "bin/input_corrections/DataMC_ReReco_HT_Correction_lowM_QCD_HT.root").c_str(),"read");
	hCorrections1D_["hHtWeight"] =  (TH1D*) fCorrections_["fHtWeight"] -> Get("h_clone") ;
	//Declare weights:
	pWeight_ = std::unique_ptr<Weights>(new Weights(lowM_));
	// BTag Reader
	std::vector<std::string> csvv2 = {defaultInputDir_ + "bin/input_corrections/CSVv2_Moriond17_B_F.csv",
									  defaultInputDir_ + "bin/input_corrections/CSVv2_Moriond17_G_H.csv"};//{defaultInputDir_ + "bin/input_corrections/SFbLib.csv"};//
	if(isMC()) BTagLib_.Setup("csvv2",csvv2);

	//TEST weights for bbx
}

const bool JetAnalysisBase::leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets){
//	if(TEST) std::cout<<"I'm in JetAnalysisBase::leadingJetSelection::shared_ptr"<<std::endl;
	//Selection of good Leading Jets:
	//Only jets that pass Loose identification will be considered
	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);

	if(jet1.pt() < pt1_ || jet2.pt() < pt2_) return false;
	if(jet1.eta() < eta1_ || jet2.eta() < eta2_) return false;

	return true;
}

const bool JetAnalysisBase::OnlineSelection(const analysis::tools::Jet &fLeadOfflineJet,
							  const analysis::tools::Jet &sLeadOfflineJet){

	for(const auto & triggerObject : triggerObjectName_){
		const Candidate *onlineJet1 = fLeadOfflineJet.matched(triggerObject);
		const Candidate *onlineJet2 = sLeadOfflineJet.matched(triggerObject);

		if(!onlineJet1 || !onlineJet2 || onlineJet1 == onlineJet2) return false;
		//Check dEta condition
		if(triggerObject.find("MaxDeta1p6") != std::string::npos){
			if(std::abs(onlineJet1->eta() - onlineJet2-> eta()) > 1.6) return false;
		}
	}

	return true;
}

void JetAnalysisBase::makeHistograms(const int & size){
	histo_.Make(size,lowM_);
}

void JetAnalysisBase::SetupStandardOutputFile(const std::string & outputFileName){

		//get the full file name and path from the Tree
		std::string fullName = this->tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
		std::string outputName = baseOutputName_;
		if(outputFileName != "") outputName+="_"+outputFileName;
		if(lowM_) outputName += "_lowM_";
		else {outputName += "_highM_";}

		//TODO: find better solution to get the file name
		size_t found;
		found=fullName.find_last_of("/");
		fullName = fullName.substr(0,found);
		found = fullName.find_last_of("/");
		fullName = fullName.substr(0,found);
		found = fullName.find_last_of("/");
		fullName = fullName.substr(0,found);
		found = fullName.find_last_of("/");

		if(is_mc_ && findStrings(fullName,"asymptotic")){
			fullName = fullName.substr(0,found);
			found = fullName.find_last_of("/");
			fullName = fullName.substr(found+1);
		}
		else{
			fullName = fullName.substr(found+1);
		}

		outputName += fullName;// + ".root";
		this->createOutputFile(outputName);
//	else{
//		this->createOutputFile(outputFileName);
//	}

}

void JetAnalysisBase::createOutputFile(const std::string &name){

	std::string finale_name = name;
	const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));

	if(name.find(".root") == std::string::npos) finale_name += ".root";
	if(finale_name.find("/") != std::string::npos){
		outputFile_ = std::make_unique<TFile>(finale_name.c_str(),"RECREATE");
	}
	else {
		finale_name = cmsswBase + "/src/Analysis/MssmHbb/output/" + finale_name;
//		finale_name = "/nfs/dust/cms/user/" + std::string(getenv("USERNAME")) + "/output/" + finale_name;
		outputFile_ = std::make_unique<TFile>(finale_name.c_str(),"RECREATE");
	}
	std::cout<<"File: "<<finale_name<<" was created"<<std::endl;
}

void JetAnalysisBase::writeHistograms(){

	for(const auto & h : histo_.getHisto()){
		if(h.second->GetEntries() != 0) h.second->Write();
	}
	outputFile_->Close();
}

void JetAnalysisBase::fillHistograms(const std::shared_ptr<Collection<Jet> > &offlineJets, const double & weight){

//	if(TEST) std::cout<<"I'm in JetAnalysisBase::fillHistograms"<<std::endl;

	Jet jet1 = offlineJets->at(0);
	Jet jet2 = offlineJets->at(1);

	(histo_.getHisto())["jet_pt1"]->Fill(jet1.pt(),weight);
	(histo_.getHisto())["jet_pt2"]->Fill(jet2.pt(),weight);

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

	TLorentzVector obj12;
	obj12 = jet1.p4() + jet2.p4();
	(histo_.getHisto())["diJet_pt"]->Fill(obj12.Pt(),weight);
	(histo_.getHisto())["diJet_eta"]->Fill(obj12.Eta(),weight);
	(histo_.getHisto())["diJet_phi"]->Fill(obj12.Phi(),weight);

}

void JetAnalysisBase::combineBTagSFs(const std::array<BTagScaleFactor::ScaleFactor,3>& sf){
//	Method to combine SFs from different jets(virtual)
	for(int i = 0; i<nJets_;++i){
		if(sf[i].flavour != 0){ //only for b/c
			weight_["SFb_central"] 	*= sf[i].central;
			weight_["SFb_up"]			*= (2.*(sf[i].up - sf[i].central) + sf[i].central);
			weight_["SFb_down"]		*= (2.*(sf[i].down - sf[i].central) + sf[i].central);
		}
		else {
			weight_["SFl_central"] 	*= sf[i].central;
			weight_["SFl_up"]			*= (2.*(sf[i].up - sf[i].central) + sf[i].central);
			weight_["SFl_down"]		*= (2.*(sf[i].down - sf[i].central) + sf[i].central);
		}
	}
}

void JetAnalysisBase::combineBTagSFs(const std::array< std::vector<BTagScaleFactor::ScaleFactor> , 3 >& sfs, const std::vector<double>& lumis, const double& lumi_tot){
//	Method to combine SFs from different jets(virtual)

	for(int i = 0; i<nJets_;++i){
		//Loop over the b-tag jets
		double SF_central = 0, eSF_up = 0, eSF_down = 0;
		for(auto j = 0; j != int(sfs[i].size()); ++j){
			SF_central 	+= sfs[i].at(j).central * lumis.at(j) / lumi_tot;
			eSF_up 		+= (sfs[i].at(j).up) * lumis.at(j) / lumi_tot;
			eSF_down	+= (sfs[i].at(j).down) * lumis.at(j) / lumi_tot;
		}
		if(sfs[i].at(0).flavour != 0){ //only for b/c
			weight_["SFb_central"] 	*= SF_central;
			weight_["SFb_up"]		*= (2.*(eSF_up - SF_central) + SF_central);
			weight_["SFb_down"]		*= (2.*(eSF_down - SF_central) + SF_central);
		}
		else {
			weight_["SFl_central"] 	*= SF_central;
			weight_["SFl_up"]		*= (2.*(eSF_up - SF_central) + SF_central);
			weight_["SFl_down"]		*= (2.*(eSF_down - SF_central) + SF_central);
		}
	}

}

void JetAnalysisBase::runAnalysis(const std::string &json, const std::string &output, const int &size){
	this->setupAnalysis(json);
	this->SetupStandardOutputFile(output);
	this->makeHistograms(size);
	this->applySelection();
	this->writeHistograms();
}

const double JetAnalysisBase::assignWeight(){
	double weight = 1;
	if(isMC()) {
		weight = weight_["dEta"] * weight_["Lumi"] * weight_["2DPt"] * weight_["BTag"] * weight_["PU_central"];
	}
	return weight;
}

std::shared_ptr<tools::Collection<tools::Jet> > JetAnalysisBase::modifyJetCollection(tools::Jet & jet,
						std::shared_ptr<tools::Collection<tools::Jet> > & initialJets
						){
//	if(TEST) std::cout<<"I'm in JetAnalysisBase::modifyJetCollection"<<std::endl;

	initialJets->add(jet);
	return initialJets;
}

int JetAnalysisBase::returnMassPoint() const {
	int Mpoint = 0;
	if(!signalMC_){
		return 0;
	}
	std::string MassPos = "_M-";
	auto p1 = inputFilelist_.find(MassPos) + 3;
	if(p1 == std::string::npos) {
		std::cerr<<"FileNames were cahnged!!!!"<<std::endl;
		exit(1);
	}
	auto p2 = inputFilelist_.find("_",p1);
	std::string MpointString = inputFilelist_.substr(p1,size_t(p2-p1));
	Mpoint = std::stoi(MpointString);
	return Mpoint;
}

void JetAnalysisBase::setupXSections(){
    xsection_[100] = 1.;//5.23;	//5.22917222713
    xsection_[120] = 1.;//41.79;	//41.7887314999
    xsection_[160] = 1.;//75.32;	//75.3185132566
    xsection_[200] = 1.;//35.42;	//35.4180244078
    xsection_[250] = 1.;//15.55;	//15.5527782822
    xsection_[300] = 1.;//7.64;	//7.63976194847
    xsection_[350] = 1.;//4.1; 	//4.09571617045
    xsection_[400] = 1.;//2.34;	//2.34044663823
    xsection_[500] = 1.;//0.88;	//0.877385402091
    xsection_[600] = 1.;//0.377;	//0.376766971183
    xsection_[700] = 1.;//0.18;	//0.17940535717
    xsection_[750] = 1.;//0.18;	//0.17940535717
    xsection_[800] = 1.;//0.18;	//0.17940535717
    xsection_[900] = 1.;//0.051;	//0.0507664661841
    xsection_[1100] = 1.;//0.018;//0.0180548176636
    xsection_[1300] = 1.;//0.008;//0.00802737010696
}

const double JetAnalysisBase::mHat(){
	double p_prot = 13000. /2.;
	double p1 = this->pdf().x.first * p_prot;
	double p2 = this->pdf().x.second * p_prot;
	double mHat_ = std::sqrt((p1+p2)*(p1+p2) - (p1-p2)*(p1-p2));
	return mHat_;
}


const void JetAnalysisBase::pdfUncertainties(const double& genScale, const analysis::tools::PDF& partons) {

	const LHAPDF::PDFSet set("PDF4LHC15_nlo_mc");
	const size_t nmem = set.size()-1;
  const std::vector<LHAPDF::PDF*> pdfs = set.mkPDFs();

  std::vector<double> pdf1Vec, pdf2Vec;
  for (size_t imem = 0; imem <= nmem; imem++) {
    pdf1Vec.push_back(pdfs[imem]->xfxQ(partons.id.first, partons.x.first, genScale));
    pdf2Vec.push_back(pdfs[imem]->xfxQ(partons.id.second, partons.x.second, genScale));
  }
	std::cout << pdf1Vec[0] << " - " << pdf2Vec[0] << std::endl;
	double meanValXsec = 0;
	for (size_t imem = 0; imem <= nmem; imem++) {
		meanValXsec += pdf1Vec[imem]*pdf2Vec[imem];
	}
	meanValXsec /= nmem;
	double stdDeviation = 0;
	for (size_t imem = 0; imem <= nmem; imem++) {
		double diff = (pdf1Vec[imem]*pdf2Vec[imem] - meanValXsec);
		stdDeviation += diff*diff;
	}
	stdDeviation = std::sqrt(stdDeviation/(nmem-1));


}
