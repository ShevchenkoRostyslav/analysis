/**\class MssmHbb MssmHbb.cc Analysis/Tools/src/MssmHbb.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <iostream>
// 
// user include files
#include "Analysis/MssmHbb/interface/MssmHbb.h"

//Delcare functions

//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;

//
// constructors and destructor
//

MssmHbb::MssmHbb(const std::string & inputFilelist, const std::string & evtinfo) : Analysis(inputFilelist,evtinfo), BasicTree(this->isMC())
{
	if(this->isMC()){
		// Add MC information:
		this->crossSections("MssmHbb/Metadata/CrossSections");
		this->generatorFilter("MssmHbb/Metadata/GeneratorFilter");
		this->ShowMCInformation();
	}

	// Tree for Jets
	// For compatibility between 7.4.X and 7.6.X
	if(CMSSW.find("7_4_") != std::string::npos){
		this->addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsPuppi");
	}
	else {
		this->addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");
	}
	this->triggerResults("MssmHbb/Events/TriggerResults");
	// Tree for Vertices
	this->addTree<Vertex> ("Vertices","MssmHbb/Events/offlineSlimmedPrimaryVertices");

        if(findStrings(inputFilelist,"susy")) signalMC_ = true;
        else signalMC_ = false;


}

MssmHbb::~MssmHbb()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

void MssmHbb::SetupStandardOutputFile(const std::string & outputFileName){
	//get the full file name and path from the Tree
	std::string outputFileName__ = this->tree<Jet>("Jets")->PhysicsObjectTree<Jet>::PhysicsObjectTreeBase::TreeBase::tree()->GetFile()->GetName();
	std::string outputName = outputFileName;
	if(lowMSelection_) outputName = outputFileName + "_lowMTrigger_";
	else if(!lowMSelection_) outputName = outputFileName + "_highMTrigger_";
	BasicTree::SetupStandardOutputFile(outputName,outputFileName__);
//	this->setBranches();

}

void MssmHbb::addTriggerObjects(const std::vector<std::string> &triggerObjectName, const std::string & path)
{

	triggerObjectName_ = triggerObjectName;
	if(triggerObjectName_.size() == 0)
	{
		std::cerr<<"Error: Empty vector of triggerObjectNames were sepcified. Interupt!"<<std::endl;
		exit(0);
	}
	lowMSelection_ = false;

	//Add trees with different Trigger Objects and specify Trigger Logic name
	for(const auto & triggerObject : triggerObjectName_)
	{
		this->addTree<TriggerObject>(triggerObject,(path + triggerObject).c_str());

		if(triggerObject.find("MaxDeta1p6") != std::string::npos) {
			triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV0p9_DoublePFJetsC100MaxDeta1p6_v";
			lowMSelection_ = true;
		}

	}

	if(!lowMSelection_){
		triggerLogicName_ = "HLT_DoubleJetsC100_DoubleBTagCSV0p85_DoublePFJetsC160_v";
	}
	this->SetupConstants();

}

void MssmHbb::SetupConstants(){
	if(lowMSelection_){
		pt1_cut__ = 100;
		pt2_cut__ = 100;
		pt3_cut__ = 40;

		deta_cut__ = 1.55;
		btag1_cut__ = 0.80; //CSVv2T
		btag2_cut__ = 0.80; //CSVv2T
		btag3_cut__ = 0.46; 	 //CSVv2M
	}
	else{
		pt1_cut__ = 160;
		pt2_cut__ = 160;
		pt3_cut__ = 40;

		deta_cut__ = 100;
		btag1_cut__ = 0.80;
		btag2_cut__ = 0.80;
		btag3_cut__ = 0.46;
	}
}

void MssmHbb::setBranches(){
	//set special branches for MssmHbb analysis
	//Vertices
	OutTree_->Branch("XPrimaryVTX",XPrimaryVTX_,"XPrimaryVTX[50]/D");
	OutTree_->Branch("YPrimaryVTX",YPrimaryVTX_,"YPrimaryVTX[50]/D");
	OutTree_->Branch("ZPrimaryVTX",ZPrimaryVTX_,"ZPrimaryVTX[50]/D");
	OutTree_->Branch("NPrimaryVTX",&NPrimaryVTX_,"NPrimaryVTX/I");


	if(this->isMC()){
		//Reweighting variables
		OutTree_->Branch("BTagWeight",&BTagWeight_,"BTagWeight/D");
		OutTree_->Branch("WeightHt",&WeightHt_,"WeightHt/D");

		//BTag SF variables:
		OutTree_->Branch("BTagSFcentral",btagSFcentral_,"BTagSFcentral[20]/D");
		OutTree_->Branch("BTagSFup",btagSFup_,"BTagSFup[20]/D");
		OutTree_->Branch("BTagSFdown",btagSFdown_,"BTagSFdown[20]/D");

		//PileUp
		OutTree_->Branch("WeightPileUp",WeightPileUp_,"WeightPileUp[3]/D");

		//Flavour composition
		OutTree_->Branch("cc",&cc_,"cc/I");
		OutTree_->Branch("bb",&bb_,"bb/I");
		OutTree_->Branch("qq",&qq_,"qq/I");
		OutTree_->Branch("bc",&bc_,"bc/I");
		OutTree_->Branch("bq",&bq_,"bq/I");
		OutTree_->Branch("qc",&qc_,"qc/I");

                if(this->isSignalMc()){
                        OutTree_->Branch("mHat",&mHat_,"mHat/D");
                }
	}

	//Set basic branches
	BasicTree::setBranches();
}

void MssmHbb::cleanVariables(){

	//Clean Basic variables
	BasicTree::cleanVariables();

	//Clean Specific for MssmHbb analysis variables
	NPrimaryVTX_ = 0;
	std::fill_n(XPrimaryVTX_,50,-100);
	std::fill_n(YPrimaryVTX_,50,-100);
	std::fill_n(ZPrimaryVTX_,50,-100);

	if(isMC()){
		std::fill_n(btagSFcentral_,20,-100.);
		std::fill_n(btagSFdown_,20,-100);
		std::fill_n(btagSFup_,20,-100);

		std::fill_n(WeightPileUp_,3,-100);

		BTagWeight_ = -100;
		qq_ = 0;
		bb_ = 0;
		cc_ = 0;
		bq_ = 0;
		bc_ = 0;
		qc_ = 0;

                mHat_ = -100;
//                if(isSignalMc()) {
//                double p_prot = 13000. /2.;
//                double p1 = this->pdf().x.first * p_prot;
//                double p2 = this->pdf().x.second * p_prot;
//                mHat_ = std::sqrt((p1+p2)*(p1+p2) - (p1-p2)*(p1-p2));
//                }
	}
}

//For double BTag Selection
bool MssmHbb::OnlineSelection(const analysis::tools::Jet &fLeadOfflineJet,
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

void MssmHbb::ShowMCInformation(){
	//Show MC information
	this->listCrossSections();
	this->listGeneratorFilter();
}

//Methods for calculation
//Calculate BTag SF
void MssmHbb::calculateBTagSF(const BTagCalibrationReader & reader, const BTagCalibrationReader &reader_up, const BTagCalibrationReader &reader_down){

	double btagJetPt = 0;
	btagJetPt = LeadPt_[jetCounter_];
	bool DoubleUncertainty = false;
	if(btagJetPt >= maxBJetPt__){
		btagJetPt = maxBJetPt__ - 1; //Doesn't work at the limit, so that I add "-1"
		DoubleUncertainty = true;
	}

	btagSFcentral_[jetCounter_] = reader.eval(BTagEntry::FLAV_B, LeadEta_[jetCounter_], btagJetPt);
	btagSFup_[jetCounter_] = reader_up.eval(BTagEntry::FLAV_B, LeadEta_[jetCounter_], btagJetPt);
	btagSFdown_[jetCounter_] = reader_down.eval(BTagEntry::FLAV_B, LeadEta_[jetCounter_], btagJetPt);

	// If Pt > Maximum - uncertinties should be doubled and calculated with the extremum value
	if(DoubleUncertainty){
		btagSFup_[jetCounter_] = 2*(btagSFup_[jetCounter_] - btagSFcentral_[jetCounter_]) + btagSFcentral_[jetCounter_];
		btagSFdown_[jetCounter_] = 2*(btagSFdown_[jetCounter_] - btagSFcentral_[jetCounter_]) + btagSFcentral_[jetCounter_];
	}

	//For 2 Sigma variation:
	btagSFup_[jetCounter_] = 2*(btagSFup_[jetCounter_] - btagSFcentral_[jetCounter_]) + btagSFcentral_[jetCounter_];
	btagSFdown_[jetCounter_] = 2*(btagSFdown_[jetCounter_] - btagSFcentral_[jetCounter_]) + btagSFcentral_[jetCounter_];

}

//Compute Flavour composition
void MssmHbb::calculateFlavourComposition(){

	if (LeadJet_[0].flavour("Hadron") == 0 && LeadJet_[1].flavour("Hadron") == 0	)   qq_ = 1; //qq
	else if (LeadJet_[0].flavour("Hadron") == 4 && LeadJet_[1].flavour("Hadron") == 4 ) cc_ = 1; //cc
	else if (LeadJet_[0].flavour("Hadron") == 5 && LeadJet_[1].flavour("Hadron") == 5 ) bb_ = 1;
	else if ((LeadJet_[0].flavour("Hadron") == 0 || LeadJet_[1].flavour("Hadron") == 0) && ( LeadJet_[0].flavour("Hadron") == 4 || LeadJet_[1].flavour("Hadron") == 4 )) qc_ = 1;
	else if ((LeadJet_[0].flavour("Hadron") == 5 || LeadJet_[1].flavour("Hadron") == 5) && ( LeadJet_[0].flavour("Hadron") == 0 || LeadJet_[1].flavour("Hadron") == 0 )) bq_ = 1;
	else if ((LeadJet_[0].flavour("Hadron") == 5 || LeadJet_[1].flavour("Hadron") == 5) && ( LeadJet_[0].flavour("Hadron") == 4 || LeadJet_[1].flavour("Hadron") == 4 )) bc_ = 1;
	else {
		std::cerr<<"Error: Undefined color. Interrupt!"<<std::endl;
		exit(5);
	}
}


void MssmHbb::calculateWeights(TH1F *btag_low_eta,TH1F * btag_mid_eta, TH1F * btag_high_eta, TH2F *pt, const double & dataLumi){
	//TODO: Find better solution
	//Pt trigger efficiency weight based on factorisation apporach
	FactorizationPtWeight_ = factorizationPtWeight1D(LeadJet_[0].pt()) * factorizationPtWeight1D(LeadJet_[1].pt());
	//dEta trigger efficiency weight
	dEtaWeight_ = dEtaWeight(dEtaFS_);
	//Online BTag trigger efficiency
	BTagWeight_ = BTagWeight(btag_low_eta,btag_mid_eta,btag_high_eta,LeadJet_[0].pt(),LeadJet_[0].eta()) *
				  BTagWeight(btag_low_eta,btag_mid_eta,btag_high_eta,LeadJet_[1].pt(),LeadJet_[1].eta());
	//Pt trigger efficiency weight based on 2D distribution
	TwoDPtWeight_ = twoDPtWeight(pt,LeadJet_[0].pt(),LeadJet_[1].pt());
	lumiWeight_ = dataLumi / this->luminosity();

}

void MssmHbb::addVertexInfo(analysis::tools::Vertex & vert){

		XPrimaryVTX_[NPrimaryVTX_] = vert.x();
		YPrimaryVTX_[NPrimaryVTX_] = vert.y();
		ZPrimaryVTX_[NPrimaryVTX_] = vert.z();
		NPrimaryVTX_ ++;
}

double MssmHbb::factorizationPtWeight1D(const double &pt)
{
	if(lowMSelection_){
		return 1./(1.+exp( -8.87968e-02 * (pt - 1.03879e+02) )) * 1./ (1.+exp( -1.88620e-01 * (pt - 1.02891e+02) ));
	}
	else {
		return 1./(1.+exp( -1.60938e-01 * (pt - 1.64206e+02) )) * 1./ (1.+exp( -1.53883e-01 * (pt + 3.50936e+01) ));
	}
}
double MssmHbb::dEtaWeight(const double & dEta)
{
	if(lowMSelection_){
		return 1./(1.+exp( 1.41132e+02 * (dEta - 1.60010e+00) )) * 1./ (1.+exp( -2.06352e-02 * (dEta + 2.27738e+02) ));
	}
	else return 1.;
}
double MssmHbb::BTagWeight(TH1F*region1,TH1F* region2,TH1F* region3, const double &pt, const double &eta){
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

double MssmHbb::twoDPtWeight(TH2F *histo, const double &pt1, const double &pt2){

	if(histo == 0) std::cerr<<"ERROR"<<std::endl;

	if (pt1 > 500) {
		if (pt2 > 500){
			return histo->Interpolate(499,499);
		}
		else return histo->Interpolate(pt2,499);
	}
	else return histo -> Interpolate(pt2,pt1);
}
//
// member functions
//
// ------------ method called for each event  ------------

int MssmHbb::returnMassPoint() const {
        int Mpoint = 0;
        if(!isSignalMc()){
                return 0;
        }
        std::string MassPos = "_M-";
        auto p1 = inputFilelist_.find(MassPos) + 3;
        if(p1 == std::string::npos) {
                std::cerr<<"FileNames were changed!!!!"<<std::endl;
                exit(1);
        }
        auto p2 = inputFilelist_.find("_",p1);
        std::string MpointString = inputFilelist_.substr(p1,size_t(p2-p1));
        Mpoint = std::stoi(MpointString);
        return Mpoint;
}

