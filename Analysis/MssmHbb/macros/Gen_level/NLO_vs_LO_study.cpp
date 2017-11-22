/*
 * mass_dependence.cpp
 *
 *  Created on: 31 May 2017
 *      Author: shevchen
 *
 *      Macro to plot mass dependence of the
 *      gen distributions: dR12; pT_phi; pt_3,4; dR34
 */

#include <vector>
#include <string>

#include "TH1.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace analysis;
using namespace analysis::tools;
using namespace mssmhbb;

std::map<std::string,TH1D*> GetHistograms(const int& mass, const bool& nlo);
void SelectNLOvsLOshapes(const std::vector<int>& mass_points, const std::string& LO_samples_path, const std::string& NLO_samples_path, const std::string& output_folder);
std::map<std::string,TH1D*> PerformGenEventSelection(const int& mass, const std::string& sample_path, const std::string& output_folder, const bool& nlo);
void MakeComparisonPlots(std::map<std::string,TH1D*> & h_NLO,std::map<std::string,TH1D*> & h_LO, const std::string& output_folder);
double GetMaximum(TH1 * h1, TH1 * h2);

template <typename T> struct less_than_pt
{
    inline bool operator() (const T& struct1, const T& struct2)
    {
        return (struct1.pt() < struct2.pt());
    }
};

//std::vector<std::string> dEtas = {"higgs","bquark_Hdaughter","bbarquark_Hdaughter","bquark_0","bquark_1","bquark_2","bquark_3","bquark_4","bquark_5","b_leading_spectator"};
//std::vector<std::string> dPhis = {"dPhi_gen_BBar_higgs_daughters","dPhi_gen_higgs_daughters","dPhi_gen12","dPhi_gen13","dPhi_gen14","dPhi_gen23","dPhi_gen24","dPhi_gen34","dPhi12","dPhi13","dPhi14","dPhi23","dPhi24","dPhi34","dPhi_recoJet_higgs"};
//std::vector<std::string> dRs = {"dR_gen_BBar_higgs_daughters","dR_gen_higgs_daughters","dR_gen12","dR_gen13","dR_gen14","dR_gen23","dR_gen24","dR_gen34","dR12","dR13","dR14","dR23","dR24","dR34","dR_recoJet_higgs"};
std::vector<std::string> pts = {"higgs","bquark_Hdaughter","bbarquark_Hdaughter","bquark_0","bquark_1","bquark_2","bquark_3","bquark_4","bquark_5","b_leading_spectator",
								"genJet_0","genJet_1","genJet_2","genJet_3","genJet_4","genJet_5","genJet_6","genJet_7",
								"recoJet_0","recoJet_1","recoJet_2","recoJet_3","recoJet_4","recoJet_5","recoJet_6","recoJet_7"};
std::vector<std::string> etas = {"higgs","bquark_Hdaughter","bbarquark_Hdaughter","bquark_0","bquark_1","bquark_2","bquark_3","bquark_4","bquark_5","b_leading_spectator",
								"genJet_0","genJet_1","genJet_2","genJet_3","genJet_4","genJet_5","genJet_6","genJet_7",
								"recoJet_0","recoJet_1","recoJet_2","recoJet_3","recoJet_4","recoJet_5","recoJet_6","recoJet_7"};
std::vector<std::string> phis = {"higgs","bquark_Hdaughter","bbarquark_Hdaughter","bquark_0","bquark_1","bquark_2","bquark_3","bquark_4","bquark_5","b_leading_spectator",
								"genJet_0","genJet_1","genJet_2","genJet_3","genJet_4","genJet_5","genJet_6","genJet_7",
								"recoJet_0","recoJet_1","recoJet_2","recoJet_3","recoJet_4","recoJet_5","recoJet_6","recoJet_7"};
std::vector<std::string> flavour = {"recoJet_0","recoJet_1","recoJet_2","recoJet_3","recoJet_4","recoJet_5","recoJet_6","recoJet_7"};
std::vector<std::string> btag = {"recoJet_0","recoJet_1","recoJet_2","recoJet_3","recoJet_4","recoJet_5","recoJet_6","recoJet_7"};
//std::vector<std::string> Ms = {"M_true_higgs","M_genBBar_higgs","M_genJet_higgs","M_gen_12","M_12","M_recoJet_higgs"};

int main(){
	HbbStyle style;
	style.set(PRELIMINARY_SIMULATION);

	std::string output_folder = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/";
	std::string NLO_samples_path = "/afs/desy.de/user/c/chayanit/cms/forRostyslav/21072016/";
	std::string LO_samples_path  = cmsswBase +"/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/";
	std::vector<int> mass_points = {350,600,900};

	SelectNLOvsLOshapes(mass_points,LO_samples_path, NLO_samples_path, output_folder);
}

void SelectNLOvsLOshapes(const std::vector<int>& mass_points, const std::string& LO_samples_path, const std::string& NLO_samples_path, const std::string& output_folder){
	/*
	 * Function to perform comparison of NLO vs
	 */
	for(const auto& mass : mass_points){
		std::string mass_str = std::to_string(int(mass));
		std::string full_LO_sample_path 	= LO_samples_path 	+ "SUSYGluGluToBBHToBB_NarrowWidth_M-" + mass_str + "_TuneCUETP8M1_13TeV-pythia8.txt";
		std::string full_NLO_sample_path 	= NLO_samples_path 	+ "SUSYGluGluToBBHToBB_M-" + mass_str + "_cfg_ntuple.txt";
		auto h_NLO = PerformGenEventSelection(mass,full_NLO_sample_path,output_folder,true);
		auto h_LO = PerformGenEventSelection(mass,full_LO_sample_path,output_folder,false);
		MakeComparisonPlots(h_NLO,h_LO,output_folder);
	}
}

void MakeComparisonPlots(std::map<std::string,TH1D*> & h_NLO,std::map<std::string,TH1D*> & h_LO, const std::string& output_folder){
	/*
	 * Draw comparison
	 */
	std::cout<<"WTFFF: "<<std::endl;
	for(const auto& h : h_NLO){

		auto hLO = h_LO[h.first];
		auto hNLO= h_NLO[h.first];
		hLO->SetLineColor(kRed);
		hNLO->SetLineColor(kBlue);
		hLO->Scale(1./hLO->GetSumOfWeights());
		hNLO->Scale(1./hNLO->GetSumOfWeights());
		double min = hLO->GetXaxis()->GetXmin();
		double max = hLO->GetXaxis()->GetXmax();
		double ymax = 1.1 * GetMaximum(hLO, hNLO);

		TCanvas can;
		TH2D *topFrame = new TH2D("topFrame","",1,min,max,1,0,ymax);
		HbbStyle::setRatioTopFrame(topFrame,can.YtoPixel(can.GetY1()) );
		auto *topPad = HbbStyle::getRatioTopPad();
		topPad->Draw();
		topPad->cd();
		auto *leg = HbbStyle::legend("top,right",3);//new TLegend(0.6,0.6,0.95,0.95);
		leg->SetY1(leg->GetY1()*0.8);
		leg->SetY2(leg->GetY2()*0.95);
		leg->AddEntry(hLO,"Pythia8","l");
		leg->AddEntry(hNLO,"aMC@NLO","l");
//		hLO->SetMaximum(1.1 * GetMaximum(hLO, hNLO));
		topFrame->Draw();
		hLO->Draw("HIST same");
		hNLO->Draw("same HIST");
		leg->Draw();

		can.cd();
		auto *botPad = HbbStyle::getRatioBottomPad();
		botPad->Draw();
		botPad->cd();
		TH2D *frame2 = new TH2D("frame2","",1,min,max,1,0,3);
		HbbStyle::setRatioBottomFrame(frame2->GetXaxis(),frame2->GetYaxis(),can.YtoPixel(can.GetY1()),topPad->YtoPixel(topPad->GetY1()));
		frame2->SetTitle( ("; " + h.first + "; #frac{aMC@NLO}{Pythia8}").c_str());
		frame2->Draw();

		auto ratio = static_cast<TH1D*>(hNLO->Clone("ratio"));
		ratio->SetLineColor(kBlack);
		ratio->Divide(hLO);
		ratio->Draw("HISTO same");
		can.Print( (output_folder + h.first + ".pdf").c_str());
	}
}

std::map<std::string,TH1D*> PerformGenEventSelection(const int& mass, const std::string& sample_path, const std::string& output_folder, const bool& nlo){
	/*
	 * Event selection
	 */
	Analysis analysis(sample_path);
	auto histos = GetHistograms(mass,nlo);

	// Physics Objects Collections
	analysis.addTree<GenParticle> ("GenParticles","MssmHbb/Events/prunedGenParticles");
	analysis.addTree<GenJet> ("GenJets","MssmHbb/Events/slimmedGenJets");
	analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");
	std::string additional_histo_name = "";
//	if(nlo) additional_histo_name

	std::string mass_str = std::to_string(mass);
	//Event loop
	int nevents = analysis.size();
//	nevents = 1000;
	for (int i = 0 ; i < nevents; ++i){
		analysis.event(i);

		double gen_weight = std::abs(analysis.genWeight())/analysis.genWeight();

		auto genParticles = analysis.collection<GenParticle>("GenParticles");
		GenParticle higgs, b_from_higgs, bbar_from_higgs, spectator_b, spectator_bbar, leading_spectator_b;
		std::vector<GenParticle> bquarks;

		int n_higgs = 0, n_bquark = 0;
		double pt_min = -999;
		for ( int j = 0 ; j < genParticles->size() ; ++j )
		{
			GenParticle gp = genParticles->at(j);
			//For the LO
			if(!nlo){
				//Higgs boson
				if ( gp.pdgId() == 36 && gp.status() > 60 ) {
					higgs = gp;
					++n_higgs;
				}
				// b-quarks from Higgs decay
				if (fabs(gp.pdgId()) == 5 && gp.higgsDaughter())
				{
					if ( gp.pdgId() > 0 )  b_from_higgs = gp;
					if ( gp.pdgId() < 0 )  bbar_from_higgs = gp;
					++n_bquark;
				}
				//spectator partons
				if(fabs(gp.pdgId()) == 5 && (j==1 || j ==2)){
					if ( gp.pdgId() > 0 ) spectator_b = gp;
					if ( gp.pdgId() < 0 ) spectator_bbar = gp;
					if(gp.pt() > pt_min) {
						//Selecting leading pt
						leading_spectator_b = gp; pt_min = gp.pt();
					}
				}
			}
			//For the NLO
			else{
				if ( gp.pdgId() == 25 && gp.status() > 60 ) {
					higgs = gp;
					++n_higgs;
				}
				// b-quarks from Higgs decay
				if (fabs(gp.pdgId()) == 5 && (j!=1 && j!=2) && gp.status() == 23){
					if ( gp.pdgId() > 0 )  b_from_higgs = gp;
					if ( gp.pdgId() < 0 )  bbar_from_higgs = gp;
				}
				//spectator partons
				if(fabs(gp.pdgId()) == 5 && (j==1 || j ==2)){
					if ( gp.pdgId() > 0 ) spectator_b = gp;
					if ( gp.pdgId() < 0 ) spectator_bbar = gp;
					if(gp.pt() > pt_min) {
						//Selecting leading pt
						leading_spectator_b = gp; pt_min = gp.pt();
					}
				}
			}
			//All b-quarks
			if(fabs(gp.pdgId()) == 5 && gp.status()<40) bquarks.push_back(gp);
	    }

		//Order b-quarks by pT
//		std::sort(bquarks.begin(),bquarks.end(),less_than_pt<GenParticle>());
		//Fill histograms with gen-particles
		auto n_b_quarks = 6; if(bquarks.size() < 6) n_b_quarks =bquarks.size();
		for(auto i = 0 ; i < n_b_quarks;++i){
			histos["bquark_" + std::to_string(int(i)) + "_" + mass_str + "_pt"]->Fill(bquarks.at(i).pt(),gen_weight);
			histos["bquark_" + std::to_string(int(i)) + "_" + mass_str + "_eta"]->Fill(bquarks.at(i).eta(),gen_weight);
			histos["bquark_" + std::to_string(int(i)) + "_" + mass_str + "_phi"]->Fill(bquarks.at(i).phi(),gen_weight);
		}
		histos["higgs_" + mass_str + "_pt"]	->Fill(higgs.pt(),gen_weight);
		histos["higgs_" + mass_str + "_eta"]	->Fill(higgs.eta(),gen_weight);
		histos["higgs_" + mass_str + "_phi"]	->Fill(higgs.phi(),gen_weight);

		histos["bquark_Hdaughter_" + mass_str + "_pt"]->Fill(b_from_higgs.pt(),gen_weight);
		histos["bquark_Hdaughter_" + mass_str + "_eta"]->Fill(b_from_higgs.eta(),gen_weight);
		histos["bquark_Hdaughter_" + mass_str + "_phi"]->Fill(b_from_higgs.phi(),gen_weight);
		histos["bbarquark_Hdaughter_" + mass_str + "_pt"]->Fill(bbar_from_higgs.pt(),gen_weight);
		histos["bbarquark_Hdaughter_" + mass_str + "_eta"]->Fill(bbar_from_higgs.eta(),gen_weight);
		histos["bbarquark_Hdaughter_" + mass_str + "_phi"]->Fill(bbar_from_higgs.phi(),gen_weight);
		histos["b_leading_spectator_" + mass_str + "_pt"]->Fill(leading_spectator_b.pt(),gen_weight);
		histos["b_leading_spectator_" + mass_str + "_eta"]->Fill(leading_spectator_b.eta(),gen_weight);
		histos["b_leading_spectator_" + mass_str + "_phi"]->Fill(leading_spectator_b.phi(),gen_weight);
		/*
		 * Generated jets
		 */
		auto genJets = analysis.collection<GenJet>("GenJets");
		std::vector<GenJet> genJetVec;
		for ( int j = 0 ; j < genJets->size() ; ++j ){
			GenJet gjet = genJets->at(j);
			genJetVec.push_back(gjet);
		}

		//Fill gen-jet histograms
		auto n_gen_jets = 4; if(genJetVec.size() < 4) n_gen_jets =genJetVec.size();
		if(n_gen_jets > 7) n_gen_jets = 7;
		for(auto i = 0 ; i < n_gen_jets;++i){
			histos["genJet_" + std::to_string(int(i)) + "_" + mass_str + "_pt"]->Fill(genJetVec.at(i).pt(),gen_weight);
			histos["genJet_" + std::to_string(int(i)) + "_" + mass_str + "_eta"]->Fill(genJetVec.at(i).eta(),gen_weight);
			histos["genJet_" + std::to_string(int(i)) + "_" + mass_str + "_phi"]->Fill(genJetVec.at(i).phi(),gen_weight);
//			histos["genJet_" + std::to_string(int(i)) + "_" + mass_str + "_flavour"]->Fill(genJetVec.at(i).flavour(),gen_weight);
		}

		/*
		 * Reco jets
		 */
		auto recoJets = analysis.collection<Jet>("Jets");
		if(recoJets->size() < 1) continue;
		std::vector<Jet> recoJetVec;
		for ( int j = 0 ; j < recoJets->size() ; ++j ){
			Jet jet = recoJets->at(j);
			recoJetVec.push_back(jet);
		}

		//Fill reco-jet histograms
		int n_reco_jets = recoJetVec.size();
		if(n_reco_jets > 7) n_reco_jets = 7;
		for(int i = 0 ; i < n_reco_jets;++i){
			histos["recoJet_" + std::to_string(int(i)) + "_" + mass_str + "_pt"]->Fill(recoJetVec.at(i).pt(),gen_weight);
			histos["recoJet_" + std::to_string(int(i)) + "_" + mass_str + "_eta"]->Fill(recoJetVec.at(i).eta(),gen_weight);
			histos["recoJet_" + std::to_string(int(i)) + "_" + mass_str + "_phi"]->Fill(recoJetVec.at(i).phi(),gen_weight);
			histos["recoJet_" + std::to_string(int(i)) + "_" + mass_str + "_flavour"]->Fill(recoJetVec.at(i).flavour(),gen_weight);
			histos["recoJet_" + std::to_string(int(i)) + "_" + mass_str + "_btag"]->Fill(recoJetVec.at(i).btag(),gen_weight);
		}

	}
	return histos;
}
/*
int main(){

	auto graphs = GetGraphs(ntuples.size());
	map<string,double> dR, dEta, dPhi;
	map<string,double> pT, M;
	for(const auto& dr : dRs) dR[dr] = -100;
	for(const auto& deta : dEtas) dEta[deta] = -100;
	for(const auto& dphi : dPhis) dPhi[dphi] = -100;
	for(const auto& pt : pts) pT[pt] = -100;
	for(const auto& m : Ms) M[m] = -100;
	int counter = 0;

	for(const auto& nt : ntuples){

		string mass_string = to_string(nt.first);
			if(n_higgs == 0 || n_bquark < 2) continue;
			TLorentzVector bbar_higgs = bquark.p4() + bbarquark.p4();

			// Generated jets
			// Let's match with the Higgs daughters
			GenJet b_higgs_genjet;
			GenJet bbar_higgs_genjet;
			int nb = 0;
			int nbbar = 0;
			for ( int j = 0 ; j < genJets->size() ; ++j )
			{
				GenJet gjet = genJets->at(j);
				if ( gjet.deltaR(bquark) < 0.3 && nb == 0)
				{
					b_higgs_genjet = gjet ;
					++nb;
				}
				if ( gjet.deltaR(bbarquark) < 0.3 && nbbar == 0 )
				{
					bbar_higgs_genjet = gjet ;
		            ++nbbar;
				}
				if ( nb > 0 && nbbar > 0 ) break;
			}
			if(nb < 1 || nbbar < 1) continue;

			TLorentzVector bbar_jet_higgs = b_higgs_genjet.p4() + bbar_higgs_genjet.p4();

			GenJet leading_genJet1 = genJets->at(0);
			GenJet leading_genJet2 = genJets->at(1);
			GenJet leading_genJet3 = genJets->at(2);
			GenJet leading_genJet4 = genJets->at(3);

			TLorentzVector gen12 = leading_genJet1.p4() + leading_genJet2.p4();

			auto Jets = analysis.collection<Jet>("Jets");
			if(Jets->size() < 4) continue;
			//Reco jets
			Jet b_higgs_recojet;
			Jet bbar_higgs_recojet;
			nb = 0;
			nbbar = 0;
			for ( int j = 0 ; j < Jets->size() ; ++j )
			{
				Jet jet = Jets->at(j);
				if ( jet.deltaR(bquark) < 0.3 && nb == 0)
				{
					b_higgs_recojet = jet ;
					++nb;
				}
				if ( jet.deltaR(bbarquark) < 0.3 && nbbar == 0 )
				{
					bbar_higgs_recojet = jet ;
		            ++nbbar;
				}
				if ( nb > 0 && nbbar > 0 ) break;
			}
			if(nb < 1 || nbbar < 1) continue;

			Jet leadingJet1 = Jets->at(0);
			Jet leadingJet2 = Jets->at(1);
			Jet leadingJet3 = Jets->at(2);
			Jet leadingJet4 = Jets->at(3);

			TLorentzVector jet12 = leadingJet1.p4() + leadingJet2.p4();
			TLorentzVector jet34 = leadingJet3.p4() + leadingJet4.p4();

			TLorentzVector bbar_recojet_higgs = b_higgs_recojet.p4() + bbar_higgs_recojet.p4();

			//Delta Rs
			dR["dR_gen_BBar_higgs_daughters"] = bquark.deltaR(bbarquark);
			dR["dR_gen_higgs_daughters"] = b_higgs_genjet.deltaR(bbar_higgs_genjet);
			dR["dR_gen12"] = leading_genJet1.deltaR(leading_genJet2);
			dR["dR_gen13"] = leading_genJet1.deltaR(leading_genJet3);
			dR["dR_gen14"] = leading_genJet1.deltaR(leading_genJet4);

			dR["dR_gen23"] = leading_genJet2.deltaR(leading_genJet3);
			dR["dR_gen24"] = leading_genJet2.deltaR(leading_genJet4);

			dR["dR_gen34"] = leading_genJet3.deltaR(leading_genJet4);

			dR["dR_recoJet_higgs"] = b_higgs_recojet.deltaR(bbar_higgs_recojet);
			dR["dR12"] = leadingJet1.deltaR(leadingJet2);
			dR["dR13"] = leadingJet1.deltaR(leadingJet3);
			dR["dR14"] = leadingJet1.deltaR(leadingJet4);

			dR["dR23"] = leadingJet2.deltaR(leadingJet3);
			dR["dR24"] = leadingJet2.deltaR(leadingJet4);

			dR["dR34"] = leadingJet3.deltaR(leadingJet4);

			//Delta Etas
			dEta["dEta_gen_BBar_higgs_daughters"] = bquark.eta() - bbarquark.eta();
			dEta["dEta_gen_higgs_daughters"] = b_higgs_genjet.eta() - bbar_higgs_genjet.eta();
			dEta["dEta_gen12"] = leading_genJet1.eta() - leading_genJet2.eta();
			dEta["dEta_gen13"] = leading_genJet1.eta() - leading_genJet3.eta();
			dEta["dEta_gen14"] = leading_genJet1.eta() - leading_genJet4.eta();

			dEta["dEta_gen23"] = leading_genJet2.eta() - leading_genJet3.eta();
			dEta["dEta_gen24"] = leading_genJet2.eta() - leading_genJet4.eta();

			dEta["dEta_gen34"] = leading_genJet3.eta() - leading_genJet4.eta();

			dEta["dEta_recoJet_higgs"] = b_higgs_recojet.eta() - bbar_higgs_recojet.eta();
			dEta["dEta12"] = leadingJet1.eta() - leadingJet2.eta();
			dEta["dEta13"] = leadingJet1.eta() - leadingJet3.eta();
			dEta["dEta14"] = leadingJet1.eta() - leadingJet4.eta();
			dEta["dEta23"] = leadingJet2.eta() - leadingJet3.eta();
			dEta["dEta24"] = leadingJet2.eta() - leadingJet4.eta();

			dEta["dEta34"] = leadingJet3.eta() - leadingJet4.eta();

			//Delta Phis
			dPhi["dPhi_gen_BBar_higgs_daughters"] = TVector2::Phi_mpi_pi(bquark.phi() - bbarquark.phi());
			dPhi["dPhi_gen_higgs_daughters"] = TVector2::Phi_mpi_pi(b_higgs_genjet.phi() - bbar_higgs_genjet.phi());
			dPhi["dPhi_gen12"] = TVector2::Phi_mpi_pi(leading_genJet1.phi() - leading_genJet2.phi());
			dPhi["dPhi_gen13"] = TVector2::Phi_mpi_pi(leading_genJet1.phi() - leading_genJet3.phi());
			dPhi["dPhi_gen14"] = TVector2::Phi_mpi_pi(leading_genJet1.phi() - leading_genJet4.phi());

			dPhi["dPhi_gen23"] = TVector2::Phi_mpi_pi(leading_genJet2.phi() - leading_genJet3.phi());
			dPhi["dPhi_gen24"] = TVector2::Phi_mpi_pi(leading_genJet2.phi() - leading_genJet4.phi());

			dPhi["dPhi_gen34"] = TVector2::Phi_mpi_pi(leading_genJet3.phi() - leading_genJet4.phi());

			dPhi["dPhi_recoJet_higgs"] = TVector2::Phi_mpi_pi(b_higgs_recojet.phi() - bbar_higgs_recojet.phi());
			dPhi["dPhi12"] = TVector2::Phi_mpi_pi(leadingJet1.phi() - leadingJet2.phi());
			dPhi["dPhi13"] = TVector2::Phi_mpi_pi(leadingJet1.phi() - leadingJet3.phi());
			dPhi["dPhi14"] = TVector2::Phi_mpi_pi(leadingJet1.phi() - leadingJet4.phi());
			dPhi["dPhi23"] = TVector2::Phi_mpi_pi(leadingJet2.phi() - leadingJet3.phi());
			dPhi["dPhi24"] = TVector2::Phi_mpi_pi(leadingJet2.phi() - leadingJet4.phi());

			dPhi["dPhi34"] = TVector2::Phi_mpi_pi(leadingJet3.phi() - leadingJet4.phi());

			//PTs
			pT["pT_true_higgs"] = higgs.pt();
			pT["pT_genBBar_higgs"] = bbar_higgs.Pt();
			pT["pT_genJet_higgs"] =  bbar_jet_higgs.Pt();
			pT["pT_gen_12"] = gen12.Pt();
			pT["pT_12"] = jet12.Pt();
			pT["pT_34"] = jet34.Pt();
			pT["pT_1"] = leadingJet1.pt();
			pT["pT_2"] = leadingJet2.pt();
			pT["pT_3"] = leadingJet3.pt();
			pT["pT_4"] = leadingJet4.pt();
			pT["pT_recoJet_higgs"] =  bbar_recojet_higgs.Pt();

			//Masses
			M["M_true_higgs"] = higgs.m();
			M["M_genBBar_higgs"] = bbar_higgs.M();
			M["M_genJet_higgs"] = bbar_jet_higgs.M();
			M["M_gen_12"] = gen12.M();
			M["M_12"] = jet12.M();
			M["M_recoJet_higgs"] =  bbar_recojet_higgs.M();

			for(const auto &dr : dR) 	histos[dr.first + "_" + mass_string]->Fill(dr.second);
			for(const auto &pt : pT) 	histos[pt.first + "_" + mass_string]->Fill(pt.second);
			for(const auto &m : M) 		histos[m.first + "_" + mass_string]->Fill(m.second);
			for(const auto &deta : dEta)histos[deta.first + "_" + mass_string]->Fill(deta.second);
			for(const auto &dphi : dPhi)histos[dphi.first + "_" + mass_string]->Fill(dphi.second);
			
			//std::cout<<" WTF "<<pT["pT_true_higgs"]<<std::endl;
		}

		for(const auto &dr : dR){
			TCanvas c("c","c",600,600);
			string name = dr.first + "_" + mass_string;
			histos[name]->Draw("E");
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + name + ".pdf").c_str());
			graphs[dr.first] -> SetPoint(counter, nt.first, histos[name]->GetMean());
			graphs[dr.first] ->SetPointError(counter,0,histos[name]->GetMeanError());
		}

		for(const auto &pt : pT){
			TCanvas c("c","c",600,600);
			string name = pt.first + "_" + mass_string;
			histos[name]->Draw("E");
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + name + ".pdf").c_str());
			graphs[pt.first] -> SetPoint(counter, nt.first, histos[name]->GetMean());
			graphs[pt.first] ->SetPointError(counter,0,histos[name]->GetMeanError());
		}

		for(const auto &pt : M){
			TCanvas c("c","c",600,600);
			string name = pt.first + "_" + mass_string;
			histos[name]->Draw("E");
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + name + ".pdf").c_str());
			graphs[pt.first] -> SetPoint(counter, nt.first, histos[name]->GetMean());
			graphs[pt.first] ->SetPointError(counter,0,histos[name]->GetMeanError());
		}

		for(const auto &v : dEta){
			TCanvas c("c","c",600,600);
			string name = v.first + "_" + mass_string;
			histos[name]->Draw("E");
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + name + ".pdf").c_str());
			graphs[v.first] -> SetPoint(counter, nt.first, histos[name]->GetMean());
			graphs[v.first] ->SetPointError(counter,0,histos[name]->GetMeanError());
		}

		for(const auto &v : dPhi){
			TCanvas c("c","c",600,600);
			string name = v.first + "_" + mass_string;
			histos[name]->Draw("E");
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + name + ".pdf").c_str());
			graphs[v.first] -> SetPoint(counter, nt.first, histos[name]->GetMean());
			graphs[v.first] ->SetPointError(counter,0,histos[name]->GetMeanError());
		}
		++counter;

//		auto *gen_jets = static_cast<TTree*>(f.Get("MssmHbb/Events/slimmedGenJets"));
//		gen_jets->Draw(" sqrt( pow((eta[0] - eta[1]),2) + pow((phi[0] - phi[1]),2) ) ");
	}
	for(const auto &dr : dR) {
		TCanvas c2("c2","c2",600,600);
		graphs[dr.first] -> Draw("AP");
		c2.Print((cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + dr.first + ".pdf").c_str());
	}

	for(const auto &dr : pT) {
		TCanvas c2("c2","c2",600,600);
		graphs[dr.first] -> Draw("AP");
		c2.Print((cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + dr.first + ".pdf").c_str());
	}

	for(const auto &dr : M) {
		TCanvas c2("c2","c2",600,600);
		graphs[dr.first] -> Draw("AP");
		c2.Print((cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + dr.first + ".pdf").c_str());
	}

	for(const auto &dr : dEta) {
		TCanvas c2("c2","c2",600,600);
		graphs[dr.first] -> Draw("AP");
		c2.Print((cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + dr.first + ".pdf").c_str());
	}

	for(const auto &dr : dPhi) {
		TCanvas c2("c2","c2",600,600);
		graphs[dr.first] -> Draw("AP");
		c2.Print((cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/" + dr.first + ".pdf").c_str());
	}

	return 0;
}
*/

std::map<std::string,TH1D*> GetHistograms(const int& mass, const bool& nlo){

	std::string mass_string = std::to_string(mass);
	std::map<std::string,TH1D*> histos;
	std::string name;
//	for(const auto& dR : dRs) {
//		name = dR + "_" + mass_string;
//		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,0,7);
//		histos[name]->SetStats(kFALSE);
//		histos[name]->GetXaxis()->SetTitle(name.c_str());
//		histos[name]->SetMarkerStyle(20);
//		histos[name]->SetMarkerSize(1.2);
//	}
	for(const auto& pt : pts) {
		name = pt + "_" + mass_string + "_pt";
//		if(nlo) name += "_NLO";
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,0,mass);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
		std::cout<<"Histo name: "<<name<<std::endl;
	}
	for(const auto& v : etas) {
		name = v + "_" + mass_string + "_eta";
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,-7,7);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& v : phis) {
		name = v + "_" + mass_string + "_phi";
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,-1.1*TMath::Pi(),1.1*TMath::Pi());
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& v : flavour) {
		name = v + "_" + mass_string + "_flavour";
		histos[name]	= new TH1D( name.c_str(),name.c_str(),12,-6,6);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& v : btag) {
		name = v + "_" + mass_string + "_btag";
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,0,1);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}

//	for(const auto& pt : Ms) {
//		name = pt + "_" + mass_string;
//		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,0,mass*1.5);
//		histos[name]->SetStats(kFALSE);
//		histos[name]->GetXaxis()->SetTitle(name.c_str());
//		histos[name]->SetMarkerStyle(20);
//		histos[name]->SetMarkerSize(1.2);
//	}
//	for(const auto& dEta : dEtas) {
//		name = dEta + "_" + mass_string;
//		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,-7,7);
//		histos[name]->SetStats(kFALSE);
//		histos[name]->GetXaxis()->SetTitle(name.c_str());
//		histos[name]->SetMarkerStyle(20);
//		histos[name]->SetMarkerSize(1.2);
//	}
//	for(const auto& dPhi : dPhis) {
//		name = dPhi + "_" + mass_string;
//		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,-1.1*TMath::Pi(),1.1*TMath::Pi());
//		histos[name]->SetStats(kFALSE);
//		histos[name]->GetXaxis()->SetTitle(name.c_str());
//		histos[name]->SetMarkerStyle(20);
//		histos[name]->SetMarkerSize(1.2);
//	}

	return histos;
}

double GetMaximum(TH1 * h1, TH1 * h2){
	if(h1->GetMaximum() > h2->GetMaximum()) return h1->GetMaximum();
	else return h2->GetMaximum();
}

