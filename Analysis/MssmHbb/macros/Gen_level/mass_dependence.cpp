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

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace analysis;
using namespace analysis::tools;
using namespace mssmhbb;

map<string,TH1D*> GetHistograms(const int& mass);
map<string,TGraphErrors*> GetGraphs(const int& npoints);

std::vector<std::string> dEtas = {"dEta_gen_BBar_higgs_daughters","dEta_gen_higgs_daughters","dEta_gen12","dEta_gen13","dEta_gen14","dEta_gen23","dEta_gen24","dEta_gen34","dEta12","dEta13","dEta14","dEta23","dEta24","dEta34","dEta_recoJet_higgs"};
std::vector<std::string> dPhis = {"dPhi_gen_BBar_higgs_daughters","dPhi_gen_higgs_daughters","dPhi_gen12","dPhi_gen13","dPhi_gen14","dPhi_gen23","dPhi_gen24","dPhi_gen34","dPhi12","dPhi13","dPhi14","dPhi23","dPhi24","dPhi34","dPhi_recoJet_higgs"};
std::vector<std::string> dRs = {"dR_gen_BBar_higgs_daughters","dR_gen_higgs_daughters","dR_gen12","dR_gen13","dR_gen14","dR_gen23","dR_gen24","dR_gen34","dR12","dR13","dR14","dR23","dR24","dR34","dR_recoJet_higgs"};
std::vector<std::string> pts = {"pT_true_higgs","pT_genBBar_higgs","pT_genJet_higgs","pT_gen_12","pT_12","pT_34","pT_1","pT_2","pT_3","pT_4","pT_recoJet_higgs"};
std::vector<std::string> Ms = {"M_true_higgs","M_genBBar_higgs","M_genJet_higgs","M_gen_12","M_12","M_recoJet_higgs"};

int main(){

	HbbStyle style;
	style.set(PRELIMINARY);

	map<int,string> ntuples = {
			{300,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.txt"},
			{350,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.txt"},
			{400,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.txt"},
			{500,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.txt"},
			{600,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.txt"},
			{700,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.txt"},
			{900,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.txt"},
			{1100,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.txt"},
			{1300,cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.txt"}
//			{300,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151441/0000/ntuple_1.root"},
//			{350,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170106_103107/0000/ntuple_1.root"},
//			{400,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151506/0000/ntuple_1.root"},
//			{500,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151518/0000/ntuple_1.root"},
//			{600,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151529/0000/ntuple_1.root"},
//			{700,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151538/0000/ntuple_1.root"},
//			{900,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151547/0000/ntuple_1.root"},
//			{1100,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151556/0000/ntuple_1.root"},
//			{1300,"/pnfs/desy.de/cms/tier2/store/user/rwalsh/Analysis/Ntuples/MC/Moriond17/80x_moriond_2017_tranchiv_v1/SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/170105_151604/0000/ntuple_1.root"}
	};

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

		Analysis analysis(nt.second);

		auto histos = GetHistograms(nt.first);

		// Physics Objects Collections
		analysis.addTree<GenParticle> ("GenParticles","MssmHbb/Events/prunedGenParticles");
		analysis.addTree<GenJet> ("GenJets","MssmHbb/Events/slimmedGenJets");
		analysis.addTree<Jet> ("Jets","MssmHbb/Events/slimmedJetsReapplyJEC");

		string mass_string = to_string(nt.first);

		int nevents = analysis.size();
		std::cout<<"Number of events: "<<nevents<<std::endl;
		for (int i = 0 ; i < nevents; ++i){
			analysis.event(i);

			auto genParticles = analysis.collection<GenParticle>("GenParticles");
			GenParticle bquark;
			GenParticle bbarquark;
			GenParticle higgs;
			int n_higgs = 0, n_bquark = 0;
			for ( int j = 0 ; j < genParticles->size() ; ++j )
			{
				GenParticle gp = genParticles->at(j);
				if ( gp.pdgId() == 36 && gp.status() > 60 ) {
					higgs = gp;
					++n_higgs;
				}
		         if (fabs(gp.pdgId()) == 5 && gp.higgsDaughter())
		         {
		            if ( gp.pdgId() > 0 )  bquark = gp;
		            if ( gp.pdgId() < 0 )  bbarquark = gp;
		            ++n_bquark;
		         }

		    }
			if(n_higgs == 0 || n_bquark < 2) continue;
			TLorentzVector bbar_higgs = bquark.p4() + bbarquark.p4();

			// Generated jets
			// Let's match with the Higgs daughters
			auto genJets = analysis.collection<GenJet>("GenJets");
			if(genJets->size()<4) continue;

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

map<string,TH1D*> GetHistograms(const int& mass){

	string mass_string = to_string(mass);
	map<string,TH1D*> histos;
	string name;
	for(const auto& dR : dRs) {
		name = dR + "_" + mass_string;
		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,0,7);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& pt : pts) {
		name = pt + "_" + mass_string;
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,0,mass);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}

	for(const auto& pt : Ms) {
		name = pt + "_" + mass_string;
		histos[name]	= new TH1D( name.c_str(),name.c_str(),100,0,mass*1.5);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& dEta : dEtas) {
		name = dEta + "_" + mass_string;
		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,-7,7);
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}
	for(const auto& dPhi : dPhis) {
		name = dPhi + "_" + mass_string;
		histos[name] 	= new TH1D( name.c_str(),name.c_str(),100,-1.1*TMath::Pi(),1.1*TMath::Pi());
		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}

	return histos;
}

map<string,TGraphErrors*> GetGraphs(const int& npoints){
	map<string,TGraphErrors*> histos;
	//VS mass
	for(const auto& dR : dRs) {
		histos[dR] = new TGraphErrors(npoints);
		histos[dR]->SetTitle(("<" + dR + ">").c_str());
		histos[dR]->GetXaxis()->SetTitleOffset(1.05);
		histos[dR]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
		histos[dR]->GetYaxis()->SetTitle(("<" + dR + ">").c_str());
		histos[dR]->SetMarkerStyle(20);
		histos[dR]->SetMarkerSize(1.2);
	}
	for(const auto& dEta : dEtas) {
		histos[dEta] = new TGraphErrors(npoints);
		histos[dEta]->SetTitle(("<" + dEta + ">").c_str());
		histos[dEta]->GetXaxis()->SetTitleOffset(1.05);
		histos[dEta]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
		histos[dEta]->GetYaxis()->SetTitle(("<" + dEta + ">").c_str());
		histos[dEta]->SetMarkerStyle(20);
		histos[dEta]->SetMarkerSize(1.2);
	}
	for(const auto& dPhi : dPhis) {
		histos[dPhi] = new TGraphErrors(npoints);
		histos[dPhi]->SetTitle(("<" + dPhi + ">").c_str());
		histos[dPhi]->GetXaxis()->SetTitleOffset(1.05);
		histos[dPhi]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
		histos[dPhi]->GetYaxis()->SetTitle(("<" + dPhi + ">").c_str());
		histos[dPhi]->SetMarkerStyle(20);
		histos[dPhi]->SetMarkerSize(1.2);
	}
	for(const auto& pt : pts) {
		histos[pt] = new TGraphErrors(npoints);
		histos[pt]->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
		histos[pt]->GetXaxis()->SetTitleOffset(1.05);
		histos[pt]->GetYaxis()->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->SetMarkerStyle(20);
		histos[pt]->SetMarkerSize(1.2);
	}

	for(const auto& pt : Ms) {
		histos[pt] = new TGraphErrors(npoints);
		histos[pt]->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->GetXaxis()->SetTitleOffset(1.05);
		histos[pt]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
		histos[pt]->GetYaxis()->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->SetMarkerStyle(20);
		histos[pt]->SetMarkerSize(1.2);
	}

	return histos;
}
