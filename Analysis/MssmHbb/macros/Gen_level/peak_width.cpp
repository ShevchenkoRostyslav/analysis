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
#include "TF1.h"

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace analysis;
using namespace analysis::tools;
using namespace mssmhbb;

map<string,TH1D*> GetHistograms(const int& mass);
map<string,TGraphErrors*> GetGraphs(const int& npoints);
double RelBreitWigner(double *x, double *par);
double mybw(Double_t* x, Double_t* par);

std::vector<std::string> Ms = {"M_higgs"};

int main(){
	gStyle->SetOptFit(111111);
	HbbStyle style;
	style.set(PRELIMINARY);
	gStyle->SetOptStat(2000);

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
	};

	auto graphs = GetGraphs(ntuples.size());
	map<string,double> M;
	for(const auto& m : Ms) M[m] = -100;
	int counter = 0;

	for(const auto& nt : ntuples){

		Analysis analysis(nt.second);

		auto histos = GetHistograms(nt.first);

		// Physics Objects Collections
		analysis.addTree<GenParticle> ("GenParticles","MssmHbb/Events/prunedGenParticles");

		string mass_string = to_string(nt.first);

		int nevents = analysis.size();
		std::cout<<"Number of events: "<<nevents<<std::endl;
		for (int i = 0 ; i < nevents; ++i){
			analysis.event(i);

			auto genParticles = analysis.collection<GenParticle>("GenParticles");
			GenParticle higgs;
			int n_higgs = 0;
			for ( int j = 0 ; j < genParticles->size() ; ++j )
			{
				GenParticle gp = genParticles->at(j);
				if ( gp.pdgId() == 36 && gp.status() > 60 ) {
					higgs = gp;
					++n_higgs;
				}
				if(n_higgs != 0) break;
		    }
			if(n_higgs == 0) continue;

			//Masses
			histos["M_higgs_" + mass_string]->Fill(higgs.m());

		}

		for(const auto &pt : M){
			TCanvas c("c","c",600,600);
			string name = pt.first + "_" + mass_string;
			histos[name]->Draw("E");
			//Fitting
			double xmin = histos[name]->GetXaxis()->GetXmin();
			double xmax = histos[name]->GetXaxis()->GetXmax();
			xmin = 0.6*nt.first;
			xmax = 1.4*nt.first;
//			TF1 fit("fit",RelBreitWigner,xmin,xmax,4);
			TF1 fit("fit",RelBreitWigner,xmin,xmax,3);
			fit.SetParameters(1e06,40.,nt.first);
			fit.SetParNames("const", "#Gamma", "M");
			histos[name]->Fit(&fit,"QR");
			std::cout<<"RMS = "<<histos[name]->GetRMS()<<", G = "<<fit.GetParameter(1)<<std::endl;

			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/Width_study_" + name + ".pdf").c_str());
			gPad->SetLogy();
			c.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/GenLevel/Width_study_" + name + "_log.pdf").c_str());
			graphs[pt.first] -> SetPoint(counter, nt.first, fit.GetParameter(1));
			graphs[pt.first] ->SetPointError(counter,0,fit.GetParError(1));
		}
		++counter;
	}

	for(const auto &dr : M) {
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

	for(const auto& pt : Ms) {
		name = pt + "_" + mass_string;
		histos[name]	= new TH1D( name.c_str(),name.c_str(),2000,0,mass*2.5);
//		histos[name]->SetStats(kFALSE);
		histos[name]->GetXaxis()->SetTitle(name.c_str());
		histos[name]->SetMarkerStyle(20);
		histos[name]->SetMarkerSize(1.2);
	}

	return histos;
}

map<string,TGraphErrors*> GetGraphs(const int& npoints){
	map<string,TGraphErrors*> histos;
	//VS mass

	for(const auto& pt : Ms) {
		histos[pt] = new TGraphErrors(npoints);
//		histos[pt]->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->GetXaxis()->SetTitleOffset(1.05);
		histos[pt]->GetXaxis()->SetTitle("m_{A/H} (GeV)");
//		histos[pt]->GetYaxis()->SetTitle(("<" + pt + ">").c_str());
		histos[pt]->SetMarkerStyle(20);
		histos[pt]->SetMarkerSize(1.2);
	}

	return histos;
}

double RelBreitWigner(double *x, double *par){

//	double pol = par[0]*(par[1]*x[0]*x[0] + par[2]*x[0] + 1);
//	double norm = par[0]*(par[1]*par[3]*par[3] + par[2]*par[3] + 1);

	double G = par[1];
	double G2 = G*G;
	double norm = par[0];
	double M = par[2];
	double M2 = M*M;
	double gamma = std::sqrt(M2 * (M2 + G2));
	double k = 2*sqrt(2)*M*G*gamma/(TMath::Pi()*sqrt(M2+gamma));
	double bw = k / (pow((x[0]*x[0] - M2),2)+M2*G2);
	return bw*norm;

}

//Breit-Wigner function
double mybw(Double_t* x, Double_t* par)
{

  return par[0]*TMath::BreitWigner(x[0],par[2],par[1]);
}
