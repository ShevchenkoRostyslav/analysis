#include <memory>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TGraph.h"
#include "../Drawer/HbbStyle.cc"

double xsection(const int & point);
void efficiency(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda, const std::vector<double>& lumis);
void getSignalSamples(const std::string& sample_name_prefix, std::map<int,TFile* > & file_map);

const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));

HbbStyle style;

void preapproval_check_signalEfficiency()
{
   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();

   std::map<int,TFile* > inF_nominal;
   std::map<int,TFile* > inF_nobtag;
   std::map<int,TFile* > inF_nobtag_nodeta;
   std::string samples_nominal		 = "";//"ReReco_PU_test_subrange_presc";//
   std::string samples_nobtag		 = "_efficiency_test_nobtag_ptweight";
   std::string samples_nobtag_nodeta	 = "_efficiency_nobtag_nodeta";//"_signalEff_test";//

   getSignalSamples(samples_nominal,inF_nominal);
   getSignalSamples(samples_nobtag,inF_nobtag);
   getSignalSamples(samples_nobtag_nodeta,inF_nobtag_nodeta);

   style.set(PRELIMINARY);
   style.drawStandardTitle();
   gStyle->SetTitleXOffset(1.);

   std::vector<std::map<int,TFile* > > datasets {inF_nominal,inF_nobtag,inF_nobtag_nodeta};
   std::vector<std::string> legends {"default selection","no b-tag selection","no b-tag, no #Delta#eta"};
   std::vector<double> lumis {35673.,35673.,35673.};

   //Signal Efficiency calculation
   efficiency(datasets,legends,lumis);

}

double xsection(const int &point){
	switch (point){
		case 100:
			return 5.23;
		case 120:
			return 41.79;
		case 160:
			return 75.32;
		case 200:
			return 35.42;
		case 250:
			return 15.55;
		case 300:
			return 7.64;
		case 350:
			return 4.1;
		case 400:
			return 2.34;
		case 500:
			return 0.88;
		case 600:
			return 0.377;
		case 700:
			return 0.18;
		case 900:
			return 0.051;
		case 1100:
			return 0.018;
		case 1300:
			return 0.008;
		default:
			std::cerr<<"No such mass point in the list of xsections"<<std::endl;
			exit(-1);
	}
}

void efficiency(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda, const std::vector<double>& lumis){
	std::vector<short int> colours = {2,4,1,3,4,6,9,12,28,46};
	std::vector<short int> markers = {22,20,21,23,24,25,26,27};
	std::vector<TGraphErrors> Efficiencies;
	TCanvas can("can","",700,800);
	bool log = false;
	double xmin = 0.25, xmax = 0.7, ymin = 0.65, ymax = 0.9;
	if(log) {xmin = 0.4; xmax = 0.85; ymin = 0.2; ymax = 0.45;}
	TLegend leg(xmin,ymin,xmax,ymax);
	//Loop  over the datasets:
	int j =0;
	for(const auto& set : files){
		//Create efficiency plot
		TGraphErrors *eff = new TGraphErrors(set.size());
		eff->SetMarkerStyle(markers.at(j));
		eff->SetMarkerSize(1.5);
		eff->SetMarkerColor(colours.at(j));
		eff->SetLineColor(colours.at(j));
		eff->SetTitle("; di-jet M, [GeV]; #epsilon");
		int i = 0;
		for(const auto& f : set){
			double low_border = f.first*0.8;
			double up_border	 = f.first*1.2;

			TH1* denum =  (TH1*) f.second->Get("distributions/NumberOfGenEvents_afterMHat_rewPU");
			TH1 *integ = (TH1*) f.second->Get("bbH_Mbb");
			double nTot = lumis.at(j);
			double selected = integ->Integral();
			eff->SetPoint(i,f.first,selected/nTot);
			std::cout<<"M12 = "<<f.first<<" Eff. = "<<selected/nTot<<" Tot = "<<nTot<<" Sel. = "<<selected<<std::endl;
			++i;
		}
		i = 0;
		++j;
		Efficiencies.push_back(*eff);
	}

	for(unsigned int i = 0; i != Efficiencies.size(); ++i){
		TGraphErrors& eff = Efficiencies.at(i);
		if(log)gPad->SetLogy();
		double eff_max = 0.5, eff_min = 1e-03;
		if(log) {
			eff_min = 0.8*1e-05;
			eff_max = 1;
		}
//		eff.SetMaximum(0.021);
		eff.SetMinimum(eff_min);
		eff.SetMaximum(eff_max);

		if(i==0) eff.Draw("APL");
		else eff.Draw("PLsame");

		leg.AddEntry(&eff,(legenda.at(i)).c_str(),"pl");
	}
	style.drawStandardTitle();
	leg.Draw();
	std::string save_name =  "SignalEfficiency_preapp_check";
	if(log) save_name += "_log";
	can.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/" + save_name + ".pdf").c_str());
}

void getSignalSamples(const std::string& sample_name_prefix, std::map<int,TFile* > & file_map){
	/*
	 * Function to get signal sample names
	 */
	file_map[300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[350] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[400] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[500]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[600]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[700] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[900] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[1100] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root").c_str());
	file_map[1300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + sample_name_prefix + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
}
