/*
 * data_vs_mc.cpp
 *
 *  Created on: 25 Mar 2018
 *      Author: rostyslav
 */

#include <iostream>
#include <vector>
#include <string>

#include "TH1.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TExec.h"
#include "TFile.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;

void PerformDataVsMCComparison(const string& data_sample, const string& mc_sample, const vector<string> & variables, const string& output_path);
void PlotDataVsMC(TH1 *hdata,TH1 *hmc, const string& output_path);
void renormalise(TH1 * h);
string DetermineLegendPosition(const string& name);
void DivideTwoHistograms(TH1 * numerator, TH1 * denominator);
void DivideStatHistogram(TH1 * stat_err, TH1 * denominator, TH1 * data);

HbbStyle style;

int main(){
	style.setTDRstyle(PRIVATE);
	string data_sample 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/DataMC_ReReco_35673fb_lowM_Run2016.root";
	string mc_sample 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/DataMC_ReReco_35673fb_lowM_bEnriched_Tot.root";
	string output_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/";
	vector<string> variables = {"diJet_b_m","jet_b_Ht","NumberOfJets_b","jet_b_pt1","jet_b_pt2","jet_b_eta1","jet_b_eta2","jet_b_phi1","jet_b_phi2","jet12_b_assym","jet_b_btag_csv1",
						"jet_b_btag_csv2","jet_b_btag_cmva1","jet_b_btag_cmva2","jet_b_deta12","jet_b_dR12", "jet_b_dphi12",
						"diJet_b_pt","diJet_b_eta","diJet_b_phi"};
	CheckOutputDir(output_path);
	PerformDataVsMCComparison(data_sample,mc_sample,variables,output_path);
}

void PerformDataVsMCComparison(const string& data_sample, const string& mc_sample, const vector<string> & variables, const string& output_path){
	/*
	 * General method to perform data vs mc
	 * comparison
	 */
	TFile fdata(data_sample.c_str(),"READ");
	TFile fmc(mc_sample.c_str(),"READ");
	for(const auto & var : variables){
		string varname = var;
		auto *hdata = GetFromTFile<TH1D>(fdata,varname);
		auto *hmc = GetFromTFile<TH1D>(fmc,varname + "/" + varname);

		//Normalise MC to data
		hmc->Scale(hdata->Integral()/hmc->Integral());
		//Correct for the bin width
		renormalise(hdata);
		renormalise(hmc);

		PlotDataVsMC(hdata,hmc,output_path);
	}
}

void PlotDataVsMC(TH1 *hdata,TH1 *hmc, const string& output_path){
	/*
	 * Plot data and MC histograms
	 */
	auto *er_0 = new TExec("er_0","gStyle->SetErrorX(0)");
	auto *er_1 = new TExec("er_1","gStyle->SetErrorX(0.5)");

	TCanvas can;
	auto mc_stat_err = static_cast<TH1D*>(hmc->Clone("hstat"));
	mc_stat_err->SetLineColor(kBlack);
	mc_stat_err->SetFillColor(kBlack);
	mc_stat_err->SetFillStyle(3013); //3002

	//top pad and frame
	auto topPad = HbbStyle::getRatioTopPad(0.7);
	topPad->Draw();
	topPad->cd();

	string legend_pos = DetermineLegendPosition(hdata->GetName());
	auto leg = style.legend(legend_pos,4,0.4);
//	AdjustLegendPosition(leg);
	leg->AddEntry(hdata,"Data","p");
	leg->AddEntry(hmc,"QCD b-Enriched","f");
	leg->AddEntry(mc_stat_err,"MC Stat.","lf");

	if(findStrings(hdata->GetName(),"pt")) hdata->GetXaxis()->SetRangeUser(0,hdata->GetXaxis()->GetXmax());
	if(findStrings(hdata->GetName(), "deta")) hdata->SetAxisRange(hdata->GetMinimum(), hdata->GetMaximum()*1.2, "Y");
	//Hardcoded solutions for some variables
//	if(findStrings(hdata->GetName(), "pt") || findStrings(hdata->GetName(), "phi") || findStrings(hdata->GetName(), "Ht")){
//		hdata->SetAxisRange(0.01,hdata->GetMaximum(),"Y");
//		gPad->RedrawAxis();
//		gPad->SetLogy();
//	}

	auto *top_frame = new TH2D("top_frame","",1,hdata->GetXaxis()->GetXmin(),hdata->GetXaxis()->GetXmax(),1,hdata->GetYaxis()->GetXmin(),1.3*hdata->GetMaximum());
	top_frame->GetXaxis()->SetTitle(hdata->GetXaxis()->GetTitle());
	if(findStrings(hdata->GetName(), "assym")) top_frame->GetXaxis()->SetTitle("(p^{(1)}_{T} - p^{(2)}_{T})/(p^{(1)}_{T} + p^{(2)}_{T})");
	style.setRatioTopFrame(top_frame, can.YtoPixel(can.GetY1()));
	top_frame->GetYaxis()->SetTitle("Events/bin");
	top_frame->Draw();

	style.applyStandardMCHistoStyle(*hmc);
	hmc->SetFillColor(kRed-10);
	style.applyStandardDataHistoStyle(*hdata);
	top_frame->SetMinimum(0);

	top_frame->Draw();
	hmc->Draw("hist same");
	er_1->Draw();
	mc_stat_err->Draw("E2 same");
	er_0->Draw();
	hdata->Draw("E same");
	leg->Draw();

	gPad->RedrawAxis();
	can.cd();

	auto botPad = HbbStyle::getRatioBottomPad(0.3);
	botPad->Draw();
	botPad->cd();

	auto *botFrame = new TH2D("botFrame","",1,top_frame->GetXaxis()->GetXmin(),top_frame->GetXaxis()->GetXmax(),1,0.5,1.5);
	HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
	botFrame->GetXaxis()->SetTitle(top_frame->GetXaxis()->GetTitle());
	botFrame->GetYaxis()->SetTitle("Data/MC");
	botPad->cd();
	auto ratio = static_cast<TH1D*>(hdata->Clone("ratio"));
	ratio->Sumw2();
	//Divide two histograms
	DivideTwoHistograms(ratio,hmc);
	botFrame->Draw();
	ratio->Draw("E same");

	//Ratio for the MCstats
	auto ratio_mc_stat = static_cast<TH1D*>(mc_stat_err->Clone("ratio_mc_stat"));
	DivideStatHistogram(ratio_mc_stat,hmc,hdata);
	er_1->Draw();
	ratio_mc_stat->Draw("E2 same");
	er_0->Draw();

	TLine *horizLine = new TLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
	horizLine -> SetLineStyle(2);
	horizLine -> Draw();

	can.cd();
	if(findStrings(hdata->GetName(), "pt") || findStrings(hdata->GetName(), "phi") || findStrings(hdata->GetName(), "Ht")){
		style.drawStandardTitle("top,right");
	}
	else style.drawStandardTitle("top,left");

	string var_name = hdata->GetName();
	can.SaveAs((output_path + var_name + ".pdf").c_str());
}

void renormalise(TH1 * h){
	for(int i=1; i<= h->GetNbinsX();++i){
		h->SetBinContent(i,h->GetBinContent(i)/h->GetBinWidth(i));
		h->SetBinError(i,h->GetBinError(i)/h->GetBinWidth(i));
	}
}

string DetermineLegendPosition(const string& v){
	/*
	 * Returns the position of the
	 * TLegend according to the variable
	 */
	string leg_position = "top,right";
	if(findStrings(v,"NumberOfJets_b")) leg_position = "top,right";
	if(findStrings(v,"dR") || findStrings(v,"dphi")) leg_position = "bottom,left";
	if(findStrings(v,"deta")) leg_position = "top,right";
	if(findStrings(v,"btag")) leg_position = "bottom,left";
	return leg_position;
}

void DivideTwoHistograms(TH1 * numerator, TH1 * denominator){
	for(int i = 1; i<= numerator->GetNbinsX(); ++i){
		auto nmc = denominator->GetBinContent(i);
		auto ndata = numerator->GetBinContent(i);
		auto edata = numerator->GetBinError(i);
		if(nmc != 0) {
			numerator -> SetBinContent(i,ndata/nmc);
			numerator -> SetBinError(i,edata/nmc);
		}
		else {
			numerator -> SetBinContent(i,0);
			numerator -> SetBinError(i,0);
		}
	}
}

void DivideStatHistogram(TH1 * stat_err, TH1 * denominator, TH1 * data){
	for(int i = 1; i<= stat_err->GetNbinsX(); ++i){
		stat_err -> SetBinError(i,stat_err->GetBinError(i)/denominator->GetBinContent(i));
		stat_err -> SetBinContent(i,1);
		if(stat_err->GetBinContent(i) == 0) stat_err->SetBinError(i,999);
		if(stat_err->GetBinContent(i) == 0 && data->GetBinContent(i) == 0) stat_err->SetBinError(i,0);
	}
}
