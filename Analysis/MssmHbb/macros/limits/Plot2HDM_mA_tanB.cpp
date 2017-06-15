/*
 * plotTanBetaLimits.cpp
 *
 *  Created on: 9 Mar 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot 2HDM limits
 *  as a function of mA and tanBeta for particluar
 *  value of cos(beta-alpha)
 */


#include <iostream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TPad.h>
#include <TFile.h>

#include "Analysis/MssmHbb/interface/HbbLimits.h"
#include "Analysis/MssmHbb/interface/Limit.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"

void Plot_finale_1(const std::vector<Limit>& limits,
		TLegend& leg,
		const std::string& output,
		const float& yMin,
		const float& yMax,
		const float& xMin,
		const float& xMax,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const std::string& thdm_type,
		const bool& logY);
//TGraph GetAtlasZhll_flipped();

TGraph GetAtlasZhll_flipped();
TGraph GetAtlasZhll_type2();

HbbStyle style;

void AtlasMeasurementsStyle(TGraph &gr);

using namespace std;
using namespace analysis::mssmhbb;

int main(){

	HbbLimits limits(true,true);
	style.set(PRELIMINARY);
	//paths with results of the combine tool
	string path2016 = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201705/18/Asymptotic/mssm/Hbb.limits";
	//Details of the 2HDM produciton
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION";
	// type of the 2hdm: type2 or type3
	string thdm_type = "type2";
	string thdm_scans = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + thdm_type + "_mA_mH.root";
	//value of cos(beta-alpha)
	double cB_A = 0.1;

	//higgs boson: H/A/both
	string boson = "both";
	string output_dir = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170518/2hdm/";

	vector<Limit> GBR2016 = limits.ReadCombineLimits(path2016);
	vector<Limit> TanB_2HDM_Limits;
	TH3D GxBR_2hdm_3D;
	TH2D GxBR_2hdm_cB_A;

	limits.SetHiggsBoson(boson);

	GxBR_2hdm_3D = limits.Get2HDM_GxBR_3D(thdm_scans);
	GxBR_2hdm_cB_A = limits.Get2HDM_GxBR_2D(GxBR_2hdm_3D,cB_A,"z");
	TanB_2HDM_Limits = limits.Get2HDM_1D_Limits(GxBR_2hdm_cB_A,GBR2016);

	// 2HDM tanBeta vs mA limits for cos(beta-alpha) = cB_A
	TLegend leg_2HDM_tB(0.62,0.25,0.85,0.45);
	leg_2HDM_tB.SetFillColor(0);
	leg_2HDM_tB.SetTextSize(0.035);
	leg_2HDM_tB.SetBorderSize(0);
	string output_2HDM_tanB_limits = output_dir + boson + "_2HDM_" + thdm_type + "_Limits_tanB_cB_A" + to_string((double)cB_A) + "_brazil";
	Plot_finale_1(TanB_2HDM_Limits,leg_2HDM_tB,output_2HDM_tanB_limits,0.5,60,200,900,"35.7","M_{#Phi} [GeV]","tan(#beta)",thdm_type,true);
	output_2HDM_tanB_limits += "_zoomed";
	leg_2HDM_tB.Clear();
	Plot_finale_1(TanB_2HDM_Limits,leg_2HDM_tB,output_2HDM_tanB_limits,0.5,60,300,900,"35.7","M_{#Phi} [GeV]","tan(#beta)",thdm_type,true);

}

void Plot_finale_1(const std::vector<Limit>& limits,
		TLegend& leg,
		const std::string& output,
		const float& yMin,
		const float& yMax,
		const float& xMin,
		const float& xMax,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const std::string& thdm_type,
		const bool& logY){

	bool blindData_ = true;


	if(limits.size() == 0) {
		std::cerr<<"Error: No limits with this name. Please check spelling";
		exit(-1);
	}

	gStyle->SetPadTopMargin   (0.08);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);

	int nPointsX = limits.size();
	std::vector<double> X;
	std::vector<double> obs;
	std::vector<double> median;
	std::vector<double> minus1;
	std::vector<double> minus2;
	std::vector<double> plus1;
	std::vector<double> plus2;
	std::vector<double> zero;

	for(const auto& l : limits){
		X.push_back(l.getX());
		obs.push_back(l.getObserved());
		median.push_back(l.getMedian());
		minus2.push_back(l.getMedian() - l.getMinus2G());
		minus1.push_back(l.getMedian() - l.getMinus1G());
		plus2.push_back(l.getPlus2G() - l.getMedian());
		plus1.push_back(l.getPlus1G() - l.getMedian());
		zero.push_back(0);

	}

	std::vector<double> x_cosmetic = X;
	std::vector<double> median_cosmetic = median;
	x_cosmetic.push_back(300);
	median_cosmetic.push_back(60);

	TGraph * obsG = new TGraph(nPointsX, X.data(), obs.data());
	obsG->SetLineWidth(3);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(2);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(20);
	obsG->SetMarkerSize(1.4);

	TGraph * expG = new TGraph(nPointsX+1, x_cosmetic.data(),median_cosmetic.data());
	expG->SetLineWidth(3);
	expG->SetLineColor(2);
	expG->SetLineStyle(2);

	TGraphAsymmErrors * innerBand = new TGraphAsymmErrors(nPointsX, X.data(), median.data(), zero.data(), zero.data(), minus1.data(), plus1.data());
	innerBand->SetFillColor(kGreen);
	innerBand->SetLineColor(kGreen);

	TGraphAsymmErrors * outerBand = new TGraphAsymmErrors(nPointsX, X.data(), median.data(), zero.data(), zero.data(), minus2.data(), plus2.data());
	outerBand->SetFillColor(kYellow);
	outerBand->SetLineColor(kYellow);

	TH2F * frame = nullptr;

	frame = new TH2F("frame","",2,xMin,xMax,2,yMin,yMax);
	frame->GetXaxis()->SetTitle(xtitle.c_str());
	frame->GetYaxis()->SetTitle(ytitle.c_str());
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetNdivisions(206);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetXaxis()->SetTitleOffset(1.05);
	frame->GetYaxis()->SetTitleSize(0.048);
	frame->SetStats(0);

//	TGraph GetAtlasZhll_flipped();

	TGraph atlas_zhll_flipped 	= GetAtlasZhll_flipped();
	TGraph atlas_zhll_type2		= GetAtlasZhll_type2();

	TCanvas *canv = new TCanvas("canv", "histograms", 600, 600);

	frame->Draw();

	outerBand->Draw("3same");
	innerBand->Draw("3same");

	expG->SetFillStyle(3002);
	expG->SetFillColor(kRed-10);
	expG->Draw("lfsame");
	if (!blindData_)
		obsG->Draw("lpsame");

	AtlasMeasurementsStyle(atlas_zhll_flipped);
	AtlasMeasurementsStyle(atlas_zhll_type2);
	if(thdm_type == "type3"){
		atlas_zhll_flipped.Draw("CFsame");
	}
	else if (thdm_type == "type2"){
		atlas_zhll_type2.Draw("CFsame");
	}

	if (!blindData_)
		leg.AddEntry(obsG,"Observed","lp");
	leg.AddEntry(expG,"Expected","fl");
	leg.AddEntry(innerBand,"#pm1#sigma Expected","f");
	leg.AddEntry(outerBand,"#pm2#sigma Expected","f");
	leg.AddEntry(&atlas_zhll_flipped,"by A#rightarrow Zh, ATLAS","f");
	style.drawStandardTitle();

	TPad * pad = (TPad*)canv->GetPad(0);
//	Luminosity lum;
//	lum.writeExtraText = false;
//	lum.lumi_13TeV = Lumi;
//	lum.extraText = "Private Work";
//	int offs = 33;
//	if(ytitle.find("tan") != std::string::npos) offs = 11;
//	lum.CMS_lumi(pad,4,offs);

	pad->RedrawAxis();

	leg.Draw();

	if(logY) canv->SetLogy();
	canv->Update();
	canv->Print( (output+".pdf").c_str() ,"Portrait pdf");
	canv->Print( (output+".png").c_str()) ;
}

void AtlasMeasurementsStyle(TGraph &gr){
	gr.SetFillStyle(3002);
	gr.SetFillColor(kAzure-4);
}


TGraph GetAtlasZhll_flipped(){
	vector<pair<double,double>> lower_left= {
			make_pair(220,0.5),
			make_pair(225,1),
			make_pair(230,1.8),
			make_pair(235,2.),
			make_pair(240,2.3),
			make_pair(245,2.5),
			make_pair(250,2.7),
			make_pair(255,2.8),
			make_pair(260,2.9),
			make_pair(265,2.95),
			make_pair(270,3.),
			make_pair(275,3.05),
			make_pair(280,3.1),
			make_pair(285,3.15),
			make_pair(290,3.2),
			make_pair(295,3.3),
			make_pair(300,3.4),
			make_pair(305,3.5),
			make_pair(310,3.6),
			make_pair(315,3.7),
			make_pair(320,3.8),
			make_pair(325,3.9),
			make_pair(330,3.95),
			make_pair(335,4),
			make_pair(340,4),
			make_pair(345,3.7),
			make_pair(347.5,3),
			make_pair(350,0.5),

				//ADDITIONAL POINT TO DRAW!!!
//				make_pair(-0.95,0.5)
		};

		std::vector<double> cos_lower_left, tan_lower_left;
		for(const auto& val : lower_left) {
			cos_lower_left.push_back(val.first);
			tan_lower_left.push_back(val.second);
		}

		TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());

		return fl_lower_left;
}

TGraph GetAtlasZhll_type2(){
	vector<pair<double,double>> lower_left= {
			make_pair(220,0.5),
			make_pair(225,1.21),
			make_pair(230,1.68),
			make_pair(235,1.95),
			make_pair(240,2.18),
			make_pair(245,2.27),
			make_pair(250,2.36),
			make_pair(255,2.50),
			make_pair(260,2.59),
			make_pair(265,2.64),
			make_pair(270,2.68),
			make_pair(275,2.77),
			make_pair(280,2.82),
			make_pair(285,2.86),
			make_pair(290,2.91),
			make_pair(295,2.95),
			make_pair(300,3.00),
			make_pair(305,3.10),
			make_pair(310,3.12),
			make_pair(315,3.19),
			make_pair(320,3.25),
			make_pair(325,3.38),
			make_pair(330,3.50),
			make_pair(335,3.63),
			make_pair(337,3.69),
			make_pair(340,3.50),
			make_pair(345,3.13),
			make_pair(350,0.5)
	};

	std::vector<double> cos_lower_left, tan_lower_left;
	for(const auto& val : lower_left) {
		cos_lower_left.push_back(val.first);
		tan_lower_left.push_back(val.second);
	}

	TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());

	return fl_lower_left;
}
