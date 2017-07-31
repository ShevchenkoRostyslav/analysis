/*
 * plotTanBetaLimits.cpp
 *
 *  Created on: 9 Mar 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot 2HDM limits
 *  as a function of cos(beta-alpha) and tan(beta)
 *  for particluar value of mA
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

HbbStyle style;

void Plot_finale_1(const std::vector<Limit>& limits,
		TLegend& leg,
		const bool&,
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

TGraph GetAtlasZhll_flipped_lower_left();
TGraph GetAtlasZhll_flipped_upper_left();
TGraph GetAtlasZhll_flipped_lower_right();
TGraph GetAtlasZhll_flipped_upper_right();

TGraph GetAtlasZhll_type2_left();
TGraph GetAtlasZhll_type2_lower_right();
TGraph GetAtlasZhll_type2_upper_right();

void AtlasMeasurementsStyle(TGraph &gr);

using namespace std;
using namespace analysis::mssmhbb;

int main(int argc, const char **argv){

	HbbLimits limits(true,true);
	style.set(PRELIMINARY);
	//paths with results of the combine tool
	string path2016 = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201707/26/blinded/mssm/bias/Hbb.limits";
	//Details of the 2HDM produciton
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION";
	// type of the 2hdm: type2 or type3
	string thdm_type = "type";
	if(argc == 1)thdm_type = "type2";
	else {
		thdm_type += string(argv[1]);
	}
	string thdm_scans = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + thdm_type + "_mA_mH.root";
	//value of mA
	double mass =300;

	//higgs boson: H/A/both
	string boson = "both";
	string output_dir = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170726/2hdm/";
	CheckOutputDir(output_dir);

	vector<Limit> GBR2016 = limits.ReadCombineLimits(path2016);
	TH3D GxBR_2hdm_3D;
	TH2D GxBR_2hdm_mA;
	limits.SetHiggsBoson(boson);
	vector<Limit> THDM_limits;
	for(const auto& l : GBR2016){
		cout<<"M: "<<l.getX()<<" exp = "<<l.getMedian()<<" 1G = "<<l.getPlus1G()<<endl;
		if(l.getX() == mass){
			GxBR_2hdm_3D = limits.Get2HDM_GxBR_3D(thdm_scans);
			GxBR_2hdm_mA = limits.Get2HDM_GxBR_2D(GxBR_2hdm_3D,mass,"x");
			GxBR_2hdm_mA.GetXaxis()->SetRangeUser(-1.,1.);
			if(thdm_type!="type1" && thdm_type!="type4") THDM_limits = limits.Get2HDM_Limits(GxBR_2hdm_mA,l);
		}
	}

	TLegend leg(0.65,0.65,0.87,0.9); // In case of
	leg.SetFillColor(0);
	leg.SetTextSize(0.035);
	leg.SetBorderSize(0);

	std::vector<Limit> null_vc;
	double xmin = -0.8 ,xmax = 0.8,ymin = 4,ymax = 50;
	string output_brazil = output_dir + boson + "_2HDM_" + thdm_type + "_Limits_mA" + to_string((int)mass) + "_brazil_finale";
	if(thdm_type!="type1" && thdm_type!="type4"){
		// plot full cos(b-a) range
		Plot_finale_1(THDM_limits,leg,false,output_brazil,0,100,-1.,1.,"35.7","cos(#beta-#alpha)","tan#beta",thdm_type,false);
		leg.SetX1NDC(0.5); leg.SetX2NDC(0.75);
		leg.SetY1NDC(0.18); leg.SetY2NDC(0.43);
		leg.Clear();
		output_brazil += "_zoomed";
		// plot zoomed cos(b-a) range
		if(thdm_type == "type2"){
//			xmin = -0.9; xmax = 0.9; ymin = 0.5; ymax = 50;
			leg.SetX1NDC(0.45); leg.SetX2NDC(0.7);
			leg.SetY1NDC(0.2); leg.SetY2NDC(0.45);
		}
		Plot_finale_1(THDM_limits,leg,true,output_brazil,ymin,ymax,xmin,xmax,"35.7","cos(#beta-#alpha)","tan#beta",thdm_type,true);
	}

}

void Plot_finale_1(const std::vector<Limit>& limits,
		TLegend& leg,
		const bool& draw_atlas,
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

	TGraph * obsG = new TGraph(nPointsX, X.data(), obs.data());
	obsG->SetLineWidth(3);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(2);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(20);
	obsG->SetMarkerSize(1.4);

	TGraph * expG = new TGraph(nPointsX, X.data(),median.data());
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

	TGraph atlas_zhll_flipped_lower_left 	= GetAtlasZhll_flipped_lower_left();
	TGraph atlas_zhll_flipped_upper_left 	= GetAtlasZhll_flipped_upper_left();
	TGraph atlas_zhll_flipped_lower_right	= GetAtlasZhll_flipped_lower_right();
	TGraph atlas_zhll_flipped_upper_right	= GetAtlasZhll_flipped_upper_right();

	TGraph atlas_zhll_type2_lower_left 	= GetAtlasZhll_type2_left();
	TGraph atlas_zhll_type2_lower_right	= GetAtlasZhll_type2_lower_right();
	TGraph atlas_zhll_type2_upper_right	= GetAtlasZhll_type2_upper_right();

	TCanvas *canv = new TCanvas("canv", "histograms", 600, 600);

	frame->Draw();

	outerBand->Draw("3same");
	innerBand->Draw("3same");
//	expG->SetFillStyle(3002);
//	expG->SetFillColor(kRed-10);
	expG->Draw("lsame");
	if (!blindData_)
		obsG->Draw("lpsame");

	AtlasMeasurementsStyle(atlas_zhll_flipped_lower_left);
	AtlasMeasurementsStyle(atlas_zhll_flipped_upper_left);
	AtlasMeasurementsStyle(atlas_zhll_flipped_lower_right);
	AtlasMeasurementsStyle(atlas_zhll_flipped_upper_right);

	AtlasMeasurementsStyle(atlas_zhll_type2_lower_left);
	AtlasMeasurementsStyle(atlas_zhll_type2_lower_right);
	AtlasMeasurementsStyle(atlas_zhll_type2_upper_right);
	if(draw_atlas){
		if(thdm_type == "type3"){
			atlas_zhll_flipped_lower_left.Draw("LFsame");
			atlas_zhll_flipped_upper_left.Draw("LFsame");
			atlas_zhll_flipped_lower_right.Draw("LFsame");
			atlas_zhll_flipped_upper_right.Draw("LFsame");
		}
		if(thdm_type == "type2"){
			atlas_zhll_type2_lower_left.Draw("LFsame");
			atlas_zhll_type2_lower_right.Draw("LFsame");
			atlas_zhll_type2_upper_right.Draw("LFsame");
		}
	}


	if (!blindData_)
		leg.AddEntry(obsG,"Observed","lp");
	leg.AddEntry(expG,"Expected","l");
	leg.AddEntry(innerBand,"#pm1#sigma Expected","f");
	leg.AddEntry(outerBand,"#pm2#sigma Expected","f");
	if(draw_atlas) leg.AddEntry(&atlas_zhll_flipped_lower_left,"by A#rightarrow Zh, ATLAS","f");
	style.drawStandardTitle();

	TPad * pad = (TPad*)canv->GetPad(0);
//	Luminosity lum;
//	lum.writeExtraText = true;
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


TGraph GetAtlasZhll_flipped_lower_left(){
	vector<pair<double,double>> lower_left= {
			make_pair(-0.02,0.5),
			make_pair(-0.021,0.6),
			make_pair(-0.022,0.7),
			make_pair(-0.023,0.8),
				make_pair(-0.024,0.9),
				make_pair(-0.025,1.),
				make_pair(-0.026,1.3),
				make_pair(-0.05,2.),
				make_pair(-0.054,2.5),
				make_pair(-0.1,3.),
				make_pair(-0.15,3.5),
				make_pair(-0.21,4),
				make_pair(-0.3,4.5),
				make_pair(-0.45,5),
				make_pair(-0.64,5.5),
				make_pair(-0.95,6),

				//ADDITIONAL POINT TO DRAW!!!
				make_pair(-0.95,0.5)
		};

		std::vector<double> cos_lower_left, tan_lower_left;
		for(const auto& val : lower_left) {
			cos_lower_left.push_back(val.first);
			tan_lower_left.push_back(val.second);
		}

		TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());

		return fl_lower_left;
}
TGraph GetAtlasZhll_flipped_upper_left(){
	std::vector<pair<double,double> > upper_left = {
			make_pair(-0.95,8.5),
			make_pair(-0.75,9.),
			make_pair(-0.7,9.5),
			make_pair(-0.675,10.),
			make_pair(-0.65,11),
			make_pair(-0.6,12),
			make_pair(-0.55,12.5),
			make_pair(-0.5,13.5),
			make_pair(-0.45,14),
			make_pair(-0.4,18),
			make_pair(-0.375,25.),
			make_pair(-0.36,40.),
			make_pair(-0.35,50.),

			//ADDITIONAL POINT TO DRAW!!!
			make_pair(-0.95,50)
	};

	std::vector<double> cos_upper_left, tan_upper_left;
	for(const auto& val : upper_left) {
		cos_upper_left.push_back(val.first);
		tan_upper_left.push_back(val.second);
	}

	TGraph fl_upper_left(cos_upper_left.size(),cos_upper_left.data(),tan_upper_left.data());

	return fl_upper_left;
}
TGraph GetAtlasZhll_flipped_lower_right(){
	std::vector<pair<double,double> > lower_right= {
			make_pair(0.95,6),
			make_pair(0.725,5.5),
			make_pair(0.575,5),
			make_pair(0.5,4.5),
			make_pair(0.425,4),
			make_pair(0.4,3.7),
			make_pair(0.35,3.5),
			make_pair(0.3,3.5),
			make_pair(0.275,3.5),
			make_pair(0.25,3.4),
			make_pair(0.15,3.2),
			make_pair(0.12,3.),
			make_pair(0.1,2.7),
			make_pair(0.07,2.5),
			make_pair(0.05,2.),
			make_pair(0.026,1.3),
			make_pair(0.025,1.),
			make_pair(0.024,0.9),
			make_pair(0.023,0.8),
			make_pair(0.022,0.7),
			make_pair(0.021,0.6),
			make_pair(0.02,0.5),

			//ADDITIONAL POINT TO DRAW!!!
			make_pair(0.95,0.5)
	};

	std::vector<double> cos_lower_right, tan_lower_right;
	for(const auto& val : lower_right) {
		cos_lower_right.push_back(val.first);
		tan_lower_right.push_back(val.second);
	}

	TGraph fl_lower_right(cos_lower_right.size(),cos_lower_right.data(),tan_lower_right.data());

	return fl_lower_right;
}
TGraph GetAtlasZhll_flipped_upper_right(){
	std::vector<pair<double,double> > upper_right = {
			make_pair(0.35,50),
			make_pair(0.36,40.),
			make_pair(0.375,25.),
			make_pair(0.4,20),
			make_pair(0.45,14),
			make_pair(0.5,13.5),
			make_pair(0.55,12.5),
			make_pair(0.6,12),
			make_pair(0.65,11),
			make_pair(0.675,10.),
			make_pair(0.7,9.7),
			make_pair(0.75,9.),
			make_pair(0.95,8.5),

			//ADDITIONAL POINT TO DRAW!!!
			make_pair(0.95,50)
	};

	std::vector<double> cos_upper_right, tan_upper_right;
	for(const auto& val : upper_right) {
		cos_upper_right.push_back(val.first);
		tan_upper_right.push_back(val.second);
	}

	TGraph fl_upper_right(cos_upper_right.size(),cos_upper_right.data(),tan_upper_right.data());

	return fl_upper_right;
}

TGraph GetAtlasZhll_type2_left(){
	vector<pair<double,double>> curve_left= {
			make_pair(-0.9,50), // ADDITIONAL POINT TO DRAW

			make_pair(-0.35,50),
			make_pair(-0.32,32),
			make_pair(-0.35,15),
			make_pair(-0.4,10),
			make_pair(-0.42,8),
			make_pair(-0.4,7),
			make_pair(-0.35,6.2),
			make_pair(-0.32,5.9),
			make_pair(-0.3,5.7),
			make_pair(-0.25,5.2),
			make_pair(-0.2,4.8),
			make_pair(-0.15,4.1),
			make_pair(-0.1,3.5),
			make_pair(-0.05,2.5),
			make_pair(0.0,0.5),

			make_pair(-0.9,0.5),	//ADDITIONAL POINT TO DRAW
	};

	std::vector<double> cos_left, tan_left;
	for(const auto& val : curve_left) {
		cos_left.push_back(val.first);
		tan_left.push_back(val.second);
	}

	TGraph fl_left(cos_left.size(),cos_left.data(),tan_left.data());

	return fl_left;
}

TGraph GetAtlasZhll_type2_lower_right(){
	vector<pair<double,double>> curve_lright = {
			make_pair(0.0,0.5),		//ADDITIONAL POINT TO DRAW

			make_pair(0.05,2.4),
			make_pair(0.1,3.0),
			make_pair(0.15,3.1),
			make_pair(0.2,2.9),
			make_pair(0.25,2.75),
			make_pair(0.3,2.4),
			make_pair(0.35,2.1),
			make_pair(0.4,1.9),
			make_pair(0.45,1.7),
			make_pair(0.5,1.6),
			make_pair(0.55,1.5),
			make_pair(0.6,1.32),
			make_pair(0.65,1.16),
			make_pair(0.7,0.97),
			make_pair(0.75,0.85),
			make_pair(0.8,0.7),
			make_pair(0.85,0.6),
			make_pair(0.88,0.5),

	};

	std::vector<double> cos_lright, tan_lright;
	for(const auto& val : curve_lright) {
		cos_lright.push_back(val.first);
		tan_lright.push_back(val.second);
	}

	TGraph fl_lright(cos_lright.size(),cos_lright.data(),tan_lright.data());

	return fl_lright;
}

TGraph GetAtlasZhll_type2_upper_right(){
	vector<pair<double,double>> curve_hright= {
			make_pair(0.9,50),		//ADDITIONAL POINT TO DRAW

			make_pair(0.35,50),
			make_pair(0.31,32),
			make_pair(0.35,15),
			make_pair(0.4,12),
			make_pair(0.45,7.5),
			make_pair(0.4,5.7),
			make_pair(0.37,4.5),
			make_pair(0.4,3.5),
			make_pair(0.45,2.7),
			make_pair(0.5,2.3),
			make_pair(0.55,1.7),
			make_pair(0.6,1.53),
			make_pair(0.65,1.26),
			make_pair(0.7,1.11),
			make_pair(0.75,0.92),
			make_pair(0.8,0.8),
			make_pair(0.85,0.65),
			make_pair(0.9,0.5)


	};

	std::vector<double> cos_hright, tan_hright;
	for(const auto& val : curve_hright) {
		cos_hright.push_back(val.first);
		tan_hright.push_back(val.second);
	}

	TGraph fl_hright(cos_hright.size(),cos_hright.data(),tan_hright.data());

	return fl_hright;
}
