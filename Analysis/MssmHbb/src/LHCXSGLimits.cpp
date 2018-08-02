/*
 * LHCXSGLimits.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/LHCXSGLimits.h"

namespace analysis {
namespace mssmhbb {

LHCXSGLimits::LHCXSGLimits(const bool& blindData, const std::string& boson, const std::string& banchmark_path, const double& xMin, const double& xMax, const double& yMin, const double& yMax, const bool& test) :
		LimitsInterpretation(blindData,boson,xMin,xMax,yMin,yMax,test),
		xs_tool_(mssm_xs_tools(banchmark_path.c_str(),true,0)){
}

void LHCXSGLimits::ReadCombineLimits(const std::string& file_name, const std::string& uncert, const bool& UP, const std::string& benchmark_ref_path, const double& tanBref){
	/*
	 * Method to get MSSM interpretation of GxBR limits
	 */

	//	input stream:
	std::ifstream inputList(file_name);
	if(is_empty(inputList)) std::logic_error("ERROR in HbbLimits::ReadCombineLimits. WRONG INPUT FILE NAME");
	//	TFile name
	std::string tfileName;
	// Clear limit vectors for safety
	GxBr_limits_.clear();
	limits_.clear();
	//	Loop over the lines in the input file:
	while (inputList >> tfileName) {
		//GxBR limits from combine
		Limit limit = ReadCombineLimit(tfileName,blindData_);
		GxBr_limits_.push_back(limit);
		//Translate into LHCXSG limits
		Limit tan_b_limit;
		tan_b_limit.setX(limit.getX());
		//		Set limits
		tan_b_limit.setMinus2G(double(TanBeta(tan_b_limit.getX(),limit.getMinus2G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setMinus1G(double(TanBeta(tan_b_limit.getX(),limit.getMinus1G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setMedian(double(TanBeta(tan_b_limit.getX(),limit.getMedian(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setPlus1G(double(TanBeta(tan_b_limit.getX(),limit.getPlus1G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setPlus2G(double(TanBeta(tan_b_limit.getX(),limit.getPlus2G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setObserved(double(TanBeta(tan_b_limit.getX(),limit.getObserved(),uncert,UP,benchmark_ref_path,tanBref)));

		limits_.push_back(tan_b_limit);
	}

}

double LHCXSGLimits::TanBeta(double mA, double xsection, const std::string& uncert, const bool& UP, const std::string& benchmark_ref_path, const double& tanBref){
	/*
	 * Method to translate GxBR value to tanBeta for MSSM
	 */
	//for combination - check reference path:
	bool combination = false;
	mssm_xs_tools ref(benchmark_ref_path.c_str(),true,0);
	if(benchmark_ref_path != ""){
		combination = true;
	}
	double minimalDifference = 1e+10;
	bool rangeExceeded = true;
	double tanBetaTarget = -1;
	double xsecTarget    = -1;
	double maxTanBeta = 60;
	double minTanBeta = 1;

	double sigmaBBA = -100,sigmaBBH = -100;
	double tanBeta, BrAbb, BrHbb, totXSec, difference;
	int ibmax = 10*int(maxTanBeta-minTanBeta);
	for (int ib=0; ib<ibmax; ++ib) {

		tanBeta = minTanBeta + 0.1*double(ib);

		//Check systyematics:
		if(uncert == ""){
		    sigmaBBA = xs_tool_.bbHSantander_A(mA,tanBeta);
		    sigmaBBH = xs_tool_.bbHSantander_H(mA,tanBeta);
		}
		else if (uncert == "pdfas"){
		    sigmaBBA = xs_tool_.bbHSantander_A_pdfas(mA,tanBeta,UP);
		    sigmaBBH = xs_tool_.bbHSantander_H_pdfas(mA,tanBeta,UP);
		}
		else if (uncert == "scale"){
		    sigmaBBA = xs_tool_.bbHSantander_A_scale(mA,tanBeta,UP);
		    sigmaBBH = xs_tool_.bbHSantander_H_scale(mA,tanBeta,UP);
		}

	    BrAbb = xs_tool_.br_Abb(mA,tanBeta);
	    BrHbb = xs_tool_.br_Hbb(mA,tanBeta);

	    totXSec = sigmaBBA*BrAbb + sigmaBBH*BrHbb;
	    // for the combination of 7 + 8 + 13 TeV crosssection should be divided by the xsection at reference tanB
	    if(combination){
	    	totXSec /= (ref.bbH5F_A(mA, tanBref) * ref.br_Abb(mA, tanBref) + ref.bbH5F_H(mA, tanBref) *  ref.br_Hbb(mA, tanBref));
	    }
	    difference = TMath::Abs(totXSec-xsection);

	    if (difference<minimalDifference) {
	      minimalDifference = difference;
	      tanBetaTarget = tanBeta;
	      xsecTarget = totXSec;
	//	std::cout<<"Brs: A = "<<BrAbb<<" B = "<<BrHbb<<" Sigmas: A = "<<sigmaBBA<<" B = "<<sigmaBBH<<std::endl;
	    }

	    if (totXSec>xsection) {
	      rangeExceeded = false;
	      break;
	    }

	}

	if (rangeExceeded)
		tanBetaTarget = tanBetaTarget*TMath::Sqrt(xsection/xsecTarget);

	return tanBetaTarget;
}

void LHCXSGLimits::AddPlottingObjects(TH2F &frame, TLegend &leg, TGraph& obs, TGraph& exp, TGraphAsymmErrors& inner_band, TGraphAsymmErrors& outer_band, TCanvas &can){
	/*
	 * Virtual function that is actually plotting objects on the frame
	 */
	 
	/*
	 * Appearence of observed limits.
	 * For the model dependent interpretation - there is a line, shows the central value
	 * and shaded area with fully invisible border. Needed to not cover the axis on top
	 * Also - an artificial TGraph is created to go the TLegend object with both - line and area
	 * proper colored/painted
	 */

	//Modify appearence of the observed data to match 8 TeV paper
	obs.SetLineColor(kRed);

	//Shaded area with transperent borders
	TGraph *obs_shaded = new TGraph(obs);
	obs_shaded->SetLineColor(kWhite);
	obs_shaded->SetLineWidth(0);
	obs_shaded->SetFillColor(kRed-9);
	obs_shaded->SetFillColorAlpha(kRed-9,0.35);
	obs_shaded->SetPoint(obs_shaded->GetN(),xMax_,yMax_);		//Needed to make a closed area to be filled
	obs_shaded->SetPoint(obs_shaded->GetN(),xMin_,yMax_);		//Needed to make a closed area to be filled

	//Artificial TGraph to be used in the TLegend
	TGraph *obs_legend = new TGraph();
	obs_legend->SetLineColor(obs.GetLineColor());
	obs_legend->SetLineWidth(obs.GetLineWidth());
	obs_legend->SetFillColorAlpha(obs_shaded->GetFillColor(), 0.35);


	 //Modify appearence of the observed data to match 8 TeV paper
//	obs.SetLineColor(kRed);
//	obs.SetFillColorAlpha(kRed-9,0.35);
//	obs.SetPoint(obs.GetN(),xMax_,yMax_);
//	obs.SetPoint(obs.GetN(),xMin_,yMax_);

	if (!blindData_){
		obs_shaded->Draw("FLsame");
//		obs.Draw("FLsame");
		obs.Draw("Lsame");
		leg.AddEntry(obs_legend,"Observed","lf");
	}
	//Add standard legend for the limits plot
	exp.SetLineColor(kBlue);
	leg.AddEntry(&exp,"Expected","l");
	leg.AddEntry(&inner_band,"68% expected","f");
	leg.AddEntry(&outer_band,"95% expected","f");
    // add excluded areas in phase space, according to the SMH
	if(drawExcludedRegion_ != 0){
		std::vector<TGraphAsymmErrors*> excludedAreas = getH125IncompatibleAreas(xs_tool_,drawExcludedRegion_);
		for(const auto& it : excludedAreas) {
			it->SetLineWidth(0);
		    it->SetFillColor(15);
		    it->SetFillStyle(3645);
		    it->Draw("FX same");
		}
		//Add legend entry for it.
		leg.AddEntry(excludedAreas.front(),("m_{h,H} #neq 125 #pm " + std::to_string(drawExcludedRegion_) + " GeV").c_str(),"FX");
	}

	//Compare with previous results:
	if(compareWithPrevious_ != "" && scenario_->previousExists()){
		scenario_->checkResultsToCompareWith(compareWithPrevious_);
		auto previous_res = scenario_->getPreviousResults(compareWithPrevious_);
		for(auto & gr : previous_res){
			TGraph *p_gr = new TGraph(gr);
			p_gr->Draw("Lsame");
		}
//		TGraph *previous_res = new TGraph(scenario_->getPreviousResults(compareWithPrevious_));
		TText *previous_lab = new TText(scenario_->getPreviousResultsLabel(compareWithPrevious_));
//		previous_res->Draw("Lsame");
		previous_lab->Draw("same");

//		/* Temproary for the Journal review
// 		*  putting the observed limits as well
// 		*/
//		std::array<double,8> x = {{100,140,160,200,300,350,400,500}};
//		std::array<double,8> y = {{14.2,22.4,18.2,14.8,34.1,37.1,34.4,50.0}};
//	        TGraph *gr_obs = new TGraph(8,&x[0],&y[0]);
//        	gr_obs->SetLineColor(kBlack);
//		gr_obs->Draw("Lsame");
	}

	auto legendHeader = scenario_->getLabel();
	leg.SetHeader(legendHeader.c_str());
	leg.SetX1(0.8*leg.GetX1());

	leg.Draw();
	//In case of the long header like at mhmod+ scenario
	if(legendHeader.length() > 30) HbbStyle::drawLegendSplittedHeader(&leg,",");

	HbbStyle::drawStandardTitle("out");

}

} /* namespace mssmhbb */
} /* namespace analysis */
