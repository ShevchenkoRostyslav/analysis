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

void LHCXSGLimits::AddPlottingObjects(TH2F &frame, TLegend &leg, TGraph& obs, TGraph& exp, TGraphAsymmErrors& inner_band, TGraphAsymmErrors& outer_band){
	/*
	 * Virtual function that is actually plotting objects on the frame
	 */
	 
	//Modify appearence of the observed data to match 8 TeV paper
	obs.SetLineColor(kRed);
	obs.SetFillStyle(3002);
	obs.SetFillColor(kRed-9);
	obs.SetPoint(obs.GetN(),xMax_,yMax_);
	obs.SetPoint(obs.GetN(),xMin_,yMax_);
	 
	if (!blindData_){
		obs.Draw("FLsame");
		obs.Draw("Lsame");
		leg.AddEntry(&obs,"Observed","lf");
	}
	//Add standard legend for the limits plot
	exp.SetLineColor(kBlue);
	leg.AddEntry(&exp,"Expected","l");
	leg.AddEntry(&inner_band,"#pm1#sigma Expected","f");
	leg.AddEntry(&outer_band,"#pm2#sigma Expected","f");
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
//		//Add legend to the exclusion plot
//		  TLegend *leg_excl = new TLegend(0.21,.705 ,.45,.845);
//		  leg_excl->SetHeader(scenarioLabel_.c_str());
//		  leg_excl->SetBorderSize(1);
//		  leg_excl->SetLineColor(kBlack);
//		  leg_excl->SetLineWidth(2);
//		  leg_excl->SetTextSize(0.035);
//		  leg_excl->AddEntry(excludedAreas.front(),"m_{h,H} #neq 125 #pm 3 GeV","FX");
//		  leg_excl->Draw();
	}
	leg.SetHeader(scenario_->getLabel().c_str());
	leg.SetX1(0.8*leg.GetX1());
	
	//Compare with previous results:
	if(compareWithPrevious_ && scenario_->previousExists()){
		TGraph *previous_res = new TGraph(scenario_->getPreviousResults());
		TText *previous_lab = new TText(scenario_->getPreviousResultsLabel());
		previous_res->Draw("Lsame");
		previous_lab->Draw("same");
	}

//	  //Add text block describing scenario
//	  TPaveText *text = new TPaveText(0.2,.7,.4,.85,"NDC");
//	  text->SetFillColor(0);
//	  text->SetBorderSize(1);
//	  text->SetTextSize(0.035);
//	  text->SetLineWidth(2);
//	  text->AddText("m_{h}^{mod+} scenario");
//	  text->AddText("#mu=-200");
//	  text->Draw();

	  //Add legend to the exclusion plot
//	  TLegend *leg_excl = new TLegend(0.21,.705 ,.45,.845 - 2*HbbStyle::lineHeight());
//	  leg_excl->SetLineColor(kBlack);
//	  leg_excl->SetLineWidth(1);
//	  leg_excl->SetTextSize(0.035);
//	  leg_excl->AddEntry(excludedAreas.front(),"m_{h,H} #neq 125#pm3 GeV","FX");
//	  leg_excl->Draw();

}

} /* namespace mssmhbb */
} /* namespace analysis */
