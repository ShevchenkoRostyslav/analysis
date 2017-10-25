/*
 * plotGxBRLimits.cpp
 *
 *  Created on: 19 Aug 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot GxBR limits
 *  from the output of the combine tool.
 *  Comparison of 13 TeV and 7+8 TeV combination
 *  can be performed as well.
 */

#include <iostream>
#include <string>

#include "Analysis/MssmHbb/interface/HbbLimits.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

int main(){

	HbbStyle style;
	style.setTDRstyle(PRELIMINARY);

	//Prefix to the output
	string output_prefix = "13TeV_limits";
	//paths with results
	string path2016_solo 	= cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/bias/Hbb.limits";
//	string path_wo_pdf 		= cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/bias/Hbb.limits";
//	string path_wo_JES 		= cmsswBase + "/src/Analysis/MssmHbb/datacards/201709/20/unblinded/independent/bias/no_JES/Hbb.limits";
	//ouptut folder
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170823/";
	if(!mssmhbb::blinded) output += "unblinded/";
	string output_name = output + output_prefix;
	CheckOutputDir(output);

	HbbLimits limits(mssmhbb::blinded);
	limits.ReadCombineLimits(path2016_solo);
	limits.Write(output_name + "_w_PDF");

//	HbbLimits limits_wo_PDF(mssmhbb::blinded);
//	limits_wo_PDF.ReadCombineLimits(path_wo_pdf);
//	limits_wo_PDF.Write(output_name + "_wo_PDF");

//	HbbLimits limits_wo_JES(mssmhbb::blinded);
//	limits_wo_JES.ReadCombineLimits(path_wo_JES);
//	limits_wo_JES.Write(output_name + "_wo_JES_wo_PDF");

	HbbLimits::LimitsToCompare Limits_PDF;
	Limits_PDF.legend = "_PDF_unc";
	Limits_PDF.limits = limits.getLimits();

//	HbbLimits::LimitsToCompare Limits_no_JES;
//	Limits_no_JES.legend = "no_JES_PDF_unc";
//	Limits_no_JES.limits = limits_wo_JES.getLimits();

//	HbbLimits::LimitsToCompare Limits_no_PDF;
//	Limits_no_PDF.legend = "no_PDF_unc";
//	Limits_no_PDF.limits = limits_wo_PDF.getLimits();


	//if comparison needed - insert vector of CompareLimits and TLegend in front of the LimitPlotter method
//	limits.LimitPlotter(output_name,"35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%C.L. limit on #sigma x BR [pb]",true);

//	TLegend leg(0.62,0.6,0.9,0.85);
//	limits_wo_PDF.LimitPlotter(Limits_no_JES,leg,output_name + "_vs_no_JES","35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%C.L. limit on #sigma x BR [pb]",true);
//	leg.Clear();
//	limits.LimitPlotter(leg,output_name,"35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%CL limit on #sigma x BR [pb]",true);
	limits.LimitPlotter(output_name,"35.7 (2016) fb^{-1}",string(HbbStyle::axisTitleMAH().Data()),"#sigma B [pb]",true);

}
