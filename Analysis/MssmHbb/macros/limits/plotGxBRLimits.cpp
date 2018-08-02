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
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

int main(){

	HbbStyle style;
//	style.setTDRstyle(mssmhbb::publication_status);
	style.setTDRstyle(PRIVATE);

	//blinded
	string prefix = "";
	if(!mssmhbb::blinded) prefix = "unblinded"; else prefix = "blinded";
	//Prefix to the output
	string output_prefix = "13TeV_limits";
	//paths with results
	string path2016_solo 	= cmsswBase + "/src/Analysis/MssmHbb/datacards/201712/13/" + prefix + "/independent/bias/Hbb.limits";
	string path_to_compare 		= cmsswBase + "/src/Analysis/MssmHbb/datacards/201712/13/" + prefix + "/independent/adjoint_bias/500_sr1/Hbb.limits";
	string path_to_compare1100 		= cmsswBase + "/src/Analysis/MssmHbb/datacards/201712/13/" + prefix + "/independent/adjoint_bias/1100_sr2/Hbb.limits";
//	string path_to_compare 		= cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/" + prefix + "/independent/bias/Hbb.limits";
	//ouptut folder
//	string output = mssmhbb::pictures_output + "ParametricLimits/20171213/";//cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20171213/";
	string output = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20180427/";
	if(!mssmhbb::blinded) output += "unblinded/";
	string output_name = output + output_prefix;
	CheckOutputDir(output);

	HbbLimits limits(mssmhbb::blinded,290,1310,1.5e-01,78);
	limits.ReadCombineLimits(path2016_solo);
	limits.Write(output_name + "_w_NLO");

	HbbLimits limits_to_compare(mssmhbb::blinded,290,1310,1.5e-01,78);
	limits_to_compare.ReadCombineLimits(path_to_compare);
	limits_to_compare.Write(output_name + "_500SR1");

	HbbLimits limits_to_compare1100(mssmhbb::blinded,290,1310,1.5e-01,78);
	limits_to_compare1100.ReadCombineLimits(path_to_compare1100);
	limits_to_compare1100.Write(output_name + "_1100SR2");

//	HbbLimits limits_wo_JES(mssmhbb::blinded);
//	limits_wo_JES.ReadCombineLimits(path_wo_JES);
//	limits_wo_JES.Write(output_name + "_wo_JES_wo_PDF");

//	HbbLimits::LimitsToCompare limits_to_compare2;
//	limits_to_compare2.legend = "_no_NLO";
//	limits_to_compare2.limits = limits_to_compare.getLimits();

	HbbLimits::LimitsToCompare Limits_500_sr1;
	Limits_500_sr1.legend = "500 GeV, SR1";
	Limits_500_sr1.limits = limits_to_compare.getLimits();

	HbbLimits::LimitsToCompare Limits_1100_sr2;
	Limits_1100_sr2.legend = "1100 GeV, SR2";
	Limits_1100_sr2.limits = limits_to_compare1100.getLimits();

//	HbbLimits::LimitsToCompare Limits_no_PDF;
//	Limits_no_PDF.legend = "no_PDF_unc";
//	Limits_no_PDF.limits = limits_to_compare.getLimits();


	//if comparison needed - insert vector of CompareLimits and TLegend in front of the LimitPlotter method
//	limits.LimitPlotter(output_name,"35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%C.L. limit on #sigma x BR [pb]",true);

//	TLegend leg(0.62,0.6,0.9,0.85);
	auto &leg0 = *HbbStyle::legend("top,right",6);
//	leg0.SetX1(leg0.GetX1()*0.9);
	limits.LimitPlotter(Limits_500_sr1,leg0,output_name + "_vs_500SR1","35.7 (2016) fb^{-1}",string(HbbStyle::axisTitleMAH().Data()),"#sigma(b#bar{b}A/H) #bf{#it{#Beta}}(A/H#rightarrow b#bar{b}) [pb]",true);
	limits.LimitPlotter(Limits_1100_sr2,leg0,output_name + "_vs_1100SR2","35.7 (2016) fb^{-1}",string(HbbStyle::axisTitleMAH().Data()),"#sigma(b#bar{b}A/H) #bf{#it{#Beta}}(A/H#rightarrow b#bar{b}) [pb]",true);
//	leg.Clear();
//	limits.LimitPlotter(leg,output_name,"35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%CL limit on #sigma x BR [pb]",true);
	limits.LimitPlotter(output_name,"35.7 (2016) fb^{-1}",string(HbbStyle::axisTitleMAH().Data()),"#sigma(b#bar{b}A/H) #bf{#it{#Beta}}(A/H#rightarrow b#bar{b}) [pb]",true);

}
