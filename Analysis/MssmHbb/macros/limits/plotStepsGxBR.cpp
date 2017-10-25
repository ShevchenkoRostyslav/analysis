/*
 * plotStepsGxBR.cpp
 *
 *  Created on: 6 Sep 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot bias step GxBR limits
 *  from the output of the combine tool.
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
	string output_prefix = "adjoint_13TeV_limits";
	//paths with results
	string path2016_sr1 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/adjoint_bias/sr1/Hbb.limits";
	string path2016_sr2 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/adjoint_bias/sr2/Hbb.limits";
	string path2016_sr3 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/adjoint_bias/sr3/Hbb.limits"; 
	//ouptut folder
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170823/";
	if(!mssmhbb::blinded) output += "unblinded/";
	string output_name = output + output_prefix;
	CheckOutputDir(output);

	HbbLimits limits_sr1(mssmhbb::blinded);
	limits_sr1.ReadCombineLimits(path2016_sr1);
	limits_sr1.Write(output_name + "_sr1");
	
	HbbLimits limits_sr2(mssmhbb::blinded);
	limits_sr2.ReadCombineLimits(path2016_sr2);
	limits_sr2.Write(output_name + "_sr2");
	
	HbbLimits limits_sr3(mssmhbb::blinded);
	limits_sr3.ReadCombineLimits(path2016_sr3);
	limits_sr3.Write(output_name + "_sr3");

	//if comparison needed - insert vector of CompareLimits and TLegend in front of the LimitPlotter method
//	limits_sr1.PlotSubRangeSteps(limits_sr1.getLimits(),limits_sr2.getLimits(),limits_sr3.getLimits(),output_name,"35.7 (2016) fb^{-1}","M_{A/H} [GeV]","95%CL limit on #sigma x BR [pb]",true);
	limits_sr1.PlotSubRangeSteps(limits_sr1.getLimits(),limits_sr2.getLimits(),limits_sr3.getLimits(),output_name,"35.7 (2016) fb^{-1}",string(HbbStyle::axisTitleMAH().Data()),"#sigma B [pb]",true);
}
