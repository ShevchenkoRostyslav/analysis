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
	style.set(PRELIMINARY);

	//Prefix to the output
	string output_prefix = "solo_13TeV_limits";
	//paths with results
	string path2016_solo 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201707/26/unblinded/mssm/bias/Hbb.limits";
	//ouptut folder
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170818/";
	if(!mssmhbb::blinded) output += "unblinded/";
	string output_name = output + output_prefix;
	CheckOutputDir(output);

	HbbLimits limits(mssmhbb::blinded);
	limits.ReadCombineLimits(path2016_solo);
	limits.Write(output_name);

	//if comparison needed - insert vector of CompareLimits and TLegend in front of the LimitPlotter method
	limits.LimitPlotter(output_name,"35.7 (2016) fb^{-1}","M_{A} [GeV]","95%C.L. limit on #sigma x BR [pb]",true);
}
