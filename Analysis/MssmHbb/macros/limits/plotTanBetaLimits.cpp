/*
 * plotTanBetaLimits.cpp
 *
 *  Created on: 9 Mar 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot MSSM limits
 *  from the output of the combine tool.
 *  Comparison of 13 TeV and 7+8 TeV combination
 *  can be performed as well.
 */

#include <iostream>
#include <string>
#include <vector>

//ROOT includes
#include <TFile.h>

#include "Analysis/MssmHbb/interface/LHCXSGLimits.h"
#include "Analysis/MssmHbb/interface/LHCXSGScenarious.h"
#include "Analysis/MssmHbb/interface/Limit.h"
// #include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

HbbStyle style;

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

std::string getBanchmarkPath(AvailableScenarios scenario);

int main(){

	style.set(PRELIMINARY);
	//Prefix to the output
	string output_prefix = "solo_13TeV_limits";
	//paths with results
	string path2016_solo 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201707/26/unblinded/mssm/bias/Hbb.limits";
	string path2016_combined = cmsswBase + "/src/Analysis/MssmHbb/datacards/201705/15/Asymptotic/mssm/No_Bias/Hbb.limits";	//If combination is performed
	//ouptut folder
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170818/";
	if(!mssmhbb::blinded) output += "unblinded/";
	AvailableScenarios scenario = LIGHT_STOP; //MHMODP_200,LIGHT_STOP,LIGHT_STAU,HMSSM,TAU_PHOBIC
	output += AvailableScenariosToString(scenario) + "/";
	CheckOutputDir(output);
	//benchmark scenario path
	string boson = "both";
	LHCXSGLimits limits(mssmhbb::blinded,boson,getBanchmarkPath(scenario));
	limits.compareWithPrevious(true);
	limits.setScenario(scenario);
	limits.drawExcludedRegion(3);	// Region NOT compatible with H(125)
	
	string benchmark_name = AvailableScenariosToString(scenario);
	string output_name = output + boson + "_" + benchmark_name + "_" + output_prefix;

	// combined results UNCOMMENT if combination is performed
//	vector<Limit> GxBR_combined = limits.ReadCombineLimits(path2016_combined);
//	string benchmark_ref = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/signal/newmhmax_mu200_13TeV.root";
//	vector<Limit> mssm_limits_combined = limits.GetMSSMLimits(GxBR_combined,benchmark,"",false,benchmark_ref,30);
//	limits.Write(mssm_limits_combined, "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170309/combined_limits.txt");

	//solo 2016 13 TeV
	limits.ReadCombineLimits(path2016_solo);
	limits.Write(output_name + ".txt");

	//solo 7+8 TeV
	vector<Limit> tanBeta2012 = {
			Limit(300,24.8,34.1,31.1,20.5,39.2,17.6),
			Limit(350,30.0,37.1,37.5,24.8,47.1,21.4),
			Limit(400,36.5,34.4,46.9,29.8,00.0,25.5),
			Limit(500,53.5,50.0,00.0,43.2,00.0,36.8)
	};
	HbbLimits::LimitsToCompare limits_to_comp;
	limits_to_comp.limits = tanBeta2012;
	limits_to_comp.legend = "7+8 TeV exp.";

	//If one wants to comapare two different results

	//if comparison is not needed - insert null_vec as thr second argument
	limits.LimitPlotter(output_name,"24.6 vs 35.7(2016)","M_{A/H} [GeV]","tan(#beta)",false);
}

std::string getBanchmarkPath(AvailableScenarios scenario){
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/signal/";
	switch (scenario){
		case MHMODP_200:
			output += "mhmodp_mu200_13TeV.root";
			break;
		case LIGHT_STOP:
			output += "lightstopmod_13TeV.root";
			break;
		case LIGHT_STAU:
			output += "lightstau1_13TeV.root";
			break;
		case HMSSM:
			output += "hMSSM_13TeV.root";
			break;
		case TAU_PHOBIC:
			output += "tauphobic_13TeV.root";
			break;
		default: 
			break;
}
	return output;
}
