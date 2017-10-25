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
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

HbbStyle style;

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

std::string getBanchmarkPath(AvailableScenarios scenario);
void setFolderNamesAccordingToBlinding(string& input, string & output);

int main(){

	style.setTDRstyle(PRELIMINARY);
	//Prefix to the output
	string output_prefix = "solo_13TeV_limits";
	//paths with results
	string path2016_solo 	 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/";
	string path2016_combined = cmsswBase + "/src/Analysis/MssmHbb/datacards/201705/15/Asymptotic/mssm/No_Bias/Hbb.limits";	//If combination is performed
	//ouptut folder
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170823/";
	setFolderNamesAccordingToBlinding(path2016_solo,output);
	string boson = "both";
	std::vector<AvailableScenarios> scenarious = {MHMODP_200,LIGHT_STOP,LIGHT_STAU,HMSSM,TAU_PHOBIC};//{MHMODP_200};//
	for(const auto& scenario : scenarious){
		std::cout<<"Scenario: "<<AvailableScenariosToString(scenario)<<std::endl;
		string output_mod 			= output + AvailableScenariosToString(scenario) + "/";
		string path2016_solo_mod 	= path2016_solo + AvailableScenariosToString(scenario) + "/bias/Hbb.limits";
		CheckOutputDir(output);
		//benchmark scenario path
		LHCXSGLimits limits(mssmhbb::blinded,boson,getBanchmarkPath(scenario));
		limits.compareWithPrevious("CMS-AN-2013-229v10");
		limits.setScenario(scenario);
        limits.setXMax(1000);
//        limits.setXMin(90);
		if(scenario != HMSSM)limits.drawExcludedRegion(3);	// Region NOT compatible with H(125)
		else limits.drawExcludedRegion(0);
	
		string benchmark_name = AvailableScenariosToString(scenario);
		string output_name = output_mod + boson + "_" + benchmark_name + "_" + output_prefix;

	// combined results UNCOMMENT if combination is performed
//	vector<Limit> GxBR_combined = limits.ReadCombineLimits(path2016_combined);
//	string benchmark_ref = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/signal/newmhmax_mu200_13TeV.root";
//	vector<Limit> mssm_limits_combined = limits.GetMSSMLimits(GxBR_combined,benchmark,"",false,benchmark_ref,30);
//	limits.Write(mssm_limits_combined, "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170309/combined_limits.txt");

		//solo 2016 13 TeV
		limits.ReadCombineLimits(path2016_solo_mod);
		limits.Write(output_name + ".txt");

		limits.LimitPlotter(output_name,"24.6 vs 35.7(2016)",string(HbbStyle::axisTitleMA().Data()),"tan#beta",false);
	}
	

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

void setFolderNamesAccordingToBlinding(string& input, string & output){
	if(!mssmhbb::blinded) {
		input += "unblinded/";
		output += "unblinded/";
	}
	else {
		input += "blinded/";
		output += "blinded/";
	}
}
