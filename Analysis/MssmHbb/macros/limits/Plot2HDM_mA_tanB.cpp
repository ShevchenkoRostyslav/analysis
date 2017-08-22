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

#include "Analysis/MssmHbb/interface/THDMLimits.h"
#include "Analysis/MssmHbb/interface/Limit.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

HbbStyle style;

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

int main(int argc, const char** argv){

	style.set(PRELIMINARY);
	//Prefix to the output
	string output_prefix = "13TeV_limits";
	//paths with results of the combine tool
	string path2016 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201707/26/unblinded/mssm/bias/Hbb.limits";
	//value of cos(beta-alpha)
	double cB_A = 0.1;
	//Details of the 2HDM produciton
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION";
	// type of the 2hdm: type2 or type3
	string thdm_type = "type2";
	if(argc != 1) {
		thdm_type = string(argv[1]);
	}
	string banch_name = "type2";
	if(thdm_type == "flipped") banch_name = "type3";
	else if (thdm_type == "lepton_specific") banch_name = "type4";
	AvailableScenarios scenario = AvailableScenariosFromString(thdm_type);
	string thdm_scans = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + banch_name + "_mA_mH.root";

	//higgs boson: H/A/both
	string boson = "both";
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170818/";
	if(!mssmhbb::blinded) output += "unblinded/";
	output += "2hdm/" + thdm_type + "/";
	CheckOutputDir(output);

	THDMLimits limits(mssmhbb::blinded,boson,200,900,0.5,60);
	limits.setScenario(scenario);
	limits.SetHiggsBoson(boson);
	limits.ReadCombineLimits(path2016);
	limits.Get2HDM_1D_Limits(thdm_scans,cB_A,"z");
	limits.compareWithPrevious(true);

	// 2HDM tanBeta vs mA limits for cos(beta-alpha) = cB_A
	TLegend leg_2HDM_tB(0.55,0.17,0.93,0.5);
	leg_2HDM_tB.SetFillColor(0);
	leg_2HDM_tB.SetTextSize(0.035);
	leg_2HDM_tB.SetBorderSize(0);
	output += boson + "_bosons_tanB_cB_A_" + output_prefix;
	limits.LimitPlotter(leg_2HDM_tB,output,"35.7 fb^{-1}","M_{#Phi} [GeV]","tan(#beta)",true);
	output += "_zoomed";
	leg_2HDM_tB.Clear();
	limits.setXMin(300);
	std::cout<<"WTF"<<std::endl;
	limits.LimitPlotter(leg_2HDM_tB,output,"35.7 fb^{-1}","M_{#Phi} [GeV]","tan(#beta)",true);

}