/*
 * plotTanBetaLimits.cpp
 *
 *  Created on: 10 Mar 2017
 *  Author: shevchen
 *
 *  Macro to calculate and plot 2D 2HDM limits
 *  for particular mA
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

//#include "Analysis/MssmHbb/interface/HbbLimits.h"
//#include "Analysis/MssmHbb/interface/Limit.h"
//#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/THDMLimits.h"
#include "Analysis/MssmHbb/interface/Limit.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

HbbStyle style;

using namespace std;
using namespace analysis::mssmhbb;
using namespace mssmhbb;

void Draw2HDM_2D_Limits(const std::string& output, THDMLimits::THDMScan& scan);

int main(int argc,const char ** argv){

	style.setTDRstyle(PRELIMINARY);
	gStyle->SetPadRightMargin(0.15);
	gStyle->SetPadTopMargin   (0.09);
	//paths with results of the combine tool
	string path2016 = cmsswBase + "/src/Analysis/MssmHbb/datacards/201712/13/unblinded/mhmodp_200/bias/Hbb.limits";//"/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201707/26/blinded/mssm/bias/Hbb.limits";
	//Details of the 2HDM produciton
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION";
	// type of the 2hdm: type1/2/3/4
	string thdm_type = "type2";
	if(argc != 1) {
		thdm_type = string(argv[1]);
	}
	string banch_name = thdm_type;
	if(thdm_type == "flipped") banch_name = "type3";
	else if (thdm_type == "lepton_specific") banch_name = "type4";
	string thdm_scans = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + banch_name + "_mA_mH.root";
	AvailableScenarios scenario = AvailableScenariosFromString(thdm_type);
	//value of mA
	double mass =300;

	//higgs boson: H/A/both
	string boson = "both";
	string output = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20171213/";
	if(!mssmhbb::blinded) output += "unblinded/";
	output += "2hdm/" + thdm_type + "/";
	CheckOutputDir(output);

	THDMLimits limits(mssmhbb::blinded,boson,-1,1,0,100);
	limits.setScenario(scenario);
	limits.SetHiggsBoson(boson);
	limits.ReadCombineLimits(path2016);
	vector<Limit> GBR2016 = limits.getGxBrLimits();
	THDMLimits::THDMScan scan;
	TH3D GxBR_2hdm_3D;
	TH2D GxBR_2hdm_mA;
	vector<Limit> THDM_limits;
	for(const auto& l : GBR2016){
		cout<<"M: "<<l.getX()<<" exp = "<<l.getMedian()<<" 1G = "<<l.getPlus1G()<<endl;
		if(l.getX() == mass){
			GxBR_2hdm_3D = limits.Get2HDM_GxBR_3D(thdm_scans);
			GxBR_2hdm_mA = limits.Get2HDM_GxBR_2D(GxBR_2hdm_3D,mass,"x");
			scan = limits.Get2HDMmuScan(GxBR_2hdm_mA,l);
		}
	}
	string output_2hdm_2D = output + thdm_type + "_" + boson + "_bosons_mA" + to_string((int)mass) + "_2D";
	Draw2HDM_2D_Limits(output_2hdm_2D,scan);

}

void Draw2HDM_2D_Limits(const std::string& output, THDMLimits::THDMScan& scan){
	/*
	 * Draw 2D 2HDM limits
	 */
	TCanvas can("can","can",800,600);
	gPad->SetLogz();
	gPad->SetLogy();
	double mu = 1;
	scan.expected.SetTitle("(#sigma x B)_{theory} / (#sigma x B)_{95#%C.L.})");
	scan.expected.GetXaxis()->SetRangeUser(-1.,1.);
	scan.expected.SetStats(kFALSE);
	scan.expected.DrawCopy("COLZ");
	scan.expected.SetContour(1,&mu);
	scan.expected.SetLineColor(2);
	scan.expected.SetLineStyle(2);
	scan.expected.DrawCopy("cont3 same");
	style.drawStandardTitle();
	can.Print( (output + ".pdf").c_str());
}
