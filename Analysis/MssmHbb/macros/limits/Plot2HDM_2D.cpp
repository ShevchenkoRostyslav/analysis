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

#include "Analysis/MssmHbb/interface/HbbLimits.h"
#include "Analysis/MssmHbb/interface/Limit.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"

HbbStyle style;

using namespace std;
using namespace analysis::mssmhbb;

void Draw2HDM_2D_Limits(const std::string& output, HbbLimits::THDMScan& scan);

int main(){

	HbbLimits limits(true,true);
	style.set(PRELIMINARY);
	//paths with results of the combine tool
	string path2016 = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/datacards/201705/18/Asymptotic/mssm/Hbb.limits";
	//Details of the 2HDM produciton
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION";
	// type of the 2hdm: type1/2/3/4
	string thdm_type = "type4";
	string thdm_scans = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + thdm_type + "_mA_mH.root";
	//value of mA
	double mass =300;

	//higgs boson: H/A/both
	string boson = "both";
	string output_dir = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170518/2hdm/";

	vector<Limit> GBR2016 = limits.ReadCombineLimits(path2016);
	HbbLimits::THDMScan scan;
	TH3D GxBR_2hdm_3D;
	TH2D GxBR_2hdm_mA;
	limits.SetHiggsBoson(boson);
	vector<Limit> THDM_limits;
	for(const auto& l : GBR2016){
		cout<<"M: "<<l.getX()<<" exp = "<<l.getMedian()<<" 1G = "<<l.getPlus1G()<<endl;
		if(l.getX() == mass){
			GxBR_2hdm_3D = limits.Get2HDM_GxBR_3D(thdm_scans);
			GxBR_2hdm_mA = limits.Get2HDM_GxBR_2D(GxBR_2hdm_3D,mass,"x");
			scan = limits.Get2HDMmuScan(GxBR_2hdm_mA,l);
		}
	}
	string output_2hdm_2D = output_dir + boson + "_2HDM_" + thdm_type + "_Limits2D_mA" + to_string((int)mass);
	Draw2HDM_2D_Limits(output_2hdm_2D,scan);

}

void Draw2HDM_2D_Limits(const std::string& output, HbbLimits::THDMScan& scan){
	/*
	 * Draw 2D 2HDM limits
	 */
	TCanvas can("can","can",800,600);
	gPad->SetLogz();
	gPad->SetLogy();
	double mu = 1;
	scan.expected.GetXaxis()->SetRangeUser(-1.,1.);
	scan.expected.SetStats(kFALSE);
	scan.expected.DrawCopy("COLZ");
	scan.expected.SetContour(1,&mu);
	scan.expected.SetLineColor(2);
	scan.expected.SetLineStyle(2);
	scan.expected.DrawCopy("cont3 same");
	can.Print( (output + ".pdf").c_str());
}
