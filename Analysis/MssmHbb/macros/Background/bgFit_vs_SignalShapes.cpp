/*
 * bgFit_vs_SignalShapes.cpp
 *
 *  Created on: 30 May 2017
 *      Author: shevchen
 *
 *  Macro is used to plot background fits to the data
 *  in each sub-range together with a corresponded
 *  mass points
 */
#include <string>
#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "TAxis.h"

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace mssmhbb;
using namespace analysis;

void Draw(const vector<int>& signal_sr1, const vector<int>& signal_sr2, const vector<int>& signal_sr3, const string& output_path, const int& draw_single_mass = -1);
void Draw(const vector<int>& signal_points,const string& output_path, const int& sub_range);
void Draw(const int& signal_point,const string& output_path, const int& sub_range);
int DefineSubRange(const int& signal_point,const vector<int>& signal_sr1, const vector<int>& signal_sr2, const vector<int>& signal_sr3);
int GetSignalXSection(const int& signal_point);

int main(){

	//Output path
	string output_path = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Background/";

	//If want to draw single mass point instead of all
	int single_mass = -1;

	//List of the signal mass points in each sub-range
	vector<int> signal_sr1 = {300,350,400,500,600};
	vector<int> signal_sr2 = {500,600,700,900,1100};
	vector<int> signal_sr3 = {900,1100,1300};

	Draw(signal_sr1,signal_sr2,signal_sr3,output_path,single_mass);

	return 0;
}

void Draw(const vector<int>& signal_sr1, const vector<int>& signal_sr2, const vector<int>& signal_sr3, const string& output_path, const int& draw_single_mass){
	/*
	 * function to draw bgf + signal shapes
	 */
	if(draw_single_mass != -1) {
		int sub_range = DefineSubRange(draw_single_mass,signal_sr1,signal_sr2,signal_sr3);
		Draw(draw_single_mass, output_path,sub_range);
	}
	else {
		Draw(signal_sr1,output_path,1);
		Draw(signal_sr2,output_path,2);
		Draw(signal_sr3,output_path,3);
	}
}

int DefineSubRange(const int& signal_point,const vector<int>& signal_sr1, const vector<int>& signal_sr2, const vector<int>& signal_sr3){
	/*
	 * function to define to which sub-range single point corresponds to
	 */
	if( find(signal_sr1.begin(),signal_sr1.end(),signal_point) != signal_sr1.end() ) return 1;
	else if( find(signal_sr2.begin(),signal_sr2.end(),signal_point) != signal_sr2.end() ) return 2;
	else if( find(signal_sr3.begin(),signal_sr3.end(),signal_point) != signal_sr3.end() ) return 3;
	else {
		throw invalid_argument("bgFit_vs_SignalShapes::DefineSubRange - invalid mass point, not in the sub-ranges: m = " + to_string(signal_point));
	}
}

void Draw(const vector<int>& signal_points,const string& output_path, const int& sub_range){
	/*
	 * Overloaded function to draw bg vs signal for the vector of signal points
	 */
	for(const auto& signal_point : signal_points){
		Draw(signal_point, output_path,sub_range);
	}
}

void Draw(const int& signal_point,const string& output_path, const int& sub_range){
	/*
	 * Overloaded function to draw bg vs signal shape
	 */

	HbbStyle style;
	style.set(PRELIMINARY);

	CheckOutputDir(output_path);
	TFile fBg(path_bg_fits.at(sub_range).c_str(),"READ");
	TFile fSg(signal_workspaces.at(signal_point).c_str(),"READ");

	auto& wBg = *GetFromTFile<RooWorkspace>(fBg, "workspace");
	auto& wSg = *GetFromTFile<RooWorkspace>(fSg, "workspace");

	auto& x = *GetFromRooWorkspace<RooRealVar>(wBg,"mbb");
	auto& data_obs 	= *GetFromRooWorkspace<RooDataHist>(wBg,"data_obs");

	auto& pdf_sgn = *GetFromRooWorkspace<RooAbsPdf>(wSg,"signal");
	auto& pdf_bkg = *GetFromRooWorkspace<RooAbsPdf>(wBg,"background");

	TCanvas c1("c1", "canvas", 600, 600);
    auto& frame1 = *x.frame(RooFit::Name("frame1"));
    frame1.SetTitle("");
    frame1.GetYaxis()->SetTitle(( ("Events / " + to_string_with_precision(x.getBinWidth(1),2) + " GeV").c_str() ));
    frame1.GetYaxis()->SetTitleSize(24);
    frame1.GetYaxis()->SetTitleFont(43);
    frame1.GetYaxis()->SetTitleOffset(1.18);
    frame1.GetYaxis()->SetLabelFont(43);
    frame1.GetYaxis()->SetLabelSize(20);
    frame1.GetXaxis()->SetTitle("m_{b#bar{b}} (GeV)");
    frame1.GetXaxis()->SetTitleSize(25);
    frame1.GetXaxis()->SetTitleFont(43);
//    frame1.GetXaxis()->SetTitleOffset(1.5);
    frame1.GetXaxis()->SetLabelFont(43);
    frame1.GetXaxis()->SetLabelSize(20);

    TLegend leg(0.65,0.7,0.9,0.92,"","brNDC");
    leg.SetHeader( ("m_{#phi} = " + to_string(signal_point) + " GeV").c_str() );
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.02);
    leg.SetFillColor(10);

    data_obs.plotOn(&frame1, RooFit::Name("data_obs"));
    pdf_bkg.plotOn(&frame1,
    		RooFit::LineWidth(3),
    		RooFit::LineColor(kBlue),
			RooFit::LineStyle(kSolid),
			RooFit::Name(pdf_bkg.GetName()),
			RooFit::Normalization(data_obs.sumEntries(), RooAbsReal::NumEvent)
    );

    int visible_xsec = GetSignalXSection(signal_point);
    pdf_sgn.plotOn(&frame1,
    		RooFit::LineWidth(3),
			RooFit::LineColor(46),
			RooFit::LineStyle(kDotted),
			RooFit::Name(pdf_sgn.GetName()),
			RooFit::Normalization(visible_xsec, RooAbsReal::NumEvent)
    );
    data_obs.plotOn(&frame1, RooFit::Name("data_obs"));

    leg.AddEntry(frame1.findObject("data_obs"),"Data","PE");
    leg.AddEntry(frame1.getCurve(pdf_bkg.GetName()),"Bkg.","L");
    leg.AddEntry(frame1.getCurve(pdf_sgn.GetName()), ("Signal, #sigma= " + to_string(visible_xsec) + " pb").c_str(),"L");
    frame1.SetMaximum(frame1.GetMaximum()*2.5);
    frame1.SetMinimum(1);
    gPad->SetLogy();
    frame1.Draw();
    leg.Draw();

    c1.Modified();

    c1.Print( (output_path + "M" + to_string(signal_point) + "_sr" + to_string(sub_range) + ".pdf").c_str() );




}

int GetSignalXSection(const int& signal_point){
	/*
	 * function to return crossection
	 * that allows visulaze signal better
	 */
	switch (signal_point){
	case 300: return 1000;
	case 350: return 800;
	case 400: return 700;
	case 500: return 500;
	case 600: return 400;
	case 700: return 300;
	case 900: return 200;
	case 1100: return 250;
	case 1300: return 250;
	default: return 1;
	}
}
