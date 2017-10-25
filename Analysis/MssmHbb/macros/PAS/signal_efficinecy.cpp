/*
 * signal_efficinecy.cpp
 *
 *  Created on: 12 Oct 2017
 *      Author: rostyslav
 *
 *
 */
#include "TCanvas.h"
#include "TGraphErrors.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
HbbStyle style;

using namespace std;

void drawSignalEfficiency(const string& output, const float& window);
TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax);
void drawTextBox();

int main(int argc, char** argv){
	style.setTDRstyle(PRELIMINARY_SIMULATION);
	string output_folder 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/";
	string output_name 		= "PAS_signal_efficiency";
	drawSignalEfficiency(output_folder + output_name,0.4);
	return 0;
}

void drawSignalEfficiency(const string& output, const float& window){
//	gStyle->SetPadRightMargin(0.04);
//	gStyle->SetPadLeftMargin(0.14);
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	TH2D *frame = new TH2D("frame","frame",1,200,1399.99,1,0.0001,0.02);
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMAH());
	frame->GetYaxis()->SetTitle("#epsilon");
	frame->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X")*1.02);
	frame->SetMinimum(1e-06);
	frame->GetXaxis()->SetRangeUser(200,1399.99);
	//Create TGraph
	TGraphErrors *gr = new TGraphErrors(mssmhbb::signal_templates.size());
	gr->SetMarkerColor(kRed);
	//Iterate through the signal samples
	int i = 0;
	for(const auto& sample : mssmhbb::signal_templates){
//		double lower_border = sample.first * ( 1 - window/2.);
//		double upper_border = sample.first * ( 1 + window/2.);
		TFile *f = new TFile(sample.second.c_str());
//		auto* h_denum 	= GetFromTFile(sample.second,"distributions/NumberOfGenEvents_afterMHat_rewPU");
		auto* h_num		= GetFromTFile<TH1>(*f,"templates/bbH_Mbb_VIS");
		h_num->Sumw2(true);
		double ntot 		= 35673; //From luminosity
		double e_selected = 0;
		double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
		gr->SetPoint(i,sample.first,selected/ntot);
		gr->SetPointError(i,0,0);
		std::cout<<"M_{A/H} = "<<sample.first<<" eff = "<<selected/ntot<<std::endl;
		++i;
	}
//	gPad->SetLogy();
	frame->Draw();
	gr->Draw("PL same");
	HbbStyle::drawStandardTitle("out");
	drawTextBox();
	can.Print( (output + ".pdf").c_str());
}

TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax){
	TLegend *leg = new TLegend(xmin,ymin,xmax,ymax);
	style.setLegendStyle(leg);
	return leg;
}

void drawTextBox(){
	/*
	 * Draw TLatex textbox with a sign
	 */
	  TLatex latex;
	  latex.SetTextColor(kBlue+2);
	  latex.DrawLatex(1100, 0.018,
	                     ("A/H #rightarrow b#bar{b}"));
//	  latex.DrawLatexNDC(0.56, 0.9,
//	                     ("A/H #rightarrow b#bar{b}"));
}
