/*
 * mass_overlap_shape.cpp
 *
 *  Created on: 6 Jul 2017
 *      Author: shevchen
 *      description: Macro to check
 *      the shape of a background functions
 *      in the overlap region between different
 *      sub-ranges
 *
 */


#include "TCanvas.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

void compare_srs(RooWorkspace & w1,
		RooWorkspace & w2,
		const double& overlap_lo,
		const double& overlap_hi,
		const string output_name,
		const double& y_max,
		const string legend_pdf1,
		const string legend_pdf2,
		const string legend_data,
		const string function1 = "background",
		const string function2 = "background");

using namespace std;

int mass_overlap_shape(){

	TFile *fsr1 = new TFile("../../output/ReReco_bg_fit/sr1/FitContainer_workspace_SR1.root");
	TFile *fsr2 = new TFile("../../output/ReReco_bg_fit/sr2/FitContainer_workspace_SR2.root");
	TFile *fsr3 = new TFile("../../output/ReReco_bg_fit/sr3/FitContainer_workspace_SR3.root");

	auto &wsr1 = *GetFromTFile<RooWorkspace>(*fsr1,"workspace");
	auto &wsr2 = *GetFromTFile<RooWorkspace>(*fsr2,"workspace");
	auto &wsr3 = *GetFromTFile<RooWorkspace>(*fsr3,"workspace");
	compare_srs(wsr1,wsr2,350,650,"sr1_sr2",20000,"SuperNovosibirsk-1*Turn-on, sr1","Novosibirsk, sr2","data, sr1 #cap sr2");
	compare_srs(wsr2,wsr3,500,1190,"sr2_sr3",3000,"Novosibirsk, sr2","Novosibirsk, sr3","data, sr2 #cap sr3");

	return 0;
}

void compare_srs(RooWorkspace & w1,
		RooWorkspace & w2,
		const double& overlap_lo,
		const double& overlap_hi,
		const string output_name,
		const double& max,
		const string legend_pdf1,
		const string legend_pdf2,
		const string legend_data,
		const string function1,
		const string function2){
	auto &pdf1 	= *analysis::GetFromRooWorkspace<RooAbsPdf>(w1,function1);
	auto &pdf2 	= *analysis::GetFromRooWorkspace<RooAbsPdf>(w2,function2);
	auto &data1 	= *analysis::GetFromRooWorkspace<RooDataHist>(w1,"data_container");
	auto &data2 	= *analysis::GetFromRooWorkspace<RooDataHist>(w2,"data_container");
	auto &x		= *analysis::GetFromRooWorkspace<RooRealVar>(w1,"mbb");
	auto &x2		= *analysis::GetFromRooWorkspace<RooRealVar>(w2,"mbb");

	//Make frame
	RooPlot &frame = *x.frame(overlap_lo,overlap_hi);
	RooPlot &frame2 = *x2.frame(overlap_lo,overlap_hi);
	frame.SetTitle(legend_data.c_str());
	int nbins2 = (x2.getMax() - x2.getMin()) / 10;
	int nbins1 = (x.getMax() - x.getMin()) / 10;
	data1.plotOn(&frame,RooFit::MarkerSize(0.8),RooFit::LineColor(kRed),RooFit::Name("data_points"),RooFit::Binning(nbins1));
	data2.plotOn(&frame2,RooFit::MarkerSize(0.8),RooFit::LineColor(kBlue),RooFit::Name("data_points2"),RooFit::Binning(nbins2));
	pdf1.plotOn(&frame,RooFit::LineColor(kRed),RooFit::Name("pdf1"),RooFit::Normalization(data1.sumEntries("1","fit_range"),RooAbsReal::NumEvent));
	pdf2.plotOn(&frame2,RooFit::LineColor(kBlue),RooFit::Name("pdf2"),RooFit::Normalization(data2.sumEntries("1","fit_range"),RooAbsReal::NumEvent));

	TLegend leg(0.5,0.5,0.85,0.85);
	leg.SetBorderSize(0);
	leg.SetFillColor(0);
	leg.SetFillStyle(0);
	leg.AddEntry(frame.findObject("data_points"),legend_data.c_str(),"p");
	leg.AddEntry(frame.findObject("pdf1"),legend_pdf1.c_str(),"l");
	leg.AddEntry(frame2.findObject("pdf2"),legend_pdf2.c_str(),"l");
	TCanvas can("can","can",800,600);
	frame.SetMinimum(1);
	frame.SetMaximum(max);
	frame.Draw();
	frame2.Draw("same");
	leg.Draw();
	can.Print((output_name + ".pdf").c_str());
	gPad->SetLogy();
	leg.SetY1NDC(0.2);
	leg.SetY2NDC(0.5);
	can.Print((output_name + "_log.pdf").c_str());
}
