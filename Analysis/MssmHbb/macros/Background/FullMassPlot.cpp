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
#include "TFile.h"
#include "TColor.h"

#include <string>
#include "RooDataHist.h"

#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "../Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace analysis;
using namespace RooFit;

void setup_plotting(const int& nbins, const std::map<std::string,int>& signals,const std::map<int,double>& signal_sfs,const std::string& output_folder);

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

int main(){

	int nbins = 100;
	string output_folder = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Background/";
	map<string,int> signals = {
			{"sr1" , 300},
			{"sr2" , 600},
			{"sr3" , 1300}
	};
	map<int,double> signal_sfs = {
			{300 , 0.02},
			{600 , 0.02},
			{1300, 0.02}
	};

	setup_plotting(nbins,signals,signal_sfs,output_folder);

	return 0;
}

void setup_plotting(const int& nbins, const std::map<std::string,int>& signals,const std::map<int,double>& signal_sfs,const std::string& output_folder){
	/*
	 *
	 */
	HbbStyle style;
	style.set(PRELIMINARY);

	TFile fbg_sr1(mssmhbb::path_bg_sr1.c_str(),"read");
	TFile fbg_sr2(mssmhbb::path_bg_sr2.c_str(),"read");
	TFile fbg_sr3(mssmhbb::path_bg_sr3.c_str(),"read");

	auto &bg_pdf1 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr1,"background");
	auto &bg_pdf2 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr2,"background");
	auto &bg_pdf3 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr3,"background");

	auto &data1 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr1,"data_container");
	auto &data2 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr2,"data_container");
	auto &data3 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr3,"data_container");

	auto &x1 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr1,"mbb");
	auto &x2 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr2,"mbb");
	auto &x3 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr3,"mbb");

	RooRealVar x("mbb","mbb",200,1700);
	x.setRange("sr1",200,650);
	x.setRange("sr2",350,1190);
	x.setRange("sr3",500,1700);

	TCanvas can("can","can",800,600);
	RooPlot *frame_sr1 = x1.frame(200,650);
	RooPlot *frame_sr2 = x2.frame(350,1190);
	RooPlot *frame_sr3 = x3.frame(500,1700);
	RooPlot *frame = x.frame();
	data1.plotOn(frame,Binning(150),Name("data"));
	bg_pdf1.plotOn(frame_sr1,LineColor(kRed),Normalization(data1.sumEntries("1","fit_range")),Name("bg_pdf1"));
	bg_pdf2.plotOn(frame_sr2,LineColor(kBlue),Normalization(data2.sumEntries("1","fit_range")/2),Name("bg_pdf2"));
	bg_pdf3.plotOn(frame_sr3,LineColor(kGreen),Normalization(data3.sumEntries("1","fit_range")/2.5,RooAbsReal::NumEvent),Name("bg_pdf3"));

	TLegend leg(0.6,0.6,0.9,0.9);
	leg.AddEntry(frame->findObject("data"),"Data","p");
	leg.AddEntry(frame_sr1->findObject("bg_pdf1"),"SuperNovosibirsk-1, sub-range1","l");
	leg.AddEntry(frame_sr2->findObject("bg_pdf2"),"SuperNovosibirsk-0, sub-range2","l");
	leg.AddEntry(frame_sr3->findObject("bg_pdf3"),"SuperNovosibirsk-0, sub-range3","l");
	leg.SetBorderSize(0);
	leg.SetFillColor(0);
	//Draw signals
	map<string,int> colours = {
			{"sr1",2},
			{"sr2",4},
			{"sr3",3}
	};
	for(const auto& s : signals){
		TFile fsg(mssmhbb::signal_fits.at(s.second).c_str(),"read");
		auto &sg_pdf = *GetRooObjectFromTFile<RooAbsPdf>(fsg,"signal");
		sg_pdf.plotOn(frame,LineColor(colours.at(s.first)),LineStyle(2),Normalization(signal_sfs.at(s.second)),Range(s.first.c_str()),Name( ("sg_" + to_string(s.second)).c_str() ));
		leg.AddEntry(frame->findObject(("sg_" + to_string(s.second)).c_str()),("bb#phi(" + to_string(s.second) + " GeV)").c_str(),"l");
	}
	frame->SetMinimum(1);
	frame->SetTitle(";M_{12} [GeV]");
	frame->Draw();
	frame_sr1->Draw("same");
	frame_sr2->Draw("same");
	frame_sr3->Draw("same");

	leg.Draw();

	style.drawStandardTitle();

	can.Print( (output_folder + "plot.pdf").c_str()	);
	gPad->SetLogy();
	frame->SetMaximum(5*1e+05);
	can.Print( (output_folder + "plot_log.pdf").c_str()	);

}

//void compare_srs(RooWorkspace & w1,
//		RooWorkspace & w2,
//		const double& overlap_lo,
//		const double& overlap_hi,
//		const string output_name,
//		const double& max,
//		const string legend_pdf1,
//		const string legend_pdf2,
//		const string legend_data,
//		const string function1,
//		const string function2){
//	auto &pdf1 	= *analysis::GetFromRooWorkspace<RooAbsPdf>(w1,function1);
//	auto &pdf2 	= *analysis::GetFromRooWorkspace<RooAbsPdf>(w2,function2);
//	auto &data1 	= *analysis::GetFromRooWorkspace<RooDataHist>(w1,"data_container");
//	auto &data2 	= *analysis::GetFromRooWorkspace<RooDataHist>(w2,"data_container");
//	auto &x		= *analysis::GetFromRooWorkspace<RooRealVar>(w1,"mbb");
//	auto &x2		= *analysis::GetFromRooWorkspace<RooRealVar>(w2,"mbb");
//
//	//Make frame
//	RooPlot &frame = *x.frame(overlap_lo,overlap_hi);
//	RooPlot &frame2 = *x2.frame(overlap_lo,overlap_hi);
//	frame.SetTitle(legend_data.c_str());
//	int nbins2 = (x2.getMax() - x2.getMin()) / 10;
//	int nbins1 = (x.getMax() - x.getMin()) / 10;
//	data1.plotOn(&frame,RooFit::MarkerSize(0.8),RooFit::LineColor(kRed),RooFit::Name("data_points"),RooFit::Binning(nbins1));
//	data2.plotOn(&frame2,RooFit::MarkerSize(0.8),RooFit::LineColor(kBlue),RooFit::Name("data_points2"),RooFit::Binning(nbins2));
//	pdf1.plotOn(&frame,RooFit::LineColor(kRed),RooFit::Name("pdf1"),RooFit::Normalization(data1.sumEntries("1","fit_range"),RooAbsReal::NumEvent));
//	pdf2.plotOn(&frame2,RooFit::LineColor(kBlue),RooFit::Name("pdf2"),RooFit::Normalization(data2.sumEntries("1","fit_range"),RooAbsReal::NumEvent));
//
//	TLegend leg(0.5,0.5,0.85,0.85);
//	leg.SetBorderSize(0);
//	leg.SetFillColor(0);
//	leg.SetFillStyle(0);
//	leg.AddEntry(frame.findObject("data_points"),legend_data.c_str(),"p");
//	leg.AddEntry(frame.findObject("pdf1"),legend_pdf1.c_str(),"l");
//	leg.AddEntry(frame2.findObject("pdf2"),legend_pdf2.c_str(),"l");
//	TCanvas can("can","can",800,600);
//	frame.SetMinimum(1);
//	frame.SetMaximum(max);
//	frame.Draw();
//	frame2.Draw("same");
//	leg.Draw();
//	can.Print((output_name + ".pdf").c_str());
//	gPad->SetLogy();
//	leg.SetY1NDC(0.2);
//	leg.SetY2NDC(0.5);
//	can.Print((output_name + "_log.pdf").c_str());
//}
