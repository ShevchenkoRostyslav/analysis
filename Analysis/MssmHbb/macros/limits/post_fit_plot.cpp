#include <string>

#include <TFile.h>
#include <TH1F.h>

#include <RooFit.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooPlot.h>

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"

using namespace std;
using namespace RooFit;

HbbStyle  style;
const string cmsswBase = gSystem->Getenv("CMSSW_BASE");

TH1* getDataHist(RooWorkspace& w, TH1& binning);
void draw(TH1F* data, TH1F *total, TH1F *signal, TH1F *bg);

//int main(int argc, char **argv) {
int post_fit_plot(){
	/*
	 * Macro to plot post-fit M_{12} fit
	 */
	style.set(PRIVATE);
	string mass_point = "500";

	TFile *fIn = new TFile( (cmsswBase + "/src/Analysis/MssmHbb/datacards/201705/18/mll/independent/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	//201702/13/1xBins/mlfit_1.root
	TFile *fBg = new TFile( (cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_bg_fit/sr1/background_workspace_TurnOnFix_5000bins.root").c_str(),"READ");

	//Get Histograms
	TH1F *bg 		= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total_background");	// post-fit BG
	TH1F *signal 	= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total_signal"); 	// post-fit Signal
	TH1F *total 	= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total");			// post-fit Total
	TH1F *data 		= (TH1F*) getDataHist(*GetFromTFile<RooWorkspace>(*fBg,"workspace"),*bg);

	//Rebin histograms:
	total->Scale(total->GetBinWidth(10));
	signal->Scale(signal->GetBinWidth(10));
	bg->Scale(9);
//	bg->Scale(bg->GetBinWidth(10));
	total->Rebin(100);
	signal->Rebin(100);
	data->Rebin(100);
//	bg->Rebin(100);

	draw(data,total,signal,bg);

//	bg->Draw("E");
//	RooWorkspace *wBg = (RooWorkspace*) fBg->Get("workspace");
//	RooAbsPdf* dpf_bg = wBg->pdf("background");
//	RooPlot *frame = wBg->var("mbb")->frame();
//	RooRealVar *mbb = wBg->var("mbb");
//	mbb->setBins(50);
//	RooDataHist *hBg = new RooDataHist("hBg","hBg",*mbb,bg);
//	hBg->plotOn(frame,Bins(50));
//	dpf_bg->plotOn(frame,LineColor(kBlue));
//	frame->Draw();
//	bg->Draw("Esame");
//	cout<<"Wtf:"<< bg->Integral()<<" "<<wBg->var("background_norm")->getValV()<<endl;

//	total-
//
//	total->Draw("E");
//	signal->Scale(10);
//	signal->Draw("sameh");




//	for(int i = 0; i < total->GetNbinsX(); ++i){
//
//	}

	/*
	//Get RooPlot obtained by combine tool
	RooPlot *postfit = (RooPlot*) fIn->Get("bbHTo4b_fit_s");
	//Get histo
	RooHist *h = (RooHist*) postfit->getHist("h_bbHTo4b");
	//Get curve
	RooCurve *signal = (RooCurve*) postfit->getCurve("pdf_binbbHTo4b_Norm[mbb]_Comp[shapeSig*]");
	//Scale signal by factor of 10
	const double scale_factor = 20;
	for(int i = 0; i < signal->GetN();++i) signal->GetY()[i] *= scale_factor;
	postfit->Draw();
	*/


	return 0;
}

void draw(TH1F* data, TH1F *total, TH1F *signal, TH1F *bg){
	/*
	 * function to draw post-fit plot from histograms
	 */
	const double scale_factor = 200;
	TCanvas *can = new TCanvas("can","can",800,600);
//	style.InitHist(total,"M_{12}, [GeV]","Events",1,0);
//	style.InitData(total);

	TPad *pad1 = new TPad("pad1","",0,0.1,1,1);
    pad1->SetBottomMargin(0.2);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.16);
	pad1->Draw();
	pad1->cd();

	double xMin = total->GetXaxis()->GetXmin();
	double xMax = total->GetXaxis()->GetXmax();
	double yMin = 0;//total->GetMinimum();
	double yMax = 1.2 * total->GetMaximum();

	TH2F *frame = new TH2F("frame","",2,xMin,xMax,2,yMin,yMax);
	frame->GetXaxis()->SetTitle(total->GetXaxis()->GetTitle());
	frame->GetYaxis()->SetTitle(total->GetYaxis()->GetTitle());
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetNdivisions(206);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetXaxis()->SetTitleOffset(999);
	frame->GetXaxis()->SetLabelOffset(999);
	frame->GetYaxis()->SetTitleSize(0.048);
	frame->SetStats(0);
	frame->Draw();

	total->Draw("Esame");
	signal->Scale(scale_factor);
	signal->SetFillColor(2);
	signal->SetLineColor(2);

//	bg->SetFillColor(kBlue-4);
	bg->SetLineWidth( (Width_t) 1.1);
	bg->SetLineColor(kBlue-4);
	bg->Draw("hsame");
	signal->Draw("hsame");
	total->Draw("Esame");

	TLegend *leg = new TLegend(0.6,0.6,0.9,0.8);
	leg->SetBorderSize(0);
	leg->AddEntry(total,"Data","p");
	leg->AddEntry(signal,"bb#phi(500 GeV)x20","fl");
	leg->AddEntry(bg,"Background","fl");
	leg->Draw();

	//Pad2 with ratios
    can->cd();
    TPad *pad2 = new TPad("pad2","",0,0.0,1,0.265);
    pad2->SetTopMargin(1);
    pad2->SetBottomMargin(0.33);
    pad2->SetLeftMargin(pad1->GetLeftMargin());
    pad2->SetRightMargin(pad1->GetRightMargin());
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH2F *frame2 = new TH2F("frame2","",2,xMin,xMax,2,0.93,1.07);
    frame2->SetTitle("");
    frame2->GetXaxis()->SetTitle(frame->GetXaxis()->GetTitle());
    frame2->GetXaxis()->SetTitleSize(0.15);
    frame2->GetXaxis()->SetTitleOffset(1.06);
    frame2->GetXaxis()->SetLabelSize(0.215);
    frame2->GetXaxis()->SetLabelOffset(0.010);
    frame2->GetXaxis()->SetNdivisions(505);
    frame2->GetXaxis()->SetTickLength(frame->GetXaxis()->GetTickLength()*3);
    frame2->SetYTitle("Data/Bkg");
    frame2->GetYaxis()->CenterTitle(kTRUE);
    frame2->GetYaxis()->SetTitleSize(0.14);
    frame2->GetYaxis()->SetTitleOffset(0.4);
    frame2->GetYaxis()->SetNdivisions(206);
    frame2->GetYaxis()->SetLabelSize(0.115);
    frame2->GetYaxis()->SetLabelOffset(0.011);
    frame2->Draw();

    TH1F *ratio_tot = (TH1F*) total->Clone("ratio_tot");
    ratio_tot->Divide(bg);
    ratio_tot->Draw("same");

    TH1F *ratio_sg = (TH1F*) signal->Clone("ratio_sg");
    ratio_sg->Scale(1./20.);
    ratio_sg->Divide(bg);
    for(int i =0;i<ratio_sg->GetNbinsX();++i) ratio_sg->SetBinContent(i+1,ratio_sg->GetBinContent(i+1)+1);
    ratio_sg->SetFillColor(0);
    ratio_sg->Draw("hsame");

    can->Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/post_fit_500GeV.pdf").c_str() );
}


TH1* getDataHist(RooWorkspace& w, TH1& binning){
	auto &dataHist 		= *analysis::GetFromRooWorkspace<RooDataHist>(w,"data_obs");
	auto &mbb			= *analysis::GetFromRooWorkspace<RooRealVar>(w,"mbb");
	auto *data 			=  static_cast<TH1F*>(dataHist.createHistogram("data_obs",mbb,RooFit::Binning(binning.GetNbinsX(), binning.GetXaxis()->GetXmin(), binning.GetXaxis()->GetXmax())));
	return data;
}
