#include "TExec.h"
#include "TF1.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;

void DrawCorrectionFunction(const string& bbbFilePath, const string& bbnbFilePath, const std::string& plot_name = "general/diJet_m");
double fitfunction(double *x, double *par);

int main(){
	HbbStyle style;
	style.setTDRstyle(PRIVATE);

	string bbbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/DataMC_3b_ReReco_35673fb_lowM_QCD.root";
	string bbnbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/bbx_ReReco_35673fb_lowM_QCD.root";
	DrawCorrectionFunction(bbbFilePath,bbnbFilePath);
}

void DrawCorrectionFunction(const string& bbbFilePath, const string& bbnbFilePath, const string& plot_name){
	/*
	 * Function to draw bbb/bbnb correction function
	 */
	gStyle->SetErrorX(0);
	TH1::SetDefaultSumw2();
	TFile fBBB(bbbFilePath.c_str());
	TFile fBBnB(bbnbFilePath.c_str());

	//Histograms TH1
	auto hBBB 		= GetFromTFile<TH1D>(fBBB,plot_name);
	HbbStyle::applyStandardDataHistoStyle(*hBBB);
	hBBB->Scale(1./hBBB->Integral());
	auto hBBnB 		= GetFromTFile<TH1D>(fBBnB,plot_name);
	HbbStyle::applyStandardMCHistoStyle(*hBBnB);
	hBBnB->Scale(1./hBBnB->Integral());

	TCanvas can;
	TH2D *frame = new TH2D("frame","",1,60,420,1,0,2);
	HbbStyle::setFrameStyle(frame);
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMass());
	frame->GetYaxis()->SetTitle("N^{bbb}/N^{bbnb}");
	auto ratio = static_cast<TH1D*>(hBBB->Clone("ratio"));
	ratio->Sumw2();
	ratio->Divide(hBBnB);
	frame->Draw();
	ratio->Draw("E same");

	//Fit
	TF1 *fit = new TF1("fit",fitfunction,140,420,2);
	fit->SetLineWidth(3);
	ratio->Fit(fit,"R");

	// TLegend
	auto leg = HbbStyle::legend("right,top",4,0.6);
	leg->SetHeader("#font[12]{bbnb to bbb correction}");
	leg->AddEntry(ratio,"QCD b-Enriched","p");
	leg->AddEntry(fit,"Fit","l");
	leg->Draw();
//	leg->AddEntry(hBBnB_stat,"QCD reverse b-tag stat.","fl");
//	leg->SetY1(leg->GetY1()*0.8);

	HbbStyle::drawStandardTitle();
	can.Print((mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/bbbvsbbnb_correction.pdf").c_str());

}

double fitfunction(double *x, double *par){
	return 1./(1. + exp(-par[0] * (x[0] - par[1])));
}
