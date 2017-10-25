/*
 * massShape_bbb_vs_bbnb.cpp
 *
 *  Created on: 14 Oct 2017
 *      Author: rostyslav
 */
#include "TExec.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;

void DrawMCM12shapeComparison(const string& bbbFilePath, const string& bbnbFilePath, const std::string& plot_name = "diJet_b_m/diJet_b_m");
void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator);

int main(){
	HbbStyle style;
	style.setTDRstyle(PRELIMINARY_SIMULATION);

	string bbbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/DataMC_3b_ReReco_35673fb_lowM_QCD.root";
	string bbnbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/bbx_ReReco_35673fb_lowM_QCD.root";
	DrawMCM12shapeComparison(bbbFilePath,bbnbFilePath);
}

void DrawMCM12shapeComparison(const string& bbbFilePath, const string& bbnbFilePath, const string& plot_name){
	/*
	 * Main function to draw the MC M12 comparison plot.
	 */
	gStyle->SetErrorX(0);
	TH1::SetDefaultSumw2();
	TCanvas can;
	TFile fBBB(bbbFilePath.c_str());
	TFile fBBnB(bbnbFilePath.c_str());

	//Histograms TH1
	auto hBBB 		= GetFromTFile<TH1D>(fBBB,plot_name);
	HbbStyle::applyStandardDataHistoStyle(*hBBB);
	hBBB->Scale(1./hBBB->Integral());
	auto hBBnB 		= GetFromTFile<TH1D>(fBBnB,plot_name);
	HbbStyle::applyStandardMCHistoStyle(*hBBnB);
	hBBnB->Scale(1./hBBnB->Integral());
	hBBnB->SetFillColor(kRed-10);

	// correct for the bin width
	HbbStyle::correctForTheBinWidth(*hBBB);
	HbbStyle::correctForTheBinWidth(*hBBnB);

	//histo with stat. uncertainty of bbnb
	auto hBBnB_stat = static_cast<TH1D*>(hBBnB->Clone("stat"));
	hBBnB_stat->SetLineColor(kBlack);
	hBBnB_stat->SetFillColor(kBlack);
	hBBnB_stat->SetFillStyle(3013); //3002
	hBBnB_stat->SetMarkerSize(0);
	TExec *er_0 = new TExec("er_0","gStyle->SetErrorX(0)");
	TExec *er_1 = new TExec("er_1","gStyle->SetErrorX(0.5)");

	//top pad and frame
	auto topPad = HbbStyle::getRatioTopPad(0.7);
	topPad->Draw();
	topPad->cd();

	// TLegend
	auto leg = HbbStyle::legend("right,top",4,0.6);
	leg->AddEntry(hBBB,"QCD triple b-tag","p");
	leg->AddEntry(hBBnB,"QCD reverse b-tag","f");
	leg->AddEntry(hBBnB_stat,"QCD reverse b-tag stat.","fl");

	TH2D *topFrame = new TH2D("topFrame","",1,200,1150,1,0,hBBB->GetMaximum()*1.2);
	HbbStyle::setRatioTopFrame(topFrame, can.YtoPixel(can.GetY1()));
	topFrame->GetYaxis()->SetTitle("a.u.");

	//Start drawing process
	topPad->cd();
	topFrame->Draw();
	hBBnB->Draw("hist same");
	hBBB->Draw("E1 same");
	er_1->Draw();
	hBBnB_stat->Draw("E2 same");
	er_0->Draw();
	leg->Draw();
	can.cd();

	//Bottom pad and frame
	auto botPad = HbbStyle::getRatioBottomPad(0.3);
	botPad->Draw();
	botPad->cd();
	TH2D *botFrame = new TH2D("botFrame","",1,topFrame->GetXaxis()->GetXmin(),topFrame->GetXaxis()	->GetXmax(),1,-3.9999,3.9999);
	HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
	botFrame->GetYaxis()->SetTitle("#frac{bbb-bbnb}{#sigma_{bbb}}");

	botPad->cd();
	auto ratio = static_cast<TH1D*>(hBBB->Clone("ratio"));
	ratio->Sumw2();
	//Minus and Devision
	devideTwoHistograms(ratio,hBBB,hBBnB);
	botFrame->Draw();
	ratio->Draw("same");

	TLine *horizLine = new TLine(topFrame->GetXaxis()->GetXmin(),0,topFrame->GetXaxis()->GetXmax(),0);
	horizLine -> SetLineStyle(2);
	horizLine -> Draw();

	can.cd();
	HbbStyle::drawStandardTitle("out");

	can.Print( (mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/PAS_bbb_vs_bbnb.pdf").c_str());
}

void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator){
	/*
	 * Function to perform:
	 * f = (a - b) / sigma(b)
	 * assume full correlations
	 */
	double difference, edifference, devis, edevis, bbb, ebbb, bbnb, ebbnb;
	for(int i = 1; i < nominator->GetNbinsX(); ++i){
		bbb  = nominator->GetBinContent(i);  ebbb  = nominator->GetBinError(i);
		bbnb = denominator->GetBinContent(i); ebbnb = denominator->GetBinError(i);
		difference  = bbb - bbnb;
		//Assume the full correlations
//		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb - 2 * ebbb * ebbnb );
		//Assume 0 correlation
		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb);
		devis 	= difference / ebbb;
		//assume the full correlations
		edevis 	= devis * sqrt ( pow(edifference/difference,2) + pow(ebbb/bbb,2) - 2 * edifference * ebbb / (difference * bbb) );
		//Assume the full correlations
		ratio->SetBinContent(i,devis);
		ratio->SetBinError(i,edevis);
	}
}
