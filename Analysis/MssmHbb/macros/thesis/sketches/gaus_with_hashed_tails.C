#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;

void gaus_with_hashed_tails(){

	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.16);
	gStyle->SetPadRightMargin(0.02);

	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetCanvasDefH(600); //Height of canvas
	gStyle->SetCanvasDefW(600); //Width of canvas
	gStyle->SetCanvasDefX(0);   //POsition on screen
	gStyle->SetCanvasDefY(0);

	// For the Pad:
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(kWhite);
	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);

	// For the frame:
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameFillStyle(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameLineStyle(1);
	gStyle->SetFrameLineWidth(1);

	gStyle->SetOptTitle(0);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.06, "XYZ");
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(1.25);

	TCanvas *can = new TCanvas();

	//Full gaus function
	TF1 *full_gaus = new TF1("full_gaus","TMath::Gaus(x, 5, 0.8, 1)",2,8);
	full_gaus->SetTitle("");
	full_gaus->SetLineColor(1);
	full_gaus->GetXaxis()->SetTitle("#tilde{#theta}");
	full_gaus->GetYaxis()->SetTitle("f(#tilde{#theta},#theta)");
	full_gaus->SetMaximum(0.7);
	full_gaus->Draw();

	//Right area to be hashed
	TF1 *right_gaus_part = new TF1("right_gaus_part","TMath::Gaus(x, 5, 0.8, 1)",6.2,8);
	right_gaus_part->SetLineColor(1);
	right_gaus_part->SetLineStyle(2);
	right_gaus_part->SetFillStyle(3004);
	right_gaus_part->SetFillColor(kRed-4);
	right_gaus_part->Draw("same");

	TLine *l_right = new TLine(right_gaus_part->GetXmin(),0,right_gaus_part->GetXmin(),0.4);
	l_right->SetLineStyle(2);
	l_right->Draw("same");

	TLatex *txt_falpha = new TLatex();
//	txt_falpha->SetTextColor(kBlue+2);
	txt_falpha->SetTextSize(txt_falpha->GetTextSize()*1.2);
	txt_falpha->DrawLatex(right_gaus_part->GetXmin()*1.02, 0.4,
	                     ("f_{#alpha}(#theta)"));

	//Left area to be hashed
	TF1 *left_gaus_part = new TF1("left_gaus_part","TMath::Gaus(x, 5, 0.8, 1)",2,4.2);
	left_gaus_part->SetLineColor(1);
	left_gaus_part->SetLineStyle(2);
	left_gaus_part->SetFillStyle(3005);
	left_gaus_part->SetFillColor(kBlue-4);
	left_gaus_part->Draw("same");

	TLine *l_left = new TLine(left_gaus_part->GetXmax(),0,left_gaus_part->GetXmax(),0.4);
	l_left->SetLineStyle(2);
	l_left->Draw("same");

	TLatex *txt_fbeta = new TLatex();
//	txt_fbeta->SetTextColor(kBlue+2);
	txt_fbeta->SetTextSize(txt_fbeta->GetTextSize()*1.2);
	txt_fbeta->DrawLatex(left_gaus_part->GetXmax()*0.8, 0.4,
	                     ("f_{#beta}(#theta)"));

	TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
	leg->SetBorderSize(0);
	leg->AddEntry(left_gaus_part,"#beta = P(#tilde{#theta} #leq f_{#beta}(#theta))","f");
	leg->AddEntry(right_gaus_part,"#alpha = P(#tilde{#theta} #geq f_{#alpha}(#theta))","f");
	leg->Draw("same");

	can->Print("../../pictures/Thesis/gaus_with_porbs.pdf");


}
