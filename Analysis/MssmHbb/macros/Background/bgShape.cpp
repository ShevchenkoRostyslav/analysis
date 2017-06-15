/*
 * bgShape.cpp
 *
 *  Created on: 2 Jun 2017
 *      Author: shevchen
 *      Macro to plot Bg shape with some parameters
 */
#include "TH1.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TPad.h"
#include <string>
#include "RooEffProd.h"
#include "RooFormulaVar.h"

#include "Analysis/BackgroundModel/interface/RooSuperNovosibirsk.h"
#include "Analysis/MssmHbb/src/utilLib.cpp"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"
#include "Analysis/MssmHbb/macros/Drawer/RatioPlots.cpp"

RatioPlots style;

using namespace std;
using namespace analysis;
using namespace analysis::backgroundmodel;

void AbsComparisonOfPdfs(RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooRealVar &x, const bool& logY = true, const std::string& signal_subr = "sr1");
void DrawSignalPDFs(RooPlot &frame, const std::string& sr);

int main(){

	RooRealVar x = RooRealVar("mbb","x",0,200,650,"GeV");
	x.setBins(45);

	RooRealVar peak_supernovo1  = RooRealVar( "peak_supernovo1", "peak", 267.9681,0,1000);
	RooRealVar width_supernovo1 = RooRealVar( "width_supernovo1", "width", 63.5465, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo1  = RooRealVar( "tail_supernovo1", "tail", -0.4573, -10.0, 10.0);
	RooRealVar par0_supernovo1  = RooRealVar( "par0_supernovo1", "par0_supernovo1", -0.0006, -10.0, 10.0);
	RooArgList arg_list_supernovo1;
	arg_list_supernovo1.add(par0_supernovo1);
	RooSuperNovosibirsk pdf_supernovo1("supernovo1","supernovo1",x,peak_supernovo1,width_supernovo1,tail_supernovo1,arg_list_supernovo1);
	RooRealVar slope_novoeff_supernovo1 ( "slope_novoeff_supernovo1", "slope_novoeff", 0.0160, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo1( "turnon_novoeff_supernovo1", "turnon_novoeff",223.0267, 0, 1000);
	RooFormulaVar eff_supernovo1("eff_supernovo1","0.5*(TMath::Erf(slope_novoeff_supernovo1*(mbb - turnon_novoeff_supernovo1)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo1, turnon_novoeff_supernovo1));
	RooEffProd pdf_supernovo1effprod("supernovo1effprod","supernovo1effprod",pdf_supernovo1,eff_supernovo1);

	RooRealVar peak_supernovo2  = RooRealVar( "peak_supernovo2", "peak", 271.6864503918327,0,1000);
	RooRealVar width_supernovo2 = RooRealVar( "width_supernovo2", "width", 61.8238974706976, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo2  = RooRealVar( "tail_supernovo2", "tail", -0.3388502969062, -10.0, 10.0);
	RooRealVar par0_supernovo2  = RooRealVar( "par0_supernovo2", "par0_supernovo2", -0.0012712605914, -10.0, 10.0);
	RooRealVar par1_supernovo2  = RooRealVar( "par1_supernovo2", "par1_supernovo2", 0.0000009371010, -10.0, 10.0);
	RooArgList arg_list_supernovo2;
	arg_list_supernovo2.add(par0_supernovo2);
	arg_list_supernovo2.add(par1_supernovo2);
	RooSuperNovosibirsk pdf_supernovo2("supernovo2","supernovo2",x,peak_supernovo2,width_supernovo2,tail_supernovo2,arg_list_supernovo2);
	RooRealVar slope_novoeff_supernovo2( "slope_novoeff_supernovo2", "slope_novoeff", 0.0156, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo2( "turnon_novoeff_supernovo2", "turnon_novoeff",216.7410, 0, 1000);
	RooFormulaVar eff_supernovo2("eff_supernovo2","0.5*(TMath::Erf(slope_novoeff_supernovo2*(mbb - turnon_novoeff_supernovo2)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo2, turnon_novoeff_supernovo2));
	RooEffProd pdf_supernovo2effprod("supernovo2effprod","supernovo2effprod",pdf_supernovo2,eff_supernovo2);

	RooRealVar peak_supernovo3  = RooRealVar( "peak_supernovo3", "peak", 256.7718167270958,0,1000);
	RooRealVar width_supernovo3 = RooRealVar( "width_supernovo3", "width", 43.0482018127821, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo3  = RooRealVar( "tail_supernovo3", "tail", -0.1789714056358, -10.0, 10.0);
	RooRealVar par0_supernovo3  = RooRealVar( "par0_supernovo3", "par0_supernovo3", -0.0020902569342, -10.0, 10.0);
	RooRealVar par1_supernovo3  = RooRealVar( "par1_supernovo3", "par1_supernovo3", 0.0000028633721, -10.0, 10.0);
	RooRealVar par2_supernovo3  = RooRealVar( "par2_supernovo3", "par2_supernovo3", -0.0000000018054, -10.0, 10.0);
	RooArgList arg_list_supernovo3;
	arg_list_supernovo3.add(par0_supernovo3);
	arg_list_supernovo3.add(par1_supernovo3);
	arg_list_supernovo3.add(par2_supernovo3);
	RooSuperNovosibirsk pdf_supernovo3("supernovo3","supernovo3",x,peak_supernovo3,width_supernovo3,tail_supernovo3,arg_list_supernovo3);
	RooRealVar slope_novoeff_supernovo3( "slope_novoeff_supernovo3", "slope_novoeff",0.0040, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo3( "turnon_novoeff_supernovo3", "turnon_novoeff",556.6610, 0, 1000);
	RooFormulaVar eff_supernovo3("eff_supernovo3","0.5*(TMath::Erf(slope_novoeff_supernovo3*(mbb - turnon_novoeff_supernovo3)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo3, turnon_novoeff_supernovo3));
	RooEffProd pdf_supernovo3effprod("supernovo3effprod","supernovo3effprod",pdf_supernovo3,eff_supernovo3);

	//Shapes from B-only fit
	RooRealVar peak_supernovo3_bgonly  = RooRealVar( "peak_supernovo3_bgonly", "peak", 257.1481034426199,0,1000);
	RooRealVar width_supernovo3_bgonly = RooRealVar( "width_supernovo3_bgonly", "width", 42.6450024846353, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo3_bgonly  = RooRealVar( "tail_supernovo3_bgonly", "tail", -0.1347213233628, -10.0, 10.0);
	RooRealVar par0_supernovo3_bgonly  = RooRealVar( "par0_supernovo3_bgonly", "par0_supernovo3_bgonly", -0.0024909838351, -10.0, 10.0);
	RooRealVar par1_supernovo3_bgonly  = RooRealVar( "par1_supernovo3_bgonly", "par1_supernovo3_bgonly", 0.0000040981376, -10.0, 10.0);
	RooRealVar par2_supernovo3_bgonly  = RooRealVar( "par2_supernovo3_bgonly", "par2_supernovo3_bgonly", -0.0000000030484, -10.0, 10.0);
	RooArgList arg_list_supernovo3_bgonly;
	arg_list_supernovo3_bgonly.add(par0_supernovo3_bgonly);
	arg_list_supernovo3_bgonly.add(par1_supernovo3_bgonly);
	arg_list_supernovo3_bgonly.add(par2_supernovo3_bgonly);
	RooSuperNovosibirsk pdf_supernovo3_bgonly("supernovo3_bgonly","supernovo3_bgonly",x,peak_supernovo3_bgonly,width_supernovo3_bgonly,tail_supernovo3_bgonly,arg_list_supernovo3_bgonly);
	RooRealVar slope_novoeff_supernovo3_bgonly( "slope_novoeff_supernovo3_bgonly", "slope_novoeff",0.0040, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo3_bgonly( "turnon_novoeff_supernovo3_bgonly", "turnon_novoeff",556.6610, 0, 1000);
	RooFormulaVar eff_supernovo3_bgonly("eff_supernovo3_bgonly","0.5*(TMath::Erf(slope_novoeff_supernovo3_bgonly*(mbb - turnon_novoeff_supernovo3_bgonly)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo3_bgonly, turnon_novoeff_supernovo3_bgonly));
	RooEffProd pdf_supernovo3effprod_bgonly("supernovo3effprod_bgonly","supernovo3effprod_bgonly",pdf_supernovo3_bgonly,eff_supernovo3_bgonly);

	RooRealVar peak_supernovo2_bgonly  = RooRealVar( "peak_supernovo2_bgonly", "peak", 271.6696023686517,0,1000);
	RooRealVar width_supernovo2_bgonly = RooRealVar( "width_supernovo2_bgonly", "width", 62.0892263518770, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo2_bgonly  = RooRealVar( "tail_supernovo2_bgonly", "tail", -0.3483460372342, -10.0, 10.0);
	RooRealVar par0_supernovo2_bgonly  = RooRealVar( "par0_supernovo2_bgonly", "par0_supernovo2_bgonly", -0.0012108746186, -10.0, 10.0);
	RooRealVar par1_supernovo2_bgonly  = RooRealVar( "par1_supernovo2_bgonly", "par1_supernovo2_bgonly", 0.0000008448336, -10.0, 10.0);
	RooArgList arg_list_supernovo2_bgonly;
	arg_list_supernovo2_bgonly.add(par0_supernovo2_bgonly);
	arg_list_supernovo2_bgonly.add(par1_supernovo2_bgonly);
	RooSuperNovosibirsk pdf_supernovo2_bgonly("supernovo2_bgonly","supernovo2_bgonly",x,peak_supernovo2_bgonly,width_supernovo2_bgonly,tail_supernovo2_bgonly,arg_list_supernovo2_bgonly);
	RooRealVar slope_novoeff_supernovo2_bgonly( "slope_novoeff_supernovo2_bgonly", "slope_novoeff", 0.0156, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo2_bgonly( "turnon_novoeff_supernovo2_bgonly", "turnon_novoeff",216.7410, 0, 1000);
	RooFormulaVar eff_supernovo2_bgonly("eff_supernovo2_bgonly","0.5*(TMath::Erf(slope_novoeff_supernovo2_bgonly*(mbb - turnon_novoeff_supernovo2_bgonly)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo2_bgonly, turnon_novoeff_supernovo2_bgonly));
	RooEffProd pdf_supernovo2effprod_bgonly("supernovo2effprod_bgonly","supernovo2effprod_bgonly",pdf_supernovo2_bgonly,eff_supernovo2_bgonly);

	RooRealVar peak_supernovo1_bgonly  = RooRealVar( "peak_supernovo1_bgonly", "peak", 267.9443936335616,0,1000);
	RooRealVar width_supernovo1_bgonly = RooRealVar( "width_supernovo1_bgonly", "width", 63.4007056709499, 5.0, 1000, "GeV");
	RooRealVar tail_supernovo1_bgonly  = RooRealVar( "tail_supernovo1_bgonly", "tail", -0.4566147625246, -10.0, 10.0);
	RooRealVar par0_supernovo1_bgonly  = RooRealVar( "par0_supernovo1_bgonly", "par0_supernovo1_bgonly", -0.0006159916146, -10.0, 10.0);
	RooArgList arg_list_supernovo1_bgonly;
	arg_list_supernovo1_bgonly.add(par0_supernovo1_bgonly);
	RooSuperNovosibirsk pdf_supernovo1_bgonly("supernovo1_bgonly","supernovo1_bgonly",x,peak_supernovo1_bgonly,width_supernovo1_bgonly,tail_supernovo1_bgonly,arg_list_supernovo1_bgonly);
	RooRealVar slope_novoeff_supernovo1_bgonly ( "slope_novoeff_supernovo1_bgonly", "slope_novoeff", 0.0160, 0.0, 0.1);
	RooRealVar turnon_novoeff_supernovo1_bgonly( "turnon_novoeff_supernovo1_bgonly", "turnon_novoeff",223.0267, 0, 1000);
	RooFormulaVar eff_supernovo1_bgonly("eff_supernovo1_bgonly","0.5*(TMath::Erf(slope_novoeff_supernovo1_bgonly*(mbb - turnon_novoeff_supernovo1_bgonly)) + 1)",
		                    RooArgSet(x,slope_novoeff_supernovo1_bgonly, turnon_novoeff_supernovo1_bgonly));
	RooEffProd pdf_supernovo1effprod_bgonly("supernovo1effprod_bgonly","supernovo1effprod_bgonly",pdf_supernovo1_bgonly,eff_supernovo1_bgonly);

	//Generator functions from original workspaces
	//SuperDijet2 workspace
	auto &wSuperDiJet2_sr1 = *GetRooWorkspace(mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/Prescale_v5/superdijeteffprod2_200_to_650_10GeV_G4/workspace/FitContainer_workspace.root");
	auto &pdfSuperDiJet2_sr1 = *GetFromRooWorkspace<RooAbsPdf>(wSuperDiJet2_sr1, "background");
	pdfSuperDiJet2_sr1.SetTitle("SuperDijet2EffProd");
	//Bernstein-8 x Efficiency
	auto &wBernstein8 = *GetRooWorkspace(mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/Prescale_v5/berneff8_200_to_650_10GeV/workspace/FitContainer_workspace.root");
	auto &pdfBernstein8 = *GetFromRooWorkspace<RooAbsPdf>(wBernstein8, "background");
	pdfBernstein8.SetTitle("Bernstein8EffProd");

	AbsComparisonOfPdfs(pdf_supernovo1effprod,pdf_supernovo2effprod,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod,pdf_supernovo3effprod,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod,pdf_supernovo3effprod,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod,pdfBernstein8,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod,pdfBernstein8,x);
	AbsComparisonOfPdfs(pdf_supernovo3effprod,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo3effprod,pdfBernstein8,x);

	AbsComparisonOfPdfs(pdf_supernovo1effprod_bgonly,pdf_supernovo2effprod_bgonly,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod_bgonly,pdf_supernovo3effprod_bgonly,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod_bgonly,pdf_supernovo3effprod_bgonly,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod_bgonly,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo1effprod_bgonly,pdfBernstein8,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod_bgonly,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo2effprod_bgonly,pdfBernstein8,x);
	AbsComparisonOfPdfs(pdf_supernovo3effprod_bgonly,pdfSuperDiJet2_sr1,x);
	AbsComparisonOfPdfs(pdf_supernovo3effprod_bgonly,pdfBernstein8,x);

	return 0;
}

void AbsComparisonOfPdfs(RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooRealVar &x, const bool& logY, const std::string& signal_subr){
	/*
	 * Function to compare pdfs in absolute scale
	 */
	x.setBins(1000);
	RooPlot &frame = *x.frame();
	TH1D &h1 = *static_cast<TH1D*>(pdf1.createHistogram("h1",x));
	TH1D &h2 = *static_cast<TH1D*>(pdf2.createHistogram("h2",x));

	pdf1.plotOn(&frame,RooFit::LineColor(2),RooFit::LineWidth(2));
	pdf2.plotOn(&frame,RooFit::LineColor(kBlue),RooFit::LineWidth(2));

	TCanvas can("can","can",800,600);
	TPad pad1_("pad1","pad1",0,0.1,1,1);
	pad1_ . SetBottomMargin(0.2);
	pad1_ . SetRightMargin(0.05);
	pad1_ . SetLeftMargin(0.16);
	pad1_ . Draw();
	pad1_ . cd();

	frame.GetXaxis()->SetNdivisions(505);
	frame.GetYaxis()->SetNdivisions(206);
	frame.GetYaxis()->SetTitleOffset(1.1);
	frame.GetXaxis()->SetTitleOffset(999);
	frame.GetXaxis()->SetLabelOffset(999);
	frame.GetYaxis()->SetTitleSize(0.048);
	frame.SetStats(0);
	string title = static_cast<string>(pdf1.GetTitle()) + "_vs_" + static_cast<string>(pdf2.GetTitle());
	frame.SetTitle( (title + "; M_{A} [GeV]; a.u.").c_str() );
	gPad->SetLogy(logY);
	DrawSignalPDFs(frame,signal_subr);
	frame.Draw();

	can.cd();
	TPad pad2_("pad2","pad2",0,0.0,1,0.265);
	pad2_ . SetTopMargin(0);
	pad2_ . SetLeftMargin(pad1_.GetLeftMargin());
	pad2_ . SetRightMargin(pad1_.GetRightMargin());
	pad2_ . SetBottomMargin(0.28);
	pad2_ . Draw();
	pad2_ . cd();

	h1.Add(&h2,-1);
//	h1.Divide(&h2);
	h1.Scale(1000);
//	h1.SetTitle((";M_{A} [GeV];#frac{(" + static_cast<string>(pdf1.GetTitle()) + " - " + static_cast<string>(pdf2.GetTitle()) + ")}{" + static_cast<string>(pdf2.getTitle()) + "}, %").c_str() );
	h1.SetTitle((";M_{A} [GeV];#splitline{(" + static_cast<string>(pdf1.GetTitle()) + " -}{- " + static_cast<string>(pdf2.GetTitle()) + ")  x 10^{3} }").c_str() );
	style.SetBottomStyle(&h1);
	h1.SetStats(0);
	h1.SetLineWidth(1.5);
	h1.Draw();
	TLine l(x.getMin(),0,x.getMax(),0);
	l.SetLineStyle(9);
	l.Draw();

	h1.GetYaxis()->SetTitleSize(0.07);
	h1.GetYaxis()->SetTitleOffset(0.9);
	can.Update();
	can.Print(("../macros/pictures/Background/" + title + ".pdf").c_str());
}

void DrawSignalPDFs(RooPlot &frame, const std::string& sr){
	/*
	 * Draw signal pdfs on the frame according to the sub-range
	 */
	if(sr == "sr1"){
		auto &p300 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-300/workspace/FitContainer_workspace.root","signal");
		auto &p350 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-350/workspace/FitContainer_workspace.root","signal");
		auto &p400 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-400/workspace/FitContainer_workspace.root","signal");
		auto &p500 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-500/workspace/FitContainer_workspace.root","signal");

		p300.plotOn(&frame,RooFit::LineColor(kRed-4),RooFit::LineStyle(1),RooFit::Normalization(0.01));
		p350.plotOn(&frame,RooFit::LineColor(kGreen-4),RooFit::LineStyle(2),RooFit::Normalization(0.01));
		p400.plotOn(&frame,RooFit::LineColor(kCyan-4),RooFit::LineStyle(6),RooFit::Normalization(0.01));
		p500.plotOn(&frame,RooFit::LineColor(kBlue-4),RooFit::LineStyle(10),RooFit::Normalization(0.01));
	}
	else if(sr  == "sr2"){
		auto &p600 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-600/workspace/FitContainer_workspace.root","signal");
		auto &p700 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-700/workspace/FitContainer_workspace.root","signal");
		auto &p900 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-900/workspace/FitContainer_workspace.root","signal");

		p600.plotOn(&frame,RooFit::LineColor(kRed-4),RooFit::LineStyle(1),RooFit::Normalization(0.01));
		p700.plotOn(&frame,RooFit::LineColor(kGreen-4),RooFit::LineStyle(2),RooFit::Normalization(0.01));
		p900.plotOn(&frame,RooFit::LineColor(kCyan-4),RooFit::LineStyle(6),RooFit::Normalization(0.01));
	}
	else if(sr == "sr3"){
		auto &p1100 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-1100/workspace/FitContainer_workspace.root","signal");
		auto &p1300 = *GetRooObjectFromTFile<RooAbsPdf>("../output/ReReco_signal_M-1300/workspace/FitContainer_workspace.root","signal");

		p1100.plotOn(&frame,RooFit::LineColor(kRed-4),RooFit::LineStyle(1),RooFit::Normalization(0.01));
		p1300.plotOn(&frame,RooFit::LineColor(kGreen-4),RooFit::LineStyle(2),RooFit::Normalization(0.01));
	}

}
