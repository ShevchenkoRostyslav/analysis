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
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"

#include "TH1.h"

#include <string>
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include <RooFit.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include "RooAddPdf.h"
#include "RooWorkspace.h"

#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace analysis;
using namespace RooFit;
bool TEST_ = true;
double xsec_vis = 1;

struct SubRange{
	SubRange(){};
	SubRange(const SubRange& s){ //copy constructor
		workspace = s.workspace;
		h2sigmaU = s.h2sigmaU;
		h2sigmaD = s.h2sigmaD;
		h1sigmaU = s.h1sigmaU;
		h1sigmaD = s.h1sigmaD;
		hbkg = s.hbkg;
		sub_range = s.sub_range;
		xsec_vis = s.xsec_vis;
		fMLL_bg_nobias = s.fMLL_bg_nobias;
		fMLL = s.fMLL;
	}
	SubRange& operator=(const SubRange& s){
			workspace = s.workspace;
			h2sigmaU = s.h2sigmaU;
			h2sigmaD = s.h2sigmaD;
			h1sigmaU = s.h1sigmaU;
			h1sigmaD = s.h1sigmaD;
			hbkg = s.hbkg;
			sub_range = s.sub_range;
			xsec_vis = s.xsec_vis;
			fMLL_bg_nobias = s.fMLL_bg_nobias;
			fMLL = s.fMLL;
			return *this;
	}

	RooWorkspace *workspace;
	TH1F h2sigmaU, h2sigmaD, h1sigmaU, h1sigmaD, hbkg;
	string fMLL_bg_nobias, fMLL;
	std::string sub_range;
	double xsec_vis;
};

void SetupTopFrame(RooPlot *frame, const double& pad_w, const double& pad_h);
void SetupBottomFrame(RooPlot *frame2, const double& pad_w, const double& pad_h);
void SetupMiddleFrame(RooPlot *frame2, const double& pad_w, const double& pad_h);
void SetupBottomPad(TPad * pad);

SubRange setup_subranges(const double& bin_width, const int& mass_point, const std::string& combine_input, const std::string& output_path);
void CollectSubRangesData(std::map<std::string,SubRange>& data,const double& bin_width,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path);

void PlotFullMassFit(const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path);


int main(){

	int nbins = 100;
	string output_folder = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170823/";
	string input_folder_combine = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/";
	if(!mssmhbb::blinded) { input_folder_combine += "unblinded/"; output_folder += "unblinded/";} 
	output_folder += "M12PostFits/test";
	CheckOutputDir(output_folder);
	input_folder_combine += "independent/mll/";

	map<string,int> signals = {
			{"sr1" , 400},
			{"sr2" , 600},
			{"sr3" , 1300}
	};
	//Modify for log/lin scale
	xsec_vis = 1;

	PlotFullMassFit(nbins,signals,input_folder_combine,output_folder);

	return 0;
}

void PlotFullMassFit(const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path){
	/*
	 * Method to plot the full mass fit
	 */

	HbbStyle style;
	style.setTDRstyle(PRELIMINARY);
	gStyle->SetErrorX(0);
	gROOT->ForceStyle();

	std::map<std::string,SubRange> data;
	double bin_width = (1700.0-200.0)/nbins;//GeV
	CollectSubRangesData(data,bin_width,signals,combine_input,output_path);

	std::map<string,RooRealVar*> x_i, x_sgn, n_sgn_i,n_bkg_i,n_bkg_nobias_i;
	std::map<string,RooDataHist*> data_obs_i, data_rebinned_i;
	std::map<string,RooAbsPdf*> pdf_bkg_i,pdf_bkg_nobias_i,pdf_comb_ext_i,pdf_sgn_i;
	std::map<string,RooPlot*> frame1_i,frame2_i,frame_sgn_draw;
	std::map<string,TH1F*> h2sigmaU_i,h2sigmaD_i,h1sigmaU_i,h1sigmaD_i,hbkg_i,h_pdf_comb_ext_i;
	std::map<string,TH1F> h_pdf_bkg_i, h_pdf_spbkg_i;
	std::map<string,RooFitResult*> roofitresult_nobias_i, roofitresult_i;
	std::map<string,RooAbsReal*> pdf_sgn_norm_i;
	std::map<string,TFile*> fMLL_bg_nobias_i, fMLL_i;
	std::map<string,TLine*> white_line_dn, white_line_up;

	RooRealVar x("mbb","mbb",200,1700);
	//TOP FRAME
	TCanvas c1("c1", "canvas", 600, 600);
	//Size of the main pad in pixels:
	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());
	//TPad
	c1.cd();
	double pad1_ymin = 0.3;
	TPad pad1("pad1", "pad1", 0.0, pad1_ymin, 1, 1);
	pad1.SetTopMargin(gStyle->GetPadTopMargin() * 1./(1 - pad1_ymin));
	pad1.SetRightMargin(gStyle->GetPadRightMargin());
    pad1.SetBottomMargin(0.02);
    pad1.Draw();

    c1.cd();
    TPad pad2("pad2", "pad2", 0.0, 0.02, 1, 0.28);
    pad2.SetTopMargin(0.02);
    pad2.SetBottomMargin(0);
    pad2.SetGridx();
    pad2.SetGridy();
    pad2.Draw();
    c1.cd();

	//Frame
    auto *frame1 = x.frame(Name("frame"),Range(200,1700));
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle(( ("Events / ( " + to_string(int(bin_width)) + " GeV )").c_str() ));
    frame1->SetMinimum(0.5);
    frame1->SetMaximum(7e+04);
    pad1.cd();
    SetupTopFrame(frame1,pad_wPixel,pad_hPixel);
    frame1->Draw();
	double pad1_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad1_hPixel = gPad->YtoPixel(gPad->GetY1());

    auto *frame2 = x.frame(Name("frame2"),Range(200,1700));
    frame2->GetXaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetTickLength(0);
    frame2->GetYaxis()->SetRangeUser(-3.9999,3.9999);
    pad2.cd();
    SetupBottomFrame(frame2,pad1_wPixel,pad1_hPixel);
    frame2->Draw();

    auto *leg = new TLegend(0.65,0.35,0.95,0.88,"","brNDC");
//    auto *leg = new TLegend(0.6,0.4,0.95,0.88,"","brNDC");
//    leg->SetNColumns(2);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
//     leg->SetTextSize(0.02);
    leg->SetFillColor(10);

    std::array<int,3> colours {{820,861,901}};
    
    //Full data from the backgorund file
     auto &data_full = *GetRooObjectFromTFile<RooDataHist>(mssmhbb::path_bg_sr1.c_str(),("data_container"));
     data_full.plotOn(frame1,Binning(nbins),Name("data"));
     double xsec_modified = 1;

	for(int i=1;i<=3;++i){
		if(TEST_) std::cout<<"In SR loop"<<std::endl;
		pad1.cd();
		string sr = "sr" + std::to_string(i);
		x_i[sr] 				= GetFromRooWorkspace<RooRealVar>(*data[sr].workspace,("mbb"));
		x_sgn[sr]				= GetRooObjectFromTFile<RooRealVar>(mssmhbb::signal_fits.at(signals.at(sr)),("mbb"));

		n_sgn_i[sr] 			= GetFromRooWorkspace<RooRealVar>(*data[sr].workspace,"n_sgn");
		n_bkg_i[sr] 			= GetFromRooWorkspace<RooRealVar>(*data[sr].workspace,"n_bkg");
		n_bkg_nobias_i[sr] 		= GetFromRooWorkspace<RooRealVar>(*data[sr].workspace,"n_bkg_nobias");

		data_obs_i[sr] 			= GetFromRooWorkspace<RooDataHist>(*data[sr].workspace,"data_obs");
		data_rebinned_i[sr] 	= GetFromRooWorkspace<RooDataHist>(*data[sr].workspace,"data_rebinned");

		pdf_bkg_i[sr] 			= GetFromRooWorkspace<RooAbsPdf>(*data[sr].workspace,"pdf_bkg");
		pdf_bkg_nobias_i[sr] 	= GetFromRooWorkspace<RooAbsPdf>(*data[sr].workspace,"pdf_bkg_nobias");
		pdf_comb_ext_i[sr] 		= GetFromRooWorkspace<RooAbsPdf>(*data[sr].workspace,"pdf_comb_ext");
		
		pdf_sgn_i[sr]			= GetRooObjectFromTFile<RooAbsPdf>(mssmhbb::signal_fits.at(signals.at(sr)),"signal");
//		pdf_sgn_i[sr] 			= GetFromRooWorkspace<RooAbsPdf>(mssmhbb::signal_fits.at(signals[sr]),"signal");

		pdf_sgn_norm_i[sr]		= GetFromRooWorkspace<RooAbsReal>(*data[sr].workspace,"pdf_sgn_norm");

		h2sigmaU_i[sr] 			= &data[sr].h2sigmaU;
		h2sigmaD_i[sr] 			= &data[sr].h2sigmaD;
		h1sigmaU_i[sr] 			= &data[sr].h1sigmaU;
		h1sigmaD_i[sr] 			= &data[sr].h1sigmaD;
		hbkg_i[sr] 				= &data[sr].hbkg;

		if(TEST_) std::cout<<"In SR loop: post variables"<<std::endl;

		fMLL_bg_nobias_i[sr]	= new TFile(data[sr].fMLL_bg_nobias.c_str());
		fMLL_i[sr] = new TFile(data[sr].fMLL.c_str(),"READ");	//S+Bg mll
		roofitresult_nobias_i[sr]=GetFromTFile<RooFitResult>(*fMLL_bg_nobias_i[sr],"fit_b");
		roofitresult_i[sr]			= GetFromTFile<RooFitResult>(*fMLL_i[sr],"fit_s");

		frame1_i[sr] = x_i[sr]->frame(Name("frame1"),Range(x_i[sr]->getMin(),x_i[sr]->getMax()));
		double min = x_i[sr]->getMin(), max =  x_i[sr]->getMax();
		if(i==2)min -= 200;
		x_sgn[sr]->setRange(min,max);
		frame_sgn_draw[sr] = x_sgn[sr]->frame(Name("frame_sgn"),Range(min,max));

		//Start plotting
		int nbins_per_sr = (x_i[sr]->getMax() - x_i[sr]->getMin())/bin_width;
		data_rebinned_i[sr]->plotOn(frame1_i[sr],Binning(nbins_per_sr), Name( ("data_obs_" + sr).c_str()), MarkerColor(kWhite),MarkerSize(0.0008),LineColor(kWhite),XErrorSize(0),LineWidth(0),DrawOption("PX"));

		pdf_bkg_i[sr]->plotOn(frame1_i[sr],
	    		VisualizeError(*roofitresult_i[sr], 2, true),
				LineColor(kWhite),
				LineStyle(kSolid),
				FillColor(kOrange), Name( (string(pdf_bkg_i[sr]->GetName()) +"_2sigma_" + sr).c_str()), Normalization(n_bkg_i[sr]->getVal() , RooAbsReal::NumEvent), MoveToBack() );

	    pdf_bkg_i[sr]->plotOn(frame1_i[sr],
	    		VisualizeError(*roofitresult_i[sr], 1, true),
				LineColor(kWhite),
				LineStyle(kSolid),
				FillColor(kGreen+1), Name( (string(pdf_bkg_i[sr]->GetName())+"_1sigma_" + sr).c_str()) , Normalization(n_bkg_i[sr]->getVal() , RooAbsReal::NumEvent) );

	    pdf_comb_ext_i[sr]->plotOn(frame1_i[sr],
	    			LineWidth(Width_t(1)),
				LineColor(kRed),
				LineStyle(3),
				Name(pdf_comb_ext_i[sr]->GetName()) );
	    if(TEST_) std::cout<<"In SR loop: post_1/2 sigma bands"<<std::endl;

//	    pdf_bkg_nobias_i[sr]->plotOn(frame1_i[sr],
//	    		LineWidth(3),
//				LineColor(kBlue),
//				LineStyle(kSolid),
//				Name((string(pdf_bkg_nobias_i[sr]->GetName()) + "_" + sr).c_str()),
//				Normalization(n_bkg_nobias_i[sr]->getVal(), RooAbsReal::NumEvent) );
	    //generate histo to represent S+Bg
//	    h_pdf_spbkg_i[sr] = *static_cast<TH1F*>(pdf_comb_ext_i->createHistogram(
//	    		("h_pdf_spbkg_i_" + sr).c_str(),*x_i[sr],Binning(10000,x_i[sr]->getMin(),x_i[sr]->getMax())
//	    		)	);

	    //Generate HISTO from this pdf to get read of annoying lines in the middle of the pad
	    h_pdf_bkg_i[sr] = *static_cast<TH1F*>(pdf_bkg_i[sr]->createHistogram(
	    		("h_pdf_bkg_i_" + sr).c_str(),*x_i[sr],Binning(10000,x_i[sr]->getMin(),x_i[sr]->getMax())
				) );
	    h_pdf_bkg_i[sr].Scale(n_bkg_i[sr]->getVal() * 10000/nbins_per_sr);
	    h_pdf_bkg_i[sr].SetMarkerColor(kBlue);
	    h_pdf_bkg_i[sr].SetLineColor(kBlue);
	    h_pdf_bkg_i[sr].SetMarkerSize(0.075);

	    if(TEST_) std::cout<<"In SR loop: before signal "<<std::endl;
	    
	    if(sr == "sr1") xsec_modified = xsec_vis * 500;
 	   	else if (sr == "sr2") xsec_modified = xsec_vis * 120;
		else  xsec_modified = xsec_vis * 90;
		data[sr].xsec_vis = xsec_modified;
		cout<<"xsec = "<<xsec_modified<<endl;

	    pdf_sgn_i[sr]->plotOn(frame_sgn_draw[sr],
	    		LineWidth(3),
				LineColor(colours.at(i-1)),
				LineStyle(kDotted),
				Name((string(pdf_sgn_i[sr]->GetName()) + "_" + sr ).c_str()),
//				Normalization(pdf_sgn_norm_i[sr]->getVal()*data[sr].xsec_vis , RooAbsReal::NumEvent) );
				Normalization(xsec_modified) );

	    if(TEST_) std::cout<<"In SR loop: post drawing"<<std::endl;

	    white_line_dn[sr] = new TLine(	h_pdf_bkg_i[sr].GetXaxis()->GetXmin(), 0, h_pdf_bkg_i[sr].GetXaxis()->GetXmin(), h_pdf_bkg_i[sr].GetBinContent(1));
		white_line_dn[sr]->SetLineColor(kWhite);
		white_line_dn[sr]->SetLineWidth(1.4);

	    white_line_up[sr] = new TLine(	h_pdf_bkg_i[sr].GetXaxis()->GetXmax(), 0, h_pdf_bkg_i[sr].GetXaxis()->GetXmax(), h_pdf_bkg_i[sr].GetBinContent(h_pdf_bkg_i[sr].GetNbinsX()));
		white_line_up[sr]->SetLineColor(kWhite);
		white_line_up[sr]->SetLineWidth(1.4);

		frame1->Draw();

	    frame1_i[sr]->addTH1(&h_pdf_bkg_i[sr],"P HIST SAME");
	    frame1_i[sr]->Draw("same");
//	    frame1_i[sr]->SetAxisRange(3.8, 3.8, "X");
	    frame_sgn_draw[sr]->Draw("same");

	    if(TEST_) std::cout<<"In SR loop: done"<<std::endl;

	    pad2.cd();

	    frame2_i[sr] = x_i[sr]->frame(Name("frame2"),Range(x_i[sr]->getMin(),x_i[sr]->getMax()));
	    //Special rools to draw pools:
	    if(i==1){
	    	data[sr].h2sigmaU.GetXaxis()->SetRangeUser(200,650);
	    	data[sr].h2sigmaD.GetXaxis()->SetRangeUser(200,650);
	    	data[sr].h1sigmaU.GetXaxis()->SetRangeUser(200,650);
	    	data[sr].h1sigmaD.GetXaxis()->SetRangeUser(200,650);
	    	data[sr].hbkg.GetXaxis()->SetRangeUser(200,650);
	    	data[sr].hbkg.SetFillColor(kSpring);
	    	data[sr].hbkg.SetLineColor(kSpring);
	    }
	    if(i==2){
	    	data[sr].h2sigmaU.GetXaxis()->SetRangeUser(650,1190);
	    	data[sr].h2sigmaD.GetXaxis()->SetRangeUser(650,1190);
	    	data[sr].h1sigmaU.GetXaxis()->SetRangeUser(650,1190);
	    	data[sr].h1sigmaD.GetXaxis()->SetRangeUser(650,1190);
	    	data[sr].hbkg.GetXaxis()->SetRangeUser(650,1190);
	    }
	    if(i==3){
	    	data[sr].h2sigmaU.GetXaxis()->SetRangeUser(1190,1700);
	    	data[sr].h2sigmaD.GetXaxis()->SetRangeUser(1190,1700);
	    	data[sr].h1sigmaU.GetXaxis()->SetRangeUser(1190,1700);
	    	data[sr].h1sigmaD.GetXaxis()->SetRangeUser(1190,1700);
	    	data[sr].hbkg.GetXaxis()->SetRangeUser(1190,1700);
	    	data[sr].hbkg.SetFillColor(kPink+1);
	    	data[sr].hbkg.SetLineColor(kPink+1);
	    }
	    data[sr].hbkg.SetLineColor(kBlack);
	    data[sr].hbkg.SetMarkerSize(0.75);
	    frame2_i[sr]->addTH1(&data[sr].h2sigmaU,"HIST");
	    frame2_i[sr]->addTH1(&data[sr].h2sigmaD,"HIST");
	    frame2_i[sr]->addTH1(&data[sr].h1sigmaU,"HIST");
	    frame2_i[sr]->addTH1(&data[sr].h1sigmaD,"HIST");
	    frame2_i[sr]->addTH1(&data[sr].hbkg,"E1");
	    frame2_i[sr]->Draw("same");
	    
	    pad1.cd();

		c1.cd();
//		c1.Draw();
		data[sr].h2sigmaU.GetXaxis()->SetRangeUser(x_i[sr]->getMin(),x_i[sr]->getMax());		
                data[sr].h2sigmaD.GetXaxis()->SetRangeUser(x_i[sr]->getMin(),x_i[sr]->getMax());
                data[sr].h1sigmaD.GetXaxis()->SetRangeUser(x_i[sr]->getMin(),x_i[sr]->getMax());
                data[sr].h1sigmaU.GetXaxis()->SetRangeUser(x_i[sr]->getMin(),x_i[sr]->getMax());
		data[sr].hbkg.GetXaxis()->SetRangeUser(x_i[sr]->getMin(),x_i[sr]->getMax());
	}
	//Fill the legend
	leg->AddEntry(frame1->findObject("data"),"Data","PE");
	leg->AddEntry(frame1_i["sr1"]->findObject((h_pdf_bkg_i["sr1"].GetName())),"Background","L");
    leg->AddEntry(frame2_i["sr1"]->findObject("h1sigmaU"), "#pm1 std. deviation", "F");
    leg->AddEntry(frame2_i["sr1"]->findObject("h2sigmaU"), "#pm2 std. deviation", "F");
//     //leg->AddEntry(h_pdf_comb_ext_i["sr1"], ("Bkg. + signal"), "L");
    leg->AddEntry(frame_sgn_draw["sr1"]->getCurve((string(pdf_sgn_i["sr1"]->GetName())+"_sr1").c_str()), ("m_{A/H} =   " + to_string(signals.at("sr1")) + " GeV, #sigma = " + to_string(int(data["sr1"].xsec_vis)) + " pb" ).c_str(), "L");
    leg->AddEntry(frame_sgn_draw["sr2"]->getCurve((string(pdf_sgn_i["sr2"]->GetName())+"_sr2").c_str()), ("m_{A/H} =   " + to_string(signals.at("sr2")) + " GeV, #sigma = " + to_string(int(data["sr2"].xsec_vis)) + " pb" ).c_str(), "L");
    leg->AddEntry(frame_sgn_draw["sr3"]->getCurve((string(pdf_sgn_i["sr3"]->GetName())+"_sr3").c_str()), ("m_{A/H} = " + to_string(signals.at("sr3")) + " GeV, #sigma = " + to_string(int(data["sr3"].xsec_vis)) + " pb" ).c_str(), "L");
    leg->Draw();
    style.drawStandardTitle("out");
//    "m_{A/H} = " + mass_point + " GeV"

	frame1->SetMaximum(3.3e+04);
	c1.Print("something.pdf");
	frame1->SetMaximum(7e+04);
	frame1->SetMinimum(0.5);
	pad1.SetLogy();
	c1.Print("something_log.pdf");

	//Draw 3 bottom frames

	TCanvas can_4pads("can_4pads","",800,600);
	can_4pads.cd();
//	TPad pad1_0("pad1_0", "pad1", 0.02, 0.50, 1, 0.96);
	pad1_ymin = 0.5;
	TPad pad1_0("pad1_0", "pad1", 0.0, 0.50, 1, 1.);
	pad1_0.SetTopMargin(gStyle->GetPadTopMargin() * 1./(1 - pad1_ymin));
	pad1_0.SetBottomMargin(0.0);
    pad1_0.Draw();

	pad1_wPixel = gPad->XtoPixel(gPad->GetX2());
	pad1_hPixel = gPad->YtoPixel(gPad->GetY1());

    TPad pad2_0("pad2_0","",0.0,0.37,1,0.5);
	SetupBottomPad(&pad2_0); 
	pad2_0.SetTopMargin(0.01);
	pad2_0.SetBottomMargin(0.0);
	can_4pads.cd();
	TPad pad2_1("pad2_1","",0.0,0.24,1,0.37);
	SetupBottomPad(&pad2_1); 	
	pad2_1.SetTopMargin(pad2_0.GetTopMargin());
	pad2_1.SetBottomMargin(0.0);
	can_4pads.cd();
	TPad pad2_2("pad2_2","",0.0,0.0,1,0.24);
	SetupBottomPad(&pad2_2); 
	pad2_2.SetTopMargin(pad2_0.GetTopMargin() * pad2_0.GetAbsHNDC() / pad2_2.GetAbsHNDC());
	pad2_2.SetBottomMargin(gStyle->GetPadBottomMargin() * 0.8 / (0.24 - 0));
	can_4pads.cd();

	pad1_0.cd();
	SetupTopFrame(frame1,pad_wPixel,pad_hPixel);
	gPad->SetLogy();
	frame1->GetXaxis()->SetTitleOffset(0);
	frame1->Draw();
	frame1_i["sr1"]->Draw("same");
	frame_sgn_draw["sr1"]->Draw("same");
	frame1_i["sr2"]->Draw("same");
	frame_sgn_draw["sr2"]->Draw("same");
	frame1_i["sr3"]->Draw("same");
	frame_sgn_draw["sr3"]->Draw("same");
	leg->Draw();

	white_line_up["sr1"]->Draw();
	white_line_up["sr2"]->Draw();
	white_line_dn["sr3"]->Draw();
	white_line_dn["sr2"]->Draw();
	can_4pads.cd();
	style.drawStandardTitle("out");
	pad2_0.cd();
	TLine line(200,0,1700,0);
	line.SetLineStyle(2);
	frame2->SetYTitle("");
	SetupMiddleFrame(frame2,pad1_wPixel,pad1_hPixel);
	frame2->Draw();
	SetupMiddleFrame(frame2_i["sr1"],pad1_wPixel,pad1_hPixel);
	frame2_i["sr1"]->Draw("same");
	line.Draw();
	can_4pads.cd();
	pad2_1.cd();
	frame2->Draw();
	SetupMiddleFrame(frame2_i["sr2"],pad1_wPixel,pad1_hPixel);
	frame2_i["sr2"]->Draw("same");
	line.Draw();
	can_4pads.cd();
	pad2_2.cd();
	SetupBottomFrame(frame2,pad1_wPixel,pad1_hPixel);
 	frame2->Draw();
 	SetupBottomFrame(frame2_i["sr3"],pad1_wPixel,pad1_hPixel);
 	frame2_i["sr3"]->GetXaxis()->SetTickLength(frame2_i["sr3"]->GetXaxis()->GetTickLength() * pad2_0.GetAbsHNDC() / pad2_2.GetAbsHNDC());
	frame2_i["sr3"]->GetYaxis()->SetTickLength(frame2_i["sr3"]->GetYaxis()->GetTickLength() * pad2_2.GetAbsHNDC() / pad2_0.GetAbsHNDC() );
	frame2_i["sr3"]->Draw("same");
	line.Draw();
	can_4pads.cd();
	//Make HardCoded Y-title for all 3 bottom pads
 	TLatex y_axis;
 	y_axis.SetTextAlign(12);
 	y_axis.SetTextAngle(90);
 	y_axis.SetTextFont(frame1->GetYaxis()->GetTitleFont());
 	y_axis.SetTextSize(frame1->GetYaxis()->GetTitleSize());
 	y_axis.DrawLatex(0.06,0.2,"#frac{Data-Bkg.}{#sqrt{Bkg.}}");
	
 	can_4pads.SaveAs("wtf.root");
 	string out_name = output_path + "SpBg_AllInOne_" + to_string(signals.at("sr1")) + "_" + to_string(signals.at("sr2")) + "_" + to_string(signals.at("sr3"));

	can_4pads.Print( (out_name  + ".pdf").c_str());
}

SubRange setup_subranges(const double& bin_width, const int& mass, const std::string& combine_input, const std::string& plots_output){
	/*
	 * Function to get individual combine results for the particular mass point (and therefore sub-range)
	 */
	int nbins;
	const string mass_point = to_string(mass);
	string sr;
	if(find(mssmhbb::sr1.begin(), mssmhbb::sr1.end(),mass) != mssmhbb::sr1.end()) sr = "sr1";
	else if (find(mssmhbb::sr2.begin(), mssmhbb::sr2.end(),mass) != mssmhbb::sr2.end()) sr = "sr2";
	else sr = "sr3";
	double xsec_modified = 0;
	//Files with mll fit results
	TFile fMLL((combine_input + "SpBg/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	//S+Bg mll
	TFile fMLL_bg_nobias((combine_input + "bg_only/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	// Bg only, no-bias
	//Files with post-fit workspaces
	TFile fWorkspace((combine_input + "SpBg/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");
	TFile fWorkspace_bg_nobias((combine_input + "bg_only/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");

	//Get RooFitResults and workspaces
	auto *roofitresult 			= GetFromTFile<RooFitResult>(fMLL,"fit_s");
	auto *roofitresult_nobias 	= GetFromTFile<RooFitResult>(fMLL_bg_nobias,"fit_b");
	auto *workspace				= GetFromTFile<RooWorkspace>(fWorkspace,"MaxLikelihoodFitResult");
	auto *workspace_bg_nobias	= GetFromTFile<RooWorkspace>(fWorkspace_bg_nobias,"MaxLikelihoodFitResult");

	//Get data, pdfs, vars etc
	auto *x 		= GetFromRooWorkspace<RooRealVar>(*workspace,"mbb");
	auto *data_obs 	= GetFromRooWorkspace<RooDataHist>(*workspace,"data_obs");
	//pdfs
	auto *pdf_sgn			= GetFromRooWorkspace<RooAbsPdf>(*workspace, "shapeSig_bbH" + mass_point + "_bbHTo4b");
	auto *pdf_bkg			= GetFromRooWorkspace<RooAbsPdf>(*workspace, "shapeBkg_QCD_Mbb_bbHTo4b");
	auto *pdf_bkg_nobias	= GetFromRooWorkspace<RooAbsPdf>(*workspace_bg_nobias, "shapeBkg_QCD_Mbb_bbHTo4b");
	//pdfs normalisations
	auto *pdf_sgn_norm		= GetFromRooWorkspace<RooFormulaVar>(*workspace, "shapeSig_bbH" + mass_point + "_bbHTo4b__norm");
	//finale norms
	auto *n_sgn_fit			= GetFromRooWorkspace<RooFormulaVar>(*workspace,"n_exp_final_binbbHTo4b_proc_bbH"  + mass_point);
	auto *n_bkg_fit			= GetFromRooWorkspace<RooFormulaVar>(*workspace,"n_exp_final_binbbHTo4b_proc_QCD_Mbb");
	auto *bkg_sys_fit_nobias	= static_cast<RooRealVar*>(roofitresult_nobias->floatParsFinal().find("CMS_bkgd_qcd_13TeV"));
	auto *n_bkg_fit_nobias	= static_cast<RooRealVar*>(roofitresult_nobias->constPars().find("shapeBkg_QCD_Mbb_bbHTo4b__norm"));

	//Derive numbers
	RooRealVar n_sgn 		("n_sgn","",n_sgn_fit->getVal());
	RooRealVar n_bkg 		("n_bkg","",n_bkg_fit->getVal());
	RooRealVar n_bkg_nobias ("n_bkg_nobias","",n_bkg_fit_nobias->getVal() * bkg_sys_fit_nobias->getVal());

	nbins = (x->getMax() - x->getMin()) / bin_width;

	RooAddPdf pdf_comb_ext("pdf_comb_ext","",RooArgList(*pdf_sgn,*pdf_bkg),RooArgList(n_sgn,n_bkg));
	auto &h_rebinned = *data_obs->createHistogram("h_data_obs_rebinned", *x, RooFit::Binning( nbins , x->getMin(), x->getMax()) );
	RooDataHist data_rebinned("data_obs_rebinned","", RooArgList(*x), &h_rebinned, 1.0);

	//TCanvas
	TCanvas c1("c1", "canvas", 600, 600);
	//Size of the main pad in pixels:
	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());
	//TPad
	c1.cd();
	TPad pad1("pad1", "pad1", 0.02, 0.30, 1, 1.0);
    pad1.SetBottomMargin(0.02);
    pad1.SetLeftMargin(0.13);
    pad1.Draw();
    pad1.cd();

	//Frame
    auto *frame1 = x->frame(Name("frame1"));
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle(( ("Events / " + to_string(int(h_rebinned.GetBinWidth(1))) + " GeV").c_str() ));
    SetupTopFrame(frame1,pad_wPixel,pad_hPixel);

    TH1F hdummy("hdummy", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());

    TH1F hbkg("hbkg", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    hbkg.SetFillColor(kAzure+1);
    hbkg.SetLineColor(kAzure+1);
//     hbkg.SetLineWidth(0);
    hbkg.SetFillStyle(3001);

    TH1F h1sigmaU("h1sigmaU", "", x->numBins(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h1sigmaU.SetFillColor(kGreen+1);
    h1sigmaU.SetLineColor(kGreen+1);
    h1sigmaU.SetFillStyle(1001);

    TH1F h1sigmaD("h1sigmaD", "", x->numBins(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h1sigmaD.SetFillColor(kGreen+1);
    h1sigmaD.SetLineColor(kGreen+1);
    h1sigmaD.SetFillStyle(1001);

    TH1F h2sigmaU("h2sigmaU", "", x->numBins(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h2sigmaU.SetFillColor(kOrange);
    h2sigmaU.SetLineColor(kOrange);
    h2sigmaU.SetFillStyle(1001);

    TH1F h2sigmaD("h2sigmaD", "", x->numBins(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h2sigmaD.SetFillColor(kOrange);
    h2sigmaD.SetLineColor(kOrange);
    h2sigmaD.SetFillStyle(1001);

    TLegend leg(0.55,0.48,0.88,0.88,"","brNDC");
    leg.SetHeader( ("m_{A/H} = " + mass_point + " GeV").c_str() );
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.05);
    leg.SetFillColor(10);

    cout<<"START plotting"<<endl;
    cout<<"Plot data_obs..."<<endl;
    data_rebinned.plotOn(frame1, Name("data_obs"));
    cout<<"Plot bkg +/- 2 sigma..."<<endl;

    pdf_bkg->plotOn(frame1,
    		VisualizeError(*roofitresult, 2, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kOrange), Name( (string(pdf_bkg->GetName()) +"_2sigma").c_str()), Normalization(n_bkg.getVal() , RooAbsReal::NumEvent), MoveToBack() );

    cout<<"Plot bkg +/- 1 sigma..."<<endl;
    pdf_bkg->plotOn(frame1,
    		VisualizeError(*roofitresult, 1, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kGreen+1), Name( (string(pdf_bkg->GetName())+"_1sigma").c_str()) , Normalization(n_bkg.getVal() , RooAbsReal::NumEvent) );

    cout<<"Plot bkg..."<<endl;
    pdf_bkg->plotOn(frame1,
    		LineWidth(3),
			LineColor(kBlue),
			LineStyle(kSolid),
			Name(pdf_bkg->GetName()),
			Normalization(n_bkg.getVal(), RooAbsReal::NumEvent) );

    cout<<"Plot sgn+bkg..."<<endl;
    pdf_comb_ext.plotOn(frame1,
    		LineWidth(0.25),
			LineColor(kRed),
			LineStyle(kDashed),
			Name(pdf_comb_ext.GetName()) );

    cout<<"Plot sgn..."<<endl;
    if(sr == "sr1") xsec_modified = xsec_vis * 3000 / pdf_sgn->createIntegral(RooArgSet(*x))->getVal();
    else if (sr == "sr2") xsec_modified = xsec_vis * 160 / pdf_sgn->createIntegral(RooArgSet(*x))->getVal();
	else  xsec_modified = xsec_vis * 90 / pdf_sgn->createIntegral(RooArgSet(*x))->getVal();
	
    pdf_sgn->plotOn(frame1,
    		LineWidth(3),
			LineColor(46),
			LineStyle(kDotted),
            Name(pdf_sgn->GetName()),
			Normalization(xsec_modified) );

    cout<<"Plot data_obs..."<<endl;
    data_rebinned.plotOn(frame1, Name("data_obs"));

    cout<<"END plotting"<<endl;
//    auto chi2_s = frame1->chiSquare(pdf_comb_ext.GetName(), "data_obs", 2+1 if bkg_pdf=='dijet' else 3+1 )
//    chi2_b = frame1.chiSquare(pdf_bkg_b.GetName(), "data_obs", 2 if bkg_pdf=='dijet' else 3 )
//    ndof = (data_rebinned.numEntries()-(2 if bkg_pdf=='dijet' else 3))
//    print ROOT.TMath.Prob(chi2_b*ndof, ndof)
    leg.AddEntry(frame1->findObject("data_obs"), "Data", "PE");
    leg.AddEntry(frame1->getCurve(pdf_bkg->GetName()),  ("Bkg."), "L");
    leg.AddEntry(&h1sigmaU, "#pm1 std. deviation", "F");
    leg.AddEntry(&h2sigmaU, "#pm2 std. deviation", "F");
    leg.AddEntry(frame1->getCurve(pdf_comb_ext.GetName()), ("Bkg. + signal"), "L");
    leg.AddEntry(frame1->getCurve(pdf_sgn->GetName()), ("Signal, #sigma= " + to_string_with_precision(xsec_vis,2) + " pb" ).c_str(), "L");
    frame1->SetMaximum( h_rebinned.GetMaximum()*2.0 );
    frame1->SetMinimum( h_rebinned.GetMinimum()*0.8 );
    pad1.SetLogy();
    frame1->Draw();
    leg.Draw();

    c1.Modified();

    c1.cd();
    TPad pad2("pad2", "pad2", 0.02, 0.02, 1, 0.28);
    pad2.SetTopMargin(0.02);
    pad2.SetLeftMargin(0.13);
    pad2.SetBottomMargin(0.35);
    pad2.SetGridx();
    pad2.SetGridy();
    pad2.Draw();
    pad2.cd();

    auto *frame2 = x->frame(Name("frame2"));
//     SetupBottomFrame(frame2);

    auto *hresid = frame1->pullHist("data_obs", pdf_bkg->GetName());
    hresid->GetYaxis()->SetRangeUser(-4.,4.);

    auto *sigma2 	= frame1->getCurve( (string(pdf_bkg->GetName())+"_2sigma").c_str());
    auto *sigma1 	= frame1->getCurve( (string(pdf_bkg->GetName())+"_1sigma").c_str());
    auto *nominal 	= frame1->getCurve(pdf_bkg->GetName());

    TGraph up1Bound(nominal->GetN());
    TGraph lo1Bound(nominal->GetN());
    TGraph up2Bound(nominal->GetN());
    TGraph lo2Bound(nominal->GetN());
    double err1U,err1D,err2U,err2D;
    bool divide_by_sqrt_bg = true;
    std::cout<<"NUMBER OF BINS: "<<nominal->GetN()<<std::endl;

    for( int j = 0; j < sigma2->GetN(); ++j ){
    	if( j < nominal->GetN() ){
    		up1Bound.SetPoint(j, sigma1->GetX()[j], sigma1->GetY()[j]);
    		up2Bound.SetPoint(j, sigma2->GetX()[j], sigma2->GetY()[j]);
    	}
    	else{
    		lo1Bound.SetPoint(j, sigma1->GetX()[j], sigma1->GetY()[j]);
    		lo2Bound.SetPoint(j, sigma2->GetX()[j], sigma2->GetY()[j]);
    	}
    }

    int k = 0;
    std::cout<<"SUB-RANGE: "<<sr<<std::endl;
    for( int j = 0; j < nominal->GetN(); ++j ){
    		double x = nominal->GetX()[j];
    	//        std::cout<<"WTFFF: x = "<<x<<" y = "<<nominal->Eval(x)<<" up1 = "<<up1Bound.Eval(x)<<" lo1 = "<<lo1Bound.Eval(x)<<" up2 = "<<up2Bound.Eval(x)<<" lo2 = "<<lo2Bound.Eval(x)<<std::endl;
    		auto bin_i = h1sigmaU.FindBin( x );
    		if(bin_i <= 0 || bin_i > h1sigmaU.GetNbinsX()) continue;
    		++k;

    		auto n_i = nominal->Eval(x);
    		auto n_1up_i = up1Bound.Eval(x);
    		auto n_1down_i = lo1Bound.Eval(x);
    		auto n_2up_i = up2Bound.Eval(x);
    		auto n_2down_i = lo2Bound.Eval(x);

    		std::cout<<"j = "<<j<<" k = "<<k<<"x = "<<x<<" y = "<<n_i<<" up1 = "<<n_1up_i<<" up2 = "<<n_2up_i<<std::endl;

    		//Fill band histograms
    		err1U = (-n_i+n_1up_i)  ;
    		err1D = (-n_i+n_1down_i);
    		err2U = (-n_i+n_2up_i)  ;
    		err2D = (-n_i+n_2down_i);

    		if(divide_by_sqrt_bg){
    			err1U /= sqrt(n_i);
    			err1D /= sqrt(n_i);
    			err2U /= sqrt(n_i);
    			err2D /= sqrt(n_i);
    		}

    		h1sigmaU.SetBinContent( k, err1U);
    		h1sigmaD.SetBinContent( k, err1D);
    		h2sigmaU.SetBinContent( k, err2U);
    		h2sigmaD.SetBinContent( k, err2D);
    }

    //redefine pulls here:
    auto bkg_int = pdf_bkg->createIntegral(RooArgSet(*x))->getVal();
    for(int bin_i = 1; bin_i < h_rebinned.GetNbinsX()+1; ++bin_i){
        auto n_data_i = h_rebinned.GetBinContent( bin_i );
        auto e_data_i = h_rebinned.GetBinError( bin_i );
        auto bin_width_i = h_rebinned.GetBinWidth( bin_i );
        x->setRange( ("Bin" + to_string(bin_i)).c_str(), h_rebinned.GetBinLowEdge(bin_i), h_rebinned.GetBinLowEdge(bin_i)+bin_width_i);
        auto bkg_int_i = pdf_bkg->createIntegral(RooArgSet(*x), ("Bin" + to_string(bin_i)).c_str())->getVal() / bkg_int * n_bkg.getVal();
        double divide_by = sqrt(bkg_int_i);
        auto pull_i = (n_data_i-bkg_int_i);

        cout<<"Data: "<<n_data_i<<" ---- Bkg: "<<bkg_int_i<<endl;
        if(divide_by_sqrt_bg){
        		pull_i /= divide_by;
        		e_data_i /= divide_by;
        }
        hbkg.SetBinContent( bin_i, pull_i);
		hbkg.SetBinError( bin_i, e_data_i);
        hdummy.SetBinContent( bin_i, 0.);
        hdummy.SetBinError( bin_i, 0.);
    }


    frame2->addTH1(&h2sigmaU, "HIST");
    frame2->addTH1(&h2sigmaD, "HIST");
    frame2->addTH1(&h1sigmaU, "HIST");
    frame2->addTH1(&h1sigmaD, "HIST");
	frame2->addTH1(&hbkg, "HIST");
	frame2->addTH1(&hdummy, "HIST");
	frame2->Draw();
	frame2->GetYaxis()->SetRangeUser(-4,4);
	c1.cd();
	c1.Draw();

// 	c1.Print( (plots_output + "Post_fit_" + sr + "_m_" + mass_point + ".pdf").c_str());

	//Adjust names of the pdfs for simplification
	SubRange subrange;
	x->SetName("mbb");
	pdf_bkg->SetName("pdf_bkg");
	pdf_bkg_nobias->SetName("pdf_bkg_nobias");
	pdf_comb_ext.SetName("pdf_comb_ext");
	pdf_sgn->SetName("pdf_sgn");
	data_obs->SetName("data_obs");
	pdf_sgn_norm->SetName("pdf_sgn_norm");
	data_rebinned.SetName("data_rebinned");
	RooWorkspace *output_wp = new RooWorkspace();
	output_wp->import(*x);
	output_wp->import(n_sgn);
	output_wp->import(n_bkg);
	output_wp->import(n_bkg_nobias);
	output_wp->import(*pdf_bkg);
	output_wp->import(*pdf_bkg_nobias,RenameConflictNodes("nobias"));
	output_wp->import(pdf_comb_ext,RenameConflictNodes("ext"));
	output_wp->import(*pdf_sgn,RenameConflictNodes("sgn"));
	output_wp->import(*data_obs);
	output_wp->import(data_rebinned);
	output_wp->import(*pdf_sgn_norm);
	subrange.h1sigmaD 		= h1sigmaD;
	subrange.h1sigmaU 		= h1sigmaU;
	subrange.h2sigmaD 		= h2sigmaD;
	subrange.h2sigmaU 		= h2sigmaU;
	subrange.hbkg     		= hbkg;
	subrange.sub_range		= sr;
	subrange.workspace		= output_wp;
	subrange.fMLL			= combine_input + "SpBg/mll_M-" + mass_point + "/mlfit.root";
	subrange.fMLL_bg_nobias	= combine_input + "bg_only/mll_M-" + mass_point + "/mlfit.root";
	subrange.xsec_vis = 1;
	if(TEST_) std::cout<<"Out of setup_subranges"<<std::endl;

	return subrange;
}

void CollectSubRangesData(std::map<std::string,SubRange>& data,const double& bin_width,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path){
	/*
	 * Method to collect data from sub-ranges
	 */
	for(const auto & s : signals){
		//Loop over the subranges to collect the data
		data[s.first] = setup_subranges(bin_width, s.second, combine_input, output_path);
	}
	if(TEST_) std::cout<<"Out of CollectSubRangesData"<<std::endl;
}

void SetupTopFrame(RooPlot *frame1, const double& pad_w, const double& pad_h){

	//Pixel sizes of the current pad
	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());
	pad_wPixel += pad_w; pad_hPixel += 0;

	//Use relative values
//	frame1->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y"));
//	frame1->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y"));
//	frame1->GetYaxis()->SetRangeUser(frame1->GetMinimum(), frame1->GetMaximum()+200);
//
//	frame1->GetXaxis()->SetTickLength(gStyle->GetTickLength("X") * pad_h / pad_hPixel);
//	frame1->GetXaxis()->SetTitleSize(0);
//	frame1->GetXaxis()->SetLabelSize(0);
//	frame1->GetXaxis()->SetTitleOffset(999); //Effectively turn off x axis title on main plot
//	frame1->GetXaxis()->SetLabelOffset(999); //Effectively turn off x axis label on main plot

	frame1->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y")*0.25);
	frame1->GetYaxis()->SetLabelFont(43);
	frame1->GetYaxis()->SetTitleFont(43);
//    frame1->GetYaxis()->SetTitleSize(24);
    frame1->GetYaxis()->SetTitleSize(30);
    frame1->GetYaxis()->SetTitleOffset(1.25);
//    frame1->GetYaxis()->SetLabelSize(20);
    frame1->GetYaxis()->SetLabelSize(24);
    frame1->GetYaxis()->SetLabelOffset(0.007);
    frame1->GetXaxis()->SetTitleSize(0);
    frame1->GetXaxis()->SetLabelSize(0);
    frame1->GetXaxis()->SetTickLength(gStyle->GetTickLength("X"));
}

void SetupMiddleFrame(RooPlot *frame2, const double& pad_w, const double& pad_h){
    frame2->SetTitle("");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetTitleSize(30);
    frame2->GetYaxis()->SetTitleFont(43);
    frame2->GetYaxis()->SetLabelFont(43);
    frame2->GetYaxis()->SetTitleOffset(1.25);
    frame2->GetYaxis()->SetLabelOffset(0.007);
    frame2->GetYaxis()->SetNdivisions(105);
    frame2->GetYaxis()->SetLabelSize(24);
    frame2->GetYaxis()->SetTickLength(0.015);

    frame2->GetXaxis()->SetTitle("");
    frame2->GetXaxis()->SetTitleFont(43);
    frame2->GetXaxis()->SetLabelFont(43);
//    frame2->GetXaxis()->SetTitleSize(25);
    frame2->GetXaxis()->SetTitleSize(30);
    frame2->GetXaxis()->SetTitleOffset(3.9);
    frame2->GetXaxis()->SetLabelOffset(0.04);
    frame2->GetXaxis()->SetLabelSize(24);
    frame2->GetXaxis()->SetTickLength(0.1);


//	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
//	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());

	//Y-axis options
//	frame2->GetYaxis()->SetTitleFont(43);			// This gives the sizes in pixels!!!!
//	frame2->GetYaxis()->SetLabelFont(43);			// This gives the sizes in pixels!!!!
//	frame2->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y") * pad_w / pad_wPixel);
//	frame2->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y") * pad_w / pad_wPixel );
//	frame2->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y") * 1.15);
//	frame2->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y") * pad_h);
//	frame2->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y") * pad_h);
//    frame2->GetYaxis()->CenterTitle();
//    frame2->GetYaxis()->SetNdivisions(105);
//    frame2->SetMinimum(-3.9999);
//    frame2->SetMaximum(+3.9999);
//
//    //X-axis options
//	frame2->GetXaxis()->SetTitleFont(43);			// This gives the sizes in pixels!!!!
//	frame2->GetXaxis()->SetLabelFont(43);			// This gives the sizes in pixels!!!!
//    frame2->GetXaxis()->SetTickLength(gStyle->GetTickLength("X") * pad_h / pad_hPixel);
//    frame2->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X") * pad_h);
//    frame2->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X") * pad_h);
//    	frame2->GetXaxis()->SetTitleOffset(3.5);
//    	frame2->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X") * pad_h / pad_hPixel );
}

void SetupBottomFrame(RooPlot *frame2, const double& pad_w, const double& pad_h){
	SetupMiddleFrame(frame2,pad_w,pad_h);
    frame2->GetXaxis()->SetTitle(HbbStyle::axisTitleMass());
}

void SetupBottomPad(TPad * pad){
	pad->SetTopMargin(0.01);
//	pad->SetLeftMargin(0.13);
	pad->SetBottomMargin(0.0);
// 	pad->SetGridx();
// 	pad->SetGridy();
	pad->Draw();
}
