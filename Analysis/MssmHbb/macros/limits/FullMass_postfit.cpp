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

//void setup_plotting(const int& nbins, const std::map<std::string,int>& signals,const std::map<int,double>& signal_sfs,const std::string& output_folder);

struct SubRange{
	RooWorkspace workspace;
	TH1F h2sigmaU, h2sigmaD, h1sigmaU, h1sigmaD, hbkg;
	std::string sub_range;
};

//typedef std::pair<SubRange,RooWorkspace> SubRangeData;

void SetupTopFrame(RooPlot *frame);
void SetupBottomFrame(RooPlot *frame2);

SubRange setup_subranges(const int& nbins, const int& mass_point, const std::string& combine_input, const std::string& output_path);
void CollectSubRangesData(std::map<std::string,SubRange>& data,const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path);

void PlotFullMassFit(const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path);


int main(){

	int nbins = 30;
	string output_folder = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170818/PostFits/";
	CheckOutputDir(output_folder);
	string input_folder_combine = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/201707/26/";
	if(!mssmhbb::blinded) input_folder_combine += "unblinded/";
	input_folder_combine += "independent/mll/";

	map<string,int> signals = {
			{"sr1" , 300},
			{"sr2" , 600},
			{"sr3" , 1300}
	};

	PlotFullMassFit(nbins,signals,input_folder_combine,output_folder);

//	SubRange data = setup_subranges(nbins, 300, input_folder_combine, output_folder);
//	setup_plotting(nbins,signals,signal_sfs,output_folder);

	return 0;
}

void PlotFullMassFit(const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path){
	/*
	 * Method to plot the full mass fit
	 */
	std::map<std::string,SubRange> data;
	CollectSubRangesData(data,nbins,signals,combine_input,output_path);

	RooRealVar x("mbb","mbb",200,1700);
	auto &x1 = *data["sr1"].workspace.var("x");
	auto &x2 = *data["sr2"].workspace.var("x");
	auto &x3 = *data["sr3"].workspace.var("x");

	auto &data_obs1 = *data["sr1"].workspace.data("data_obs");
	auto &data_obs2 = *data["sr2"].workspace.data("data_obs");
	auto &data_obs3 = *data["sr3"].workspace.data("data_obs");

	auto &pdf_bkg1 = *data["sr1"].workspace.pdf("pdf_bkg");
	auto &pdf_bkg2 = *data["sr2"].workspace.pdf("pdf_bkg");
	auto &pdf_bkg3 = *data["sr3"].workspace.pdf("pdf_bkg");

	auto &pdf_bkg_nobias1 = *data["sr1"].workspace.pdf("pdf_bkg_nobias");
	auto &pdf_bkg_nobias2 = *data["sr2"].workspace.pdf("pdf_bkg_nobias");
	auto &pdf_bkg_nobias3 = *data["sr3"].workspace.pdf("pdf_bkg_nobias");

	auto &pdf_comb_ext1 = *data["sr1"].workspace.pdf("pdf_comb_ext");
	auto &pdf_comb_ext2 = *data["sr2"].workspace.pdf("pdf_comb_ext");
	auto &pdf_comb_ext3 = *data["sr3"].workspace.pdf("pdf_comb_ext");

	auto &pdf_sgn1 = *data["sr1"].workspace.pdf("pdf_sgn");
	auto &pdf_sgn2 = *data["sr2"].workspace.pdf("pdf_sgn");
	auto &pdf_sgn3 = *data["sr3"].workspace.pdf("pdf_sgn");

	//TOP FRAME
	//TCanvas
	TCanvas c1("c1", "canvas", 600, 600);
	//TPad
	c1.cd();
	TPad pad1("pad1", "pad1", 0.02, 0.30, 1, 1.0);
    pad1.SetBottomMargin(0.02);
    pad1.SetLeftMargin(0.13);
    pad1.Draw();
    pad1.cd();

	//Frame
    double bin_width = (1700.0-20.0)/nbins;//GeV
    auto *frame = x.frame(Name("frame"),200,1700);
    frame->SetTitle("");
    frame->GetYaxis()->SetTitle(( ("Events / " + to_string_with_precision(bin_width,1) + " GeV").c_str() ));
    SetupTopFrame(frame);
    auto *frame1 = x1.frame(Name("frame1"),Range(200,650));
    auto *frame2 = x2.frame(Name("frame2"),Range(350,1190));
    auto *frame3 = x3.frame(Name("frame3"),Range(500,1700));

    pdf_bkg_nobias1->plotOn(frame1,
//    		VisualizeError(*roofitresult_nobias, 2, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kOrange), Name( (string(pdf_bkg_nobias->GetName()) +"_2sigma").c_str()), Normalization(n_bkg_nobias->getVal() , RooAbsReal::NumEvent), MoveToBack() );

}

SubRangeData setup_subranges(const int& nbins, const int& mass, const std::string& combine_input, const std::string& plots_output){
	/*
	 * Function to get individual combine results for the particular mass point (and therefore sub-range)
	 */
	double xsec_vis = 1;
	const string mass_point = to_string(mass);
	//Files with mll fit results
	TFile fMLL((combine_input + "SpBg/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	//S+Bg mll
	TFile fMLL_bg_nobias((combine_input + "bg_only/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	// Bg only, no-bias
	//Files with post-fit workspaces
	TFile fWorkspace((combine_input + "SpBg/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");
	TFile fWorkspace_bg_nobias((combine_input + "bg_only/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");

	//Get RooFitResults and workspaces
//	auto *roofitresult 			= GetFromTFile<RooFitResult>(fMLL,"fit_s");
	auto *roofitresult_nobias 	= GetFromTFile<RooFitResult>(fMLL_bg_nobias,"fit_b");
	auto *workspace				= GetFromTFile<RooWorkspace>(fWorkspace,"MaxLikelihoodFitResult");
	auto *workspace_bg_nobias	= GetFromTFile<RooWorkspace>(fWorkspace_bg_nobias,"MaxLikelihoodFitResult");

	//Get data, pdfs, vars etc
	auto *x 		= GetFromRooWorkspace<RooRealVar>(*workspace,"mbb");
//	auto *r 		= GetFromRooWorkspace<RooRealVar>(*workspace,"r");
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
	auto *bkg_sys_fit_nobias= static_cast<RooRealVar*>(roofitresult_nobias->floatParsFinal().find("CMS_bkgd_qcd_13TeV"));
	auto *n_bkg_fit_nobias	= static_cast<RooRealVar*>(roofitresult_nobias->constPars().find("shapeBkg_QCD_Mbb_bbHTo4b__norm"));

	//Derive numbers
	RooRealVar n_sgn 		("n_sgn","",n_sgn_fit->getVal());
	RooRealVar n_bkg 		("n_bkg","",n_bkg_fit->getVal());
	RooRealVar n_bkg_nobias ("n_bkg_nobias","",n_bkg_fit_nobias->getVal() * bkg_sys_fit_nobias->getVal());

	RooAddPdf pdf_comb_ext("pdf_comb_ext","",RooArgList(*pdf_sgn,*pdf_bkg),RooArgList(n_sgn,n_bkg));
	auto &h_rebinned = *data_obs->createHistogram("h_data_obs_rebinned", *x, RooFit::Binning( nbins , x->getMin(), x->getMax()) );
	RooDataHist data_rebinned("data_obs_rebinned","", RooArgList(*x), &h_rebinned, 1.0);

	//TCanvas
	TCanvas c1("c1", "canvas", 600, 600);
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
    frame1->GetYaxis()->SetTitle(( ("Events / " + to_string_with_precision(h_rebinned.GetBinWidth(1),1) + " GeV").c_str() ));
    SetupTopFrame(frame1);

    TH1F hdummy("hdummy", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());

    TH1F hbkg("hbkg", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    hbkg.SetFillColor(kAzure+1);
    hbkg.SetLineColor(kAzure+1);
    hbkg.SetLineWidth(0);
    hbkg.SetFillStyle(3001);

    TH1F h1sigmaU("h1sigmaU", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h1sigmaU.SetFillColor(kGreen+1);
    h1sigmaU.SetFillStyle(1001);

    TH1F h1sigmaD("h1sigmaD", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h1sigmaD.SetFillColor(kGreen+1);
    h1sigmaD.SetFillStyle(1001);

    TH1F h2sigmaU("h2sigmaU", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h2sigmaU.SetFillColor(kOrange);
    h2sigmaU.SetFillStyle(1001);

    TH1F h2sigmaD("h2sigmaD", "", h_rebinned.GetNbinsX(), h_rebinned.GetXaxis()->GetXmin(), h_rebinned.GetXaxis()->GetXmax());
    h2sigmaD.SetFillColor(kOrange);
    h2sigmaD.SetFillStyle(1001);

    TLegend leg(0.55,0.48,0.88,0.88,"","brNDC");
    leg.SetHeader( ("m_{#phi} = " + mass_point + " GeV").c_str() );
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetTextSize(0.05);
    leg.SetFillColor(10);

    cout<<"START plotting"<<endl;
    cout<<"Plot data_obs..."<<endl;
    data_rebinned.plotOn(frame1, Name("data_obs"));
    cout<<"Plot bkg +/- 2 sigma..."<<endl;

    pdf_bkg_nobias->plotOn(frame1,
    		VisualizeError(*roofitresult_nobias, 2, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kOrange), Name( (string(pdf_bkg_nobias->GetName()) +"_2sigma").c_str()), Normalization(n_bkg_nobias.getVal() , RooAbsReal::NumEvent), MoveToBack() );

    cout<<"Plot bkg +/- 1 sigma..."<<endl;
    pdf_bkg_nobias->plotOn(frame1,
    		VisualizeError(*roofitresult_nobias, 1, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kGreen+1), Name( (string(pdf_bkg_nobias->GetName())+"_1sigma").c_str()) , Normalization(n_bkg_nobias.getVal() , RooAbsReal::NumEvent) );

    cout<<"Plot bkg..."<<endl;
    pdf_bkg_nobias->plotOn(frame1,
    		LineWidth(3),
			LineColor(kBlue),
			LineStyle(kSolid),
			Name(pdf_bkg_nobias->GetName()),
			Normalization(n_bkg_nobias.getVal(), RooAbsReal::NumEvent) );

    cout<<"Plot sgn+bkg..."<<endl;
    pdf_comb_ext.plotOn(frame1,
    		LineWidth(3),
			LineColor(kRed),
			LineStyle(kDashed),
			Name(pdf_comb_ext.GetName()) );

    cout<<"Plot sgn..."<<endl;
    pdf_sgn->plotOn(frame1,
    		LineWidth(3),
			LineColor(46),
			LineStyle(kDotted),
            Name(pdf_sgn->GetName()),
			Normalization(pdf_sgn_norm->getVal()*xsec_vis , RooAbsReal::NumEvent) );

    cout<<"Plot data_obs..."<<endl;
    data_rebinned.plotOn(frame1, Name("data_obs"));

    cout<<"END plotting"<<endl;
//    auto chi2_s = frame1->chiSquare(pdf_comb_ext.GetName(), "data_obs", 2+1 if bkg_pdf=='dijet' else 3+1 )
//    chi2_b = frame1.chiSquare(pdf_bkg_b.GetName(), "data_obs", 2 if bkg_pdf=='dijet' else 3 )
//    ndof = (data_rebinned.numEntries()-(2 if bkg_pdf=='dijet' else 3))
//    print ROOT.TMath.Prob(chi2_b*ndof, ndof)
    leg.AddEntry(frame1->findObject("data_obs"), "Data", "PE");
    leg.AddEntry(frame1->getCurve(pdf_bkg_nobias->GetName()),  ("Bkg."), "L");
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
    SetupBottomFrame(frame2);

    auto *hresid = frame1->pullHist("data_obs", pdf_bkg_nobias->GetName());
    hresid->GetYaxis()->SetRangeUser(-4,4);

    auto *sigma2 = frame1->getCurve( (string(pdf_bkg_nobias->GetName())+"_2sigma").c_str());
    auto *sigma1 = frame1->getCurve( (string(pdf_bkg_nobias->GetName())+"_1sigma").c_str());
    auto *nominal = frame1->getCurve(pdf_bkg_nobias->GetName());
    auto *nominalX = nominal->GetX();
    auto *nominalY = nominal->GetY();
    auto *n1sigmaY = sigma1->GetY();
    auto *n2sigmaY = sigma2->GetY();

    for(int i = 0; i < nominal->GetN(); ++i){
        auto x_i = nominalX[i];
        auto n_i = nominalY[i];
        auto n_1up_i = n1sigmaY[2*nominal->GetN()-i-1];
        auto n_1down_i = n1sigmaY[i];
        auto n_2up_i = n2sigmaY[2*nominal->GetN()-i-1];
        auto n_2down_i = n2sigmaY[i];
        auto bin_i = h1sigmaU.FindBin( x_i );
        if(bin_i <= 0 || bin_i > h1sigmaU.GetNbinsX()) continue;
//        auto bin_width_i = h_rebinned.GetBinWidth( bin_i );
//        auto n_data_i = h_rebinned.GetBinContent( bin_i );
        h1sigmaU.SetBinContent( bin_i, (-n_i+n_1up_i));
        h1sigmaD.SetBinContent( bin_i, (-n_i+n_1down_i));
        h2sigmaU.SetBinContent( bin_i, (-n_i+n_2up_i));
        h2sigmaD.SetBinContent( bin_i, (-n_i+n_2down_i));
        if(x_i < 600) cout<<"bin"<<i<<" x = "<<x_i<<" y = "<<n_i<<" +-1 = "<<n_1up_i<<" +-2 = "<<n_2up_i<<endl;
    }

    //redefine pulls here:
    auto bkg_int = pdf_bkg_nobias->createIntegral(RooArgSet(*x))->getVal();
    for(int bin_i = 1; bin_i < h_rebinned.GetNbinsX()+1; ++bin_i){
        auto n_data_i = h_rebinned.GetBinContent( bin_i );
        auto bin_width_i = h_rebinned.GetBinWidth( bin_i );
        x->setRange( ("Bin" + to_string(bin_i)).c_str(), h_rebinned.GetBinLowEdge(bin_i), h_rebinned.GetBinLowEdge(bin_i)+bin_width_i);
        auto bkg_int_i = pdf_bkg_nobias->createIntegral(RooArgSet(*x), ("Bin" + to_string(bin_i)).c_str())->getVal() / bkg_int * n_bkg_nobias.getVal();
        auto pull_i = (n_data_i-bkg_int_i)/sqrt(bkg_int_i);
        cout<<"Data: "<<n_data_i<<" ---- Bkg: "<<bkg_int_i<<endl;
        hbkg.SetBinContent( bin_i, pull_i);
        hbkg.SetBinError( bin_i, 0.);
        hdummy.SetBinContent( bin_i, 0.);
        hdummy.SetBinError( bin_i, 0.);
        h1sigmaU.SetBinContent( bin_i, h1sigmaU.GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h1sigmaD.SetBinContent( bin_i, h1sigmaD.GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h2sigmaU.SetBinContent( bin_i, h2sigmaU.GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h2sigmaD.SetBinContent( bin_i, h2sigmaD.GetBinContent( bin_i )/sqrt( bkg_int_i ) );
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

	string sr;
	if(find(mssmhbb::sr1.begin(), mssmhbb::sr1.end(),mass) != mssmhbb::sr1.end()) sr = "sr1";
	else if (find(mssmhbb::sr2.begin(), mssmhbb::sr2.end(),mass) != mssmhbb::sr2.end()) sr = "sr2";
	else sr = "sr3";
	c1.Print( (plots_output + "Post_fit_" + sr + "_m_" + mass_point + ".pdf").c_str());

	//Adjust names of the pdfs for simplification
	x->SetName("x");
	pdf_bkg->SetName("pdf_bkg");
	pdf_bkg_nobias->SetName("pdf_bkg_nobias");
	pdf_comb_ext->SetName("pdf_comb_ext");
	pdf_sgn->SetName("pdf_sgn");
	data_obs->SetName("data_obs");
	RooWorkspace output_wp;
	output_wp.import(*x);
	output_wp.import(*pdf_bkg);
	output_wp.import(*pdf_bkg_nobias);
	output_wp.import(pdf_comb_ext);
	output_wp.import(*pdf_sgn);
	output_wp.import(*data_obs);
	SubRange subrange;
	subrange.h1sigmaD 		= h1sigmaD;
	subrange.h1sigmaU 		= h1sigmaU;
	subrange.h2sigmaD 		= h2sigmaD;
	subrange.h2sigmaU 		= h2sigmaU;
	subrange.hbkg     		= hbkg;
	subrange.sub_range		= sr;
	subrange.workspace		= output_wp;
//	SubRangeData data = std::make_pair(subrange, output_wp);

	return subrange;
}

void CollectSubRangesData(std::map<std::string,SubRange>& data,const int& nbins,const std::map<std::string,int>& signals, const std::string& combine_input, const std::string& output_path){
	/*
	 * Method to collect data from sub-ranges
	 */
	for(const auto & s : signals){
		//Loop over the subranges to collect the data
		data[s.first] = setup_subranges(nbins, s.second, combine_input, output_path);
	}
}

//void setup_plotting(const int& nbins, const std::map<std::string,int>& signals,const std::map<int,double>& signal_sfs,const std::string& output_folder){
//	/*
//	 *
//	 */
//	HbbStyle style;
//	style.set(PRELIMINARY);
//
//	TFile fbg_sr1(mssmhbb::path_bg_sr1.c_str(),"read");
//	TFile fbg_sr2(mssmhbb::path_bg_sr2.c_str(),"read");
//	TFile fbg_sr3(mssmhbb::path_bg_sr3.c_str(),"read");
//
//	auto &bg_pdf1 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr1,"background");
//	auto &bg_pdf2 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr2,"background");
//	auto &bg_pdf3 = *GetRooObjectFromTFile<RooAbsPdf>(fbg_sr3,"background");
//
//	auto &data1 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr1,"data_container");
//	auto &data2 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr2,"data_container");
//	auto &data3 = *GetRooObjectFromTFile<RooDataHist>(fbg_sr3,"data_container");
//
//	auto &x1 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr1,"mbb");
//	auto &x2 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr2,"mbb");
//	auto &x3 = *GetRooObjectFromTFile<RooRealVar>(fbg_sr3,"mbb");
//
//	RooRealVar x("mbb","mbb",200,1700);
//	x.setRange("sr1",200,650);
//	x.setRange("sr2",350,1190);
//	x.setRange("sr3",500,1700);
//
//	TCanvas can("can","can",800,600);
//	RooPlot *frame_sr1 = x1.frame(200,650);
//	RooPlot *frame_sr2 = x2.frame(350,1190);
//	RooPlot *frame_sr3 = x3.frame(500,1700);
//	RooPlot *frame = x.frame();
//	data1.plotOn(frame,Binning(150),Name("data"));
//	bg_pdf1.plotOn(frame_sr1,LineColor(kRed),Normalization(data1.sumEntries("1","fit_range")),Name("bg_pdf1"));
//	bg_pdf2.plotOn(frame_sr2,LineColor(kBlue),Normalization(data2.sumEntries("1","fit_range")/2),Name("bg_pdf2"));
//	bg_pdf3.plotOn(frame_sr3,LineColor(kGreen),Normalization(data3.sumEntries("1","fit_range")/2.5,RooAbsReal::NumEvent),Name("bg_pdf3"));
//
//	TLegend leg(0.6,0.6,0.9,0.9);
//	leg.AddEntry(frame->findObject("data"),"Data","p");
//	leg.AddEntry(frame_sr1->findObject("bg_pdf1"),"SuperNovosibirsk-1, sub-range1","l");
//	leg.AddEntry(frame_sr2->findObject("bg_pdf2"),"SuperNovosibirsk-0, sub-range2","l");
//	leg.AddEntry(frame_sr3->findObject("bg_pdf3"),"SuperNovosibirsk-0, sub-range3","l");
//	leg.SetBorderSize(0);
//	leg.SetFillColor(0);
//	//Draw signals
//	map<string,int> colours = {
//			{"sr1",2},
//			{"sr2",4},
//			{"sr3",3}
//	};
//	for(const auto& s : signals){
//		TFile fsg(mssmhbb::signal_fits.at(s.second).c_str(),"read");
//		auto &sg_pdf = *GetRooObjectFromTFile<RooAbsPdf>(fsg,"signal");
//		sg_pdf.plotOn(frame,LineColor(colours.at(s.first)),LineStyle(2),Normalization(signal_sfs.at(s.second)),Range(s.first.c_str()),Name( ("sg_" + to_string(s.second)).c_str() ));
//		leg.AddEntry(frame->findObject(("sg_" + to_string(s.second)).c_str()),("bb#phi(" + to_string(s.second) + " GeV)").c_str(),"l");
//	}
//	frame->SetMinimum(1);
//	frame->SetTitle(";M_{12} [GeV]");
//	frame->Draw();
//	frame_sr1->Draw("same");
//	frame_sr2->Draw("same");
//	frame_sr3->Draw("same");
//
//	leg.Draw();
//
//	style.drawStandardTitle();
//
//	can.Print( (output_folder + "plot.pdf").c_str()	);
//	gPad->SetLogy();
//	frame->SetMaximum(5*1e+05);
//	can.Print( (output_folder + "plot_log.pdf").c_str()	);
//
//}

void SetupTopFrame(RooPlot *frame1){
    frame1->GetYaxis()->SetTitleSize(24);
    frame1->GetYaxis()->SetTitleFont(43);
    frame1->GetYaxis()->SetTitleOffset(1.18);
    frame1->GetYaxis()->SetLabelFont(43);
    frame1->GetYaxis()->SetLabelSize(20);
    frame1->GetXaxis()->SetTitleSize(0);
    frame1->GetXaxis()->SetTitleFont(43);
    frame1->GetXaxis()->SetLabelFont(43);
    frame1->GetXaxis()->SetLabelSize(0);
}

void SetupBottomFrame(RooPlot *frame2){
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle("#frac{Data-Bkg.}{#sqrt{Bkg.}}");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetYaxis()->SetTitleSize(24);
    frame2->GetYaxis()->SetTitleFont(43);
    frame2->GetYaxis()->SetTitleOffset(1.18);
    frame2->GetYaxis()->SetLabelFont(43);
    frame2->GetYaxis()->SetLabelSize(20);
    frame2->GetXaxis()->SetTitle("m_{b#bar{b}} (GeV)");
    frame2->GetXaxis()->SetTitleSize(25);
    frame2->GetXaxis()->SetTitleFont(43);
    frame2->GetXaxis()->SetTitleOffset(3.5);
    frame2->GetXaxis()->SetLabelFont(43);
    frame2->GetXaxis()->SetLabelSize(20);
}

