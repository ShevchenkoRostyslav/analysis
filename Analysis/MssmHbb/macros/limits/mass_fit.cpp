/*
 * Macro to plot post-fit mass distribution
 * author: shevchen
 * notes:WORKS ONLY with the same CMSSW realease
 * as used to run combine tool. otherwise - couldn't
 * load libraries
 */

#include <string>

#include <TFile.h>
#include <TH1F.h>

#include <RooFit.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooAbsPdf.h>
#include <RooFormulaVar.h>

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"

using namespace std;
using namespace RooFit;
using namespace analysis;

HbbStyle  style;
const string cmsswBase = gSystem->Getenv("CMSSW_BASE");//"/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1";//

TH1* getDataHist(RooWorkspace& w, TH1& binning);
void draw(TH1F* data, TH1F *total, TH1F *signal, TH1F *bg);

//int main(int argc, char **argv) {
int mass_fit(){
	/*
	 * Macro to plot post-fit M_{12} fit
	 */
	gSystem->Load("libHiggsAnalysisCombinedLimit");
	style.set(PRIVATE);

	string mass_point = "1100";
	int nbins = 50;
	double xsec_vis = 10;
	string output = "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/20170518/";
	string output_append = "";
	string campaign = "201705/18/mll/independent/";

	//Files with mll fit results
	TFile *fMLL 		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "SpBg/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	//S+Bg mll
	TFile *fMLL_bg_nobias 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "bg_only/test_no_bias/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	// Bg only, no-bias
	//Files with post-fit workspaces
	TFile *fWorkspace			= new TFile((cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "SpBg/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");
	TFile *fWorkspace_bg_nobias 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "bg_only/test_no_bias/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");

	//Get RooFitResults and workspaces
	auto *roofitresult 			= GetFromTFile<RooFitResult>(*fMLL,"fit_s");
	auto *roofitresult_nobias 	= GetFromTFile<RooFitResult>(*fMLL_bg_nobias,"fit_b");
	auto *workspace				= GetFromTFile<RooWorkspace>(*fWorkspace,"MaxLikelihoodFitResult");
	auto *workspace_bg_nobias	= GetFromTFile<RooWorkspace>(*fWorkspace_bg_nobias,"MaxLikelihoodFitResult");

	//Get data, pdfs, vars etc
	auto *x 		= GetFromRooWorkspace<RooRealVar>(*workspace,"mbb");
	auto *r 		= GetFromRooWorkspace<RooRealVar>(*workspace,"r");
	auto *data_obs 	= GetFromRooWorkspace<RooDataHist>(*workspace,"data_obs");
	//pdfs
	auto *pdf_sgn			= GetFromRooWorkspace<RooAbsPdf>(*workspace, "shapeSig_bbH_Mbb_bbHTo4b");
	auto *pdf_bkg			= GetFromRooWorkspace<RooAbsPdf>(*workspace, "shapeBkg_QCD_Mbb_bbHTo4b");
	auto *pdf_bkg_nobias	= GetFromRooWorkspace<RooAbsPdf>(*workspace_bg_nobias, "shapeBkg_QCD_Mbb_bbHTo4b");
	//pdfs normalisations
	auto *pdf_sgn_norm		= GetFromRooWorkspace<RooFormulaVar>(*workspace, "shapeSig_bbH_Mbb_bbHTo4b__norm");
	//finale norms
	auto *n_sgn_fit			= GetFromRooWorkspace<RooFormulaVar>(*workspace,"n_exp_final_binbbHTo4b_proc_bbH_Mbb");
	auto *n_bkg_fit			= GetFromRooWorkspace<RooFormulaVar>(*workspace,"n_exp_final_binbbHTo4b_proc_QCD_Mbb");
	auto *bkg_sys_fit_nobias= static_cast<RooRealVar*>(roofitresult_nobias->floatParsFinal().find("CMS_bkgd_qcd_13TeV"));
	auto *n_bkg_fit_nobias	= static_cast<RooRealVar*>(roofitresult_nobias->constPars().find("shapeBkg_QCD_Mbb_bbHTo4b__norm"));

	//Derive numbers
	auto *n_sgn 		= new RooRealVar("n_sgn","",n_sgn_fit->getVal());
	auto *n_bkg 		= new RooRealVar("n_bkg","",n_bkg_fit->getVal());
	auto *n_bkg_nobias 	= new RooRealVar("n_bkg_nobias","",n_bkg_fit_nobias->getVal() * bkg_sys_fit_nobias->getVal());

	auto *pdf_comb_ext = new RooAddPdf("pdf_comb_ext","",RooArgList(*pdf_sgn,*pdf_bkg),RooArgList(*n_sgn,*n_bkg));
	auto *h_rebinned = data_obs->createHistogram("h_data_obs_rebinned", *x, RooFit::Binning( nbins , x->getMin(), x->getMax()) );
	auto *data_rebinned = new RooDataHist("data_obs_rebinned","", RooArgList(*x), h_rebinned, 1.0);

	//TCanvas
	auto *c1 = new TCanvas("c1", "canvas", 600, 600);
	//TPad
	c1->cd();
    auto *pad1 = new TPad("pad1", "pad1", 0.02, 0.30, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.13);
    pad1->Draw();
    pad1->cd();

	//Frame
    auto *frame1 = x->frame(Name("frame1"));
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle(( ("Events / " + to_string_with_precision(h_rebinned->GetBinWidth(1),1) + " GeV").c_str() ));
    frame1->GetYaxis()->SetTitleSize(24);
    frame1->GetYaxis()->SetTitleFont(43);
    frame1->GetYaxis()->SetTitleOffset(1.18);
    frame1->GetYaxis()->SetLabelFont(43);
    frame1->GetYaxis()->SetLabelSize(20);
    frame1->GetXaxis()->SetTitleSize(0);
    frame1->GetXaxis()->SetTitleFont(43);
    frame1->GetXaxis()->SetLabelFont(43);
    frame1->GetXaxis()->SetLabelSize(0);

    auto *hdummy = new TH1F("hdummy", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());

    auto *hbkg = new TH1F("hbkg", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    hbkg->SetFillColor(kAzure+1);
    hbkg->SetLineColor(kAzure+1);
    hbkg->SetLineWidth(0);
    hbkg->SetFillStyle(3001);

    auto *h1sigmaU = new TH1F("h1sigmaU", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h1sigmaU->SetFillColor(kGreen+1);
    h1sigmaU->SetFillStyle(1001);

    auto *h1sigmaD = new TH1F("h1sigmaD", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h1sigmaD->SetFillColor(kGreen+1);
    h1sigmaD->SetFillStyle(1001);

    auto *h2sigmaU = new TH1F("h2sigmaU", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h2sigmaU->SetFillColor(kOrange);
    h2sigmaU->SetFillStyle(1001);

    auto *h2sigmaD = new TH1F("h2sigmaD", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h2sigmaD->SetFillColor(kOrange);
    h2sigmaD->SetFillStyle(1001);

    auto *leg = new TLegend(0.55,0.48,0.88,0.88,"","brNDC");
    leg->SetHeader( ("m_{#phi} = " + mass_point + " GeV").c_str() );
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->SetFillColor(10);

    cout<<"START plotting"<<endl;
    cout<<"Plot data_obs..."<<endl;
    data_rebinned->plotOn(frame1, Name("data_obs"));
    cout<<"Plot bkg +/- 2 sigma..."<<endl;

    pdf_bkg_nobias->plotOn(frame1,
    		VisualizeError(*roofitresult_nobias, 2, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kOrange), Name( (string(pdf_bkg_nobias->GetName()) +"_2sigma").c_str()), Normalization(n_bkg_nobias->getVal() , RooAbsReal::NumEvent), MoveToBack() );

    cout<<"Plot bkg +/- 1 sigma..."<<endl;
    pdf_bkg_nobias->plotOn(frame1,
    		VisualizeError(*roofitresult_nobias, 1, true),
			LineColor(kBlue),
			LineStyle(kSolid),
			FillColor(kGreen+1), Name( (string(pdf_bkg_nobias->GetName())+"_1sigma").c_str()) , Normalization(n_bkg_nobias->getVal() , RooAbsReal::NumEvent) );

    cout<<"Plot bkg..."<<endl;
    pdf_bkg_nobias->plotOn(frame1,
    		LineWidth(3),
			LineColor(kBlue),
			LineStyle(kSolid),
			Name(pdf_bkg_nobias->GetName()),
			Normalization(n_bkg_nobias->getVal(), RooAbsReal::NumEvent) );

    cout<<"Plot sgn+bkg..."<<endl;
    pdf_comb_ext->plotOn(frame1,
    		LineWidth(3),
			LineColor(kRed),
			LineStyle(kDashed),
			Name(pdf_comb_ext->GetName()) );

    cout<<"Plot sgn..."<<endl;
    pdf_sgn->plotOn(frame1,
    		LineWidth(3),
			LineColor(46),
			LineStyle(kDotted),
            Name(pdf_sgn->GetName()),
			Normalization(pdf_sgn_norm->getVal()*xsec_vis , RooAbsReal::NumEvent) );

    cout<<"Plot data_obs..."<<endl;
    data_rebinned->plotOn(frame1, Name("data_obs"));

    cout<<"END plotting"<<endl;
//    auto chi2_s = frame1->chiSquare(pdf_comb_ext->GetName(), "data_obs", 2+1 if bkg_pdf=='dijet' else 3+1 )
//    chi2_b = frame1.chiSquare(pdf_bkg_b.GetName(), "data_obs", 2 if bkg_pdf=='dijet' else 3 )
//    ndof = (data_rebinned.numEntries()-(2 if bkg_pdf=='dijet' else 3))
//    print ROOT.TMath.Prob(chi2_b*ndof, ndof)
    leg->AddEntry(frame1->findObject("data_obs"), "Data", "PE");
    leg->AddEntry(frame1->getCurve(pdf_bkg_nobias->GetName()),  ("Bkg."), "L");
    leg->AddEntry(h1sigmaU, "#pm1 std. deviation", "F");
    leg->AddEntry(h2sigmaU, "#pm2 std. deviation", "F");
    leg->AddEntry(frame1->getCurve(pdf_comb_ext->GetName()), ("Bkg. + signal"), "L");
    leg->AddEntry(frame1->getCurve(pdf_sgn->GetName()), ("Signal, #sigma= " + to_string_with_precision(xsec_vis,2) + " pb" ).c_str(), "L");
    frame1->SetMaximum( h_rebinned->GetMaximum()*2.0 );
    frame1->SetMinimum( h_rebinned->GetMinimum()*0.8 );
    pad1->SetLogy();
    frame1->Draw();
    leg->Draw();

    c1->Modified();

    c1->cd();
    auto *pad2 = new TPad("pad2", "pad2", 0.02, 0.02, 1, 0.28);
    pad2->SetTopMargin(0.02);
    pad2->SetLeftMargin(0.13);
    pad2->SetBottomMargin(0.35);
    pad2->SetGridx();
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    auto *frame2 = x->frame(Name("frame2"));
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
        auto bin_i = h1sigmaU->FindBin( x_i );
        if(bin_i <= 0 || bin_i > h1sigmaU->GetNbinsX()) continue;
        auto bin_width_i = h_rebinned->GetBinWidth( bin_i );
        auto n_data_i = h_rebinned->GetBinContent( bin_i );
        h1sigmaU->SetBinContent( bin_i, (-n_i+n_1up_i));
        h1sigmaD->SetBinContent( bin_i, (-n_i+n_1down_i));
        h2sigmaU->SetBinContent( bin_i, (-n_i+n_2up_i));
        h2sigmaD->SetBinContent( bin_i, (-n_i+n_2down_i));
        if(x_i < 600) cout<<"bin"<<i<<" x = "<<x_i<<" y = "<<n_i<<" +-1 = "<<n_1up_i<<" +-2 = "<<n_2up_i<<endl;
    }

    //redefine pulls here:
    auto bkg_int = pdf_bkg_nobias->createIntegral(RooArgSet(*x))->getVal();
    for(int bin_i = 1; bin_i < h_rebinned->GetNbinsX()+1; ++bin_i){
        auto n_data_i = h_rebinned->GetBinContent( bin_i );
        auto bin_width_i = h_rebinned->GetBinWidth( bin_i );
        x->setRange( ("Bin" + to_string(bin_i)).c_str(), h_rebinned->GetBinLowEdge(bin_i), h_rebinned->GetBinLowEdge(bin_i)+bin_width_i);
        auto bkg_int_i = pdf_bkg_nobias->createIntegral(RooArgSet(*x), ("Bin" + to_string(bin_i)).c_str())->getVal() / bkg_int * n_bkg_nobias->getVal();
        auto pull_i = (n_data_i-bkg_int_i)/sqrt(bkg_int_i);
        cout<<"Data: "<<n_data_i<<" ---- Bkg: "<<bkg_int_i<<endl;
        hbkg->SetBinContent( bin_i, pull_i);
        hbkg->SetBinError( bin_i, 0.);
        hdummy->SetBinContent( bin_i, 0.);
        hdummy->SetBinError( bin_i, 0.);
        h1sigmaU->SetBinContent( bin_i, h1sigmaU->GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h1sigmaD->SetBinContent( bin_i, h1sigmaD->GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h2sigmaU->SetBinContent( bin_i, h2sigmaU->GetBinContent( bin_i )/sqrt( bkg_int_i ) );
        h2sigmaD->SetBinContent( bin_i, h2sigmaD->GetBinContent( bin_i )/sqrt( bkg_int_i ) );
    }


    frame2->addTH1(h2sigmaU, "HIST");
    frame2->addTH1(h2sigmaD, "HIST");
    frame2->addTH1(h1sigmaU, "HIST");
    frame2->addTH1(h1sigmaD, "HIST");
	frame2->addTH1(hbkg, "HIST");
	frame2->addTH1(hdummy, "HIST");
	frame2->Draw();
	frame2->GetYaxis()->SetRangeUser(-4,4);
	c1->cd();
	c1->Draw();

	c1->Print( (cmsswBase + output + "Post_fit_m_" + mass_point + output_append + ".pdf").c_str());

//	TH1F *bg 		= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total_background");	// post-fit BG
//	TH1F *signal 	= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total_signal"); 	// post-fit Signal
//	TH1F *total 	= GetFromTFile<TH1F>(*fIn,"shapes_fit_s/bbHTo4b/total");			// post-fit Total
//	TH1F *data 		= (TH1F*) getDataHist(*GetFromTFile<RooWorkspace>(*fBg,"workspace"),*bg);
//
//	//Rebin histograms:
//	total->Scale(total->GetBinWidth(10));
//	signal->Scale(signal->GetBinWidth(10));
//	bg->Scale(9);
////	bg->Scale(bg->GetBinWidth(10));
//	total->Rebin(100);
//	signal->Rebin(100);
//	data->Rebin(100);
////	bg->Rebin(100);
//
//	draw(data,total,signal,bg);

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
