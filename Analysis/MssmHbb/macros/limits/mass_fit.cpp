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
#include <TSystem.h>

#include <RooFit.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooPlot.h>
#include <RooFitResult.h>
#include <RooAbsPdf.h>
#include <RooFormulaVar.h>
#include <RooAddPdf.h>
#include <RooDataHist.h>

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"

using namespace std;
using namespace RooFit;
using namespace analysis;

HbbStyle  style;

TH1* getDataHist(RooWorkspace& w, TH1& binning);
void draw(TH1F* data, TH1F *total, TH1F *signal, TH1F *bg);
int getBackgorundFitNBins(const int& mass_point);
int getMultiplicativeCrossection(const int& mass_point);
TH1 * getHackedDataObsTFile(const int& mass_point);
int getBackgroundFitUpEdge(const int& mass_point);
int getBackgroundFitLowEdge(const int& mass_point);
double findAverageErrorWithinTheBin(const vector<float>& v);

int main(int argc, char **argv) {
	/*
	 * Macro to plot post-fit M_{12} fit
	 */
	gSystem->Load("libHiggsAnalysisCombinedLimit");
	style.setTDRstyle(PRELIMINARY);
	gStyle->SetErrorX(0);

	string mass_point = argv[1];
	int mass = stoi(mass_point);
	int nbins = getBackgorundFitNBins(mass);
	double xsec_vis = getMultiplicativeCrossection(mass);
	if(mass < 500) gStyle->SetPadRightMargin(gStyle->GetPadRightMargin()*1.22);
	string output = "/src/Analysis/MssmHbb/macros/pictures/";
	string output_append = "test";
	string campaign = "201708/23/unblinded/independent/mll/";
	int n_fit_bins = 5000;
	output_append += "_nFitBins_" + to_string(n_fit_bins);

	//Files with mll fit results
//	TFile *fMLL 		= new TFile((mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "SpBg/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	//S+Bg mll
	TFile *fMLL_bg_nobias 	= new TFile((mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "bg_only/mll_M-" + mass_point + "/mlfit.root").c_str(),"READ");	// Bg only, no-bias
	//Files with post-fit workspaces
	TFile *fWorkspace			= new TFile((mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "SpBg/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");
	TFile *fWorkspace_bg_nobias 	= new TFile((mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/" + campaign + "bg_only/mll_M-" + mass_point + "/MaxLikelihoodFitResult.root").c_str(),"READ");

	//Get RooFitResults and workspaces
//	auto *roofitresult 			= GetFromTFile<RooFitResult>(*fMLL,"fit_s");
	auto *roofitresult_nobias 	= GetFromTFile<RooFitResult>(*fMLL_bg_nobias,"fit_b");
	auto *workspace				= GetFromTFile<RooWorkspace>(*fWorkspace,"MaxLikelihoodFitResult");
	auto *workspace_bg_nobias	= GetFromTFile<RooWorkspace>(*fWorkspace_bg_nobias,"MaxLikelihoodFitResult");

	//Get data, pdfs, vars etc
	auto *x 		= GetFromRooWorkspace<RooRealVar>(*workspace,"mbb");
	x->setBins(n_fit_bins);
//	auto *r 		= GetFromRooWorkspace<RooRealVar>(*workspace,"r");
	auto hacked_data_obs_p = getHackedDataObsTFile(stoi(mass_point));
//	auto *data_obs 	= GetFromRooWorkspace<RooDataHist>(*workspace,"data_obs");
//	auto *data_to_check = static_cast<TH1D*>(data_obs->createHistogram("mbb",5000));
//	cout<<"CHECK of the DATA content: "<<endl;
//	double sum = 0;
//	for(int i = 1; i<data_to_check->GetNbinsX();++i){
//		sum += data_to_check->GetBinContent(i);
//		std::cout<<"Bin "<<i<<" ["<<data_to_check->GetBinLowEdge(i)<<","<<data_to_check->GetBinLowEdge(i) + data_to_check->GetBinWidth(i)<<"] : "<<
//				"data = "<<data_to_check->GetBinContent(i)<<" sum: "<<sum<<endl;
//	}
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
	auto *n_sgn 		= new RooRealVar("n_sgn","",n_sgn_fit->getVal());
	auto *n_bkg 		= new RooRealVar("n_bkg","",n_bkg_fit->getVal());
	auto *n_bkg_nobias 	= new RooRealVar("n_bkg_nobias","",n_bkg_fit_nobias->getVal() * bkg_sys_fit_nobias->getVal());

	auto *pdf_comb_ext = new RooAddPdf("pdf_comb_ext","",RooArgList(*pdf_sgn,*pdf_bkg),RooArgList(*n_sgn,*n_bkg));
	auto *h_rebinned = static_cast<TH1D*>(hacked_data_obs_p->Clone("h_data_obs_rebinned"));
//	auto *h_rebinned = GetRooObjectFromTFile<RooDataSet>(*hacked_data_obs_p,"data_obs")->createHistogram("h_data_obs_rebinned", *x_hacked );
//	auto *h_rebinned = data_rebinned->createHistogram("h_data_obs_rebinned",*x,nbins);
//	auto *h_rebinned = data_obs->createHistogram("h_data_obs_rebinned", *x, RooFit::Binning( nbins , x->getMin(), x->getMax()) );
	auto *data_rebinned = new RooDataHist("data_obs_rebinned","", RooArgList(*x), h_rebinned, 1.0);
//	auto *data_rebinned = GetRooObjectFromTFile<RooDataSet>(hacked_data_obs_p,"data_obs")->;
cout<<nbins<<endl;
	cout<<"CHECK of the DATA REBINNED content: "<<endl;
	double sum = 0;
	for(int i = 1; i<h_rebinned->GetNbinsX();++i){
		sum += h_rebinned->GetBinContent(i);
		std::cout<<"Bin "<<i<<" ["<<h_rebinned->GetBinLowEdge(i)<<","<<h_rebinned->GetBinLowEdge(i) + h_rebinned->GetBinWidth(i)<<"] : "<<
				"data = "<<h_rebinned->GetBinContent(i)<<" sum: "<<sum<<endl;
	}

	//TCanvas
	auto c1 = new TCanvas();
	double canva_pixel_h = c1->YtoPixel(c1->GetY1());
	//TPad
	c1->cd();
	auto pad1 = HbbStyle::getRatioTopPad(0.7);
    pad1->Draw();
    pad1->cd();
	double topPad_pixel_h = pad1->YtoPixel(pad1->GetY1());

	//Frame
    auto *frame1 = static_cast<RooPlot*>(x->frame(Name("frame1")));
    frame1->SetTitle("");
    frame1->GetYaxis()->SetTitle(( ("Events / ( " + std::to_string(int(h_rebinned->GetBinWidth(1))) + " GeV )").c_str() ));
    HbbStyle::setRatioTopFrame(frame1->GetXaxis(), frame1->GetYaxis(), canva_pixel_h);

    auto *hdummy = new TH1F("hdummy", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());

    auto *hbkg = new TH1F("hbkg", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    hbkg->SetMarkerStyle(gStyle->GetMarkerStyle());
    hbkg->SetMarkerSize(gStyle->GetMarkerSize());

//    auto *h1sigmaU = new TH1F("h1sigmaU", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
//    h1sigmaU->SetFillColor(kGreen+1);
//    h1sigmaU->SetFillStyle(1001);
//
//    auto *h1sigmaD = new TH1F("h1sigmaD", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
//    h1sigmaD->SetFillColor(kGreen+1);
//    h1sigmaD->SetFillStyle(1001);
//
//    auto *h2sigmaU = new TH1F("h2sigmaU", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
//    h2sigmaU->SetFillColor(kOrange);
//    h2sigmaU->SetFillStyle(1001);
//
//    auto *h2sigmaD = new TH1F("h2sigmaD", "", h_rebinned->GetNbinsX(), h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
//    h2sigmaD->SetFillColor(kOrange);
//    h2sigmaD->SetFillStyle(1001);

    TH1F* h1sigmaU = new TH1F("h1sigmaU", "", n_fit_bins, h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h1sigmaU->SetFillColor(kGreen+1);
    h1sigmaU->SetLineColor(kGreen+1);
    h1sigmaU->SetFillStyle(1001);

    TH1F *h1sigmaD = new TH1F("h1sigmaD", "", n_fit_bins, h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h1sigmaD->SetFillColor(kGreen+1);
    h1sigmaD->SetLineColor(kGreen+1);
    h1sigmaD->SetFillStyle(1001);

    TH1F *h2sigmaU = new TH1F("h2sigmaU", "", n_fit_bins, h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h2sigmaU->SetFillColor(kOrange);
    h2sigmaU->SetLineColor(kOrange);
    h2sigmaU->SetFillStyle(1001);

    TH1F *h2sigmaD = new TH1F("h2sigmaD", "", n_fit_bins, h_rebinned->GetXaxis()->GetXmin(), h_rebinned->GetXaxis()->GetXmax());
    h2sigmaD->SetFillColor(kOrange);
    h2sigmaD->SetLineColor(kOrange);
    h2sigmaD->SetFillStyle(1001);

    auto *leg = new TLegend(0.63,0.5,0.85,0.88,"","brNDC");
    leg->SetHeader( ("m_{A/H} = " + mass_point + " GeV").c_str() );
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
    leg->AddEntry(frame1->getCurve(pdf_bkg_nobias->GetName()),  ("Background"), "L");
    leg->AddEntry(h1sigmaU, "#pm1 std. deviation", "F");
    leg->AddEntry(h2sigmaU, "#pm2 std. deviation", "F");
    leg->AddEntry(frame1->getCurve(pdf_comb_ext->GetName()), ("Bkg. + signal"), "L");
    leg->AddEntry(frame1->getCurve(pdf_sgn->GetName()), (std::to_string(int(xsec_vis)) + " x bbA/H").c_str(), "L");
//    leg->AddEntry(frame1->getCurve(pdf_sgn->GetName()), ("m_{A/H} = " + mass_point + "GeV, #sigma = " + std::to_string(int(xsec_vis)) + " pb" ).c_str(), "L");
    frame1->SetMaximum( h_rebinned->GetMaximum()*2.0 );
    frame1->SetMinimum( h_rebinned->GetMinimum()*0.8 );
    pad1->SetLogy();
    frame1->Draw();
    leg->Draw();

    c1->Modified();

    c1->cd();
    auto pad2 = HbbStyle::getRatioBottomPad(0.3);
    pad2->Draw();
    pad2->cd();

    auto *frame2 = static_cast<RooPlot*>(x->frame(Name("frame2")));
    frame2->SetTitle("");
    frame2->GetYaxis()->SetTitle("Data-Bkg.");//SetTitle("#frac{Data-Bkg.}{#sigma_{Bkg.}}");//SetTitle("#frac{Data-Bkg.}{#sqrt{Bkg.}}");
    frame2->GetYaxis()->CenterTitle();
    HbbStyle::setRatioBottomFrame(frame2->GetXaxis(), frame2->GetYaxis(), canva_pixel_h, topPad_pixel_h);

    //Retrieve information from the RooPlot
    auto *sigma2 	= frame1->getCurve( (string(pdf_bkg_nobias->GetName())+"_2sigma").c_str());
    auto *sigma1 	= frame1->getCurve( (string(pdf_bkg_nobias->GetName())+"_1sigma").c_str());
    auto *nominal 	= frame1->getCurve(pdf_bkg_nobias->GetName());

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

    for( int bin_i = 1; bin_i <= h1sigmaU->GetNbinsX(); ++bin_i ){
    		double x = h1sigmaU->GetBinCenter(bin_i);
    	//        std::cout<<"WTFFF: x = "<<x<<" y = "<<nominal->Eval(x)<<" up1 = "<<up1Bound.Eval(x)<<" lo1 = "<<lo1Bound.Eval(x)<<" up2 = "<<up2Bound.Eval(x)<<" lo2 = "<<lo2Bound.Eval(x)<<std::endl;

    		auto n_i = nominal->Eval(x);
    		auto n_1up_i = up1Bound.Eval(x);
    		auto n_1down_i = lo1Bound.Eval(x);
    		auto n_2up_i = up2Bound.Eval(x);
    		auto n_2down_i = lo2Bound.Eval(x);

    		//Fill band histograms
    		err1U = (-n_i+n_1up_i)  ;
    		err1D = (-n_i+n_1down_i);
    		err2U = (-n_i+n_2up_i)  ;
    		err2D = (-n_i+n_2down_i);

    		double divide_by = n_i;
    		if(divide_by_sqrt_bg){
    			err1U /= sqrt(divide_by);
    			err1D /= sqrt(divide_by);
    			err2U /= sqrt(divide_by);
    			err2D /= sqrt(divide_by);
    		}

    		h1sigmaU->SetBinContent( bin_i, err1U);
    		h1sigmaD->SetBinContent( bin_i, err1D);
    		h2sigmaU->SetBinContent( bin_i, err2U);
    		h2sigmaD->SetBinContent( bin_i, err2D);
    }

    for(int bin_i = 1; bin_i < h_rebinned->GetNbinsX()+1; ++bin_i){
            auto n_data_i = h_rebinned->GetBinContent( bin_i );
            auto e_data_i = h_rebinned->GetBinError( bin_i );

            //Evaluate value of the Bg-function
            double bkg_int_i = nominal->Eval(h_rebinned->GetBinCenter(bin_i));
            double divide_by = sqrt(bkg_int_i);//ssqrt(n_data_i);//
            auto pull_i = (n_data_i-bkg_int_i)/divide_by;
            e_data_i /= divide_by;

            hbkg->SetBinContent( bin_i, pull_i);
            hbkg->SetBinError( bin_i, e_data_i);
            hdummy->SetBinContent( bin_i, 0.);
            hdummy->SetBinError( bin_i, 0.);
        }

    //redefine pulls here:
//    auto bkg_int = pdf_bkg_nobias->createIntegral(RooArgSet(*x))->getVal();
//    for(int bin_i = 1; bin_i < h_rebinned->GetNbinsX()+1; ++bin_i){
//        auto n_data_i = h_rebinned->GetBinContent( bin_i );
//        auto e_data_i = h_rebinned->GetBinError( bin_i );
//        auto bin_width_i = h_rebinned->GetBinWidth( bin_i );
//        x->setRange( ("Bin" + to_string(bin_i)).c_str(), h_rebinned->GetBinLowEdge(bin_i), h_rebinned->GetBinLowEdge(bin_i)+bin_width_i);
//        auto bkg_int_i = pdf_bkg_nobias->createIntegral(RooArgSet(*x), ("Bin" + to_string(bin_i)).c_str())->getVal() / bkg_int * n_bkg_nobias->getVal();
//        double divide_by = sqrt(bkg_int_i);
//        auto pull_i = (n_data_i-bkg_int_i);
//
//        cout<<"Data: "<<n_data_i<<" ---- Bkg: "<<bkg_int_i<<endl;
//        if(divide_by_sqrt_bg){
//        		pull_i /= divide_by;
//        		e_data_i /= divide_by;
//        }
//        hbkg->SetBinContent( bin_i, pull_i);
//		hbkg->SetBinError( bin_i, e_data_i);
//        hdummy->SetBinContent( bin_i, 0.);
//        hdummy->SetBinError( bin_i, 0.);
//    }

    frame2->addTH1(h2sigmaU, "HIST");
    frame2->addTH1(h2sigmaD, "HIST");
    frame2->addTH1(h1sigmaU, "HIST");
    frame2->addTH1(h1sigmaD, "HIST");
	frame2->addTH1(hbkg, "E1");
	frame2->addTH1(hdummy, "HIST");
	frame2->Draw();
    frame2->SetMinimum(-3.9999);
    frame2->SetMaximum(+3.9999);
//    frame2->SetMinimum(-180);
//    frame2->SetMaximum(+180);
//    if(mass > 900) {
//    		frame2->SetMinimum(-100);
//    		frame2->SetMaximum(+100);
//    }
	c1->cd();
	HbbStyle::drawStandardTitle("out");
	c1->Draw();

	c1->Print( (mssmhbb::cmsswBase + output + "Post_fit_m_" + mass_point + output_append + ".pdf").c_str());
	c1->SaveAs("bla.root");

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

    can->Print( (mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/ParametricLimits/post_fit_500GeV.pdf").c_str() );
}


TH1* getDataHist(RooWorkspace& w, TH1& binning){
	auto &dataHist 		= *analysis::GetFromRooWorkspace<RooDataHist>(w,"data_obs");
	auto &mbb			= *analysis::GetFromRooWorkspace<RooRealVar>(w,"mbb");
	auto *data 			=  static_cast<TH1F*>(dataHist.createHistogram("data_obs",mbb,RooFit::Binning(binning.GetNbinsX(), binning.GetXaxis()->GetXmin(), binning.GetXaxis()->GetXmax())));
	return data;
}

int getBackgorundFitNBins(const int& mass_point){
	/*
	 * Get the sub-range number from the signal mass point
	 */
	if(mass_point < 500) return 45;
	else if (mass_point < 1100) return 42;
	else return 48;
}

int getMultiplicativeCrossection(const int& mass_point){
	/*
	 * Get the cross-section
	 */
	if(mass_point < 500) return 50;
	else if (mass_point < 1100) return 4;
	else return 2;
}

int getBackgroundFitLowEdge(const int& mass_point){
	if(mass_point < 500) return 200;
	else if (mass_point < 1100) return 350;
	else return 500;
}

int getBackgroundFitUpEdge(const int& mass_point){
	if(mass_point < 500) return 650;
	else if (mass_point < 1100) return 1190;
	else return 1700;
}

TH1 * getHackedDataObsTFile(const int& mass_point){
	/*
	 * Get the data obs with 45 bins, as it's in the PAS
	 */
	TFile *f = nullptr;
	f = new TFile( (mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/data/2016DataRereco_05_01_2017/TripleBTagSelectionBtoH2016_13TeV.root").c_str());
//	if(mass_point < 500) {
//		f = new TFile( (mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/finale_files/SR/extnovoefffixprod_200to650_10GeV_SR/workspace/FitContainer_workspace.root").c_str() );
//	}
//	else if (mass_point < 1100){
//		f = new TFile( (mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/finale_files/SR/novosibirsk_350to1190/workspace/FitContainer_workspace.root").c_str() );
//	}
//	else {
//		f = new TFile( (mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/finale_files/SR/novosibirsk_500to1700/workspace/FitContainer_workspace.root").c_str() );
//	}

	TTree *t;// = GetFromTFile<TTree>(*f,"MssmHbb_13TeV");
	f->GetObject("MssmHbb_13TeV",t);
	TH1D *h = new TH1D("histo","histo",getBackgorundFitNBins(mass_point),getBackgroundFitLowEdge(mass_point),getBackgroundFitUpEdge(mass_point));
	t->Draw("mbb>>histo");
	return h;
}

double findAverageErrorWithinTheBin(const vector<float>& v){
	return 1.0 * std::accumulate(v.begin(),v.end(),0.0)/v.size();
}

