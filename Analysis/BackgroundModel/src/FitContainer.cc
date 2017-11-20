#include <iostream>
#include <sstream>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooFormulaVar.h"
#include "RooEffProd.h"
#include "RooGenericPdf.h"
#include "RooExponential.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooBukinPdf.h"
#include "RooProdPdf.h"
#include "Analysis/BackgroundModel/interface/RooDoubleCB.h"
#include "Analysis/BackgroundModel/interface/RooExpGausExp.h"
#include "Analysis/BackgroundModel/interface/RooGausExp.h"
#include "Analysis/BackgroundModel/interface/RooExpBWExp.h"
#include "Analysis/BackgroundModel/interface/RooPhaseSpace.h"
#include "Analysis/BackgroundModel/interface/RooPhaseSpacePol4.h"
#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/BackgroundModel/interface/Tools.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"

using namespace analysis::backgroundmodel;

//private constructor, for members initialisation
FitContainer::FitContainer(const std::string& outputDir) :
		initialized_(false),
		written_(false),
		splitrange_(true),
		outputDir_( outputDir.back() == '/' ?  outputDir : outputDir + "/"),
		plotDir_(getOutputPath_("plots")),
		workspaceDir_(getOutputPath_("workspace")),
		fullRangeId_("full_range"),
		fitRangeId_("fit_range"),
		fitRangeLowId_("fit_range_low"),
		fitRangeHighId_("fit_range_high"),
		fitSplRangeId_("fit_range_low,fit_range_high"),
		blind_lowEdge_(650.),
		blind_highEdge_(650.),
		verbosity_(1),
		workspace_(RooWorkspace("workspace")),
		outRootFileName_(getOutputPath_("workspace")+"workspace.root"),
		mbb_("mbb"),
		weight_("weight"),
		data_("data_container"),
		signal_("signal_container"),
		bkg_("background_container"),
		bkgOnlyFit_("fit_b", "fit_b"),
		chi2_lowEdge_(-10000),
		chi2_highEdge_(10000),
		chi2BkgOnly_(-10000.0),
		normChi2BkgOnly_(-10000.0),
		ndfBkgOnly_(-10000),
		nbins_(73), //73
		nBinsToPlot_(48),
		lumi_(35.7), //2.69,12.9,36.62
		obs_(26760.) //SR1-259399, SR2-105053, SR3-26760
{}

//FitContainer::FitContainer(const FitContainer& cont){
//	initialized_ 	= cont.initialized_;
//	written_		= cont.written_;
//	splitrange_		= cont.splitrange_;
//	outputDir_		= cont.outputDir_;
//	plotDir_		= cont.plotDir_;
//	workspaceDir_	= cont.workspaceDir_;
//	fullRangeId_	= cont.fullRangeId_;
//	fitRangeId_		= cont.fitRangeId_;
//	fitRangeLowId_	= cont.fitRangeLowId_;
//	fitRangeHighId_ = cont.fitRangeHighId_;
//	fitSplRangeId_ 	= cont.fitSplRangeId_;
//	fitRangeMin_	= cont.fitRangeMin_;
//	fitRangeMax_	= cont.fitRangeMax_;
//	blind_lowEdge_	= cont.blind_lowEdge_;
//	blind_highEdge_	= cont.blind_highEdge_;
//	verbosity_		= cont.verbosity_;
//	workspace_		= cont.workspace_;
//	outRootFileName_= cont.outRootFileName_;
//	mbb_			= cont.mbb_;
//	weight_			= cont.weight_;
//	data_			= cont.data_;
//	signal_			= cont.signal_;
//	bkg_			= cont.bkg_;
////	Workaround to copy TTree
////	Original idea by Gregor Mittag
////	TODO: Implement it. DOesn't work out of the box
////	bkgOnlyFit_(((TTree&) cont.bkgOnlyFit_).CloneTree(0));
////	bkgOnlyFit_.SetDirectory(0);
////	bkgOnlyFit_.CopyEntries(cont.bkgOnlyFit_);
////	bkgOnlyFit_		= cont.bkgOnlyFit_;
//	chi2BkgOnly_	= cont.chi2BkgOnly_;
//	normChi2BkgOnly_= cont.normChi2BkgOnly_;
//	ndfBkgOnly_		= cont.ndfBkgOnly_;
//	nbins_			= cont.nbins_;
//}


FitContainer::FitContainer(const TH1* data, const std::string& outputDir, const std::string & type) : FitContainer(outputDir)
{
	RooRealVar mbb(mbb_.c_str(), "M_{12}",
                 data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(), "GeV");
	fitRangeMin_ = mbb.getMin();
	fitRangeMax_ = mbb.getMax();
	workspace_.import(mbb);
	nbins_ = data->GetNbinsX();
	if(type == "background") {
		data_ = bkg_;
//		RooDataHist bkgContainer(bkg_.c_str(), bkg_.c_str(), mbb, data);
//		workspace_.import(bkgContainer);
		RooDataHist dataContainer(data_.c_str(), data_.c_str(), mbb, data);
		workspace_.import(dataContainer);
	}
	else if (type == "signal") {
		data_ = signal_;
		RooDataHist signalContainer(signal_.c_str(), signal_.c_str(), mbb, data);
		workspace_.import(signalContainer);
		RooDataHist dataContainer(signal_.c_str(), signal_.c_str(), mbb, data);
		workspace_.import(dataContainer);
	}
	else if (type == "data") {
		RooDataHist dataContainer(data_.c_str(), data_.c_str(), mbb, data);
		workspace_.import(dataContainer);
	}
	else {
		std::cerr<<"Wrong type were provided to FitContainer::FitContainer. Possible types are: signal / background / data"<<std::endl;
		exit(-1);
	}
}

FitContainer::FitContainer(const TH1* data, const TH1* signal, const TH1* bkg,
			   const std::string& outputDir) : FitContainer(outputDir)
 {
	double xmin, xmax;
	if(data) {
		xmin = data->GetXaxis()->GetXmin();
		xmax = data->GetXaxis()->GetXmax();
		nbins_ = data->GetNbinsX();
	}
	else if (signal){
		xmin = signal->GetXaxis()->GetXmin();
		xmax = signal->GetXaxis()->GetXmax();
		nbins_ = signal->GetNbinsX();
	}
	else if (bkg){
		xmin = bkg->GetXaxis()->GetXmin();
		xmax = bkg->GetXaxis()->GetXmax();
		nbins_ = bkg->GetNbinsX();
	}
	else {
		std::cerr<<"Empty HistContainer was provided to FitContainer::FitContainer"<<std::endl;
		exit(-1);
	}
	RooRealVar mbb(mbb_.c_str(), "M_{12}",
                 xmin, xmax, "GeV");
	fitRangeMin_ = mbb.getMin();
	fitRangeMax_ = mbb.getMax();
	workspace_.import(mbb);

	// Name and title of the dataset MUST be identical (see initialize() method).
	if (data && signal && bkg){
		RooDataHist dataContainer(data_.c_str(), data_.c_str(), mbb, data);
		workspace_.import(dataContainer);
		RooDataHist signalContainer(signal_.c_str(), signal_.c_str(), mbb, signal);
		workspace_.import(signalContainer);
		RooDataHist bkgContainer(bkg_.c_str(), bkg_.c_str(), mbb, bkg);
		workspace_.import(bkgContainer);
	}
	else if(data) 	{
		RooDataHist dataContainer(data_.c_str(), data_.c_str(), mbb, data);
		workspace_.import(dataContainer);
	}
	else if(signal) 	{
		data_ = signal_;
		RooDataHist signalContainer(signal_.c_str(), signal_.c_str(), mbb, signal);
		workspace_.import(signalContainer);
		RooDataHist dataContainer(signal_.c_str(), signal_.c_str(), mbb, signal);
		workspace_.import(dataContainer);
	}
	else if(bkg) 	{
		data_ = bkg_;
		RooDataHist bkgContainer(bkg_.c_str(), bkg_.c_str(), mbb, bkg);
		workspace_.import(bkgContainer);
		RooDataHist dataContainer(bkg_.c_str(), bkg_.c_str(), mbb, bkg);
		workspace_.import(dataContainer);
	}




}

FitContainer::FitContainer(const HistContainer& container,
			   const std::string& outputDir) : FitContainer(container.data().get(), container.bbH().get(), container.background().get(),
				       outputDir)   {}

FitContainer::FitContainer(TTree& data, const std::string& outputDir) : FitContainer(outputDir)
 {
	bkg_ = "";
	signal_ = "";
	RooRealVar mbb(mbb_.c_str(), "M_{12}",
                 0.0, data.GetMaximum(mbb_.c_str()), "GeV");
	fitRangeMin_ = mbb.getMin();
	fitRangeMax_ = mbb.getMax();
	nbins_	     = 100;
	RooRealVar weight(weight_.c_str(), "weight", 0.0, 1000.0);
	workspace_.import(mbb);
	workspace_.import(weight);

	// Name and title of the dataset MUST be identical (see initialize() method).
	RooDataSet dataContainer(data_.c_str(), data_.c_str(), RooArgSet(mbb, weight),
                           RooFit::Import(data),
                           RooFit::WeightVar(weight_.c_str()));
	workspace_.import(dataContainer);
}


FitContainer::FitContainer(const TreeContainer& container, const std::string& outputDir) 
			: FitContainer(*container.data(), outputDir) {
}


FitContainer::~FitContainer() {
  workspace_.Print("v");
  if(!written_) workspace_.writeToFile(outRootFileName_.c_str());
  TFile out(outRootFileName_.c_str(), "update");
  bkgOnlyFit_.SetDirectory(&out);
  bkgOnlyFit_.Write();
  out.Close();
}


void FitContainer::initialize() {
  // Get back the name of the imported datasets. For some reason these are
  // deleted in the workspace after the constructor ends.
  // To get this hack here working, name and title of the dataset MUST be set
  // identical (see constructor methods).
  for (const auto& d: workspace_.allData()) d->SetName(d->GetTitle());

  // clean up possible pre-existing output:
  gSystem->Exec((std::string("rm -f "+plotDir_+"*").c_str()));
  gSystem->Exec((std::string("rm -f "+workspaceDir_+"*").c_str()));

  // set range used for normalization of the pdf and a default fit range:
  auto& mbb = *GetFromRooWorkspace<RooRealVar>(workspace_, mbb_);
  mbb.setRange(fullRangeId_.c_str(), mbb.getMin(), mbb.getMax());
  mbb.setRange(fitRangeId_.c_str(), fitRangeMin_, fitRangeMax_);
  mbb.setRange("chi2_range",chi2_lowEdge_,chi2_highEdge_);

  // perform split range simultaneously fit to blinded data by CA
  mbb.setRange(fitRangeLowId_.c_str(), fitRangeMin_, blind_lowEdge_);  //always have to give input of --fit_min
  //mbb.setRange(fitRangeMedId_.c_str(), blind_lowEdge_, blind_highEdge_);
  mbb.setRange(fitRangeHighId_.c_str(), blind_highEdge_, fitRangeMax_);
  
  // set fit bins
  mbb.setBins(nbins_);
  // plot the input data:
  auto& data = *GetFromRooWorkspace<RooAbsData>(workspace_,data_);
  std::unique_ptr<RooPlot> frame(mbb.frame());
  data.plotOn(frame.get(),RooFit::DataError(RooAbsData::Auto));
  TCanvas canvas("canvas", "", 600, 600);
  canvas.cd();
  prepareCanvas_(canvas);
  prepareFrame_(*frame);
  frame->Draw();
  canvas.SaveAs((plotDir_+"input_data.pdf").c_str());
  canvas.SetLogy();
  if(frame->GetMinimum() == 0) frame->SetMinimum(0.01);
  frame->Draw();
  canvas.SaveAs((plotDir_+"input_data_log.pdf").c_str());
  // initialize background-only fit result tree:
  //TODO: why TTree and TTree branches should be members of the class? Cna this be avoid?
  bkgOnlyFit_.Branch("chi2", &chi2BkgOnly_, "chi2/F");
  bkgOnlyFit_.Branch("normChi2", &normChi2BkgOnly_, "normChi2/F");
  bkgOnlyFit_.Branch("ndf", &ndfBkgOnly_, "ndf/I");
  bkgOnlyFit_.Branch("covMatrix", covMatrix_, "covMatrix[400]/D");
  bkgOnlyFit_.Branch("eigenVector", eigenVector_, "eigenVector[400]/D");

  for(int i = 0; i < 400; i++)
  {
	  covMatrix_[i] = -100.;
	  eigenVector_[i] = -100.;
  }	
  initialized_ = true;
}

void FitContainer::setModel(const Type& type, const std::string& name) {
  const std::vector<ParamModifier> modifiers; // empty modifier list
  setModel(type, name, modifiers);
}


void FitContainer::setModel(const Type& type, const std::string& name,
                            const std::vector<ParamModifier>& modifiers) {
  if (!initialized_) initialize();

  ProbabilityDensityFunctions pdfs(workspace_,mbb_.c_str());
//  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  double peak_pos = getPeakStart_(type,500);
  pdfs.setPeakStart(peak_pos);
  pdfs.setPdf(name,toString(type));

  applyModifiers_(*(workspace_.pdf(toString(type).c_str())), modifiers);
}

std::unique_ptr<RooFitResult> FitContainer::FitSignal(const std::string & name, const bool& plot_params) {
	/*
	 * Method to fit signal templates
	 */
	if(!initialized_) initialize();
	auto &Pdf 		= *GetFromRooWorkspace<RooAbsPdf>(workspace_,toString(Type::signal));
	auto &data 	= *GetFromRooWorkspace<RooAbsData>(workspace_,signal_);
	auto &mbb		= *GetFromRooWorkspace<RooRealVar>(workspace_, mbb_);
	RooFit::SumW2Error(kTRUE);

	//class-helper
	RooFitQuality fit_quality;

	// Split Range for blinded data and perform simultaneous fit
	std::unique_ptr<RooFitResult> fitResult(Pdf.fitTo(data,
			RooFit::Save(),
			RooFit::PrintLevel(verbosity_),
			RooFit::SumW2Error(kTRUE),
			RooFit::InitialHesse(kTRUE)
	));

	//INformation about pre- post-fit parametrs
	fit_quality.PrintParametersInfo(*fitResult);
	//Create log file:
	makeLog_(*fitResult);

	  // Top frame
	  std::unique_ptr<RooPlot> frame(mbb.frame());
	  data.plotOn(frame.get(), RooFit::MarkerSize(0.8), RooFit::Name("data_points"));
	  Pdf.plotOn(frame.get(),
	               RooFit::LineColor(kRed),
	               RooFit::Name("signal_curve"),
	               RooFit::NormRange(fullRangeId_.c_str()),
	               RooFit::Range(fitRangeId_.c_str()),
	               RooFit::Normalization(data.sumEntries("1", fitRangeId_.c_str()),
	                                     RooAbsReal::NumEvent));
	  if(plot_params){
		  double par_xmin = 0.65, par_xmax = 0.9, par_ymax = 0.6;
	//	  par_xmin = 0.2; par_xmax = 0.45; par_ymax = 0.9;
		  Pdf.paramOn(frame.get(),RooFit::Layout(par_xmin,par_xmax,par_ymax));//0.98-pad1->GetRightMargin(),0.83-pad1->GetTopMargin()));
		  frame->getAttText()->SetTextSize(0.03);
	  }
	  //Compare integrals of the data and fit in a chi2 deffinition range
	  double data_integral = data.sumEntries( (std::string("mbb > " +std::to_string(chi2_lowEdge_) + " && mbb< " + std::to_string(chi2_highEdge_))).c_str() );
	  double fit_integral  = Pdf.createIntegral(mbb,RooFit::NormSet(mbb),RooFit::Range("chi2_range"))->getVal()*data.sumEntries();
	  std::cout<<"Data vs Fit integrals"<<std::endl;
	  std::cout<<"Range: ["<<chi2_lowEdge_<<";"<<chi2_highEdge_<<"]"<<std::endl;
	  std::cout<<"Data = "<<data_integral;
	  std::cout<<" Fit = "<<fit_integral<<std::endl;
	  std::cout<<"(Data-Fit)/Data*100% = "<<(data_integral-fit_integral)/data_integral * 100<<std::endl;

	    int nPars = fitResult->floatParsFinal().getSize();

	    // Get Covariance Matrix for diagonalisation
	    TMatrixDSymEigen CM = fitResult->covarianceMatrix();
	    TMatrixD covVec = CM.GetEigenVectors();
	    TMatrixD covMat = fitResult->covarianceMatrix();
	    covMat.GetMatrix2Array(covMatrix_);
	    covVec.GetMatrix2Array(eigenVector_);

	    ////////////////////////////////////////////////////////////////////////////////
	    //////////////////////// ChiSquare by Roofit ///////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////
	    ndfBkgOnly_ = getNonZeroBins_(data) - nPars;
	    normChi2BkgOnly_ = frame->chiSquare("signal_curve", "data_points", nPars); //chi^2 from RooFit (RooPlot::chiSquare())
	    std::cout<<"ROOFIT: Chi^2/Ndf = "<<normChi2BkgOnly_<<std::endl;
	    //chi2BkgOnly_ = normChi2BkgOnly_ * ndfBkgOnly_;
	    //bkgOnlyFit_.Fill();


	    ////////////////////////////////////////////////////////////////////////////////
	    ////////////////////////ChiSquare by Chayanit///////////////////////////////////
	    ////////////////////////////////////////////////////////////////////////////////
	    Chi2Ndf chi2ndf = fit_quality.chiSquare(*frame, "signal_curve", "data_points", nPars,blind_lowEdge_,blind_highEdge_,chi2_lowEdge_,chi2_highEdge_ );
	    ndfBkgOnly_ 		= chi2ndf.ndf;
	    chi2BkgOnly_ 		= chi2ndf.chi2;
	    normChi2BkgOnly_ 	= chi2BkgOnly_/ndfBkgOnly_;
	    std::cout<<"Rostyslav: Chi^2/Ndf = "<<chi2ndf.chi2<<"/"<<chi2ndf.ndf<<" = "<<chi2ndf.chi2/chi2ndf.ndf<<std::endl;
	  	bkgOnlyFit_.Fill();

	  	std::string chi2str(Form("%.1f/%d = %.1f", chi2BkgOnly_,
	  		 	   ndfBkgOnly_, normChi2BkgOnly_));
	    std::cout << "\nNormalized chi^2: " << chi2str << std::endl;

	    // TMath::Prob p-value by CA
	    double prob = TMath::Prob(chi2BkgOnly_,ndfBkgOnly_);
	    std::string probstr(Form("%.2f", prob));
	    std::cout << "\nTMath::Prob(#chi^{2}_{RooFit},ndf) : " << probstr << std::endl;

	    // Bottom frame
	    RooHist* hpull;
	    hpull = frame->pullHist();
	    hpull->SetMarkerSize(0.8);	//0.8 for lowM
	    std::unique_ptr<RooPlot> frame2(mbb.frame());
	    frame2->addPlotable(hpull,"P");
	    //frame2->addObject(frame->pullHist());

	    //TCanvas canvas("canvas", "", 600, 600);
	    TCanvas canvas;
	    canvas.SetCanvasSize(500,500);
	    canvas.cd();
	    //prepareCanvas_(canvas);
	    prepareFrame_(*frame);
	    prepareFrame_(*frame2);

	    TPad* pad1;
	    pad1 = new TPad("pad1","",0,0.1,1,1);
	    pad1->SetBottomMargin(0.2);
	    pad1->SetRightMargin(0.05);
	    pad1->SetLeftMargin(0.16);
	    pad1->Draw();
	    pad1->cd();
	    frame->GetXaxis()->SetTitleOffset(999); //Effectively turn off x axis title on main plot
	    frame->GetXaxis()->SetLabelOffset(999); //Effectively turn off x axis label on main plot
	    frame->GetYaxis()->SetTitleSize(0.038);
	    frame->GetYaxis()->SetTitleOffset(1.6);
	    frame->GetYaxis()->SetLabelSize(0.033);
	    //frame->GetYaxis()->SetRangeUser(frame->GetMinimum(), frame->GetMaximum()+200);
	    frame->Draw();
	    //Draw TLine shows chi2 calculation borders
	    if(chi2_lowEdge_ != -10000 && chi2_highEdge_ != 10000){
	    	TLine low_edge(chi2_lowEdge_,0,chi2_lowEdge_,frame->GetMaximum());
	    	low_edge.SetLineStyle(7);
	    	low_edge.Draw();
	    	TLine high_edge(chi2_highEdge_,0,chi2_highEdge_,frame->GetMaximum());
	    	high_edge.SetLineStyle(7);
	    	high_edge.Draw();
	    }
	
	    std::string lumistr(Form("%.1f", lumi_));

	    TLatex latex;
	    latex.SetTextFont(43);
	    latex.SetTextSize(17);
	    latex.SetTextAlign(11);
	    latex.DrawLatexNDC(pad1->GetLeftMargin(), 1.02-canvas.GetTopMargin(),
			       (std::string("CMS Preliminary #sqrt{s} = 13 TeV, L = ")+lumistr+std::string(" fb^{-1}")).c_str());
	    latex.SetTextSize(15);
	    latex.SetTextAlign(33);
	    latex.SetTextColor(kBlue+2);
	    latex.DrawLatexNDC(0.98-pad1->GetRightMargin(), 0.98-pad1->GetTopMargin(),
	                       (std::string("#chi^{2}_{RooFit}/ndf = ")+chi2str).c_str());
	    latex.SetTextColor(kGreen+2);
	    latex.DrawLatexNDC(0.98-pad1->GetRightMargin(), 0.93-pad1->GetTopMargin(),
	    		       (std::string("p-value = ")+probstr).c_str());
	    latex.SetTextColor(kOrange+2);
	    std::string minstr(Form("%.0f", fitRangeMin_));
	    std::string maxstr(Form("%.0f", fitRangeMax_));
	    latex.DrawLatexNDC(0.98-pad1->GetRightMargin(), 0.88-pad1->GetTopMargin(),
	    		       (minstr+std::string(" < M_{12} < ")+maxstr).c_str());

	    canvas.cd();
	    TPad *pad2 = new TPad("pad2","",0,0.0,1,0.265);
	    pad2->SetTopMargin(1);
	    pad2->SetBottomMargin(0.33);
	    pad2->SetLeftMargin(pad1->GetLeftMargin());
	    pad2->SetRightMargin(pad1->GetRightMargin());
	    pad2->SetGridy();
	    pad2->Draw();
	    pad2->cd();
	    frame2->SetTitle("");
	    frame2->GetXaxis()->SetTitleSize(0.15);
	    frame2->GetXaxis()->SetTitleOffset(0.9);
	    frame2->GetXaxis()->SetLabelSize(0.115);
	    frame2->GetXaxis()->SetLabelOffset(0.010);
	    frame2->SetYTitle("Pulls");
	    frame2->GetYaxis()->CenterTitle(kTRUE);
	    frame2->GetYaxis()->SetTitleSize(0.14);
	    frame2->GetYaxis()->SetTitleOffset(0.4);
	    frame2->GetYaxis()->SetNdivisions(3,5,0);
	    frame2->GetYaxis()->SetLabelSize(0.115);
	    frame2->GetYaxis()->SetLabelOffset(0.011);
	    frame2->SetMinimum(-5.);
	    frame2->SetMaximum(+5.);
	    frame2->Draw();

	    canvas.Modified();
	    canvas.Update();
	    canvas.SaveAs((plotDir_+name+"_lowM_linear.pdf").c_str());
	    canvas.SaveAs((plotDir_+name+"_lowM_linear.png").c_str());
	    pad1->SetLogy();
	    frame->GetYaxis()->SetRangeUser(0.001, frame->GetMaximum()*5);
	    canvas.Modified();
	    canvas.Update();
	    canvas.SaveAs((plotDir_+name+"_lowM_log.pdf").c_str());
	    canvas.SaveAs((plotDir_+name+"_lowM_log.png").c_str());

	    return fitResult;

}


std::unique_ptr<RooFitResult> FitContainer::backgroundOnlyFit(const std::string& name, const bool& plot_params, const bool& control_region) {
  if (!initialized_) initialize();

	auto &Pdf 		= *GetFromRooWorkspace<RooAbsPdf>(workspace_,toString(Type::background));
	auto &data	 	= *GetFromRooWorkspace<RooAbsData>(workspace_,data_);
	auto &mbb		= *GetFromRooWorkspace<RooRealVar>(workspace_, mbb_);
	gStyle->SetPadRightMargin(0.025);
//	RooFit::SumW2Error(kTRUE);

	//class-helper
	RooFitQuality fit_quality;

  // Split Range for blinded data and perform simultaneous fit by CA
  std::unique_ptr<RooFitResult>
    fitResult(Pdf.fitTo(data,
			RooFit::Save(),
			RooFit::PrintLevel(verbosity_),
//			RooFit::Range(fitRangeLowId_.c_str()),
			RooFit::Range(fitSplRangeId_.c_str()),
//			RooFit::SumW2Error(kTRUE),
			RooFit::SplitRange()
    ));

  fit_quality.PrintParametersInfo(*fitResult);

  // Top frame
  std::unique_ptr<RooPlot> frame(mbb.frame());
  data.plotOn(frame.get(),
		 RooFit::Binning(nBinsToPlot_),
	     RooFit::Name("data_curve"));

    mbb.setBins(nBinsToPlot_);
  //Get TH1D from data
  TH1D *hData = static_cast<TH1D*>(data.createHistogram("data",mbb));

  std::string bg_curve_name = "background_curve";

//  Uncertainty bands
//  2 Sigma
  Pdf.plotOn(frame.get(),
		  RooFit::VisualizeError(*fitResult, 2, false),
		  RooFit::LineColor(kRed),
		  RooFit::LineStyle(kSolid),
		  RooFit::FillColor(kOrange),
		  RooFit::Name( (bg_curve_name + "_2sigma").c_str()),
		  RooFit::Normalization(data.sumEntries("1", fitRangeId_.c_str()) , RooAbsReal::NumEvent),
		  RooFit::MoveToBack() );

  //  1 Sigma
  Pdf.plotOn(frame.get(),
		  RooFit::VisualizeError(*fitResult, 1, false),
		  RooFit::LineColor(kRed),
		  RooFit::LineStyle(kSolid),
		  RooFit::FillColor(kGreen+1),
		  RooFit::Name( (bg_curve_name + "_1sigma").c_str()),
		  RooFit::Normalization(data.sumEntries("1", fitRangeId_.c_str()) , RooAbsReal::NumEvent));

//  Central
  Pdf.plotOn(frame.get(),
		//RooFit::DrawOption("L"),
		RooFit::LineColor(kRed),
		RooFit::Name(bg_curve_name.c_str()),
		RooFit::NormRange(fullRangeId_.c_str()),
		RooFit::Range(fitRangeId_.c_str()),
		RooFit::Normalization(data.sumEntries("1", fitRangeId_.c_str()), RooAbsReal::NumEvent));


  if(plot_params){
	  double par_xmin = 0.65, par_xmax = 0.9, par_ymax = 0.6;
//	  par_xmin = 0.2; par_xmax = 0.45; par_ymax = 0.9;
	  Pdf.paramOn(frame.get(),RooFit::Layout(par_xmin,par_xmax,par_ymax));//0.98-pad1->GetRightMargin(),0.83-pad1->GetTopMargin()));
	  frame->getAttText()->SetTextSize(0.03);
  }

  int nPars = fitResult->floatParsFinal().getSize();

  //TODO: Remove this part!!!!!!! This is a Fitter and not a generator. ->
  //-> this class takes care about the fitting procedure but not about the
  //generation of the new data! Please don't mix different purposes in a single class.
  //
  // Creat Asimov data for combine tool
//  RooDataSet* asimov = Pdf.generate(mbb, obs_);
//  asimov->SetName("data_obs");
//  workspace_.import(*asimov);
 
  // Get Covariance Matrix for diagonalisation 
  TMatrixDSymEigen CM = fitResult->covarianceMatrix();
  TMatrixD covVec = CM.GetEigenVectors();
  TMatrixD covMat = fitResult->covarianceMatrix(); 
  covMat.GetMatrix2Array(covMatrix_);
  covVec.GetMatrix2Array(eigenVector_);

  ////////////////////////////////////////////////////////////////////////////////
  //////////////////////// ChiSquare by Roofit ///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  Chi2Ndf chi2ndf;
  if (!splitrange_) {
	  //TODO: implement method(reload existed one in RooFitQuality)
  	ndfBkgOnly_ = getNonZeroBins_(data) - nPars;
  	normChi2BkgOnly_ = frame->chiSquare("background_curve", "data_curve", nPars); //chi^2 from RooFit (RooPlot::chiSquare())
  	chi2BkgOnly_ = normChi2BkgOnly_ * ndfBkgOnly_;
  	bkgOnlyFit_.Fill();
  }

  ////////////////////////////////////////////////////////////////////////////////
  ////////////////////////ChiSquare by Chayanit///////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  else {
    chi2ndf = fit_quality.chiSquare(*frame, "background_curve", "data_curve", nPars, blind_lowEdge_, blind_highEdge_);
    ndfBkgOnly_ 		= chi2ndf.ndf;
    chi2BkgOnly_ 		= chi2ndf.chi2;
    normChi2BkgOnly_ 	= chi2BkgOnly_/ndfBkgOnly_;
  	bkgOnlyFit_.Fill();
  }
  std::string chi2str(Form("%.1f/%d = %.1f", chi2BkgOnly_,
		 	   ndfBkgOnly_, normChi2BkgOnly_));
  std::cout << "\nNormalized chi^2: " << chi2str << std::endl;

  // TMath::Prob p-value by CA
  double prob = TMath::Prob(chi2BkgOnly_,ndfBkgOnly_);
  std::string probstr(Form("%.2f", prob));
  std::cout << "\nTMath::Prob(#chi^{2},nof) : " << probstr << std::endl;
 
  // Bottom frame
  RooHist* hpull;
  hpull = frame->pullHist();
  std::unique_ptr<RooPlot> frame2(mbb.frame(hData->GetXaxis()->GetXmin(),hData->GetXaxis()->GetXmax()));

  TCanvas canvas;
  canvas.cd();
  prepareFrame_(*frame);
  prepareFrame_(*frame2);

  //Size of the canva in pixels:
  double pad_hPixel = canvas.YtoPixel(canvas.GetY1());

  double pad1_ymin = 0.3;
//  auto pad1 = new TPad("pad1","",0,pad1_ymin,1,1);
  auto pad1 = HbbStyle::getRatioTopPad(1 - pad1_ymin);
  pad1->Draw();
  pad1->cd();
  //Size of the pad1 pad in pixels:
  double pad1_hPixel = pad1->YtoPixel(pad1->GetY1());
  HbbStyle::setRatioTopFrame(frame->GetXaxis(),frame->GetYaxis(),pad_hPixel);
  frame->Draw();

  std::string lumistr(Form("%.1f", lumi_));

  TLatex latex;
  latex.SetTextFont(43);

  latex.SetTextSize(19);
  latex.SetTextAlign(33);
  latex.SetTextColor(kBlue+2);
  latex.DrawLatexNDC(0.96-pad1->GetRightMargin(), 0.96-pad1->GetTopMargin(),
                     (std::string("#chi^{2}/dof = ")+chi2str).c_str());
  latex.SetTextColor(kGreen+2);
  latex.DrawLatexNDC(0.96-pad1->GetRightMargin(), 0.91-pad1->GetTopMargin(),
  		       (std::string("p-value = ")+probstr).c_str());
  latex.SetTextColor(kOrange+2);
  std::string minstr(Form("%.0f", fitRangeMin_));
  std::string maxstr(Form("%.0f", fitRangeMax_));
  latex.DrawLatexNDC(0.96-pad1->GetRightMargin(), 0.86-pad1->GetTopMargin(),
  		       (minstr+std::string(" < M_{12} < ")+maxstr).c_str());

  //Draw TLegend
  float legendx1, legendx2, legendy1, legendy2;
  legendx1 = 0.6; legendx2 = 0.96 - pad1->GetRightMargin();
  legendy1 = 0.6; legendy2 = 0.81 - pad1->GetTopMargin();
  TLegend leg(legendx1,legendy1,legendx2,legendy2);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetTextFont(43);
//  leg.SetTextAlign(33);
  leg.SetTextSize(latex.GetTextSize());
  std::string which_data = "";
  if (control_region) which_data += " control region";
  else which_data += " signal region";
  leg.AddEntry(frame->findObject("data_curve"),("Data" + which_data).c_str(),"p");
  leg.AddEntry(frame->getCurve("background_curve"),"Fit","l");
  leg.Draw();

  canvas.Modified();
  canvas.cd();
//  TPad *pad2 = new TPad("pad2","",0,0.0,1,pad1_ymin);
  auto pad2 = HbbStyle::getRatioBottomPad(pad1_ymin);
  pad2->Draw();
  pad2->cd();

  bool drawUncertaintyBands = true;
  hBands bands;
  if(drawUncertaintyBands) {
	  bands = getPullBands_(frame.get(),bg_curve_name, hData, mbb, Pdf);
	  frame2->addTH1(bands.h2sigmaU,"HIST");
	  frame2->addTH1(bands.h2sigmaD,"HIST");
	  frame2->addTH1(bands.h1sigmaU,"HIST");
	  frame2->addTH1(bands.h1sigmaD,"HIST");
	  frame2->addTH1(bands.central,"E1");
	  frame2->SetYTitle("#frac{Data-Fit}{#sqrt{Fit}}");
  }
  else {
	  frame2->addPlotable(hpull,"P");
	  frame2->SetYTitle("#frac{Data-Fit}{#sqrt{Data}}");
  }

  HbbStyle::setRatioBottomFrame(frame2->GetXaxis(),frame2->GetYaxis(),pad_hPixel,pad1_hPixel);
  frame2->GetYaxis()->SetRangeUser(-3.9999, 3.9999);
  frame2->Draw();
  //Draw TLine for 0 pulls
  TLine line(frame2->GetXaxis()->GetXmin(),0,frame2->GetXaxis()->GetXmax(),0);
  line.SetLineStyle(2);
  line.Draw();
  canvas.cd();
  HbbStyle::drawStandardTitle("out");

  canvas.SaveAs((plotDir_+name+"_lowM_linear.pdf").c_str());
  pad1->SetLogy();
  frame->GetYaxis()->SetRangeUser(0.1, frame->GetMaximum()*5);
  canvas.Modified();
  canvas.Update();
  canvas.SaveAs((plotDir_+name+"_lowM_log.pdf").c_str());

  return fitResult;
}


void FitContainer::profileModel(const Type& type) {
  RooAbsPdf& model= *(workspace_.pdf(toString(type).c_str()));
  if (&model == nullptr) {
    std::stringstream msg;
    msg << "No model of type '" << toString(type) << "' is set, yet.";
    throw std::logic_error(msg.str());
  }

  // get the objects from the workspace:
  RooAbsData& data = *workspace_.data(data_.c_str());
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());

  std::unique_ptr<RooAbsReal> nll(model.createNLL(data));
  std::unique_ptr<RooArgSet> parameters(model.getParameters(mbb));
  std::unique_ptr<TIterator> iter(parameters->createIterator());
  // use raw pointer for 'parameter' because 'model' owns the object it points to:
  auto parameter = static_cast<RooRealVar*>(iter->Next());
  while (parameter) {
    if (!(parameter->isConstant())) {
      std::unique_ptr<RooAbsReal> profile(nll->createProfile(*parameter));
      std::unique_ptr<RooPlot> frame(parameter->frame());
      if (frame == nullptr) {
        std::stringstream msg;
        msg << "Problems creating frame for '" << parameter->GetName() << "'.";
        throw std::runtime_error(msg.str());
      }
      profile->plotOn(frame.get());
      TCanvas canvas("canvas", "", 600, 600);
      prepareCanvas_(canvas);
      prepareFrame_(*frame);
      frame->Draw();
      canvas.SaveAs((plotDir_+toString(type)+"_profile_"+
                     parameter->GetName()+".pdf").c_str());
      canvas.SetLogy();
      canvas.SaveAs((plotDir_+toString(type)+"_profile_"+
                     parameter->GetName()+"_log.pdf").c_str());
    }
    parameter = static_cast<RooRealVar*>(iter->Next());
  }
}


void FitContainer::showModels() const {
  std::cout << "\n=============================================" << std::endl;
  std::cout << "Defined Models" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  RooArgSet models(workspace_.allPdfs());
  std::unique_ptr<TIterator> itModel(models.createIterator());
  // use raw pointer for 'model' because 'models' owns the object it points to:
  auto model = static_cast<RooAbsPdf*>(itModel->Next());
  while (model) {
    model->Print();
    std::unique_ptr<RooArgSet> parameters(model->getParameters(mbb));
    std::unique_ptr<TIterator> itPar(parameters->createIterator());
    // use raw pointer for 'parameter' because 'model' owns the object it points to:
    auto parameter = static_cast<RooRealVar*>(itPar->Next());
    while (parameter) {
      parameter->Print();
      parameter = static_cast<RooRealVar*>(itPar->Next());
    }
    model = static_cast<RooAbsPdf*>(itModel->Next());
    std::cout << "---------------------------------------------" << std::endl;
  }
  std::cout << std::endl;
}

FitContainer& FitContainer::verbosity(int level) {
  verbosity_ = level;
  return *this;
}


FitContainer& FitContainer::fitRangeMin(float min) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  fitRangeMin_ = min;
  mbb.setMin(min);
  mbb.setRange(fullRangeId_.c_str(), mbb.getMin(), mbb.getMax());
  mbb.setRange(fitRangeId_.c_str(), fitRangeMin_, fitRangeMax_);
  mbb.setRange(fitRangeLowId_.c_str(), fitRangeMin_, blind_lowEdge_);
  mbb.setRange(fitRangeHighId_.c_str(), blind_highEdge_, fitRangeMax_);

  return *this;
}


FitContainer& FitContainer::fitRangeMax(float max) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  fitRangeMax_ = max;
  mbb.setMax(max);
  mbb.setRange(fullRangeId_.c_str(), mbb.getMin(), mbb.getMax());
  mbb.setRange(fitRangeId_.c_str(), fitRangeMin_, fitRangeMax_);
  mbb.setRange(fitRangeLowId_.c_str(), fitRangeMin_, blind_lowEdge_);
  mbb.setRange(fitRangeHighId_.c_str(), blind_highEdge_, fitRangeMax_);

  return *this;
}

FitContainer& FitContainer::setNBins(int nbins) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  nbins_ = nbins;
  mbb.setBins(nbins_); 

  return *this;
}

std::string FitContainer::getOutputPath_(const std::string& subdirectory) {
  std::string path = outputDir_ + subdirectory;
  gSystem->mkdir(path.c_str(), true);
  path += "/FitContainer_";
  return path;
}


void FitContainer::prepareCanvas_(TCanvas& raw) {
  raw.SetLeftMargin(0.15);
  raw.SetRightMargin(0.05);
}


void FitContainer::prepareFrame_(RooPlot& raw) {
//  raw.GetYaxis()->SetTitleOffset(2);
  raw.SetTitle("");
}

int FitContainer::getNonZeroBins_(const RooAbsData& data) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  std::unique_ptr<TH1> hist(data.createHistogram(mbb_.c_str(), mbb));
  int nonZeroBins = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double center = hist->GetBinCenter(i);
    if (hist->GetBinContent(i) > 0 &&
        center > fitRangeMin_ && center < fitRangeMax_) ++nonZeroBins;
  }
  std::cout << "number of Non-Zero bins : " << nonZeroBins << std::endl;
  return nonZeroBins;
}


int FitContainer::getBlindedBins_(const RooAbsData& data, double blind_lowEdge, double blind_highEdge) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  std::unique_ptr<TH1> hist(data.createHistogram(mbb_.c_str(), mbb));
  int blindedBins = 0;
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double center = hist->GetBinCenter(i);
    if (hist->GetBinContent(i) > 0. && center > blind_lowEdge && center < blind_highEdge) ++blindedBins;
  }
  std::cout << "number of blinded bins : " << blindedBins << std::endl;
  return blindedBins;
}

bool FitContainer::applyModifiers_(RooAbsPdf& pdf,
                                   const std::vector<ParamModifier>& modifiers) {
  bool modified = false;
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  std::unique_ptr<RooArgSet> parameters(pdf.getParameters(mbb));
  std::unique_ptr<TIterator> iter(parameters->createIterator());
  // use raw pointer for 'parameter' because 'pdf' owns the object it points to:
  RooRealVar* parameter = static_cast<RooRealVar*>(iter->Next());
  while (parameter) {
    for (const auto& m : modifiers) {
    		m.show();
      if (m.modify(*parameter)) {
        modified = true;
      }
    }
    parameter = static_cast<RooRealVar*>(iter->Next());
  }
  if (modifiers.size() > 0 && !modified) {
    std::cerr << ">>> None of the modifiers provided to '" << pdf.GetName()
              << "' pdf could be applied." << std::endl;
    std::cerr << ">>> Provided modifiers: ";
    for (auto m = modifiers.cbegin(); m != modifiers.cend(); ++m) {
      if (m != modifiers.cbegin()) std::cerr << ",";
      std::cerr << " " << m->name();
    }
    std::cerr << std::endl;
    std::cerr << ">>> Found in pdf: ";
    parameters->Print();
  }
  return modified;
}

double FitContainer::getPeakStart_(const Type& type) {
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  return getPeakStart_(type, mbb.getMax());
}


double FitContainer::getPeakStart_(const Type& type,const double& max) {
	/*
	 * Method to get estimated position of the peak.
	 * BG - max bin content OR (max + min)/2
	 * Signal - max bin content
	 */
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  double peakStart = (mbb.getMin() + max) / 2.0;
  switch (type) {
  case Type::signal:
    if (signal_ != "") {
      RooAbsData& signal = *workspace_.data(signal_.c_str());
      peakStart = getMaxPosition_(signal);
      return peakStart;
    }
    break;
  case Type::background:
    if (bkg_ != "") {
      RooAbsData& bkg = *workspace_.data(bkg_.c_str());
      peakStart = getMaxPosition_(bkg);
    }
    break;
  }
  return peakStart < max ? peakStart : max;
}


double FitContainer::getMaxPosition_(const RooAbsData& data) {
	/*
	 * Method to return x of the bin with a max content
	 */
  RooRealVar& mbb = *workspace_.var(mbb_.c_str());
  std::unique_ptr<TH1> hist(data.createHistogram(mbb_.c_str(), mbb));
  int maximumBin = hist->GetMaximumBin();
  return hist->GetBinCenter(maximumBin);
}

void FitContainer::makeLog_(const RooFitResult& fitResult){
	std::filebuf fb;
	fb.open((plotDir_ + "log.txt").c_str(),std::ios::out);
	std::ostream f(&fb);
	f<<"\n Normalized chi^2: "<<normChi2BkgOnly_<<" Probability: "<<TMath::Prob(chi2BkgOnly_,ndfBkgOnly_);
	f<<"\n constant parameters: \n";
	fitResult.constPars().printMultiline(f,1111,1);
	f<<"\n floating parameters (init): \n";
	fitResult.floatParsInit().printMultiline(f,1111,1);
	f<<"\n floating parameters (final): \n";
	fitResult.floatParsFinal().printMultiline(f,1111,1);
	f<<"\n cov.matrix: I HAVE NO IDEA HOW TO WRITE IT!!!!!!\n";
	fb.close();
}

void FitContainer::SetupTopFrame(RooPlot *frame1, const double& pad_w, const double& pad_h){

//	//Pixel sizes of the current pad
	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());

	//Use relative values
	frame1->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y") * pad_w / pad_wPixel);
	frame1->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y"));
	frame1->GetYaxis()->SetRangeUser(frame1->GetMinimum(), frame1->GetMaximum()+200);

	frame1->GetXaxis()->SetTickLength(gStyle->GetTickLength("X") * pad_h / pad_hPixel);
    frame1->GetXaxis()->SetTitleSize(0);
    frame1->GetXaxis()->SetLabelSize(0);
	frame1->GetXaxis()->SetTitleOffset(999); //Effectively turn off x axis title on main plot
	frame1->GetXaxis()->SetLabelOffset(999); //Effectively turn off x axis label on main plot
}

void FitContainer::SetupBottomFrame(RooPlot *frame2, const double& pad_w, const double& pad_h){

	double pad_wPixel = gPad->XtoPixel(gPad->GetX2());
	double pad_hPixel = gPad->YtoPixel(gPad->GetY1());

	//Y-axis options
	frame2->GetYaxis()->SetTitleFont(43);			// This gives the sizes in pixels!!!!
	frame2->GetYaxis()->SetLabelFont(43);			// This gives the sizes in pixels!!!!
	frame2->GetYaxis()->SetTickLength(gStyle->GetTickLength("Y") * pad_w / pad_wPixel);
	frame2->GetYaxis()->SetLabelOffset(gStyle->GetLabelOffset("Y") * pad_w / pad_wPixel );
	frame2->GetYaxis()->SetTitleOffset(gStyle->GetTitleOffset("Y") * 1.15);
	frame2->GetYaxis()->SetLabelSize(gStyle->GetLabelSize("Y") * pad_h);
	frame2->GetYaxis()->SetTitleSize(gStyle->GetTitleSize("Y") * pad_h);
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetNdivisions(105);
    frame2->SetMinimum(-3.9999);
    frame2->SetMaximum(+3.9999);

    //X-axis options
	frame2->GetXaxis()->SetTitleFont(43);			// This gives the sizes in pixels!!!!
	frame2->GetXaxis()->SetLabelFont(43);			// This gives the sizes in pixels!!!!
    frame2->GetXaxis()->SetTitle(HbbStyle::axisTitleMass());
    frame2->GetXaxis()->SetTickLength(gStyle->GetTickLength("X") * pad_h / pad_hPixel);
    frame2->GetXaxis()->SetLabelSize(gStyle->GetLabelSize("X") * pad_h);
    frame2->GetXaxis()->SetTitleSize(gStyle->GetTitleSize("X") * pad_h);
    	frame2->GetXaxis()->SetTitleOffset(3.5);
    	frame2->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X") * pad_h / pad_hPixel );
}

FitContainer::hBands FitContainer::getPullBands_(RooPlot * frame, const std::string& curve_name, TH1D * hData, RooRealVar & x, RooAbsPdf & fit){
	/*
	 * Method to calculate 1 and 2 sigma badns for the
	 * pull distribution
	 */
	hBands bands;
	TH1::SetDefaultSumw2(1);

    TH1D *hcentral = new TH1D("hcentral", "", hData->GetNbinsX(), hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    bands.central = hcentral;
    bands.central->SetMarkerColor(1);

    TH1D *h1sigmaU = new TH1D("h1sigmaU", "", nbins_, hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    bands.h1sigmaU = h1sigmaU;
    bands.h1sigmaU->SetFillColor(kGreen+1);
    bands.h1sigmaU->SetLineColor(kGreen+1);
    bands.h1sigmaU->SetFillStyle(1001);

    TH1D *h1sigmaD = new TH1D("h1sigmaD", "", nbins_, hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    bands.h1sigmaD = h1sigmaD;
    bands.h1sigmaD->SetFillColor(kGreen+1);
    bands.h1sigmaD->SetLineColor(kGreen+1);
    bands.h1sigmaD->SetFillStyle(1001);

    TH1D *h2sigmaU = new TH1D("h2sigmaU", "", nbins_, hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    bands.h2sigmaU = h2sigmaU;
    bands.h2sigmaU->SetFillColor(kOrange);
    bands.h2sigmaU->SetLineColor(kOrange);
    bands.h2sigmaU->SetFillStyle(1001);

    TH1D *h2sigmaD = new TH1D("h2sigmaD", "", nbins_, hData->GetXaxis()->GetXmin(), hData->GetXaxis()->GetXmax());
    bands.h2sigmaD = h2sigmaD;
    bands.h2sigmaD->SetFillColor(kOrange);
    bands.h2sigmaD->SetLineColor(kOrange);
    bands.h2sigmaD->SetFillStyle(1001);

    //Retrieve information from the RooPlot
    auto *sigma2 = frame->getCurve((curve_name + "_2sigma").c_str());
    auto *sigma1 = frame->getCurve((curve_name + "_1sigma").c_str());
    auto *nominal= frame->getCurve(curve_name.c_str());

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

        err1D = -1 * err1U;
        err2D = -1 * err2U;

//        auto n_data_i = hData->GetBinContent( bin_i );
        double divide_by = n_i;//n_data_i;//n_i;

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

        std::cout<<"j = "<<bin_i<<"WTFFF: x = "<<x<<" y = "<<n_i<<" up1 = "<<err1U<<" lo1 = "<<err1D<<" up2 = "<<err2U<<" lo2 = "<<err2D<<std::endl;
    }

    for(int bin_i = 1; bin_i < hData->GetNbinsX()+1; ++bin_i){
        auto n_data_i = hData->GetBinContent( bin_i );
        auto e_data_i = hData->GetBinError( bin_i );

        //Evaluate value of the Bg-function
        double bkg_int_i = nominal->Eval(hData->GetBinCenter(bin_i));
        double divide_by = sqrt(bkg_int_i);//ssqrt(n_data_i);//
        auto pull_i = (n_data_i-bkg_int_i)/divide_by;

        if(n_data_i == 0) {
        		bands.central->SetBinContent(bin_i,0);
        		bands.central->SetBinError( bin_i,0);
        		bands.h1sigmaU->SetBinContent( bin_i,0);
        		bands.h1sigmaD->SetBinContent( bin_i,0);
        		bands.h2sigmaU->SetBinContent( bin_i,0);
        		bands.h2sigmaD->SetBinContent( bin_i,0);
        }
        else {
            bands.central->SetBinContent( bin_i, pull_i);
            bands.central->SetBinError( bin_i, e_data_i/divide_by);
        }
    }
    x.getValV();fit.getValV();

    /*
    for(int i =0; i <nominal->GetN(); ++i){
    		int j = i + 2;
    		if( sigma1X[j] != nominalX[i] || sigma2X[j] != nominalX[i]) throw std::logic_error("SOMETHING WRONG with point in RooCurves at  FitContainer::getPullBands_");
        auto x_i = nominalX[i];
        auto n_i = nominalY[i];
        auto n_1up_i = n1sigmaY[(sigma1->GetN())-j-1];
        auto n_1down_i = n1sigmaY[j];
        auto n_2up_i = n2sigmaY[(sigma2->GetN())-j-1];
        auto n_2down_i = n2sigmaY[j];
        auto bin_i = h1sigmaU->FindBin( x_i );

        std::cout<<"bands: i = "<<i<<" x = "<<x_i<<" y = "<<n_i<<" up1 = "<<n_1up_i<<" as = "<<n1sigmaY[sigma1->GetN() - i]<<" down = "<<n1sigmaY[i]<<" lol = "<<sigma1->average(nominalX[i],nominalX[i+1])<<std::endl;
        if(bin_i <= 0 || bin_i > h1sigmaU->GetNbinsX()) continue;

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

        h1sigmaU->SetBinContent( i, err1U);
        h1sigmaD->SetBinContent( i, err1D);
        h2sigmaU->SetBinContent( i, err2U);
        h2sigmaD->SetBinContent( i, err2D);
//        bands.h1sigmaU->SetBinContent( bin_i, (-n_i+n_1up_i));
//        bands.h1sigmaD->SetBinContent( bin_i, (-n_i+n_1down_i));
//        bands.h2sigmaU->SetBinContent( bin_i, (-n_i+n_2up_i));
//        bands. h2sigmaD->SetBinContent( bin_i, (-n_i+n_2down_i));
    }
    */

    //redefine pulls here
    /*
    auto bkg_int = fit.createIntegral(RooArgSet(x))->getVal();
    for(int bin_i = 1; bin_i < hData->GetNbinsX()+1; ++bin_i){
//    		std::cout<<"WTFsss: data = "<<bin_i<<" ["<<hData->GetBinLowEdge(bin_i)<<","<<hData->GetBinLowEdge(bin_i) + hData->GetBinWidth( bin_i )<<"]"<<std::endl;
        auto n_data_i = hData->GetBinContent( bin_i );
        auto e_data_i = hData->GetBinError( bin_i );
        auto bin_width_i = hData->GetBinWidth( bin_i );
        x.setRange( ("Bin" + std::to_string(bin_i)).c_str(), hData->GetBinLowEdge(bin_i), hData->GetBinLowEdge(bin_i)+bin_width_i);
        auto bkg_int_i = fit.createIntegral(RooArgSet(x), ("Bin" + std::to_string(bin_i)).c_str())->getVal() * hData->Integral() / bkg_int;// / bkg_int * nominal->getVal();

        double divide_by = sqrt(bkg_int_i);
        auto pull_i = (n_data_i-bkg_int_i)/divide_by;

        if(n_data_i == 0) {
        		bands.central->SetBinContent(bin_i,0);
        		bands.central->SetBinError( bin_i,0);
        		bands.h1sigmaU->SetBinContent( bin_i,0);
        		bands.h1sigmaD->SetBinContent( bin_i,0);
        		bands.h2sigmaU->SetBinContent( bin_i,0);
        		bands.h2sigmaD->SetBinContent( bin_i,0);
        }
        else {
            bands.central->SetBinContent( bin_i, pull_i);
            bands.central->SetBinError( bin_i, e_data_i/divide_by);
        }
    }
    */

	return bands;
}

const std::string FitContainer::defaultOutputDir_ =
  std::string(gSystem->Getenv("CMSSW_BASE"))+"/src/Analysis/BackgroundModel/test/";
