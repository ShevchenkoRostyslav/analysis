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
		lumi_(35.7), //2.69,12.9,36.62
		obs_(259399.) //SR1-259399, SR2-105053, SR3-26760
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


std::unique_ptr<RooFitResult> FitContainer::backgroundOnlyFit(const std::string& name, const bool& plot_params) {
  if (!initialized_) initialize();

	auto &Pdf 		= *GetFromRooWorkspace<RooAbsPdf>(workspace_,toString(Type::background));
	auto &data	 	= *GetFromRooWorkspace<RooAbsData>(workspace_,data_);
	auto &mbb		= *GetFromRooWorkspace<RooRealVar>(workspace_, mbb_);
	RooFit::SumW2Error(kTRUE);

	//class-helper
	RooFitQuality fit_quality;

  // Split Range for blinded data and perform simultaneous fit by CA
  std::unique_ptr<RooFitResult>
    fitResult(Pdf.fitTo(data,
			RooFit::Save(),
			RooFit::PrintLevel(verbosity_),
			RooFit::Range(fitRangeLowId_.c_str()),
			RooFit::SumW2Error(kTRUE),
//			RooFit::Range(fitSplRangeId_.c_str()),
//			RooFit::Offset(true),
			RooFit::SplitRange()
  //RooFit::InitialHesse(kTRUE)
    ));

  fit_quality.PrintParametersInfo(*fitResult);

  // Top frame
  std::unique_ptr<RooPlot> frame(mbb.frame());
  data.plotOn(frame.get(),
	     RooFit::MarkerSize(0.8),	//0.8 for lowM
	     //RooFit::DataError(RooAbsData::Auto), 
	     RooFit::Name("data_curve"));
  Pdf.plotOn(frame.get(),
	     //RooFit::VisualizeError(*fitResult,1,kFALSE),
             //RooFit::DrawOption("L"),
             RooFit::LineColor(kRed),
             //RooFit::LineWidth(2),
	     //RooFit::LineStyle(kDashed),
	     //RooFit::FillColor(kOrange),
             RooFit::Name("background_curve"),
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

  int nPars = fitResult->floatParsFinal().getSize();

  //TODO: Remove this part!!!!!!! This is a Fitter and not a generator. -> 
  //-> this class takes care about the fitting procedure but not about the 
  //generation of the new data! Please don't mix different purposes in a single class.
  //
  // Creat Asimov data for combine tool
  RooDataSet* asimov = Pdf.generate(mbb, obs_);
  asimov->SetName("data_obs");
  workspace_.import(*asimov); 
 
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
  ////////////////////////ChiSquare by Gregor/////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////
  //double myChi2 = chiSquare_(data,
  //                           *(frame->getCurve("background_curve")),
  //				blind_lowEdge_, blind_highEdge_, nPars);
  //std::string myChi2str(Form("%.1f/%d = %.1f", myChi2, ndfBkgOnly_,
  //			     myChi2/ndfBkgOnly_));
  //std::cout << "\nMy normalized chi^2: " << myChi2str << std::endl;
  //latex.DrawLatexNDC(0.98-canvas.GetRightMargin(), 0.93-canvas.GetTopMargin(),
  //                   (std::string("#chi^{2}/ndf = ")+myChi2str).c_str());

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
  std::cout << "\nTMath::Prob(#chi^{2}_{RooFit},ndf) : " << probstr << std::endl;
 
  // Bottom frame
  RooHist* hpull;
  hpull = frame->pullHist();
  hpull->SetMarkerSize(0.8);	//0.8 for lowM
  std::unique_ptr<RooPlot> frame2(mbb.frame());
  frame2->addPlotable(hpull,"P");
  //frame2->addObject(frame->pullHist()); 
  frame2->SetMinimum(-5.); 
  frame2->SetMaximum(+5.); 

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

  std::string lumistr(Form("%.1f", lumi_));

  TLatex latex;
  latex.SetTextFont(43);
  latex.SetTextSize(17);
  latex.SetTextAlign(11);
  latex.DrawLatexNDC(0.31-pad1->GetRightMargin(), 1.02-canvas.GetTopMargin(),
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
  frame2->Draw();

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
  raw.GetYaxis()->SetTitleOffset(2);
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

const std::string FitContainer::defaultOutputDir_ =
  std::string(gSystem->Getenv("CMSSW_BASE"))+"/src/Analysis/BackgroundModel/test/";
