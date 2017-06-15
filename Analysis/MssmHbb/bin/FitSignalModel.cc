#include <vector>
#include <array>
#include <map>

#include "stdlib.h"
#include <memory>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/join.hpp>
#include "TSystem.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "TFormula.h"

#include "RooDataHist.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooWorkspace.h"

#include "TFile.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH2.h"
#include "TAxis.h"

#include "Analysis/BackgroundModel/interface/HistContainer.h"
#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/BackgroundModel/interface/ParamModifier.h"
#include "Analysis/BackgroundModel/interface/Tools.h"
#include "Analysis/BackgroundModel/interface/ProbabilityDensityFunctions.h"
#include "Analysis/BackgroundModel/interface/RooQuadGausExp.h"

//style file
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
HbbStyle style;

#include "Analysis/MssmHbb/interface/utilLib.h"

// hi my name is Rostyslav and I like bananas

namespace po = boost::program_options;
namespace ab = analysis::backgroundmodel;

typedef std::unique_ptr<TFile> pTFile;

using namespace analysis::backgroundmodel;
using namespace std;

void drawSystFits(const int & mass,const std::string &syst, const std::string &output);
void plotShapeParameters(const int & mass,const vector<string> &syst, const std::string &folder, const vector<string>& parameters);
void fixParameters(const string& name, RooWorkspace& w, const RooWorkspace& central);
void fixParameters(const string& name, RooWorkspace& w, const double& val);
void fixShiftParameter(RooWorkspace& w, const RooWorkspace& central);
void setParameters(RooWorkspace& w, const RooWorkspace& central);

int main(int argc, char* argv[]) {
  const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));
  style.set(PRELIMINARY);
  gStyle->SetTitleXOffset(1.3);
  gStyle->SetTitleYOffset(1.15);
  gStyle->SetTitleSize(0.04,"XYZ");

	std::map<int,std::string> signal;
	string selection = "";//"ReReco_36615fb_SFb";//
	string mssm_normTanB = "";//"_NormToTanB30";
	signal[300] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[350] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[400] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[500] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[600] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[700] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[900] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[1100] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
	signal[1300] = "/src/Analysis/MssmHbb/output/MssmHbbSignal" + selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8" + mssm_normTanB + ".root";
/**/
  int rebin = 1;
  const int verbosity = 1;
  double min = 0;
  double max = 1700;
  std::string histo_name = "templates/bbH_Mbb_VIS";
  std::string histo_name_240 = "bbH_Mbb";
  //List of shape uncertainties
  std::vector<std::string> shape_unc = {"CMS_scale_j_13TeV","CMS_res_j_13TeV","CMS_eff_pTonl_13TeV","CMS_eff_b_13TeV"};//{};//
  map<string,int> colors;
  colors["Up"] = 4;
  colors["Down"] = 2;
  std::string CMS = "CMS";
  std::string energy = "13TeV";
  std::array<std::string,2> sign = { {"Up","Down"} };
  std::string model = "doublegausexp";
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  for(const auto & mass : signal){

	  std::cout<<"***************************************"<<std::endl;
	  std::cout<<"*************** M = "<<mass.first<<"***************"<<std::endl;
	  std::cout<<"***************************************"<<std::endl;
	  TFile f((cmsswBase + mass.second).c_str(),"read");
	  std::string name = cmsswBase + "/src/Analysis/MssmHbb/output/" + "ReReco_signal_M-" + std::to_string(mass.first) + mssm_normTanB;
	  boost::filesystem::remove_all(name); // clean dir
	  std::string hname =  histo_name_240; //histo_name;//
	  //Central value
	  TH1D* hSignal = (TH1D*) f.Get(hname.c_str());
	  if(hSignal == nullptr) throw std::logic_error("Wrong histo name: " + hname + " has been provided");
	  if(mass.first == 200 || mass.first > 500){
		  model = "bukin";
	  }
	  else model =  "doublegausexp";//"doublegausexp";//
	  hSignal->Rebin(rebin);
	  ab::FitContainer fitter(hSignal,name,"signal");
	  min = hSignal->GetXaxis()->GetXmin();
	  max = hSignal->GetXaxis()->GetXmax();
	  fitter.fitRangeMin(min);
	  fitter.fitRangeMax(max);
	  fitter.verbosity(verbosity - 1);
	  fitter.setModel(ab::FitContainer::Type::signal,model);
	  RooWorkspace& wC = fitter.getWorkspace();
	  std::unique_ptr<RooFitResult> Signalfit = fitter.FitSignal(model);
	  TH1D* hSignal_240 = (TH1D*) f.Get(histo_name_240.c_str());
	  RooRealVar signal_norm("signal_norm","signal_norm",hSignal_240->Integral());
	  signal_norm.setConstant(true);
	  fitter.Import(signal_norm);
	  fitter.Write();
	  //Syst variation
	  for(const auto & syst : shape_unc){
		  std::cout<<"***************************************"<<std::endl;
		  std::cout<<"***************"<<syst<<"***************"<<std::endl;
		  std::cout<<"***************************************"<<std::endl;
		for(const auto & s : sign){
			  std::cout<<"***************************************"<<std::endl;
			  std::cout<<"***************"<<syst<<" : "<<s<<"***************"<<std::endl;
			  std::cout<<"***************************************"<<std::endl;

//			  hname = histo_name_240 + "_" + CMS + "_" + syst + "_" + energy + s;
			  hname = histo_name_240 + "_" + syst + s;
//			  hname = "templates/" + histo_name_240 + "_" + CMS + "_" + syst + "_VIS_" + energy + s;
			  TH1D* hSignal_syst = (TH1D*) f.Get(hname.c_str());
			  if(hSignal_syst == nullptr) throw std::logic_error("Wrong histo name: " + hname + " has been provided");
			  hSignal_syst->Rebin(rebin);
			  ab::FitContainer fit_syst(hSignal_syst,name + "/" + syst + "_" + s,"signal");
			  fit_syst.fitRangeMin(min);
			  fit_syst.fitRangeMax(max);
			  fit_syst.verbosity(verbosity - 1);
			  fit_syst.setModel(ab::FitContainer::Type::signal,model);
			  RooWorkspace& w = fit_syst.getWorkspace();
			  setParameters(w,wC);
          	  if(syst == "CMS_res_j_13TeV") {
          		  if(model == "bukin"){
          			  fixParameters("Xp",w,wC);
          			  fixParameters("rho1",w,wC);
          			  fixParameters("rho2",w,wC);
          			  fixParameters("xi",w,wC);
          		  }
          		  else {
          			fixParameters("tail_sigma",w,wC);
          			fixParameters("mean",w,wC);
          			fixParameters("tail_shift",w,wC);
          		  }
          	  }
          	  else if (syst == "CMS_scale_j_13TeV"){
          		  if(model == "bukin"){
          			  fixParameters("sigp",w,wC);
          			  fixParameters("rho1",w,wC);
          			  fixParameters("rho2",w,wC);
          			  fixParameters("xi",w,wC);
          		  }
          		  else {
          			fixParameters("tail_sigma",w,wC);
          			  fixParameters("sigmaL",w,wC);
          			  fixParameters("sigmaR",w,wC);
          			  fixShiftParameter(w,wC);
          		  }
          	  }
          	  else if (syst == "CMS_eff_pTonl_13TeV"){
          		  if(model == "bukin"){
          			fixParameters("sigp",w,wC);
          			fixParameters("Xp",w,wC);
          			fixParameters("rho2",w,wC);
          		  }
          		  else {
          			fixParameters("tail_sigma",w,wC);
          			fixParameters("mean",w,wC);
          			fixParameters("tail_shift",w,wC);
          		  }
          	  }
          	  else if (syst == "CMS_eff_b_13TeV"){
          		  if(model == "bukin"){
          			  fixParameters("Xp",w,wC);
          			fixParameters("sigp",w,wC);
          		  }
          		  else{
          			fixParameters("tail_shift",w,wC);
          			fixParameters("mean",w,wC);
          		  }
          	  }
               	std::unique_ptr<RooFitResult> Signalfit_syst = fit_syst.FitSignal(model);

               	//Write signal normalisation
               	TH1D* hSignal_240_syst = (TH1D*) f.Get((hname).c_str());
               	RooRealVar signal_norm("signal_norm","signal_norm",hSignal_240_syst->Integral());
               	signal_norm.setConstant(true);
               	fit_syst.Import(signal_norm);
               	fit_syst.Write();
	        }
		//Draw plots
		drawSystFits(mass.first,syst,name);
	  }
	  //Plot shape parameters
	  if(model == "bukin") plotShapeParameters(mass.first,shape_unc,name,{"Xp","sigp","xi","rho1","rho2"});
	  else plotShapeParameters(mass.first,shape_unc,name,{"mean","sigmaL","sigmaR","tail_shift","tail_sigma"});

/**/
  }

  return 0;

}

void drawSystFits(const int & mass,const std::string &syst, const std::string &output){
	TFile fCentral( (output + "/workspace/FitContainer_workspace.root").c_str(),"read");
	TCanvas c("c","c",1000,800);
	RooWorkspace& wCentral = *( (RooWorkspace*) fCentral.Get("workspace"));

	//Extract RooRealVar for X-axis
	RooRealVar& x = *( wCentral.var("mbb"));
	RooPlot& frame = *(x.frame());
	RooAbsPdf& PdfCentral = *(wCentral.pdf("signal"));
	RooAbsData& dataC = *wCentral.data("signal_container");
	dataC.plotOn(&frame,RooFit::MarkerSize(1.2),RooFit::MarkerStyle(20),RooFit::MarkerColor(kBlack),RooFit::LineColor(kBlack),RooFit::Name("Central"));
	PdfCentral.plotOn(&frame,RooFit::Name("central_curve"),RooFit::LineColor(kBlack),RooFit::Normalization(dataC.sumEntries("1","fit_range"),RooAbsReal::NumEvent));

	//down:
	TFile fDown( (output + "/" + syst + "_" + "Down/workspace/FitContainer_workspace.root").c_str());
	RooWorkspace& wDown = *( (RooWorkspace*) fDown.Get("workspace"));
	RooAbsPdf& PdfDown = *(wDown.pdf("signal"));
	RooAbsData& dataDown = *wDown.data("signal_container");
	dataDown.plotOn(&frame,RooFit::MarkerSize(1.2),RooFit::MarkerStyle(20),RooFit::MarkerColor(kRed),RooFit::LineColor(kRed),RooFit::Name("Down"));
	PdfDown.plotOn(&frame,RooFit::LineColor(kRed),RooFit::Name("Down_curve"),RooFit::Normalization(dataDown.sumEntries("1","fit_range"),RooAbsReal::NumEvent));

	//up
	TFile fUp( (output + "/" + syst + "_" + "Up/workspace/FitContainer_workspace.root").c_str());
	RooWorkspace& wUp = *( (RooWorkspace*) fUp.Get("workspace"));
	RooAbsPdf& PdfUp = *(wUp.pdf("signal"));
	RooAbsData& dataUp = *wUp.data("signal_container");
	dataUp.plotOn(&frame,RooFit::MarkerSize(1.2),RooFit::MarkerStyle(20),RooFit::MarkerColor(kBlue),RooFit::LineColor(kBlue),RooFit::Name("Up"));
	PdfUp.plotOn(&frame,RooFit::LineColor(kBlue),RooFit::Name("Up_curve"),RooFit::Normalization(dataUp.sumEntries("1","fit_range"),RooAbsReal::NumEvent));

//	TPad pad1("pad1","",0.,0.35,1,1);
//	pad1.SetBottomMargin(0.0);
//	pad1.SetLeftMargin(0.16);
//	pad1.SetRightMargin(0.05);
//	pad1.Draw();
//	pad1.cd();
	double xmin = 0.6,ymin = 0.6,xmax = 0.9,ymax = 0.85;
	if(mass > 900 || mass == 500) {
		xmin = 0.25; ymin = 0.6 ; xmax = 0.45 ; ymax = 0.85;
	}
	TLegend leg(xmin,ymin,xmax,ymax);
	leg.SetBorderSize(0);
	leg.SetLineColor(0);
	frame.SetTitle("");
	frame.Draw();
	leg.AddEntry(frame.findObject("Central"),"Central","pl");
	leg.AddEntry(frame.findObject("Down"),("-2#sigma variation " + syst).c_str(),"pl");
	leg.AddEntry(frame.findObject("Up"),("+2#sigma variation " + syst).c_str(),"pl");
	leg.Draw();

	style.drawStandardTitle();

//	c.cd();
//	TPad pad2("pad2","",0.,0.21,1,0.35);
//	pad2.SetBottomMargin(0.0);
//	pad2.SetTopMargin(0.);
//	pad2.SetLeftMargin(pad1.GetLeftMargin());
//	pad2.SetRightMargin(pad1.GetRightMargin());
//	pad2.Draw();
//	pad2.cd();

	c.Print( (output + "/plots/" + syst + ".pdf").c_str());
}

void plotShapeParameters(const int & mass,const vector<string> &syst, const std::string &folder, const vector<string>& parameters){
	TFile fCentral( (folder + "/workspace/FitContainer_workspace.root").c_str(),"read");
	RooWorkspace& wCn = *( (RooWorkspace*) fCentral.Get("workspace"));

	for(const auto& s: syst){
		//down
		TFile fDown( (folder + "/" + s + "_" + "Down/workspace/FitContainer_workspace.root").c_str());
		RooWorkspace& wDn = *( (RooWorkspace*) fDown.Get("workspace"));
		//up
		TFile fUp( (folder + "/" + s + "_" + "Up/workspace/FitContainer_workspace.root").c_str());
		RooWorkspace& wUp = *( (RooWorkspace*) fUp.Get("workspace"));

		TCanvas c(("c" + s).c_str(),"c",1000,800);
		TH1F hUp("hUp",(";Shape parameters;Pulls, (#theta_{up/dn} - #theta_{0})/#sigma_{#theta_{0}}"),parameters.size(),0,10);
		TH1F hDn("hDn","hDn",parameters.size(),0,10);
		hDn.SetCanExtend(TH1::kAllAxes);
		hUp.SetCanExtend(TH1::kAllAxes);
		hUp.SetStats(0);
		hUp.SetMarkerStyle(20);
		hDn.SetMarkerColor(kRed);
		hDn.SetLineColor(kRed);
		hDn.SetMarkerStyle(20);
		hUp.SetMarkerColor(kBlue);
		hUp.SetLineColor(kBlue);
		hUp.GetYaxis()->SetRangeUser(-5.,5.);
		int i = 0;
		double e_up;
		double e_down;
		for(const auto& p : parameters){
			++i;
			if(wCn.var(p.c_str())->getError() == 0) continue;
			double up = (wUp.var(p.c_str())->getValV() - wCn.var(p.c_str())->getValV()) / wCn.var(p.c_str())->getError();
			if(up != 0) {
				e_up =  (wCn.var(p.c_str())->getError() + wUp.var(p.c_str())->getError())/2 / wCn.var(p.c_str())->getError() ;
			} else e_up = 0;
			double down = (wDn.var(p.c_str())->getValV() - wCn.var(p.c_str())->getValV()) /wCn.var(p.c_str())->getError();
			if(down != 0) {
				e_down = (wCn.var(p.c_str())->getError() + wDn.var(p.c_str())->getError())/2 / wCn.var(p.c_str())->getError();
			} else e_down = 0;
			hUp.GetXaxis()->SetBinLabel(i,p.c_str());
			hUp.SetBinContent(i,up);
			hUp.SetBinError(i,e_up);
			hDn.SetBinContent(i,down);
			hDn.SetBinError(i,e_down);
		}
		hUp.Draw("E");
		hDn.Draw("Esame");
		double xmin = 0.65,ymin = 0.7,xmax = 0.85,ymax = 0.85;
		TLegend leg(xmin,ymin,xmax,ymax);
		leg.SetBorderSize(0);
		leg.SetLineColor(0);
		leg.AddEntry(&hUp,(s + " Up").c_str(),"pl");
		leg.AddEntry(&hDn,(s + " Down").c_str(),"pl");
		leg.Draw();
		style.drawStandardTitle();

		c.Print( (folder + "/plots/delta" + s + ".pdf").c_str() );
	}

}

void fixParameters(const string& name, RooWorkspace& w, const RooWorkspace& central){
	double val = central.var(name.c_str())->getValV();
	w.var(name.c_str())->setVal(val);
	w.var(name.c_str())->setConstant();
}

void fixParameters(const string& name, RooWorkspace& w, const double& val){
	w.var(name.c_str())->setVal(val);
	w.var(name.c_str())->setConstant();
}

void setParameters(RooWorkspace& w, const RooWorkspace& central){
	RooArgSet vars = central.allVars();
	auto iter = vars.createIterator();
	RooRealVar *v;
	while ((v = (RooRealVar*) iter->Next()) ){
		if(string(v->GetName()) == "signal_norm") continue;
		double val = central.var(v->GetName())->getValV();
		w.var(v->GetName())->setVal(val);
	}
}

void fixShiftParameter(RooWorkspace& w, const RooWorkspace& central){
	RooRealVar& muC = *(RooRealVar*) central.var("mean");
	RooRealVar& shiftC = *(RooRealVar*) central.var("tail_shift");
	const double d = shiftC.getValV() - muC.getValV();
	w.removeSet("tail_shift");
	RooRealVar& mu = *(RooRealVar*) w.var("mean");
	RooRealVar delta("delta","delta",d,"GeV");
	delta.setConstant();
	RooFormulaVar tail_shift("tail_shift","tail_shift","@0 + @1",RooArgSet(mu,delta));
	w.import(tail_shift);
}
