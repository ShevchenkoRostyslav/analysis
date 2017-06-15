#include <string>

#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "RooFit.h"

#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;
namespace ab = analysis::backgroundmodel;

void reweight(TH1& h);
double weight_function(const double& x);

int main() {
	/*
	 * Script to fit bbnb data (weighted and not)
	 */
	HbbStyle style;
	style.set(PRIVATE);

	const string cmsswBase = gSystem->Getenv("CMSSW_BASE");
	const string selection = "TripleBTagReverseSelectionBtoH2016_prescale_13TeV_G4.root";
	bool rew = false;

	//Preferences for the fit
//	int rebin = 1;
	const int verbosity = 1;
	int nbins = 45;
	double min = 200;
	double max = 650;
	string model = "extnovoeffprod";//"superdijeteffprod,2";//
	string output_name = cmsswBase + "/src/Analysis/BackgroundModel/test/test_FitRewBBnB/" + model + "_fit/";

	TFile f((cmsswBase + "/src/Analysis/BackgroundModel/data/2016DataRereco_05_01_2017/" + selection).c_str(),"read");
	auto &tree = *GetFromTFile<TTree>(f,"MssmHbb_13TeV");
	TH1D *histo = new TH1D("histo","histo",nbins,min,max);
	tree.Draw("mbb>>histo");
//	histo->Rebin(rebin);
	if(rew) {
		reweight(*histo);
		output_name += "Rew";
	}

	//Setup fitter
	ab::FitContainer fitter(histo,output_name,"background");
	fitter.fitRangeMax(max);
	fitter.fitRangeMin(min);
	fitter.verbosity(verbosity - 1);
	fitter.setNBins(nbins);
	fitter.setModel(ab::FitContainer::Type::background,model);
	//Set Initial values
	fitter.getWorkspace().var("slope_novoeff")->setVal(0.0159782); // 0.0159765
	fitter.getWorkspace().var("turnon_novoeff")->setVal(222.813); // 223.027
	fitter.getWorkspace().var("peak")->setVal(268.225);//350
	fitter.getWorkspace().var("width")->setVal(63.3788);
	fitter.getWorkspace().var("tail")->setVal(-0.447615);
	fitter.getWorkspace().var("par4")->setVal(-0.000641769);
//
//	fitter.getWorkspace().var("tail")->setRange(-5,5);
//	fitter.getWorkspace().var("turnon_novoeff")->setRange(200,250);
//
//	if(rew) {
//		fitter.getWorkspace().var("slope_novoeff")->setConstant();
//		fitter.getWorkspace().var("turnon_novoeff")->setConstant();
//	}

	fitter.getWorkspace().var("peak")->setConstant();
	fitter.getWorkspace().var("width")->setConstant();
	fitter.getWorkspace().var("par4")->setConstant();
	fitter.getWorkspace().var("tail")->setConstant();

	//Fit
	std::unique_ptr<RooFitResult> fit = fitter.backgroundOnlyFit(model);
}

void reweight(TH1& h){
	/*
	 * Method to perform bin-by-bin reweighting
	 */
	TRandom3 r;
	for(int i = 1;i <= h.GetNbinsX(); ++i){
		if(h.GetBinCenter(i) > 250) continue;
		double weight = weight_function(h.GetXaxis()->GetBinCenter(i));
		double n_med = h.GetBinContent(i) * weight;
		double n_rew = r.Poisson(n_med);
		h.SetBinContent(i,n_rew);
		h.SetBinError(i,sqrt(n_rew));
	}
}

double weight_function(const double& x){
	/*
	 * weight function accroding to the ratio bbb to bbnb
	 */
	double weight = 1;
	if(x < 400) weight= 1./(1.+exp( -0.0540249 * (x - 188.762) ));

	return weight;
}
