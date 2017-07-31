#include <string>

#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "RooFit.h"

#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"

using namespace std;
namespace ab = analysis::backgroundmodel;

void MixWithSample(TH1 &h1,const TH1& h2);
void prescale(TH1& h, const double& val = 9);
void reweight(TH1& h);
double weight_function(const double& x);

int main(int argc, char **argv) {
	/*
	 * Script to fit bbnb data (weighted and not)
	 */
	HbbStyle style;
	style.set(PRIVATE);

	const string cmsswBase = gSystem->Getenv("CMSSW_BASE");
	const string selection = "DataMC_3b_ReReco_35673fb_lowM_bEnriched_Tot.root";//"DataMC_3b_ReReco_35673fb_lowM_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root";//"bbx_ReReco_35673fb_lowM_Run2016.root";//
	const string selection2= "DataMC_3b_ReReco_35673fb_lowM_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root";
	string output_name = "../test/test_3b_bEnr_bg_fit";//"../test/test_tt_bg_fit";//
	bool rew = false, mix = true;
	TFile f( (cmsswBase + "/src/Analysis/MssmHbb/output/" + selection).c_str(),"read");
	TFile f2((cmsswBase + "/src/Analysis/MssmHbb/output/" + selection2).c_str(),"read");
	string histo_name = "general/diJet_m";//"diJet_b_m/diJet_m";//

	//Preferences for the fit
	int rebin = 1;
	const int verbosity = 1;
	double min = 200;
	double max = 650;//1100;//
	string model = "extnovoeffprod";

	TH1D *h = (TH1D*) f.Get(histo_name.c_str());
	TH1D *h2= (TH1D*) f2.Get(histo_name.c_str());
	h->Rebin(rebin);
	h2->Rebin(rebin);
	//Add h1 and h2 to make a mixture
	if (mix) {
		//Show mixture content
		int h1_bin_low = h->FindBin(min);
		int h1_bin_high= h->FindBin(max);
		int h2_bin_low = h2->FindBin(min);
		int h2_bin_high= h2->FindBin(max);
		double integr_h1 = h->Integral(h1_bin_low,h1_bin_high);
		double integr_h2 = h2->Integral(h2_bin_low,h2_bin_high);
		std::cout<<"Integrals for ["<<min<<";"<<max<<"]: h1 = "<<integr_h1<<" h2 = "<<integr_h2<<" sum = "<<integr_h1 + integr_h2<<endl;
		std::cout<<"Percents in sum: h1 = "<<integr_h1/(integr_h1 + integr_h2)<<" h2 = "<<integr_h2/(integr_h1 + integr_h2)<<endl;
		MixWithSample(*h, *h2);
		output_name += "_MixWithTT";
	}
//	prescale(*h);
	if(rew) {
//		reweight(*h);
		output_name += "_Rew";
	}

	//Setup fitter
	ab::FitContainer fitter(h,output_name,"background");
	fitter.fitRangeMax(max);
	fitter.fitRangeMin(min);
	fitter.verbosity(verbosity - 1);
	fitter.setModel(ab::FitContainer::Type::background,model);

	//Set Initial values
	fitter.getWorkspace().var("slope_novoeff")->setVal(0.0151088);
	fitter.getWorkspace().var("turnon_novoeff")->setVal(227.188);
	fitter.getWorkspace().var("peak")->setVal(258.566); //268.093
	fitter.getWorkspace().var("width")->setVal(57.2917); //63.3651
	fitter.getWorkspace().var("tail")->setVal(-0.219176); //-0.448627
	fitter.getWorkspace().var("par4")->setVal(-0.00118348); //-0.00063946
//
	fitter.getWorkspace().var("tail")->setRange(-2,2);
	fitter.getWorkspace().var("width")->setRange(0,100);
//	fitter.getWorkspace().var("turnon_novoeff")->setRange(200,250);

//	if(rew) {
//		fitter.getWorkspace().var("slope_novoeff")->setConstant();
//		fitter.getWorkspace().var("turnon_novoeff")->setConstant();
//	}
//	fitter.getWorkspace().var("peak")->setConstant();
//	fitter.getWorkspace().var("width")->setConstant();
//	fitter.getWorkspace().var("par4")->setConstant();
//	fitter.getWorkspace().var("tail")->setConstant();


	//Fit
	std::unique_ptr<RooFitResult> fit = fitter.backgroundOnlyFit(model,1);
	std::cout<<"Fit status: "<<fit->status()<<std::endl;
}

void prescale(TH1& h, const double& val){
//	double prob = 1. / val;
	TRandom3 r;
	for(int i = 1;i <= h.GetNbinsX(); ++i){
		int Ntot = h.GetBinContent(i);
//		int Npresc = Ntot / val;
//		double binomial_prob = TMath::Factorial(Ntot) / (TMath::Factorial(Npresc) * Tmath::Factorial(Ntot - Npresc)) * pow(prob,Npresc) * pow(1-prob,Ntot - Npresc);
		double binomial_prob = 0;
		for(int i = 1; i <= Ntot; ++i){
			if(r.Rndm() < 1./9.) binomial_prob+=1;
//			binomial_prob += TMath::BinomialI(prob,Ntot,Npresc);
//			binomial_prob += TMath::BinomialI(prob,9,1);
//			binomial_prob += TMath::Factorial(9) / (TMath::Factorial(1) * TMath::Factorial(9 - 1)) * pow(prob,1) * pow(1-prob,9 - 1);
		}
		int Nfinale = binomial_prob;
		h.SetBinContent(i,Nfinale);
		h.SetBinError(i,sqrt(Nfinale));
//		h.SetBinContent(i,(int) h.GetBinContent(i)/ val);
//		h.SetBinError(i,h.GetBinError(i) / sqrt(val));
	}
}

void reweight(TH1& h){
	/*
	 * Method to perform bin-by-bin reweighting
	 */
	for(int i = 1;i <= h.GetNbinsX(); ++i){
		double weight = weight_function(h.GetXaxis()->GetBinCenter(i));
		h.SetBinContent(i,weight * h.GetBinContent(i));
		h.SetBinError(i,h.GetBinError(i) * weight);
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

void MixWithSample(TH1& h1, const TH1& h2){
	/*
	 * Function to add h1 and h2
	 */
	h1.Add(&h2,1);
}
