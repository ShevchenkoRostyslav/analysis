#include <string>

#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "RooFit.h"

#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"

using namespace std;
namespace ab = analysis::backgroundmodel;
void PerformBackgroundQCDFit(const string& input_file, const string& name);
TH1D * GetSepcificRangeHistp(TH1* h, const double& xmin, const double& xmax, const int& nbins);

int main() {
	/*
	 * Script to run background fit for the
	 * QCD MC simulation
	 */
	string bbbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/DataMC_3b_ReReco_35673fb_lowM_QCD.root";
	string bbnbFilePath 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/bbx_ReReco_35673fb_lowM_QCD.root";

	PerformBackgroundQCDFit(bbbFilePath,"bbb");
	PerformBackgroundQCDFit(bbnbFilePath,"bbnb");
}

void PerformBackgroundQCDFit(const string& input_file, const string& name){

	HbbStyle style;
	style.setTDRstyle(PRIVATE);

	string output_name 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/BackgroundFit_QCD_" + name + "/";

	const int verbosity = 1;
	int nbins = 45;
	double min = 200;
	double max = 650;
	string model = "extnovoeffprod";

	TFile f(input_file.c_str(),"READ");
	auto h = GetFromTFile<TH1D>(f,"general/diJet_m");
	auto histo = GetSepcificRangeHistp(h,min,max,nbins);

	ab::FitContainer fitter(histo,output_name,"background");
	fitter.setBinsToPlot(nbins);
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

	fitter.getWorkspace().var("tail")->setRange(-10,10);

	std::unique_ptr<RooFitResult> fit = fitter.backgroundOnlyFit(model);
}

TH1D * GetSepcificRangeHistp(TH1* h, const double& xmin, const double& xmax, const int& nbins){
	TH1D *histo = new TH1D("histo",h->GetTitle(),nbins,xmin,xmax);
	double bin_width = (xmax - xmin) / nbins;
	int j = 1;
	for(int i = 1; i <= h->GetNbinsX();++i){
		auto xc = h->GetBinCenter(i);
		auto data = h->GetBinContent(i);
		auto edata = h->GetBinError(i);
		if((xc >= (xmin + bin_width / 2.)) && (xc <= xmax - bin_width/2. )){
			histo->SetBinContent(j,data);
			histo->SetBinError(j,edata);
			++j;
		}
	}
	return histo;
}
