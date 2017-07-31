#include "boost/program_options.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <memory>

#include <vector>

//Root includes
#include "TFile.h"
#include "TSystem.h"

//RooFit includes
#include "RooWorkspace.h"
#include "RooFormulaVar.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooCustomizer.h"

#include "boost/filesystem.hpp"

#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

namespace fs = boost::filesystem;
using namespace std;

void addSignalBias(const string& signal_wsp);

int main(int argc, char ** argv){
	//list of the inputs:
	const string cmsswBase = getenv("CMSSW_BASE");

	for(const auto& mass : mssmhbb::signal_workspaces){
		addSignalBias(mass.second);
	}

	return 0;
}

void addSignalBias(const string& signal_wsp){
	/*
	 * Functin to add bias to the signal workspace
	 */
	if(gSystem->AccessPathName(signal_wsp.c_str()) ) throw invalid_argument("Error: no signal file " + signal_wsp);
	TFile f(signal_wsp.c_str(),"read");
	RooWorkspace wOut("workspace","RECREATE");
	//Get signal workspace:
	RooWorkspace& w = (RooWorkspace&) *f.Get("workspace");
	RooAbsReal& sg_norm = *w.function("signal_norm");
	//create var for bias
	RooRealVar rBias("bias","bias",0,-10.,10);
	//Ceate Bias nomalisation parameter
	RooFormulaVar fBias("signal_bias_norm","","bias*signal_norm",RooArgList(rBias,sg_norm));

	wOut.import(*w.pdf("signal"));
	wOut.pdf("signal")->SetName("signal_bias");
	wOut.import(fBias);

	//get parent directory for signal file
	fs::path p(signal_wsp.c_str());
	string signal_dir = p.parent_path().string();
	string out_path = signal_dir + "/signal_bias_workspace.root";

	wOut.writeToFile( out_path.c_str());
	wOut.Print("v");

}
