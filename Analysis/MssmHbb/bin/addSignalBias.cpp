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

namespace fs = boost::filesystem;
using namespace std;

void addSignalBias(const string& signal_wsp);

int main(int argc, char ** argv){
	//list of the inputs:
	const string cmsswBase = getenv("CMSSW_BASE");
	map<int,string> input_signal;
//	input_signal[300] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-300/workspace/signal_workspace.root";
//	input_signal[350] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-350/workspace/signal_workspace.root";
//	input_signal[400] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-400/workspace/signal_workspace.root";
//	input_signal[500] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-500/workspace/signal_workspace_500SR2.root";
//	input_signal[600] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-600/workspace/signal_workspace_600SR1.root";
//	input_signal[700] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-700/workspace/signal_workspace.root";
//	input_signal[900] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-900/workspace/signal_workspace.root";
	input_signal[1100] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1100/workspace/signal_workspace_1100SR2.root";
//	input_signal[1300] = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1300/workspace/signal_workspace.root";

	for(const auto& mass : input_signal){
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
	string out_path = signal_dir + "/signal_bias_workspace_1100SR2.root";

	wOut.writeToFile( out_path.c_str());
	wOut.Print("v");

}
