#include <string>
#include <map>
#include <vector>

/*
 * Namespace that contains all basic variables
 * used in the analysis flow
 */

namespace mssmhbb{
	// current CMSSW release base dir
	static const std::string cmsswBase = getenv("CMSSW_BASE");
	// current signal selection
	static const std::string signal_selection = "ReReco_35673fb";
	// vector of masses:
	static const std::vector<int> masses = {300,350,400,500,600,700,900,1100,1300};
	// signal MC points and pathes
	static const std::map<int,std::string> signal_templates {
		{300,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.root"},
		{350,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.root"},
		{400,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.root"},
		{500,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root"},
		{600,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.root"},
		{700,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.root"},
		{900,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.root"},
		{1100,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root"},
		{1300,cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + signal_selection + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.root"},
	};
	// signal MC PDFs
	static const std::map<int,std::string> signal_workspaces {
		{300,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-300/workspace/signal_workspace.root"},
		{350,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-350/workspace/signal_workspace.root"},
		{400,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-400/workspace/signal_workspace.root"},
		{500,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-500/workspace/signal_workspace.root"},
		{600,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-600/workspace/signal_workspace.root"},
		{700,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-700/workspace/signal_workspace.root"},
		{900,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-900/workspace/signal_workspace.root"},
		{1100,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1100/workspace/signal_workspace.root"},
		{1300,cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1300/workspace/signal_workspace.root"},
	};
	// path to mssm_xs_tool root file
	static const std::string path_mssm_xsections = cmsswBase + "/src/Analysis/MssmHbb/macros/signal/mhmodp_mu200_13TeV.root";

	//path to bg fits
	static const std::string path_bg_sr1 = cmsswBase + "/src/Analysis/BackgroundModel/test/extnovoeffprod_200to650_10GeV_G4/workspace/FitContainer_workspace.root";
	static const std::string path_bg_sr2 = cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_350to1190_20GeV_G4/workspace/FitContainer_workspace.root";
	static const std::string path_bg_sr3 = cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_500to1700_25GeV_G4/workspace/FitContainer_workspace.root";
	static const std::map<int,std::string> path_bg_fits {
		{1,path_bg_sr1}, //sub-range 1
		{2,path_bg_sr2}, //sub-range 2
		{3,path_bg_sr3}, //sub-range 2
	};

};
