#include "boost/program_options.hpp"

#include <exception>
#include <iostream>
#include <string>
#include <memory>

#include <vector>

//Root includes
#include "TH1.h"
#include "TFile.h"
#include "TSystem.h"
#include "TDirectory.h"
#include <TKey.h>
#include "TClass.h"
#include <TROOT.h>

//RooFit includes
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooDataHist.h"
#include "RooBukinPdf.h"
#include "RooEffProd.h"
#include "RooNovosibirsk.h"
#include "RooDataSet.h"

#include "Analysis/BackgroundModel/src/classes.h"
#include "Analysis/BackgroundModel/interface/RooDoubleGausExp.h"
#include "Analysis/BackgroundModel/interface/RooQuadGausExp.h"

#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;
  const size_t ERROR_BG_WORKSPACE = 3;
  const size_t ERROR_SIGNAL_WORKSPACE = 4;
} // namespace

using namespace std;
using namespace analysis;
using namespace boost::program_options;

//defin >> stream operator to be used in boost::po
namespace std {
	istream& operator>>(istream& in, pair<string,string>& ss) {
	string s;
  	  in >> s;
  	  const size_t sep = s.find(':');
  	  if(sep ==string::npos) s.find(',');
  	  if (sep==string::npos) {
	  ss.first = string();
    	ss.second = s;
  	  } else {
	  ss.first  = s.substr(0,sep);
    	ss.second = s.substr(sep+1);
  	  }
  	  return in;
	}
}

typedef unique_ptr<TFile> pTFile;
typedef map<string, map<string, pair<string,string> > > myMap;

void setup_bg(const string& in_path,const string& out_path, TH1& data_obs);
void setup_bg(const string& in_path,const string& out_path, const bool& generate_asimov = true);
void setup_signal(const string& in_path, const vector<string>& syst, const string& out_suffix = "");
string SignalModel(const vector<string>& parameters);
string GetDirName(const string& filename);
void setMbb(RooRealVar& var, const string& path);

int main(int argc, char ** argv){
	try {
		//Additional suffix for the ouptut signal workspace name
		string name_suffix = "";

		//list of the inputs:
		vector<string> input_signal = {
//				mssmhbb::signal_folders.at(300),
//				mssmhbb::signal_folders.at(350),
//				mssmhbb::signal_folders.at(400),
//				mssmhbb::signal_folders.at(500),
//				mssmhbb::signal_folders.at(600),
//				mssmhbb::signal_folders.at(700),
//				mssmhbb::signal_folders.at(900),
				mssmhbb::signal_folders.at(1100),
//				mssmhbb::signal_folders.at(1300)
								};
		vector<string> input_background = {
				mssmhbb::path_bg_sr1,
				mssmhbb::path_bg_sr2,
				mssmhbb::path_bg_sr3,
		};
		vector<string> data_obs = {
				mssmhbb::path_data_sr1,
				mssmhbb::path_data_sr2,
				mssmhbb::path_data_sr3,
		};
		vector<string> systematics = {"CMS_scale_j_13TeV","CMS_res_j_13TeV","CMS_eff_pTonl_13TeV","CMS_eff_b_13TeV"};

		options_description desc("Options");
		desc.add_options()
				("help,h", "Print help messages")
				("background,b" , value<vector<string> >(&input_background)	->multitoken(), "Input background workspaces")
				("signal,s"     , value<vector<string> >(&input_signal)    	->multitoken(),"Input std::vector<signal folder>" )
				("data_obs,d"   , value<vector<string> >(&data_obs)			->multitoken(),"Input file with observed data")
				("syst_errors,e", value<vector<string> >(&systematics)		->multitoken(),"std::vector<systematic error>" )
				;

		//Boost container for input arguments.
		variables_map vm;
		try {
			store(parse_command_line(argc,argv,desc), vm);

			// --help option
			if(vm.count("help")){
				cout<<"Basic Command Line App to prepare signal/bg workspaces for parametric combine tool fitting"<<endl;
				cout<<desc<<endl;
				return SUCCESS;
			}

			notify(vm);
		} catch (error& e) {
			//Command line parcing can throw exception
			// should be catched
			cerr<<"ERROR: "<<e.what()<< endl << endl;
			cerr<< desc;
			return ERROR_IN_COMMAND_LINE;
		}

		//----------Main procedure----------
		//Adjust bg workspace
		for (unsigned int i = 0; i != data_obs.size(); ++i ){
			TFile temp(data_obs.at(i).c_str(),"read");
			//Read list of keys
//			TIter next(temp.GetListOfKeys());
//			TKey *key;
//			while ((key = (TKey*)next())) {
//				TClass *cl = gROOT->GetClass(key->GetClassName());
//				if (!cl->InheritsFrom("RooDataSet")) continue;
				try {
//					setup_bg(input_background.at(i),GetDirName(data_obs.at(i)),false);
					setup_bg(input_background.at(i),GetDirName(data_obs.at(i)),*((TH1D*) temp.Get("data_obs")));
				} catch (exception& e) {
					cerr<<e.what()<<endl;
					return ERROR_BG_WORKSPACE;
				}

//			}
		}
		//Create signal workspace
		for(const auto& in : input_signal){
			try {
				setup_signal(in /* + "_NormToTanB30" */,systematics,name_suffix);
			} catch (exception& e) {
				cerr<<e.what()<<endl;
				return ERROR_SIGNAL_WORKSPACE;
			}
		}

	} catch (exception& e) {
		cerr<<"Unhandled Exception reached the top of main: "<<e.what()<<", application will now exit";
		return ERROR_UNHANDLED_EXCEPTION;
	}
	return SUCCESS;
}

void setup_bg(const string& in_path,const string& out_path, const bool& generate_asimov){
	if(gSystem->AccessPathName(in_path.c_str()) ) throw invalid_argument("Error: no bg file " + in_path);
	TFile fIn(in_path.c_str(),"read");
	// Number of bins in daat_obs:
	int nbins = 5000;
	//Input workspace
	RooWorkspace& w = *( (RooWorkspace*)fIn.Get("workspace") );
	//Input RooDataSet:
	RooAbsData& roo_data = *w.data("data_obs");
	RooRealVar& mbb = *w.var("mbb");
	mbb.setBins(nbins);

	//Prepare Bg normalisation variable!!!
	string bg_pdf_name = "background";
	if(w.pdf(bg_pdf_name.c_str()) == nullptr) throw invalid_argument("Error: no <background> pdf has been found in bg workspace");
	RooRealVar bg_norm((bg_pdf_name+"_norm").c_str(),"background_norm",roo_data.sumEntries());
	bg_norm.setConstant();

	std::string add_name = "";
	//RooDataHist with data_obs
	TH1 *h;
	if(generate_asimov){
		h = w.pdf(bg_pdf_name.c_str())->createHistogram("QCD",mbb,RooFit::Binning(nbins,mbb.getMin(),mbb.getMax()));
	}
	else{
		h = roo_data.createHistogram("QCD", mbb, RooFit::Binning(nbins,mbb.getMin(),mbb.getMax()));
		add_name = "_RealBBnB";
	}

	h->Scale(roo_data.sumEntries()/h->Integral());
	h->SetTitle("data_obs");
	h->SetName("data_obs");
	RooDataHist h_data("data_obs","data_obs",mbb,RooFit::Import(*h));

	//Output workspace
	RooWorkspace wOut("workspace");
	wOut.import(*(RooAbsPdf*)w.pdf(bg_pdf_name.c_str()));
	wOut.import(h_data);
	wOut.import(bg_norm);

	/*
	 * Trying to fix Turn-on for the first sub-range
	 */
	if(out_path.find("sr1") != std::string::npos){
		RooRealVar sl("slope_novoeff","slope_novoeff",0);
		RooRealVar tn("turnon_novoeff","turnon_novoeff",0);
		wOut.import(sl);
		wOut.import(tn);
		wOut.var("slope_novoeff")->setConstant();
		wOut.var("turnon_novoeff")->setConstant();
		add_name += "_TurnOnFix";
	}
	if(out_path.find("sr2") != std::string::npos) wOut.var("tail1")->setRange(-10,10);

	add_name += "_" + std::to_string(nbins) + "bins";
	wOut.Print("v");
	wOut.writeToFile((out_path + "/background_workspace" + add_name + ".root").c_str());

}

void setup_bg(const string& in_path,const string& out_path, TH1& data_obs){
	if(gSystem->AccessPathName(in_path.c_str()) ) throw invalid_argument("Error: no bg file " + in_path);
	//Read input file
	TFile fIn(in_path.c_str(),"read");
	// Number of bins in data_obs:
	int nbins = data_obs.GetNbinsX();
	if(mssmhbb::blinded) nbins = 5000;
	//Input workspace
	auto& w = *GetFromTFile<RooWorkspace>(fIn, "workspace");
	//Prepare RooDataHist for data_obs;
	RooRealVar mbb("mbb","mbb",data_obs.GetXaxis()->GetXmin(),data_obs.GetXaxis()->GetXmax());
	mbb.setBins(nbins);

	//Prepare Bg normalisation variable
	string bg_pdf_name = "background";
	//Input pdf
	auto &pdf = *GetFromRooWorkspace<RooAbsPdf>(w,bg_pdf_name);
	RooRealVar bg_norm((bg_pdf_name+"_norm").c_str(),"background_norm",data_obs.Integral());
	bg_norm.setConstant();

	//Output workspace
	RooWorkspace wOut("workspace");
	wOut.import(pdf);

	/*
	 * Prepare data_obs to be included to the finale backgorund file
	 * Asymov data in case of a blinded analysis and real data_obs
	 * when analysis is unblinded
	 */
	TH1* h = &data_obs;
	if(mssmhbb::blinded){
		//Generate Toy data_obs
		h = w.pdf(bg_pdf_name.c_str())->createHistogram("QCD",mbb,RooFit::Binning(nbins,data_obs.GetXaxis()->GetXmin(),data_obs.GetXaxis()->GetXmax()));
		h->Scale(data_obs.Integral() / h->Integral());
	}
	h->SetTitle("data_obs");
	h->SetName("data_obs");
	RooDataHist h_data("data_obs","data_obs",mbb,RooFit::Import(*h));
	wOut.import(h_data);
	wOut.import(bg_norm);

	/*
	 * Trying to fix Turn-on for the first sub-range
	 */
	std::string add_name = "";
	if(out_path.find("sr1") != std::string::npos){
//		RooRealVar sl("slope_novoeff","slope_novoeff",0);
//		RooRealVar tn("turnon_novoeff","turnon_novoeff",0);
//		wOut.import(sl);
//		wOut.import(tn);
		GetFromRooWorkspace<RooRealVar>(wOut,"slope_novoeff")->setConstant();
		GetFromRooWorkspace<RooRealVar>(wOut,"turnon_novoeff")->setConstant();
		add_name = "_TurnOnFix";

		/*
		 * Relax BG parameters limits
		 */
//		wOut.var("par4")->setRange(-10,10);
//		wOut.var("peak")->setRange(0,1000);
//		wOut.var("tail")->setRange(-10,10);
//		wOut.var("width")->setRange(5,1000);
//		add_name += "_limitlessBG";
	}
	add_name += "_" + std::to_string(nbins) + "bins";
	if(out_path.find("sr2") != std::string::npos) GetFromRooWorkspace<RooRealVar>(wOut,"tail1")->setRange(-10,10);
//	add_name += "_inBins";
	wOut.Print("v");
	wOut.writeToFile((out_path + "/background_workspace" + add_name + ".root").c_str());
}

void setup_signal(const string& in_folder, const vector<string>& syst, const string& out_suffix){
	//Sanity check of input folder
	if(gSystem->AccessPathName(in_folder.c_str())) throw invalid_argument("Error in <prepParamWorkspace::setup_signal>: Wrong path: " + in_folder);
	//Sanity check of central WP
	string central = (in_folder + "/workspace/FitContainer_workspace.root");
	if(gSystem->AccessPathName(central.c_str())) throw invalid_argument("Error in <prepParamWorkspace::setup_signal>: file /workspace/FitContainer_workspace.root was not found in " + in_folder);
	//Read central workspace
	TFile fInCentral( (in_folder + "/workspace/FitContainer_workspace.root").c_str(),"read");
	RooWorkspace& wInCentral = *( (RooWorkspace*) fInCentral.Get("workspace") );
	RooRealVar& mbb = *( (RooRealVar*) wInCentral.var("mbb") );
	setMbb(mbb,in_folder);
	// Nbins:
//	int nbins = 100000;
//	mbb.setRange(240,1700);	//needed to define the function
//	mbb.setBins(nbins);
	//Output workspace
	RooWorkspace w("workspace");
	//Signal shape parameters list
	vector<string> parameters;
	vector<string> var_parameters; //parameters list that are not constant
	RooArgList vars = wInCentral.allVars();
	auto iter = vars.createIterator();
	RooRealVar *v;
	while ((v = (RooRealVar*) iter->Next()) ){
		if(string(v->GetName()) == string(mbb.GetName())) continue;
		parameters.push_back( string(v->GetName()) );
	}
	//Form strings for RooFormulaVars
	map<string,string> formula;
	RooArgList nuisance_list; // Needed to add nuisances
	for(const auto &par : parameters){
		RooRealVar& vCent = (RooRealVar&) *wInCentral.var(par.c_str());
		formula[par] = to_string(vCent.getValV());
	}
	//Loop over the system. sources
	for(const auto& s : syst){
		//Read Up and Down TFiles
		string dnName = (in_folder + "/" + s + "_Up/workspace/FitContainer_workspace.root");
		if(gSystem->AccessPathName(dnName.c_str()) ) throw invalid_argument("Error in <prepParamWorkspace::setup_signal>: No /" + s + "_Up/workspace/FitContainer_workspace.root in " + in_folder);
		TFile fInUp( (dnName).c_str(),"read" );
		string upName = in_folder + "/" + s + "_Down/workspace/FitContainer_workspace.root";
		if(gSystem->AccessPathName(upName.c_str()) ) throw invalid_argument("Error in <prepParamWorkspace::setup_signal> No /" + s + "_Down/workspace/FitContainer_workspace.root in " + in_folder);
		TFile fInDown( (upName).c_str(),"read");
		//Get workspaces:
		RooWorkspace& wInUp = *( (RooWorkspace*) fInUp.Get("workspace") );
		RooWorkspace& wInDown = *( (RooWorkspace*) fInDown.Get("workspace") );
		//RooRealVar for nuisance:
//		string nuisance_name = s;
		RooRealVar nuisance(s.c_str(),s.c_str(),0,-4,4);
		for(const auto &par : parameters){
			if(!wInUp.var(par.c_str())->hasError()  && par != "signal_norm" ) continue;
			RooRealVar& vUp   = (RooRealVar&) *wInUp.var(par.c_str());
			RooRealVar& vDown = (RooRealVar&) *wInDown.var(par.c_str());
			string f = "+ ("+ to_string(vUp.getValV()) + " - " + to_string(vDown.getValV()) +") * 0.25 * " + s + " ";
			formula[par] += f;
			var_parameters.push_back(par);
		}
		w.import(nuisance);
		nuisance_list.add(*w.var(nuisance.GetName()));
	}
	//Finaly make a RooFormulaVars:
	for(const auto &par : parameters){
		if(find(var_parameters.begin(),var_parameters.end(),par) != var_parameters.end()){	//If parameter is var = f(nuisances)
			RooFormulaVar fPar( (par).c_str(),(par).c_str(),formula[par].c_str(),nuisance_list);
			w.import(fPar);
		}
		else {
//			if(par == "signal_norm") continue;
			RooRealVar constant(("f_"+par).c_str(),("f_"+par).c_str(),stod(formula[par].c_str()));
			constant.setConstant();
			w.import(constant);
			RooFormulaVar fconstant(par.c_str(),par.c_str(),"@0",RooArgSet(constant));
			w.import(fconstant);
		}
	}
	/* I believe there should me more elegant way to determine a type
	 * of RooFit Pdf. But curently - this is the inly solution that I see
	 */
	//Create new RooQuadGausExp as a function of Nuisnace 1 : JES
	string model = SignalModel(parameters);
	unique_ptr<RooAbsPdf> func;
	if(model == "bukin"){
		func = make_unique<RooBukinPdf>("signal","signal",mbb,(RooFormulaVar&) *w.function("Xp"),(RooFormulaVar&) *w.function("sigp"),
												(RooFormulaVar&)*w.function("xi"),(RooFormulaVar&) *w.function("rho1"),
												(RooFormulaVar&)*w.function("rho2"));
	}
	else if (model == "doublegausexp"){
		func = make_unique<analysis::backgroundmodel::RooDoubleGausExp>("signal","signal",mbb,(RooFormulaVar&) *w.function("mean"),(RooFormulaVar&) *w.function("sigmaL"),
													(RooFormulaVar&) *w.function("sigmaR"),(RooFormulaVar&)*w.function("tail_shift"),
													(RooFormulaVar&) *w.function("tail_sigma"));
	}
	else if (model == "quadgausexp"){
		func = make_unique<analysis::backgroundmodel::RooQuadGausExp>("signal","signal",mbb,(RooFormulaVar&) *w.function("mean"),(RooFormulaVar&) *w.function("sigmaL1"),
				(RooFormulaVar&) *w.function("sigmaL2"),
				(RooFormulaVar&) *w.function("sigmaR1"),
				(RooFormulaVar&) *w.function("sigmaR2"),
				(RooFormulaVar&)*w.function("tail_shift"),
				(RooFormulaVar&) *w.function("tail_sigma"),
				(RooFormulaVar&) *w.function("norm_g1"),
				(RooFormulaVar&) *w.function("norm_g2"));
	}
	else throw logic_error("ERROR");
	w.import(*func.get());
	string out_path = in_folder + "/workspace/signal_workspace" + out_suffix + ".root";
	std::cout<<"***************************************"<<std::endl;
	std::cout<<"*************** M = "<<out_path<<"***************"<<std::endl;
	std::cout<<"***************************************"<<std::endl;
	w.Print("v");
	w.writeToFile(out_path.c_str(),1);
	cout<<"INFO: file " + out_path +" was created"<<endl;
}

string SignalModel(const vector<string>& parameters){
	map<string, vector<string> > models;
	vector<string> inserter;
	models["bukin"] = {"signal_norm","Xp","sigp","xi","rho1","rho2"};
	models["doublegausexp"] = {"signal_norm","mean","sigmaL","sigmaR","tail_shift","tail_sigma"};
	models["quadgausexp"] = {"signal_norm","mean","sigmaL1","sigmaL2","sigmaR1","sigmaR2","tail_shift","tail_sigma","norm_g1","norm_g2"};
	set<string> s1(parameters.begin(),parameters.end());
	for(const auto& model : models){
		inserter.clear();
		set<string> s2(model.second.begin(),model.second.end());
		set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),back_inserter(inserter)); //compare this vectors
		if(inserter.size() != model.second.size()) continue;
		else return model.first;
	}
	throw logic_error("Error in <SignalModel> : Couldn't match available parameters to any models specified");
}

string GetDirName(const string& filename){
	size_t found;
	found=filename.find_last_of("/\\");
	string folder = filename.substr(0,found);
	return folder;
}

void setMbb(RooRealVar& var, const string& path){
	double max;
	double min;
	int mass = returnMassPoint(path);
	if( find(mssmhbb::sr3.begin(),mssmhbb::sr3.end(),mass) != mssmhbb::sr3.end()) {
		min = 500;
		max = 1700;
	}
	else if( find(mssmhbb::sr2.begin(),mssmhbb::sr2.end(),mass) != mssmhbb::sr2.end() ){
		min = 350;
		max = 1190;
	}
	else if ( find(mssmhbb::sr1.begin(),mssmhbb::sr1.end(),mass) != mssmhbb::sr1.end() ) {
		min = 200;
		max = 650;
	}
	else {
		throw invalid_argument("AdjustMbbVariable:: No sub-range for mass: " + to_string(mass));
	}

	var.setRange(min,max);
}
