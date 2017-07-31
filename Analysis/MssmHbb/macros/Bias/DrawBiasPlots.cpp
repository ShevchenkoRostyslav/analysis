/*
 * DrawBiasPlots.cpp
 *
 *  Created on: 21 Mar 2017
 *      Author: shevchen
 *  Description:
 *  Macro to Draw plots dedicated to Bias study. Namely,
 *  Bias itself, (Nfit - Ninjected), values of shape parameters
 *  etc.
 */

#include <iostream>
#include <string>
#include <exception>
#include "boost/filesystem.hpp"

#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPad.h"

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"
#include "Analysis/MssmHbb/src/utilLib.cpp"

using namespace std;
namespace fs = boost::filesystem;

HbbStyle style;

void RunPlotBuilder(const std::map<int,std::string>& input_files, const std::string& output_path, const std::string& additional_file_name = "");
//Function to define the fitting PDF that was used
std::vector<std::string> GetBackgroundParameters(const std::string& file_name);
//Fucntions for particular plots:
void PlotGausHisto(const std::string& ttree_name, const std::string& draw_var_name, const std::string& root_input,const std::string& title,const std::string& output_name);

void PlotBias(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path); 		// Create Bias plot: (Npostfit - Nprefit)/ sigma_fit
void PlotdN(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path);			// Create (Npostfit - Nprefit) plot
void PlotSigma_fit(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path);	// Create sigma_fit plot
void PlotShapePar(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path, const std::string& shape_par); // Create plot of shape par

//Function to Draw histo from ttree
void DrawHistoFromTTree(TTree& tree, const std::string& var_name, const std::string& title, TCanvas& can, const std::string& output_name, const bool& do_fit = true);

//Adjust input names accroding to flags
map<int,string> AdjustSignalFilesName(const bool& bg_only, const int& signal_strength, const map<int,string>& my_map, const bool& fix_turnon);

int DrawBiasPlots(){

	style.set(PRIVATE);
	gStyle->SetOptStat("eMR");
	gStyle->SetOptFit(1111);
	gStyle->SetStatW(0.155);
	gStyle->SetStatH(0.12);
	gStyle->SetTitleXSize(0.04);
	gStyle->SetTitleXOffset(1.4);

	bool bg_only = false;
	auto fix_turnon = false;
//	int signal_strength = 4;
	vector<int> signal_stregths {0,1,2,4};
	string additional_file_name = "";

	map<int,string> input_sanity_files = {
//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},
//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitnovosibirsk_G4_R2_600GeV_r"},
//			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitnovosibirsk_G4_R2_700GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitnovosibirsk_G4_R2_900GeV_r"},
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitnovosibirsk_G4_R3_1100GeV_r"},
//			{1300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitnovosibirsk_G4_R3_1300GeV_r"}
	};

	map<int,string> input_bias_files = {

			/*
			 * pre-approval check: put 500 GeV to SR2, while 600 GeV to SR1
			 * 900 GeV -> SR3 and 1100 GeV -> SR2
			 */

//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitextnovoeffprod_G4_R1_600GeV_r"},
//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitextnovoeffprod_G4_R1_600GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_500GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R3_fitnovosibirsk_G4_R3_900GeV_r"},
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_1100GeV_r"},

			/*
			 * pre-approval check: check bias with more freedom in the fitting model,
			 * namely supernovosibirsk-1/2/3 has been checked
			 */

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsupernovoeffprod3_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsupernovoeffprod3_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsupernovoeffprod3_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsupernovoeffprod3_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod2_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitsupernovoeffprod3_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitsupernovoeffprod3_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitsupernovoeffprod3_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitsupernovoeffprod3_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitsupernovoeffprod3_G4_R1_500GeV_r"},

//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitsupernovosibirsk1_G4_R2_600GeV_r"},
//			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitsupernovosibirsk1_G4_R2_700GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitsupernovosibirsk1_G4_R2_900GeV_r"},
//
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R3_fitsupernovosibirsk1_G4_R3_1100GeV_r"},
//			{1300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R3_fitsupernovosibirsk1_G4_R3_1300GeV_r"},

//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitsupernovosibirsk1_G4_R2_600GeV_r"},
//			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitsupernovosibirsk1_G4_R2_700GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitsupernovosibirsk1_G4_R2_900GeV_r"},
//
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitsupernovosibirsk1_G4_R3_1100GeV_r"},
//			{1300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitsupernovosibirsk1_G4_R3_1300GeV_r"}

			/*
			 * Normal Bias tests
			 */
			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_500GeV_r"},
			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_600GeV_r"},
			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_700GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_R1_fitextnovoeffprod_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_R1_fitextnovoeffprod_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_R1_fitextnovoeffprod_R1_400GeV_r"},
//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_R1_fitextnovoeffprod_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_R1_fitextnovoeffprod_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_R1_fitextnovoeffprod_R1_400GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod8par_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toysuperdijeteffprod2_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsuperdijeteffprod2_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsuperdijeteffprod2_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsuperdijeteffprod2_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitsuperdijeteffprod2_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitberneffprod9par_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitberneffprod9par_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitberneffprod9par_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitberneffprod9par_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod9par_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod9par_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod9par_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyberneffprod9par_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetv3logprod_G4_R1_fitextnovoeffprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetv3logprod_G4_R1_fitextnovoeffprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetv3logprod_G4_R1_fitextnovoeffprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetv3logprod_G4_R1_fitextnovoeffprod_G4_R1_500GeV_r"},
//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_600GeV_r"},
//			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_700GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R2_fitnovosibirsk_G4_R2_900GeV_r"},
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R3_fitnovosibirsk_G4_R3_1100GeV_r"},
//			{1300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toydijetexp_G4_R3_fitnovosibirsk_G4_R3_1300GeV_r"},

//			{300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitdijetv3logprod_G4_R1_300GeV_r"},
//			{350,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitdijetv3logprod_G4_R1_350GeV_r"},
//			{400,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitdijetv3logprod_G4_R1_400GeV_r"},
//			{500,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toyextnovoeffprod_G4_R1_fitdijetv3logprod_G4_R1_500GeV_r"},
//			{600,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitdijetexp_G4_R2_600GeV_r"},
//			{700,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitdijetexp_G4_R2_700GeV_r"},
//			{900,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R2_fitdijetexp_G4_R2_900GeV_r"},
//			{1100,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitdijetexp_G4_R3_1100GeV_r"},
//			{1300,mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/BiasTest_toynovosibirsk_G4_R3_fitdijetexp_G4_R3_1300GeV_r"}
	};

	map<int,string> input_bias_folder_adjusted, input_sanity_files_adjusted;
	string output_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Bias_test/26072017_newSignalModels/";
	if (fs::create_directory(output_path)) cout<<"Directory: "<<output_path<<" has been created."<<endl;

	for(const auto& signal_strength : signal_stregths){
		//sanity checks
//		input_sanity_files_adjusted = AdjustSignalFilesName(bg_only,signal_strength,input_sanity_files,fix_turnon);
//		RunPlotBuilder(input_sanity_files_adjusted, output_path);
		//Bias
		input_bias_folder_adjusted = AdjustSignalFilesName(bg_only,signal_strength,input_bias_files,fix_turnon);
		RunPlotBuilder(input_bias_folder_adjusted,output_path,additional_file_name);
	}



	return 0;
}

void RunPlotBuilder(const std::map<int,std::string>& input_files, const std::string& output_path, const std::string& additional_file_name){
	/*
	 * Function to run individual plot makers.
	 *
	 * Loop over the map content is preformed
	 */
	for(const auto& file : input_files){
		//mass
		int mass = file.first;
		//construct output name and root_file name
		auto found = file.second.find_last_of("/");
		auto name = file.second.substr(found + 1) + additional_file_name;
		auto root_file_name = file.second + "/" + name + ".root";
		auto single_output_path = output_path + name + "/";
		//create directory if it doesn't exist
		if(fs::create_directory(single_output_path)) cout<<"Directory: "<<single_output_path<<" has been created."<<endl;
		//signal strength from the file name:
		auto signal_strength = file.second.substr(file.second.find("GeV_r") + 5,1);
		cout<<"\n";

		//Plot functions
		PlotBias(root_file_name,mass,signal_strength,single_output_path); 		// Create Bias plot: (Npostfit - Nprefit)/ sigma_fit
		PlotdN(root_file_name,mass,signal_strength,single_output_path);			// Create (Npostfit - Nprefit) plot
		PlotSigma_fit(root_file_name,mass,signal_strength,single_output_path);	// Create sigma_fit plot
		PlotShapePar(root_file_name,mass,signal_strength,single_output_path,"all");
	}
}

void PlotShapePar(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path, const std::string& shape_par){
	/*
	 * Plot shape parameter
	 */
	//output name
	string out_name = output_path + "ShapePar_" + shape_par + ".pdf";
	//draw variable
	string var = shape_par;
	//title:
	string title = "M_{A} = " + to_string(mass) + ", #mu = " + mu + ", " + shape_par + ";" + shape_par + ";Counts";
	PlotGausHisto("Toys", var, root_input, title, out_name);
}

void PlotBias(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path){
	/*
	 * Plot Bias distribution (Npostfit - Nprefit)/ sigma_fit
	 */
	//output name
	string out_name = output_path + "Bias.pdf";
	//draw variable
	string var = "bias";
	//title:
	string title = "M_{A} = " + to_string(mass) + ", #mu = " + mu + ", (N_{postfit} - N_{prefit})/ #sigma_{fit};#frac{(N_{postfit} - N_{prefit})}{#sigma_{fit}};Counts";
	std::cout<<"mu = "<<mu<<", M = "<<mass<<": Bias = ";
	PlotGausHisto("Toys", var, root_input, title, out_name);
}

void PlotdN(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path){
	/*
	 * Plot dN distribution (Npostfit - Nprefit)
	 */
	//output name
	string out_name = output_path + "dN_post_pre.pdf";
	//draw variable
	string var = "dN_post_pre";
	//title:
	string title = "M_{A} = " + to_string(mass) + ", #mu = " + mu + ", (N_{postfit} - N_{prefit});(N_{postfit} - N_{prefit});Counts";
	std::cout<<"mu = "<<mu<<", M = "<<mass<<": dN = ";
	PlotGausHisto("Toys", var, root_input, title, out_name);
}

void PlotSigma_fit(const std::string& root_input,const int& mass, const std::string& mu, const std::string& output_path){
	/*
	 * Plot Sigma fit distribution
	 */
	//output name
	string out_name = output_path + "ns_err_postfit.pdf";
	//draw variable
	string var = "ns_err_postfit";
	//title:
	string title = "M_{A} = " + to_string(mass) + ", #mu = " + mu + ", #sigma_{fit};#sigma_{fit};Counts";
	PlotGausHisto("Toys", var, root_input, title, out_name);
}

void PlotGausHisto(const std::string& ttree_name, const std::string& draw_var_name, const std::string& root_input,const std::string& title,const std::string& output_name){
	/*
	 * The most general function to draw gausian shape histo
	 * and fit it with gaus function
	 */
	TCanvas can("can","can",800,600);
	TFile f(root_input.c_str(),"READ");
//	try{
		CheckZombie(f);
//	}
//	catch (exception& e){
//		return;
//	}
	TTree &tree = *(TTree*) f.Get(ttree_name.c_str());
	//If draw_var_name = "all", need to retrieve names of the parameters from TTree
	if(draw_var_name == "all"){
		auto branches = *tree.GetListOfBranches();
		for(const auto& var : branches){
			string var_name = var->GetName();
			string hard_title = var_name + ";" + var_name + ";Counts";
			string hard_output_name = output_name.substr(0,output_name.find_last_of("/")) + "/ShapePar_" + var_name + ".pdf";
			bool fit = false;
			if(var_name == "Bias" || var_name == "dN_post_pre" || var_name == "minNLL" || var_name == "chi_2") fit = true;
			DrawHistoFromTTree(tree,var_name,hard_title,can,hard_output_name,fit);
		}
	}
	else{
		DrawHistoFromTTree(tree,draw_var_name,title,can,output_name);
	}

}

void DrawHistoFromTTree(TTree& tree, const std::string& var_name, const std::string& title, TCanvas& can, const std::string& output_name, const bool& do_fit){
	/*
	 * Function to draw TH1 from TTree and Fit it
	 */
	string draw_restrictions = "";
	if(var_name == "Bias" || var_name == "bias") draw_restrictions = " bias > -10 && bias < 10";
	if(var_name == "dN_post_pre") draw_restrictions = " dN_post_pre > -15000 && dN_post_pre < 15000";
	tree.Draw(var_name.c_str(),draw_restrictions.c_str());
	TH1F &h = *(TH1F*)gPad->GetPrimitive("htemp");
	h.SetTitle(title.c_str());
	h.SetStats(1);
	h.Draw();

	//Fit
	if(do_fit){
		TF1 fit("fit","gaus",h.GetXaxis()->GetXmin(),h.GetXaxis()->GetXmax());
		h.Fit(&fit, "LLq");
		std::cout<<" "<<fit.GetParName(1)<<" = "<<fit.GetParameter(1)<<" +/- "<<fit.GetParError(1)<<fit.GetParName(2)<<" = "<<fit.GetParameter(2)<<" +/- "<<fit.GetParError(2)<<" Mean/Sigma = "<<fit.GetParameter(1)/fit.GetParameter(2)<<std::endl;
	}

	//Make picture centered
	/*
	double mean = h.GetMean();
	double xmin = h.GetXaxis()->GetXmin();
	double xmax = h.GetXaxis()->GetXmax();

	if(abs(xmin - mean) > abs(xmax-mean)) h.GetXaxis()->SetRangeUser(xmin, mean + abs(xmin - mean));
	else h.GetXaxis()->SetRangeUser(mean - abs(xmax - mean),xmax);
	*/

	can.Print(output_name.c_str());
}

std::vector<std::string> GetBackgroundParameters(const std::string& file_name){
	/*
	 * Function to get background parameters according to the input root file name
	 */
	vector<string> novo_pars{"width1","peak1","tail1","width1_err","peak1_err","tail1_err"};
	vector<string> extnovo_pars{"width","peak","tail","par4","slope_novoeff","turnon_novoeff","width_err","peak_err","tail_err","par4_err","slope_novoeff_err","turnon_novoeff_err"};
	vector<string> dijetexp_pars{};
	if(file_name.find("fitnovosibirsk") != string::npos){
		return novo_pars;
	}
	else if (file_name.find("fitdijetexp") != string::npos){
		return dijetexp_pars;
	}
	else if (file_name.find("fitextnovoeffprod") != string::npos){
		return extnovo_pars;
	}
	else {
		throw logic_error("Wrong or Incomplete file name has been provided: " + file_name);
	}
}

map<int,string> AdjustSignalFilesName(const bool& bg_only, const int& signal_strength, const map<int,string>& my_map, const bool& fix_turnon){
	/*
	 * Function to adjust file names in the map
	 * accroding to the flags
	 */
	map<int,string> temp_map;
	for (auto & m : my_map){
		temp_map[m.first] = m.second + to_string(signal_strength);
		if(bg_only) temp_map[m.first] += "_BGonly";
		if(fix_turnon) temp_map[m.first] += "_fixedTurnOn";
	}
	return temp_map;
}

