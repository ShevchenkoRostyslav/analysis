/*
 * signal_shapes.cpp
 *
 *  Created on: 13 Oct 2017
 *      Author: shevchen
 */
#include "RooDataHist.h"
#include "TH1.h"

#include <string>
#include <vector>
#include <memory>
//#include "TColor.h"
#include "TLegendEntry.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/Tools/interface/Analysis.h"

#include <boost/program_options.hpp>

using namespace boost::program_options;

typedef std::unique_ptr<variables_map> pVariables_map;

std::unique_ptr<variables_map> ProcessUserInput(int argc, char** argv);
void SetupPreDefSubRangeOptions_(std::unique_ptr<variables_map> & var_map);
void DrawSignalShapes(const std::unique_ptr<variables_map> & var_map);
void DrawSignalTemplates_(TH2* frame, const std::vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation);
void SetCenteredTLegendHeader_(TLegend &leg, const char* header);
std::string PrepareOutputName(const std::unique_ptr<variables_map> & var_map);
TH1D* GetHistogramFromEventLoop(const int& mp);
TH1D* GetHistogram(const int& mp);
struct sub_range{
	sub_range() : nbins_(1), xmin_(0), xmax_(2000) {};
	sub_range(const int& nbins, const int& xmin, const int& xmax) : nbins_(nbins), xmin_(xmin), xmax_(xmax) {};
	sub_range(const sub_range& sr) {nbins_ = sr.nbins_; xmin_ = sr.xmin_ ; xmax_ = sr.xmax_;};

	int nbins_{1};
	int xmin_{0};
	int xmax_{0};
};
sub_range DefineSubrange(const int& mp);

HbbStyle style;

//*****************Test different colors*********************
TColor *col0 = new TColor(10000, 0, 0.4470588235, 0.6980392157);
TColor *col1 = new TColor(10001, 0.337254902, 0.7058823529, 0.9137254902);
TColor *col2 = new TColor(10002, 0.8, 0.4745098039, 0.6549019608);
TColor *col3 = new TColor(10003, 0, 0.6196078431, 0.4509803922);
TColor *col4 = new TColor(10004, 0.8352941176, 0.368627451, 0);
TColor *col_darkturquoise 		= new TColor(10005,25./255, 109./255, 128./255, 	"darkturquoise");
TColor *col_violet				= new TColor(10006,129./255, 31./255, 157./255, 	"violet");
TColor *col_darkblue				= new TColor(10007,34./255, 91./255, 196./255, 	"darkblue");
TColor *col_lightblue			= new TColor(10008,50./255, 160./255, 246./255, 	"lightblue");
TColor *col_pink					= new TColor(10009,247./255, 125./255, 247./255, "pink");
TColor *col_lightturquoise		= new TColor(10010,58./255, 209./255, 218./255, 	"lightturquoise");
TColor *col_darkred				= new TColor(10011,166./255, 25./255, 63./255, 	"darkred");
TColor *col_lightorange			= new TColor(10012,242./255, 85./255, 41./255,  	"lightorange");
TColor *col_green				= new TColor(10013,38./255, 178./255, 94./255,  	"green");
TColor *col_yellow				= new TColor(10014,238./255, 239./255, 78./255,  	"yellow");
TColor *col_lightgreen			= new TColor(10015,163./255, 248./255, 137./255, 	"lightgreen");
TColor *col_beige				= new TColor(10016,248./255, 230./255, 193./255, 	"beige");
//************************************************************
//vector<int> nice_colors {kRed,kBlue,kGreen+3,kOrange-3,kMagenta+1,kAzure-4};
std::vector<int> nice_colors {10011,10007,10013,10003,10004,5};
std::vector<int> nice_linestyles {1,3,2,6,4,3};
//vector<int> nice_linestyles {1,1,1,1,1,1,1};

int main(int argc,char** argv){
	style.setTDRstyle(PRIVATE);
	auto user_input = ProcessUserInput(argc, argv);
	DrawSignalShapes(user_input);
}

void DrawSignalShapes(const std::unique_ptr<variables_map> & var_map){
	/*
	 * Main plotting function
	 * Plot signal shapes: templates or pdfs
	 */
	//Read input from the boost::po::variable_map
	float xmin = var_map->at("xmin").as<float>();
	float ymin = var_map->at("ymin").as<float>();
	float xmax = var_map->at("xmax").as<float>();
	float ymax = var_map->at("ymax").as<float>();
	std::vector<int> mass_points = var_map->at("mass_points").as<std::vector<int>>();
	std::string legend_pos   = var_map->at("legend_pos").as<std::string>();
	double rebin 		= var_map->at("rebin").as<float>();
	bool normalisation = var_map->at("normalisation").as<float>();

	TCanvas can;
	//Setup plotting frame
	ymin = 0.001;
	TH2D *frame = new TH2D("frame","",1,xmin,xmax,1,ymin,ymax);
	frame->GetYaxis()->SetTitle("Arbitrary units");
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMass());
	frame->Draw();
	//Legend
	auto &leg = *HbbStyle::legend(legend_pos.c_str(),mass_points.size()+2);
	HbbStyle::setLegendStyle(&leg);

	DrawSignalTemplates_(frame, mass_points, leg, rebin, normalisation);

	leg.Draw();

	HbbStyle::drawStandardTitle();
	auto out_name = PrepareOutputName(var_map);
	gPad->RedrawAxis();
	can.RedrawAxis();
	can.Print( (out_name + ".pdf").c_str());
	TH2D *frame2 = new TH2D("frame2","",1,xmin,xmax,1,1e-03,30);
	frame2->GetYaxis()->SetTitle("Arbitrary units");
	frame2->GetXaxis()->SetTitle(HbbStyle::axisTitleMAH());
	frame2->Draw();
	leg.Clear();
	DrawSignalTemplates_(frame2, mass_points, leg, rebin, normalisation);
	gPad->RedrawAxis();
	can.RedrawAxis();
	leg.Draw();
	HbbStyle::drawStandardTitle();
	gPad->SetLogy();
	can.Print( (out_name + "_log.pdf").c_str());

}

void DrawSignalTemplates_(TH2* frame, const std::vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation){
	/*
	 * Main funtion to draw signal templates
	 */
	SetCenteredTLegendHeader_(leg,"Signal Shape");
	int i = 0;
	for(const auto& mp : mass_points){
		auto *h = GetHistogramFromEventLoop(mp);
		h->Rebin(rebin);
		h->SetLineColor(nice_colors.at(i));
		h->SetLineWidth(2);
		h->SetLineStyle(nice_linestyles.at(i));
		if(normalisation != -1) h->Scale(normalisation / h->Integral());
		leg.AddEntry(h,("m_{A/H} = " + std::to_string(mp) + " GeV").c_str(),"l");
		h->Draw("HIST same");
		++i;
	}
}

TH1D* GetHistogramFromEventLoop(const int& mp){
	/*
	 * Function to perform the event loop
	 * ofver the MC sample
	 * to fill the true Higgs mass
	 */
	std::string ntuple = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-" + std::to_string(mp) + "_TuneCUETP8M1_13TeV-pythia8.txt";
	analysis::tools::Analysis analysis(ntuple);
	analysis.addTree<analysis::tools::GenParticle> ("GenParticles","MssmHbb/Events/prunedGenParticles");
	std::string histo_name = "histo_" + std::to_string(mp);
	auto histo = GetHistogram(mp);
	int nevents = 0.2 * analysis.size();
	for (int i = 0 ; i < nevents; ++i){
		analysis.event(i);
		auto genParticles = analysis.collection<analysis::tools::GenParticle>("GenParticles");
		analysis::tools::GenParticle higgs;
		int n_higgs = 0;
		for ( int j = 0 ; j < genParticles->size() ; ++j )
		{
			analysis::tools::GenParticle gp = genParticles->at(j);
			if ( gp.pdgId() == 36 && gp.status() > 60 ) {
				higgs = gp;
				++n_higgs;
			}
			if(n_higgs != 0) break;
		}
		if(n_higgs == 0) continue;
		histo->Fill(higgs.m());
	}
	return histo;
}

std::string PrepareOutputName(const std::unique_ptr<variables_map> & var_map){
	/*
	 * Make an output name from the input
	 */
	std::string output_folder = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/";
	std::string output_name = "";
	//Templates or PDFs
	output_name = "TrueSignalTemplates_";
	//Add mass points
	auto &v = var_map->at("mass_points").as<std::vector<int>>();
	for(const auto& m : v){
		output_name+= std::to_string(m);
		if(int(&m - &v[0]) != int(v.size()) - 1) output_name += "and";
	}

	return (output_folder + output_name);
}

std::unique_ptr<variables_map> ProcessUserInput(int argc,char** argv){
	/*
	 * Function to process user input for this macro.
	 * User input:
	 * x_min - min value for X-axis
	 * x_max - mac value for X-axis
	 * y_min - min value for Y-axis
	 * y_max - max value for Y-axis
	 * Mass_points - vector of the mass points to be plotted
	 *
	 * OR
	 *
	 * SR1 / SR2 / SR3 - to plot signal shapes from
	 * corresponded sub-ranges
	 */

	// general options to show help menu
	// and control output verbosity
	options_description generalOptions("General options");
	generalOptions.add_options()
		("help,h", "produce help message")
		("verbose,v", value<int>()->default_value(0), "more verbose output")
		("rebin,r", value<float>()->default_value(1), "rebin histograms")
		("normalisation", value<float>()->default_value(1),"Total normalisation of the tempaltes/pdfs. -1 - normalise to a real area from a selection")
		("legend_pos",value<std::string>()->default_value("top,right"),"position of the legend: 'top,right', 'top,left'");

	// command line options
	options_description cmdOptions("Command line options");
	cmdOptions.add_options()
		("xmin",value<float>()->required(),"min value for X-axis")
		("xmax",value<float>()->required(),"max value for X-axis")
		("ymin",value<float>()->default_value(0),"min value for Y-axis")
		("ymax",value<float>()->default_value(1),"max value for Y-axis")
		("mass_points",value<std::vector<int>>()->multitoken(),"vector of the signal mass points to be plotted");

	//Options to be used with pre-deffined SRs input
	options_description predefSRsOptions("sub-range options");
	predefSRsOptions.add_options()
		("sub_range,s",value<std::string>(),"choice the sub-range with pre-defined settings: SR1/SR2/SR3");

    // Hidden options, will be allowed but not shown to the user.
    options_description hiddenOptions("Hidden options");
    hiddenOptions.add_options()
        ("test", value<int>()->default_value(0),"if test specify 1. Default 0.");

    /*
     * Join options together
     */
    //for CMD
    options_description allCmdOptions;
    allCmdOptions.add(generalOptions).add(cmdOptions).add(hiddenOptions);

    //for pre-def. SRs
    options_description allPreDefOptions;
    allPreDefOptions.add(generalOptions).add(predefSRsOptions).add(hiddenOptions);

    // visible for the user
    options_description visibleOptions("Allowed Options");
    visibleOptions.add(generalOptions).add(cmdOptions).add(predefSRsOptions);

    pVariables_map output_vm(new variables_map());
    //Check whether pre-deffined method is used
    store(command_line_parser(argc,argv).options(allPreDefOptions).allow_unregistered().run(),*output_vm);
    notify(*output_vm);

    //show help menu
    if (output_vm->count("help")) {
	    std::cout << visibleOptions << std::endl;
	    std::exit(0);
    }

    //if not pre-def. options
    if(!output_vm->count("sub_range")){
    		store(parse_command_line(argc,argv,allCmdOptions),*output_vm);
    		try {
    			notify(*output_vm);
    		} catch (const required_option& e) {
				throw std::exception(e);
    		}
    }
    else {
    		//Setup pre-defined options
    		SetupPreDefSubRangeOptions_(output_vm);
    }

    return output_vm;
}

void SetupPreDefSubRangeOptions_(std::unique_ptr<variables_map> & var_map){
	/*
	 * Setup default values for the sub-ranges.
	 * Default values are made for:
	 * x_min - min value for X-axis
	 * x_max - mac value for X-axis
	 * y_min - min value for Y-axis
	 * y_max - max value for Y-axis
	 * Mass_points - vector of the mass points to be plotted
	 */
	auto sub_range = (*var_map)["sub_range"].as<std::string>();
	auto normalisation = (*var_map)["normalisation"].as<float>();
	float xmin,xmax,ymin,ymax;
	float rebin = 1;
	std::vector<int> vec;
	std::string legend_pos = "right,top";
	if(findStrings(sub_range,"sr1")){
		vec = mssmhbb::sr1;
		xmin = 100.001; xmax = 750; ymin = 0; ymax = 1.;
		xmin = 100.001;
		ymax = 100;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.6;
	}
	else if (findStrings(sub_range,"sr2")){
		vec = mssmhbb::sr2;
		xmin = 100; xmax = 1199.99; ymin = 0; ymax = 1.;
		xmin = 100;
		ymax = 100;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.55;
	}
	else if (findStrings(sub_range,"sr3")){
		legend_pos = "right,top";
		vec = mssmhbb::sr3;
		xmin = 100; xmax = 1700; ymin = 0; ymax = 1.;
		ymax = 40;
		rebin = 2;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.5;

	}
	else {
		legend_pos = "right,top";
		vec = {300,500,700,900,1100,1300};
		xmin = 100; xmax = 1500;
		ymax = 0.7;
		ymin = 0.001;
	}

	var_map->insert(std::make_pair("xmin",variable_value(float(xmin),false)));
	var_map->insert(std::make_pair("xmax",variable_value(float(xmax),false)));
	var_map->insert(std::make_pair("ymin",variable_value(float(ymin),false)));
	var_map->insert(std::make_pair("ymax",variable_value(float(ymax),false)));
	var_map->insert(std::make_pair("mass_points",variable_value(std::vector<int>(vec),false)));
	var_map->at("legend_pos") = variable_value(std::string(legend_pos),false);
	var_map->insert(std::make_pair("rebin",variable_value(float(rebin),false)));
}

void SetCenteredTLegendHeader_(TLegend &leg, const char* header){
	/*
	 * Function to make TLegend header centered
	 */
	leg.SetHeader(header);
	auto head = static_cast<TLegendEntry*>(leg.GetListOfPrimitives()->First());
	head->SetTextAlign(22);
}

TH1D* GetHistogram(const int& mp){
	std::string histo_name = "histo_" + std::to_string(mp);
	auto sr = DefineSubrange(mp);
	TH1D* h = new TH1D(histo_name.c_str(),histo_name.c_str(),sr.nbins_,sr.xmin_,sr.xmax_);
	return h;
}

sub_range DefineSubrange(const int& mp){
	sub_range sr;
//	if(mp >= 300 && mp <= 400){
//		sr = sub_range(180,200,650);
//	}
//	else if (mp >=500 && mp <= 900){
//		sr = sub_range(168,350,1190);
//	}
//	else sr = sub_range(176,500,1700);
	if(mp >= 300 && mp <= 400){
		sr = sub_range(500,200,1400);
	}
	else if (mp >=500 && mp <= 900){
		sr = sub_range(500,200,1400);
	}
	else sr = sub_range(500,200,1400);

	return sr;
}
