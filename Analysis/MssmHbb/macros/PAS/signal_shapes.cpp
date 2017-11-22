/*
 * signal_shapes.cpp
 *
 *  Created on: 13 Oct 2017
 *      Author: shevchen
 */
#include "RooDataHist.h"

#include <memory>
#include "TLegendEntry.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"

#include <boost/program_options.hpp>

using namespace std;
using namespace boost::program_options;

typedef unique_ptr<variables_map> pVariables_map;

unique_ptr<variables_map> ProcessUserInput(int argc, char** argv);
void SetupPreDefSubRangeOptions_(unique_ptr<variables_map> & var_map);
void DrawSignalShapes(const unique_ptr<variables_map> & var_map);
void DrawSignalTemplates_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation);
void DrawSignalPDFs_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation);
void DrawSignalPDFsAndTemplates_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation);
void SetCenteredTLegendHeader_(TLegend &leg, const char* header);
string PrepareOutputName(const unique_ptr<variables_map> & var_map);
double GetSignalNormalisation_(const int& mp);

vector<int> nice_colors {kRed,kBlue,kGreen+3,kOrange-3,kMagenta+1,kAzure-4};
vector<int> nice_linestyles {1,2,4,6,9,10};

int main(int argc,char** argv){
	HbbStyle style;
	style.setTDRstyle(PRELIMINARY_SIMULATION);
	auto user_input = ProcessUserInput(argc, argv);
	DrawSignalShapes(user_input);
}

void DrawSignalShapes(const unique_ptr<variables_map> & var_map){
	/*
	 * Main plotting function
	 * Plot signal shapes: templates or pdfs
	 */
	//Read input from the boost::po::variable_map
	float xmin = var_map->at("xmin").as<float>();
	float ymin = var_map->at("ymin").as<float>();
	float xmax = var_map->at("xmax").as<float>();
	float ymax = var_map->at("ymax").as<float>();
	vector<int> mass_points = var_map->at("mass_points").as<vector<int>>();
	bool draw_templates = var_map->at("templates").as<bool>();
	bool draw_pdfs = var_map->at("pdfs").as<bool>();
	string legend_pos   = var_map->at("legend_pos").as<string>();
	double rebin 		= var_map->at("rebin").as<float>();
	bool normalisation = var_map->at("normalisation").as<float>();

	TCanvas can;
	//Setup plotting frame
	TH2D *frame = new TH2D("frame","",1,xmin,xmax,1,ymin,ymax);
	frame->GetYaxis()->SetTitle("a.u.");
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMass());
	frame->Draw();
	//Legend
	auto &leg = *HbbStyle::legend(legend_pos.c_str(),mass_points.size()+2);
	HbbStyle::setLegendStyle(&leg);

	if(draw_templates && !draw_pdfs) DrawSignalTemplates_(frame, mass_points, leg, rebin, normalisation);
	if(draw_pdfs && !draw_templates) DrawSignalPDFs_(frame, mass_points, leg, rebin, normalisation);
	if(draw_templates && draw_pdfs) DrawSignalPDFsAndTemplates_(frame, mass_points, leg, rebin, normalisation);
//	else DrawSignalPDFs_(frame,mass_points,leg,rebin);
//	DrawSignalTemplates_(frame,mass_points,leg,rebin);

	leg.Draw();

	HbbStyle::drawStandardTitle("out");
	auto out_name = PrepareOutputName(var_map);
	can.Print( (out_name + ".pdf").c_str());

}

void DrawSignalTemplates_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation){
	/*
	 * Main funtion to draw signal templates
	 */
	SetCenteredTLegendHeader_(leg,"Signal Shape");
	int i = 0;
	for(const auto& mp : mass_points){
		auto *f = new TFile(mssmhbb::signal_templates.at(mp).c_str());
		auto *h = GetFromTFile<TH1>(*f,"templates/bbH_Mbb_VIS");
		h->Rebin(rebin);
		h->SetLineColor(nice_colors.at(i));
		h->SetLineWidth(2);
		h->SetLineStyle(nice_linestyles.at(i));
		if(normalisation != -1) h->Scale(normalisation / h->Integral());
		leg.AddEntry(h,("m_{A/H} = " + to_string(mp) + " GeV").c_str(),"l");
		h->Draw("HIST same");
		++i;
	}
}

void DrawSignalPDFs_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation){
	/*
	 * Main function to draw RooFit PDFs
	 */
	int i = 0;
	SetCenteredTLegendHeader_(leg,"Signal Parametrisation");
	for(const auto& mp : mass_points){
		auto *x		= analysis::GetRooObjectFromTFile<RooRealVar>(mssmhbb::signal_fits.at(mp),"mbb");
		x->setRange(frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax());
		auto *pdf 	= analysis::GetRooObjectFromTFile<RooAbsPdf>(mssmhbb::signal_fits.at(mp),"signal");
		//TH1 Needed for normalisation
		auto *file = new TFile(mssmhbb::signal_templates.at(mp).c_str());
		auto *h = GetFromTFile<TH1>(*file,"templates/bbH_Mbb_VIS");
		h->Rebin(rebin);
		if(normalisation != -1) h->Scale(normalisation / h->Integral());
		RooDataHist *hist = new RooDataHist("hist","hist",RooArgList(*x),h);
		auto *fr 	= x->frame();
		hist->plotOn(fr,RooFit::Invisible());
		pdf->plotOn(fr,RooFit::LineColor(nice_colors.at(i)),RooFit::Range(frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax()),RooFit::Name((string("signal_pdf_") + to_string(i)).c_str()));//,RooFit::Normalization(norm, RooAbsReal::Relative));//,RooFit::Normalization(norm->getValV()));
		leg.AddEntry(fr->getCurve((string("signal_pdf_") + to_string(i)).c_str()),("m_{A/H} = " + to_string(mp) + " GeV").c_str(),"l");
		fr->Draw("same");
		++i;
	}
}

void DrawSignalPDFsAndTemplates_(TH2* frame, const vector<int>& mass_points, TLegend & leg, const float& rebin, const float& normalisation){
	/*
	 * Draw Both PDFs and Templates
	 */
	DrawSignalTemplates_(frame, mass_points, leg, rebin, normalisation);
	leg.Clear();
	DrawSignalPDFs_(frame, mass_points, leg, rebin, normalisation);
	leg.SetHeader("Signal");
}

string PrepareOutputName(const unique_ptr<variables_map> & var_map){
	/*
	 * Make an output name from the input
	 */
	string output_folder = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/";
	string output_name = "";
	//Templates or PDFs
	auto pdfs = var_map->at("pdfs").as<bool>();
	auto templates = var_map->at("templates").as<bool>();
	if(templates && !pdfs) output_name = "SignalTemplates_PAS_";
	if(pdfs && !templates) output_name = "SignalPDFs_PAS_";
	if(pdfs && templates) output_name =  "SignalPDFsAndTemplates_PAS_";
	//Add mass points
	auto &v = var_map->at("mass_points").as<vector<int>>();
	for(const auto& m : v){
		output_name+= to_string(m);
		if(int(&m - &v[0]) != int(v.size()) - 1) output_name += "and";
	}

	return (output_folder + output_name);
}

unique_ptr<variables_map> ProcessUserInput(int argc,char** argv){
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
		("templates", value<bool>()->default_value(false),"draw signal templates")
		("pdfs", value<bool>()->default_value(true),"draw signal pdfs")
		("normalisation", value<float>()->default_value(-1),"Total normalisation of the tempaltes/pdfs. -1 - normalise to a real area from a selection")
		("legend_pos",value<string>()->default_value("top,right"),"position of the legend: 'top,right', 'top,left'");

	// command line options
	options_description cmdOptions("Command line options");
	cmdOptions.add_options()
		("xmin",value<float>()->required(),"min value for X-axis")
		("xmax",value<float>()->required(),"max value for X-axis")
		("ymin",value<float>()->default_value(0),"min value for Y-axis")
		("ymax",value<float>()->default_value(1),"max value for Y-axis")
		("mass_points",value<vector<int>>()->multitoken(),"vector of the signal mass points to be plotted");

	//Options to be used with pre-deffined SRs input
	options_description predefSRsOptions("sub-range options");
	predefSRsOptions.add_options()
		("sub_range,s",value<string>(),"choice the sub-range with pre-defined settings: SR1/SR2/SR3");

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
	    cout << visibleOptions << endl;
	    exit(0);
    }

    //if not pre-def. options
    if(!output_vm->count("sub_range")){
    		store(parse_command_line(argc,argv,allCmdOptions),*output_vm);
    		try {
    			notify(*output_vm);
    		} catch (const required_option& e) {
				throw exception(e);
    		}
    }
    else {
    		//Setup pre-defined options
    		SetupPreDefSubRangeOptions_(output_vm);
    }

    return output_vm;
}

void SetupPreDefSubRangeOptions_(unique_ptr<variables_map> & var_map){
	/*
	 * Setup default values for the sub-ranges.
	 * Default values are made for:
	 * x_min - min value for X-axis
	 * x_max - mac value for X-axis
	 * y_min - min value for Y-axis
	 * y_max - max value for Y-axis
	 * Mass_points - vector of the mass points to be plotted
	 */
	auto sub_range = (*var_map)["sub_range"].as<string>();
	auto draw_templates = (*var_map)["templates"].as<bool>();
	auto normalisation = (*var_map)["normalisation"].as<float>();
//	auto draw_pdfs = (*var_map)["pdfs"].as<bool>();
	float xmin,xmax,ymin,ymax;
	float rebin = 1;
	vector<int> vec;
	string legend_pos = "right,top";
	if(findStrings(sub_range,"sr1")){
		vec = mssmhbb::sr1;
		xmin = 100.001; xmax = 750; ymin = 0; ymax = 0.2;
		if(draw_templates) {
			xmin = 100.001;
			ymax = 100;
		}
		ymax = 100;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.3;
	}
	else if (findStrings(sub_range,"sr2")){
		vec = mssmhbb::sr2;
		xmin = 100; xmax = 1199.99; ymin = 0; ymax = 0.2;
		if(draw_templates) {
			xmin = 100;
			ymax = 100;
		}
		ymax = 100;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.2;

//		var_map->insert(make_pair("xmin",variable_value(float(350.0),false)));
//		var_map->insert(make_pair("xmax",variable_value(float(1190.0),false)));
//		var_map->insert(make_pair("ymin",variable_value(float(0.0),false)));
//		var_map->insert(make_pair("ymax",variable_value(float(0.2),false)));
//		var_map->insert(make_pair("mass_points",variable_value(vector<int>(mssmhbb::sr2),false)));


//		var_map->at("xmin").value() = boost::any(200);
//		var_map->at("xmax").value() = boost::any(1100);
//		var_map->at("ymin").value() = boost::any(0);
//		var_map->at("ymax").value() = boost::any(1);
//		var_map->at("mass_points").value() = boost::any(mssmhbb::sr2);
	}
	else {
		legend_pos = "left,top";
		vec = mssmhbb::sr3;
		xmin = 100; xmax = 1700; ymin = 0; ymax = 0.1;
//		if(draw_templates) ymax = 20;
		ymax = 40;
		rebin = 2;

		//Adjust ymax for the normalised plots
		if(normalisation == 1) ymax = 0.1;

//		var_map->insert(make_pair("xmin",variable_value(float(500.0),false)));
//		var_map->insert(make_pair("xmax",variable_value(float(1700.0),false)));
//		var_map->insert(make_pair("ymin",variable_value(float(0.0),false)));
//		var_map->insert(make_pair("ymax",variable_value(float(0.2),false)));
//		var_map->insert(make_pair("mass_points",variable_value(vector<int>(mssmhbb::sr3),false)));
//		var_map->at("xmin").value() = boost::any(500);
//		var_map->at("xmax").value() = boost::any(1700);
//		var_map->at("ymin").value() = boost::any(0);
//		var_map->at("ymax").value() = boost::any(1);
//		var_map->at("mass_points").value() = boost::any(mssmhbb::sr3);
	}

	var_map->insert(make_pair("xmin",variable_value(float(xmin),false)));
	var_map->insert(make_pair("xmax",variable_value(float(xmax),false)));
	var_map->insert(make_pair("ymin",variable_value(float(ymin),false)));
	var_map->insert(make_pair("ymax",variable_value(float(ymax),false)));
	var_map->insert(make_pair("mass_points",variable_value(vector<int>(vec),false)));
	var_map->at("legend_pos") = variable_value(string(legend_pos),false);
	var_map->insert(make_pair("rebin",variable_value(float(rebin),false)));
//	var_map->insert(make_pair("legend_pos",variable_value(string(legend_pos),false)));

	cout<<"STANDARD VALUES: xmin = "<<var_map->at("xmin").as<float>()<<" xmax = "<<
			var_map->at("xmax").as<float>()<<" ymin = "<<var_map->at("ymin").as<float>()<<" ymax = "<<
			var_map->at("ymax").as<float>()<<" legend position: "<<var_map->at("legend_pos").as<string>()<<endl;
}

void SetCenteredTLegendHeader_(TLegend &leg, const char* header){
	/*
	 * Function to make TLegend header centered
	 */
	leg.SetHeader(header);
	auto head = static_cast<TLegendEntry*>(leg.GetListOfPrimitives()->First());
	head->SetTextAlign(22);
}

double GetSignalNormalisation_(const int& mp){
	/*
	 * Function to get the total signal normalisation from the signal template
	 */
	TFile f(mssmhbb::signal_templates.at(mp).c_str());
	auto *h_full = GetFromTFile<TH1>(f,"templates/bbH_Mbb_VIS");
//	auto *h_window = GetFromTFile<TH1>(f,"bbH_Mbb");
	double norm = h_full->Integral();
//	double norm = GetFromTFile<TH1>(f,"templates/bbH_Mbb_VIS")->Integral();

	return norm;
}
