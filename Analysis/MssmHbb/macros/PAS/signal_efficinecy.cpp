/*
 * signal_efficinecy.cpp
 *
 *  Created on: 12 Oct 2017
 *      Author: rostyslav
 *
 *
 */
#include <fstream>
#include <boost/regex.hpp>
#include <regex>
#include <memory>

#include "TCanvas.h"
#include "TGraphErrors.h"

#include <boost/program_options.hpp>

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
HbbStyle style;

using namespace std;
using namespace boost::program_options;

typedef unique_ptr<variables_map> pVariables_map;
typedef unique_ptr<TH2D> pTH2D;

struct Cut{
	Cut(const string& name) : name(name), legend_title(""), weighted_value(0.0), weight(0.0), nevents(0) {};
	Cut(const string& name, const string& legend_title) : name(name), legend_title(legend_title), weighted_value(0.0), weight(0.0), nevents(0) {};
	Cut() : name(""), legend_title(""), weighted_value(0.0), weight(0.0), nevents(0) {};
	void show() {cout<<name<<" in TLegend: "<<legend_title<<" Nevents = "<<nevents<<" weighted_events = "<<weighted_value<<" with weight = "<<weight<<endl;}
	string name = "";
	string legend_title = "";
	float weighted_value = 0.0;
	float weight = 0.0;
	int nevents = 0;
};

pVariables_map ProcessUserInput(int argc, char** argv);
void GetM12WindowInfo(const string& path, Cut & cut);
vector<Cut> GetCutFlowInformation(const std::string& cut_flow_path, const vector<pair<string,string> >& cuts_names);
void ReadCutFlowOutput(const std::string& cut_flow_path, Cut & cut);
void drawSignalEfficiency(const string& output, const std::string& cut_flow_path,  const vector<pair<string,string> >& cuts_names);
void drawSignalEfficiency(const string& output, const float& window);
void drawSignalEfficiency(const string& output);
TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax);
void drawTextBox(const string& position = "top,right");
void drawTextBox(const double& x, const double& y);

pTH2D GetFrame(const double& xmin = 200, const double& xmax = 1399.99, const double& ymin = 0.0001, const double& ymax = 0.02);

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
vector<int> nice_colors {10013,10007,10011,10003,10004,5};
vector<int> nice_linestyles {3,2,1,6,4,10};

//vector<int> nice_colors {kBlack,kRed,kGreen+1,kOrange-3,kMagenta+1,kAzure-4};
//vector<int> nice_colors {10002,10001,10000,10003,10004,5};
//vector<int> nice_linestyles {3,2,1,6,4,10};

int main(int argc, char** argv){
	PublicationStatus status = mssmhbb::publication_status;
	if(status == PRELIMINARY) status = PRELIMINARY_SIMULATION;
	else if (status == PUBLIC) status = SIMULATION;
	style.setTDRstyle(status);
	string output_folder 	= mssmhbb::pictures_output;
	string output_name 		= "PAS_signal_efficiency";

	string cut_flow_path		= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/submit_cutFlow4/";
	auto uinput = ProcessUserInput(argc,argv);
	vector<pair<string,string> > cuts = {

//			{"LoosID","LoosID"},
//
//			{"pt1","p_T^{1}"},
//			{"pt2","p_T^{2}"},
//			{"pt3","p_T^{3}"},
//
//			{"eta1 ","#eta^{1}"},
//			{"eta2","#eta^{2}"},
//			{"eta3","#eta^{3}"},
//
//			{"dR12","Offline kinematic selection"},
//			{"dR23","Offline kinematic selection"},
//			{"dR13","Offline kinematic selection"},
//
//			{"TriggerMatching","Online selection"},
//
//			{"btag1","Offline b-tag selection"},
//			{"btag2","Offline b-tag selection"},
//			{"btag3","Offline b-tag selection"},



	//		{"pt3","p_T^{3}"},
	//		{"eta3","#eta^{3}"},
			{"dR13","Offline kinematic selection"},
			{"TriggerMatching","Online selection"},
			{"btag3","Offline b tagging selection"},
//			{"M12_window","Sub-range mass window"}
//			{"LeadingJetSelection","Full"}
	};

	if( uinput->at("cut_flow").as<bool>() ) drawSignalEfficiency(output_folder + output_name + "_cutflow",cut_flow_path,cuts);
	else if ( uinput->at("standard").as<bool>() ) drawSignalEfficiency(output_folder + output_name + "_standard",0.2);
	else drawSignalEfficiency(output_folder + output_name + "_sub_ranges");
	return 0;
}

void drawSignalEfficiency(const string& output, const std::string& cut_flow_path,  const vector<pair<string,string> >& cuts_names){
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	auto frame(GetFrame(200,1399.99,0.001,7));
	frame->Draw();
	//Create TGraphs
	vector<TGraphErrors*> graphs;

	// Create TLegend
	auto *leg = HbbStyle::legend("top,right",cuts_names.size()+1,0.7);
	leg->SetX1(leg->GetX1()*0.9);
	int i = 0, j = 0;
	TGraphErrors *gr1 = nullptr, *gr2 = nullptr, *gr3 = nullptr;
	for(const auto& cut_name : cuts_names){
		TGraphErrors *gr = new TGraphErrors(mssmhbb::signal_templates.size());
		if(cut_name.first == "M12_window"){
			gr1 = new TGraphErrors(mssmhbb::sr1.size());
			gr1->SetLineColor(nice_colors.at(j));
			gr1->SetMarkerColor(nice_colors.at(j));
			gr1->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
//			gr1->SetLineColor(j+1);
//			gr1->SetMarkerColor(j+1);
//			gr1->SetMarkerStyle(24); // NOT USED
			gr2 = new TGraphErrors(mssmhbb::sr2.size());
			gr2->SetLineColor(nice_colors.at(j));
                        gr2->SetMarkerColor(nice_colors.at(j));
			gr2->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
//			gr2->SetLineColor(j+1);
//			gr2->SetMarkerColor(j+1);
			gr2->SetMarkerStyle(21);
//			gr2->SetMarkerStyle(25); // NOT USED
			gr3 = new TGraphErrors(mssmhbb::sr3.size());
			gr3->SetLineColor(nice_colors.at(j));
			gr3->SetLineStyle(nice_linestyles.at(j));
            gr3->SetMarkerColor(nice_colors.at(j));
			gr3->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
//			gr3->SetLineColor(j+1);
//			gr3->SetMarkerColor(j+1);
			gr3->SetMarkerStyle(22);
//			gr3->SetMarkerStyle(26); // NOT USED 
		}
		gr->SetLineWidth(2);
		gr->SetLineColor(nice_colors.at(j));
		gr->SetLineStyle(nice_linestyles.at(j));
        gr->SetMarkerColor(nice_colors.at(j));
		gr->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
		gr->SetMarkerStyle(20-j+int(cuts_names.size()) - 1);
//		gr->SetLineColor(j+1);
//		gr->SetMarkerColor(j+1);
                Cut cut(cut_name.first,cut_name.second);
		i=0; ++j;
		for(const auto& sample : mssmhbb::signal_templates){
//			if(sample.first != 350) continue;
			cout<<"M_A = "<<sample.first<<endl;
			string full_cut_flow_path = cut_flow_path + "SUSYGluGluToBBHToBB_NarrowWidth_M-" + to_string(sample.first) + "_TuneCUETP8M1_13TeV-pythia8.o";
//			string full_cut_flow_path = cut_flow_path + "SUSYGluGluToBBHToBB_M-350_cfg_ntuple.o";
			if(cut_name.first != "M12_window") ReadCutFlowOutput(full_cut_flow_path,cut);
			else GetM12WindowInfo(sample.second,cut);
                	cut.show();
			double ntot             = 35673; //From luminosity
			gr->SetPoint(i,sample.first,cut.weighted_value/ntot);
			gr->SetPointError(i,0,0);

			if(cut_name.first == "M12_window"){
				if(find(mssmhbb::sr1.begin(),mssmhbb::sr1.end(),sample.first) != mssmhbb::sr1.end()) gr1->SetPoint(i,sample.first,cut.weighted_value/ntot);
				else if (find(mssmhbb::sr2.begin(),mssmhbb::sr2.end(),sample.first) != mssmhbb::sr2.end()) gr2->SetPoint(i-int(mssmhbb::sr1.size()),sample.first,cut.weighted_value/ntot);
				else gr3->SetPoint(i-(int(mssmhbb::sr1.size()) + int(mssmhbb::sr2.size())),sample.first,cut.weighted_value/ntot);

/*
				//Make sub-range1
				for(auto i = 0; i < gr1->GetN(); ++i){
					//Get file :
					int mass = 0;
					string file_name = "";
					if(i == gr1->GetN() - 1) {
						file_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_SR1_NLO_correction_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root";
						mass = 500;
					}
					else {
						mass = mssmhbb::sr1[i];
						file_name = mssmhbb::signal_templates.at(mass);
					}
					TFile *f = new TFile(file_name.c_str());
					auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
					double ntot 		= 35673; //From luminosity
					double e_selected = 0;
					double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
					gr1->SetPoint(i,mass,selected/ntot);
					gr1->SetPointError(i,0,0);
					std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
				}

				//Make sub-range2
				for(auto i = 0; i < gr2->GetN(); ++i){
					//Get file :
					int mass = -1000;;
					string file_name = "";
					if(i == gr2->GetN() - 1) {
						file_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_SR2_NLO_correction_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root";
						mass = 1100;
					}
					else {
						mass = mssmhbb::sr2[i];
						file_name = mssmhbb::signal_templates.at(mass);
					}
					TFile *f = new TFile(file_name.c_str());
					auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
					double ntot 		= 35673; //From luminosity
					double e_selected = 0;
					double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
					gr2->SetPoint(i,mass,selected/ntot);
					gr2->SetPointError(i,0,0);
					std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
				}

				//Make sub-range3
				for(auto i = 0; i < gr3->GetN(); ++i){
					//Get file :
					int mass = mssmhbb::sr3[i];
					string file_name = mssmhbb::signal_templates.at(mass);
					TFile *f = new TFile(file_name.c_str());
					auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
					double ntot 		= 35673; //From luminosity
					double e_selected = 0;
					double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
					gr3->SetPoint(i,mass,selected/ntot);
					gr3->SetPointError(i,0,0);
					std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
				}

*/
			}
			cout<<"Cut: "<<cut_name.first<<" Efficiency: "<<cut.weighted_value/ntot<<endl;
			++i;
		}
		string legend_title = cut.legend_title.c_str();
		if(j != 1) legend_title = "+ " + legend_title;
		graphs.push_back(gr);

		if(cut_name.first == "M12_window"){
			gr1->Draw("PL same");
			gr2->Draw("PL same");
			gr3->Draw("PL same");
			leg->AddEntry(gr1,legend_title.c_str(),"lp");
		} else {
			gr->Draw("PL same");
			leg->AddEntry(gr,legend_title.c_str(),"lp");
		}
        }
	gPad->SetLogy();
	leg->Draw();
	drawTextBox(250,1.5e-03);
	HbbStyle::drawStandardTitle();
//	drawTextBox(250,3);
//	HbbStyle::drawStandardTitle("out");
	can.Print( (output + ".pdf").c_str());
}

void drawSignalEfficiency(const string& output, const float& window){
//	gStyle->SetPadRightMargin(0.04);
//	gStyle->SetPadLeftMargin(0.14);
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	auto frame(GetFrame());
	//Create TGraph
	TGraphErrors *gr = new TGraphErrors(mssmhbb::signal_templates.size());
	gr->SetMarkerColor(kRed);
	gr->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
	//Iterate through the signal samples
	int i = 0;
	for(const auto& sample : mssmhbb::signal_templates){
//		double lower_border = sample.first * ( 1 - window/2.);
//		double upper_border = sample.first * ( 1 + window/2.);
		TFile *f = new TFile(sample.second.c_str());
//		auto* h_denum 	= GetFromTFile(sample.second,"distributions/NumberOfGenEvents_afterMHat_rewPU");
		auto* h_num		= GetFromTFile<TH1>(*f,"templates/bbH_Mbb_VIS");
		h_num->Sumw2(true);
		double ntot 		= 35673; //From luminosity
		double e_selected = 0;
		double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
		gr->SetPoint(i,sample.first,selected/ntot);
		gr->SetPointError(i,0,0);
		std::cout<<"M_{A/H} = "<<sample.first<<" eff = "<<selected/ntot<<std::endl;
		++i;
	}
//	gPad->SetLogy();
	frame->Draw();
	gr->Draw("PL same");
	drawTextBox();
	HbbStyle::drawStandardTitle();
//	HbbStyle::drawStandardTitle("out");
//	drawTextBox();
	can.Print( (output + ".pdf").c_str());
}

void drawSignalEfficiency(const string& output){
	/*
	 * Draw signal efficiency with overlaping sub-range points
	 * but w/o cut flow
	 */
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	auto frame(GetFrame(200,1399.99,4e-03,0.023));
	//TLegend
	auto *leg = HbbStyle::legend("top,right",4);
	//Create TGraph
	TGraphErrors *gr1 = new TGraphErrors(int(mssmhbb::sr1.size()) + 1);
	TGraphErrors *gr2 = new TGraphErrors(int(mssmhbb::sr2.size()) + 1);
	gr2->SetMarkerStyle(22);
	gr2->SetMarkerSize(gStyle->GetMarkerSize()*1.1);
	gr2->SetLineColor(2);
	TGraphErrors *gr3 = new TGraphErrors(int(mssmhbb::sr3.size()));
	gr3->SetMarkerStyle(34);
	gr3->SetLineColor(kBlue);
	leg->AddEntry(gr1,"sub-range 1","pl");
	leg->AddEntry(gr2,"sub-range 2","pl");
	leg->AddEntry(gr3,"sub-range 3","pl");

	//Make sub-range1
	for(auto i = 0; i < gr1->GetN(); ++i){
		//Get file :
		int mass = 0;
		string file_name = "";
		if(i == gr1->GetN() - 1) {
			file_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_SR1_NLO_correction_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root";
			mass = 500;
		}
		else {
			mass = mssmhbb::sr1[i];
			file_name = mssmhbb::signal_templates.at(mass);
		}
		TFile *f = new TFile(file_name.c_str());
		auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
		double ntot 		= 35673; //From luminosity
		double e_selected = 0;
		double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
		gr1->SetPoint(i,mass,selected/ntot);
		gr1->SetPointError(i,0,0);
		std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
	}

	//Make sub-range2
	for(auto i = 0; i < gr2->GetN(); ++i){
		//Get file :
		int mass = -1000;;
		string file_name = "";
		if(i == gr2->GetN() - 1) {
			file_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_SR2_NLO_correction_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root";
			mass = 1100;
		}
		else {
			mass = mssmhbb::sr2[i];
			file_name = mssmhbb::signal_templates.at(mass);
		}
		TFile *f = new TFile(file_name.c_str());
		auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
		double ntot 		= 35673; //From luminosity
		double e_selected = 0;
		double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
		gr2->SetPoint(i,mass,selected/ntot);
		gr2->SetPointError(i,0,0);
		std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
	}

	//Make sub-range3
	for(auto i = 0; i < gr3->GetN(); ++i){
		//Get file :
		int mass = mssmhbb::sr3[i];
		string file_name = mssmhbb::signal_templates.at(mass);
		TFile *f = new TFile(file_name.c_str());
		auto* h_num		= GetFromTFile<TH1>(*f,"bbH_Mbb");
		double ntot 		= 35673; //From luminosity
		double e_selected = 0;
		double selected 	= h_num->IntegralAndError(1,h_num->GetNbinsX(),e_selected);
		gr3->SetPoint(i,mass,selected/ntot);
		gr3->SetPointError(i,0,0);
		std::cout<<"M_{A/H} = "<<mass<<" eff = "<<selected/ntot<<std::endl;
	}

	gPad->SetLogy();
	frame->Draw();
	gr1->Draw("PL same");
	gr2->Draw("PL same");
	gr3->Draw("PL same");
	leg->Draw();
	HbbStyle::drawStandardTitle("out");
	drawTextBox("top,left");
	can.Print( (output + ".pdf").c_str());
}

TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax){
	TLegend *leg = new TLegend(xmin,ymin,xmax,ymax);
	style.setLegendStyle(leg);
	return leg;
}

void drawTextBox(const string& position){
	/*
	 * Draw TLatex textbox with a sign
	 */
	  double x,y = 0.018;
	  if(position.find("top") != string::npos) y = 0.018;
	  if(position.find("right") != string::npos) x = 1100;
	  else x = 250;
	  drawTextBox(x, y);
}

void drawTextBox(const double& x, const double& y){
	  TLatex latex;
	  latex.SetTextColor(kBlue+2);
	  latex.SetTextSize(latex.GetTextSize()*1.2);
	  latex.DrawLatex(x, y,
	                     ("A/H #rightarrow b#bar{b}"));
}

vector<Cut> GetCutFlowInformation(const std::string& cut_flow_path, const vector<pair<string,string> >& cuts_names){
	/*
	 * Function to read info from the cut-flow .o ascii files
	 *
	 * return vector of "Cut" objects filled from the cut-flow tables
	 */
	vector<Cut> cuts;
	for(const auto& cut_name : cuts_names){
		Cut cut(cut_name.first,cut_name.second);
		ReadCutFlowOutput(cut_flow_path,cut);
		cut.show();
		cuts.push_back(cut);
	}
	return cuts;
}

void ReadCutFlowOutput(const std::string& cut_flow_path, Cut & cut){
	/*
	 * Read cut-flow ascii file into Cut object
	 */
	ifstream file(cut_flow_path.c_str());
	string line;
	float lumi_weight = 0, tot_evs = 0, weighted_evs = 0;
	while(getline(file,line)){
		//line-by-line parsing
		istringstream iss(line);
		string str_number = "", rubbish = "";
		//Get lumi:
		if(line.find("LUMI") != string::npos){
			iss >> rubbish >> str_number >> rubbish >> rubbish;
			lumi_weight = stof(str_number);
		}
		if(line.find(cut.name) != string::npos){
			//found this cut
			string substring ="";
			regex base_regex(R"(\b[0-9]+)");
			smatch base_match;
			if(regex_search(line,base_match,base_regex)) {
				substring = base_match[0].str();
			}
			cut.nevents = stof(substring);
		}
		//find selected weight and last number of events
		if(line.find("# selected events:") != string::npos){
			string substring ="";
                        regex base_regex(R"(\b[0-9]+)");
                        smatch base_match;
			if(regex_search(line,base_match,base_regex)) tot_evs = stof( base_match[0].str() );
		}
		//found weighted events
		if(line.find("# selected weighted") != string::npos){
			string substring ="";
                        regex base_regex(R"(\b[0-9]+)");
                        smatch base_match;
                        if(regex_search(line,base_match,base_regex)) weighted_evs = stof( base_match[0].str() );
			cut.weight = lumi_weight * weighted_evs / tot_evs;
			cut.weighted_value = cut.nevents * cut.weight;
			return;
		}
	}
}

pVariables_map ProcessUserInput(int argc, char** argv){
	/*
	 * Function to process user input for this macro.
	 */

	bool cut_flow = false, sub_range_devision = false, standard = false;
	// general options to show help menu
	// and control output verbosity
	options_description generalOptions("General options");
	generalOptions.add_options()
		("help,h", "produce help message")
		("verbose,v", value<int>()->default_value(0), "more verbose output")
		("cut_flow", bool_switch(&cut_flow), "Show cut flow efficiency")
		("sub_range_devision", bool_switch(&sub_range_devision), "Show an effect of the sub-range devision")
		("standard", bool_switch(&standard), "Standard efficiency curve");

    /*
     * Join options together
     */

    pVariables_map output_vm(new variables_map());
    store(command_line_parser(argc,argv).options(generalOptions).allow_unregistered().run(),*output_vm);
    notify(*output_vm);

    //show help menu
    if (output_vm->count("help")) {
	    cout << generalOptions << endl;
	    exit(0);
    }

    //if not pre-def. options
    store(parse_command_line(argc,argv,generalOptions),*output_vm);
    notify(*output_vm);
    if(!cut_flow && !sub_range_devision && !standard){
    		//At least one of the options has to be specified!
    		cout<<"\nWARNING: non of the options were specified. Go to default 'cut_flow' method\n"<<endl;
    		cut_flow = true;
    		output_vm->at("cut_flow") = variable_value(bool(cut_flow),false);
    }

    return output_vm;
}

pTH2D GetFrame(const double& xmin, const double& xmax, const double& ymin, const double& ymax ){
	/*
	 * Setup frame to be used for the efficiency
	 * calculation
	 */
//	TH2D *frame = new TH2D("frame","frame",1,200,1399.99,1,0.0001,0.02);
	unique_ptr<TH2D> frame(new TH2D("frame","frame",1,xmin,xmax,1,ymin,ymax));
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMAH());
	frame->GetYaxis()->SetTitle("Efficiency");
	frame->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X")*1.02);
	frame->GetXaxis()->SetTitleOffset(frame->GetXaxis()->GetTitleOffset()*1.05); // private preferences
	frame->SetMinimum(1e-06);
	return frame;
}

void GetM12WindowInfo(const string& path, Cut & cut){
	TFile f(path.c_str());
	auto* h_num             = GetFromTFile<TH1>(f,"bbH_Mbb");
	cut.weighted_value = h_num -> Integral();
}
