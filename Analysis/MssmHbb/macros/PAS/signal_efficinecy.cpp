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

#include "TCanvas.h"
#include "TGraphErrors.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
HbbStyle style;

using namespace std;

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
void GetM12WindowInfo(const string& path, Cut & cut);
vector<Cut> GetCutFlowInformation(const std::string& cut_flow_path, const vector<pair<string,string> >& cuts_names);
void ReadCutFlowOutput(const std::string& cut_flow_path, Cut & cut);
void drawSignalEfficiency(const string& output, const std::string& cut_flow_path,  const vector<pair<string,string> >& cuts_names);
void drawSignalEfficiency(const string& output, const float& window);
TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax);
void drawTextBox();

int main(int argc, char** argv){
	style.setTDRstyle(PRELIMINARY_SIMULATION);
	string output_folder 	= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/testLO";
	string output_name 		= "PAS_signal_efficiency";

	string cut_flow_path		= mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/bin/submit_cutFlow3/";
	vector<pair<string,string> > cuts = {

			{"LoosID","LoosID"},

			{"pt1","p_T^{1}"},
			{"pt2","p_T^{2}"},
			{"pt3","p_T^{3}"},

			{"eta1 ","#eta^{1}"},
			{"eta2","#eta^{2}"},
			{"eta3","#eta^{3}"},

			{"dR12","Offline kinematic selection"},
			{"dR23","Offline kinematic selection"},
			{"dR13","Offline kinematic selection"},

			{"TriggerMatching","Online selection"},

			{"btag1","Offline b-tag selection"},
			{"btag2","Offline b-tag selection"},
			{"btag3","Offline b-tag selection"},



	//		{"pt3","p_T^{3}"},
	//		{"eta3","#eta^{3}"},
//			{"dR13","Offline kinematic selection"},
//			{"TriggerMatching","Online selection"},
//			{"btag3","Offline b-tag selection"},
//			{"M12_window","Sub-range mass window"}
//			{"LeadingJetSelection","Full"}
	};

	drawSignalEfficiency(output_folder + output_name,cut_flow_path,cuts);
	return 0;
}

void drawSignalEfficiency(const string& output, const std::string& cut_flow_path,  const vector<pair<string,string> >& cuts_names){
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	TH2D *frame = new TH2D("frame","frame",1,200,1399.99,1,0.001,100.);
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMAH());
	frame->GetYaxis()->SetTitle("#epsilon");
	frame->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X")*1.02);
	frame->SetMinimum(1e-06);
	frame->GetXaxis()->SetRangeUser(200,1399.99);
	frame->Draw();
	//Create TGraphs
	vector<TGraphErrors*> graphs;

	// Create TLegend
	auto *leg = HbbStyle::legend("top,right",cuts_names.size(),0.7);
	int i = 0, j = 0;
	TGraphErrors *gr1 = nullptr, *gr2 = nullptr, *gr3 = nullptr;
	for(const auto& cut_name : cuts_names){
		TGraphErrors *gr = new TGraphErrors(mssmhbb::signal_templates.size());
		if(cut_name.first == "M12_window"){
			gr1 = new TGraphErrors(mssmhbb::sr1.size());
			gr1->SetLineColor(j+1); gr1->SetMarkerColor(j+1);
			gr2 = new TGraphErrors(mssmhbb::sr2.size());
			gr2->SetLineColor(j+1); gr2->SetMarkerColor(j+1); gr2->SetMarkerStyle(21);
			gr3 = new TGraphErrors(mssmhbb::sr3.size());
			gr3->SetLineColor(j+1); gr3->SetMarkerColor(j+1); gr3->SetMarkerStyle(22);
		}
		gr->SetLineColor(j+1);
		gr->SetMarkerColor(j+1);
                Cut cut(cut_name.first,cut_name.second);
		i=0; ++j;
		for(const auto& sample : mssmhbb::signal_templates){
			if(sample.first != 350) continue;
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
				else if (find(mssmhbb::sr2.begin(),mssmhbb::sr2.end(),sample.first) != mssmhbb::sr2.end()) gr2->SetPoint(i,sample.first,cut.weighted_value/ntot);
				else gr3->SetPoint(i,sample.first,cut.weighted_value/ntot);
			}
			cout<<"Cut: "<<cut_name.first<<" Efficiency: "<<cut.weighted_value/ntot<<endl;
			++i;
		}
		
		graphs.push_back(gr);
		if(cut_name.first == "M12_window"){
			gr1->Draw("PL same");
			gr2->Draw("PL same");
			gr3->Draw("PL same");
			leg->AddEntry(gr1,cut.legend_title.c_str(),"lp");
		} else {
			gr->Draw("PL same");
			leg->AddEntry(gr,cut.legend_title.c_str(),"lp");
		}
        }
	gPad->SetLogy();
	leg->Draw();
	HbbStyle::drawStandardTitle("out");
	can.Print( (output + "lol.pdf").c_str());
}

void GetM12WindowInfo(const string& path, Cut & cut){
	TFile f(path.c_str());
	auto* h_num             = GetFromTFile<TH1>(f,"bbH_Mbb");	
	cut.weighted_value = h_num -> Integral();
}

void drawSignalEfficiency(const string& output, const float& window){
//	gStyle->SetPadRightMargin(0.04);
//	gStyle->SetPadLeftMargin(0.14);
	//Create TCanvas
	TCanvas can("can","can",800,600);
	//Create TH1 as a frame
	TH2D *frame = new TH2D("frame","frame",1,200,1399.99,1,0.0001,0.02);
	frame->GetXaxis()->SetTitle(HbbStyle::axisTitleMAH());
	frame->GetYaxis()->SetTitle("#epsilon");
	frame->GetXaxis()->SetLabelOffset(gStyle->GetLabelOffset("X")*1.02);
	frame->SetMinimum(1e-06);
	frame->GetXaxis()->SetRangeUser(200,1399.99);
	//Create TGraph
	TGraphErrors *gr = new TGraphErrors(mssmhbb::signal_templates.size());
	gr->SetMarkerColor(kRed);
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
	HbbStyle::drawStandardTitle("out");
	drawTextBox();
	can.Print( (output + ".pdf").c_str());
}

TLegend * createTLegend(const double& xmin, const double& ymin, const double& xmax, const double& ymax){
	TLegend *leg = new TLegend(xmin,ymin,xmax,ymax);
	style.setLegendStyle(leg);
	return leg;
}

void drawTextBox(){
	/*
	 * Draw TLatex textbox with a sign
	 */
	  TLatex latex;
	  latex.SetTextColor(kBlue+2);
	  latex.DrawLatex(1100, 0.018,
	                     ("A/H #rightarrow b#bar{b}"));
//	  latex.DrawLatexNDC(0.56, 0.9,
//	                     ("A/H #rightarrow b#bar{b}"));
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
