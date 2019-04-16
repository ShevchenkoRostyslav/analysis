/*
 * Macro to show effect of each MC reweighting
 * correction.
 */

#include <string>
#include <vector>
#include <iostream>

#include "TH1.h"

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;

void CheckCorrections(const vector<string> &corrections, const string& histo_name, const string& input_path, const string& output_path);
void PlotTheCorrectionEffect(TH1 * hcentral, TH1 * hcorrection, const string& correction, const string& mass, const string& output_path);
string ExtractCorrectionName(const string& correction);
void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator);

int main(){
	string output_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/";
	string input_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/";
	vector<string> corrections = {"no_PU_weight","no_OnlineBTagEff_weight","no_PTEff_weight","no_SFbSFl_weight"};
	string histo_name = "bbH_Mbb";
	CheckOutputDir(output_path);
	CheckCorrections(corrections,histo_name,input_path,output_path);
}

void CheckCorrections(const vector<string> &corrections, const string& histo_name, const string& input_path, const string& output_path){
	for(const auto& mass : mssmhbb::masses){
		string smass = to_string(mass);

		auto central_file_name = mssmhbb::signal_templates.at(mass);
		TFile f_central((central_file_name).c_str(),"READ");
		auto *h_central = GetFromTFile<TH1F>(f_central,histo_name);

		for(const auto& correction : corrections){
			string correction_file_name = "MssmHbbSignal_" + correction + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-" + smass + "_TuneCUETP8M1_13TeV-pythia8.root";

			TFile f_correction((input_path + correction_file_name).c_str(),"READ");
			auto *h_correction = GetFromTFile<TH1F>(f_correction,histo_name);

			PlotTheCorrectionEffect(h_central,h_correction,correction,smass,output_path);
		}
	}
}

void PlotTheCorrectionEffect(TH1 * hcentral, TH1 * hcorrection, const string& correction, const string& mass, const string& output_path){
	HbbStyle style;
	style.setTDRstyle(PRIVATE);

	auto correction_name = ExtractCorrectionName(correction);

	TCanvas can;

	//top pad and frame
	auto topPad = HbbStyle::getRatioTopPad(0.7);
	topPad->Draw();
	topPad->cd();

	auto *top_frame = new TH2D("top_frame","",1,hcentral->GetXaxis()->GetXmin(),hcentral->GetXaxis()->GetXmax(),1,hcentral->GetYaxis()->GetXmin(),1.2*hcentral->GetMaximum());
	style.setRatioTopFrame(top_frame, can.YtoPixel(can.GetY1()));
	top_frame->SetTitle((correction_name +";" + style.axisTitleMass() + ";Events"));
	top_frame->Draw();

	style.applyStandardMCHistoStyle(*hcentral);
	hcentral->SetFillColor(kRed-10);
	style.applyStandardDataHistoStyle(*hcorrection);
	hcentral->Draw("hist same");
	hcorrection->Draw("E same");

	string leg_position = "top,right";
	bool shift_leg_toleft = false;
	if(mass == "1300" || mass == "900") {
		shift_leg_toleft = true;
		leg_position = "top,left";
	}
	TLegend &leg = *HbbStyle::legend(leg_position,3);
	leg.SetY1(leg.GetY1()*0.65);
	leg.SetY2(leg.GetY2()*0.95);
	if(!shift_leg_toleft){
		leg.SetX1(leg.GetX1()*0.9);
	}
	leg.SetHeader(("#font[12]{#splitline{" + mass + " GeV,}{" + correction_name + "}}").c_str());
	leg.AddEntry(hcentral,"Corrected","f");
	leg.AddEntry(hcorrection,"w/o Correction","p");
	leg.Draw();

	gPad->RedrawAxis();
	can.cd();

	auto botPad = HbbStyle::getRatioBottomPad(0.3);
	botPad->Draw();
	botPad->cd();

	float up_bound = 20, dn_bound = -20;
	if((mass == "700" || mass == "900" || mass == "1100" || mass == "1300") && (correction == "no_OnlineBTagEff_weight")){
		up_bound = 5; dn_bound = -35;
	}
//	auto *botFrame = new TH2D("botFrame","",1,top_frame->GetXaxis()->GetXmin(),top_frame->GetXaxis()->GetXmax(),1,-3.9999,3.9999);
	auto *botFrame = new TH2D("botFrame","",1,top_frame->GetXaxis()->GetXmin(),top_frame->GetXaxis()->GetXmax(),1,dn_bound,up_bound);
	HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
	botFrame->GetXaxis()->SetTitle(top_frame->GetXaxis()->GetTitle());
//	botFrame->GetYaxis()->SetTitle("#frac{Corr.-W/O Corr.}{#sigma_{W/O Corr.}}");
	botFrame->GetYaxis()->SetTitle("#frac{Cor.-w/o Cor.}{Cor.},%");
	botPad->cd();
	auto ratio = static_cast<TH1D*>(hcentral->Clone("ratio"));
	ratio->Sumw2();

	//Minus and Devision
	devideTwoHistograms(ratio,hcentral,hcorrection);
	botFrame->Draw();
	ratio->Draw("Esame");

	TLine *horizLine = new TLine(botFrame->GetXaxis()->GetXmin(),0,botFrame->GetXaxis()->GetXmax(),0);
	horizLine -> SetLineStyle(2);
	horizLine -> Draw();

	can.cd();
	style.drawStandardTitle("out");

	can.SaveAs((output_path + mass + "_" + correction + ".pdf").c_str());

}

string ExtractCorrectionName(const string& correction){
	/*
	 * Hardcoded names for the plots
	 */
	if(correction == "no_PU_weight") return "Pile-up correction";
	else if (correction == "no_OnlineBTagEff_weight") return "Online b-tag correction";
	else if (correction == "no_PTEff_weight") return "Kinematic trigger correction";
	else if (correction == "no_SFbSFl_weight") return " Offline b-tag correction";
	else throw logic_error("WRONG Correction name in the list");
}

void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator){
	/*
	 * Function to perform:
	 * f = (a - b) / sigma(b)
	 * assume full correlations
	 */
	double difference, edifference, devis, edevis, bbb, ebbb, bbnb, ebbnb;
	for(int i = 1; i <= nominator->GetNbinsX(); ++i){
		bbb  = nominator->GetBinContent(i);  ebbb  = nominator->GetBinError(i);
		bbnb = denominator->GetBinContent(i); ebbnb = denominator->GetBinError(i);
		difference  = bbb - bbnb;
		//Assume the full correlations
		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb - 2 * ebbb * ebbnb );
		//Assume 0 correlation
		devis 	= difference / bbb;
		//assume the full correlations
		edevis 	= abs(devis) * sqrt ( pow(edifference/difference,2) + pow(ebbb/bbb,2) - 2 * edifference * ebbb / (difference * bbb) );
		//Assume the full correlations
		if(devis != devis) { devis = 0; edevis = 0;};
		ratio->SetBinContent(i,devis*100);
		ratio->SetBinError(i,edevis*100);
	}
}
