/*
 * signal_templates_variations.cpp
 *
 *  Created on: 14 Apr 2018
 *      Author: rostyslav
 */

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TH1.h"
#include <memory>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TExec.h"
#include "TROOT.h"
#include <string>

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

HbbStyle style;

std::map<int,TFile*> GetTFiles(const std::map<int,std::string>& pathes);
TFile * GetTFile(const std::string& path);
void DrawSystematicVariations(const std::map<int,TFile*>& files);
std::string GetTLegendHeader(const std::string& s);
double correlatedDivision(const double& v1, const double& ev1,const double& v2, const double& ev2){
	// y = v1/v2
	if(v1 == 0 || v2 == 0) return 0;
	double f = v1/v2 * sqrt( pow(ev1/v1,2) + pow(ev2/v2,2) - 2*(ev1*ev2)/(v2*v1) );
	return f;
}

int main(){
	style.setTDRstyle(PRIVATE);

	std::string output_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/signal_templates_variaitons/";

	std::map<int,std::string> mass_points_pathes = mssmhbb::signal_templates;
	auto mass_points_files = GetTFiles(mass_points_pathes);

	CheckOutputDir(output_path);
	DrawSystematicVariations(mass_points_files);
}

void DrawSystematicVariations(const std::map<int,TFile*>& files){
	/*
	 * Main method to draw nice templates
	 */
	std::vector<std::string> list_of_systematics = {"CMS_scale_j","CMS_res_j","CMS_eff_pTonl","CMS_eff_b","CMS_PU","CMS_eff_bonl","CMS_eff_l"};
	for (const auto& f : files){
		//Get the entral histo
		auto central = GetFromTFile<TH1D>(*f.second,"templates/bbH_Mbb_VIS");
		style.applyStandardMCHistoStyle(*central);
		central->SetFillColor(kRed-10);

		//Prepare the frame
		TCanvas can;
		auto topPad = HbbStyle::getRatioTopPad(0.7);
		topPad->Draw();
		topPad->cd();

		auto *top_frame = new TH2D("top_frame","",1,central->GetXaxis()->GetXmin(),central->GetXaxis()->GetXmax(),1,central->GetYaxis()->GetXmin(),1.2*central->GetMaximum());
		style.setRatioTopFrame(top_frame, can.YtoPixel(can.GetY1()));
		top_frame->SetTitle(";" + style.axisTitleMass() + ";Arbitrary units");

		can.cd();
		auto botPad = HbbStyle::getRatioBottomPad(0.3);
		botPad->Draw();
		botPad->cd();

		auto *botFrame = new TH2D("botFrame","",1,top_frame->GetXaxis()->GetXmin(),top_frame->GetXaxis()->GetXmax(),1,0.5,1.5);
		HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
		botFrame->GetXaxis()->SetTitle(top_frame->GetXaxis()->GetTitle());
		botFrame->GetYaxis()->SetTitle("#frac{Down(Up)}{Central}");

		can.cd();
		topPad->cd();

		std::string leg_position = "top,right";
		if(f.first >= 900) {
			leg_position = "top,left";
		}
		auto legend = *style.legend(leg_position,10);
		legend.SetY1(legend.GetY1()*0.9);
		legend.SetY2(legend.GetY2()*0.95);
		//get up and down variations
		TH1D *up,*down;
		for(const auto & s : list_of_systematics){
			topPad->cd();

			up = GetFromTFile<TH1D>(*f.second,"templates/bbH_Mbb_" + s + "_VIS_13TeVUp");
			down = GetFromTFile<TH1D>(*f.second,"templates/bbH_Mbb_" + s + "_VIS_13TeVDown");

			style.applyStandardDataHistoStyle(*up);
			up->SetMarkerColor(kBlue);
			style.applyStandardDataHistoStyle(*down);
			down->SetMarkerColor(kRed);

			top_frame->Draw();
			central->Draw("hist same");
			up->Draw("same");
			down->Draw("same");

			std::string header_name = GetTLegendHeader(s);
			legend.Clear();
			legend.SetHeader(("#font[12]{#splitline{m_{A/H} = " + std::to_string(f.first) + " GeV}{" + header_name + "}}").c_str());
			legend.AddEntry(central,"Central","f");
			legend.AddEntry(down,"-2#sigma variation","pl");
			legend.AddEntry(up,"+2#sigma variation","pl");
			legend.Draw();

			gPad->RedrawAxis();
			can.cd();
			botPad->cd();
			auto rUp = static_cast<TH1D*>(up->Clone("rUp"));
			auto rDown = static_cast<TH1D*>(down->Clone("rDown"));

			rUp->Divide(central);
			rDown->Divide(central);
			//Calculate uncertainties
			for(int i=1;i<central->GetNbinsX();++i){
				double err_u = up->GetBinError(i);
				double u = up->GetBinContent(i);
				double err_d = down->GetBinError(i);
				double d = down->GetBinContent(i);
				double err_c = central->GetBinError(i);
				double c = central->GetBinContent(i);

				//fully correlated errors of the devision
				rUp->SetBinError(i,correlatedDivision(c, err_c, u, err_u));
				rDown->SetBinError(i,correlatedDivision(c, err_c, d, err_d));
			}
			botFrame->Draw();
			rUp->Draw("same");
			rDown->Draw("same");

			TLine *horizLine = new TLine(central->GetXaxis()->GetXmin(),1,central->GetXaxis()->GetXmax(),1);
			horizLine -> SetLineStyle(2);
			horizLine -> Draw();

			can.cd();
			style.drawStandardTitle("out");

			std::string out_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/" + std::to_string(f.first) + "_GeV_" + s + "_variation.pdf";
			can.Print(out_name.c_str());
		}
	}
}

std::map<int,TFile*> GetTFiles(const std::map<int,std::string>& pathes){
	/*
	 * Function to get the TFiles fro the pathes
	 */
	std::map<int,TFile*> out;
	for(const auto & p : pathes){
		auto f = GetTFile(p.second);
		out[p.first] = f;
	}
	return out;
}

TFile * GetTFile(const std::string& path){
	/*
	 * Get a single TFile accorindg to the path
	 */
	TFile *f = nullptr;
	if(file_exists(path)) f = new TFile(path.c_str(),"READ");
	return f;
}

std::string GetTLegendHeader(const std::string& s){
	/*
	 * Return the header for TLegend
	 * for the current systematic unc
	 */
	std::string out = "";
	if(s == "CMS_scale_j") out = "Jet Energy Scale";
	else if (s == "CMS_res_j") out = "Jet Energy Resolution";
	else if (s == "CMS_eff_pTonl") out = "p_{T} trigger eff.";
	else if (s == "CMS_eff_b") out = "Offline b-flavour";
	else if (s == "CMS_PU") out = "Pile-up";
	else if (s == "CMS_eff_bonl") out = "Online b-tag eff.";
	else if (s == "CMS_eff_l") out = "Offline udsg-flavour";

	return out;
}
