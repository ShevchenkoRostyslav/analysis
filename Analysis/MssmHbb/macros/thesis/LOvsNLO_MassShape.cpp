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

void CheckShapes( const string& histo_name, const string& output_path, const vector<int>& mass_ps = mssmhbb::masses);
void PlotTheCorrectionEffect(TH1 * hcentral, TH1 * hcorrection, const string& mass, const string& output_path);
std::string GetHardcodedNLOFilePath(const int & mp);
void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator);

int main(){
	string output_path = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/Thesis/";
//	string histo_name = "templates/bbH_Mbb_VIS";
	string histo_name = "bbH_Mbb";		
	CheckOutputDir(output_path);
	vector<int> masses = {350,600,900};
	CheckShapes(histo_name,output_path,masses);
}

void CheckShapes(const string& histo_name, const string& output_path, const vector<int>& mass_ps){
	for(const auto& mass : mass_ps){
		string smass = to_string(mass);

		auto central_file_name = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_no_LOtoNLO_weight_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-" + smass + "_TuneCUETP8M1_13TeV-pythia8.root";
		TFile f_central((central_file_name).c_str(),"READ");
		auto *h_central = GetFromTFile<TH1F>(f_central,histo_name);

		string correction_file_name = GetHardcodedNLOFilePath(mass);
 
		TFile f_correction((correction_file_name).c_str(),"READ");
		auto *h_correction = GetFromTFile<TH1F>(f_correction,histo_name);

		PlotTheCorrectionEffect(h_central,h_correction,smass,output_path);
		}
}

void PlotTheCorrectionEffect(TH1 * hcentral, TH1 * hcorrection, const string& mass, const string& output_path){
	HbbStyle style;
	style.setTDRstyle(PRIVATE);

	string correction_name = "NLO vs LO";
	string correction = "NLOvsLO";

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
	double sf = 1.;
	double sfy = 1.;
	if (mass == "900") {sf = 0.7; sfy = 1.13;}
	leg.SetY1(leg.GetY1()*0.65*sfy);
	leg.SetY2(leg.GetY2()*0.95);
	leg.SetX2(leg.GetX2()*sf);
	if(!shift_leg_toleft){
		leg.SetX1(leg.GetX1()*0.9);
	}
	leg.SetHeader(("#font[12]{#splitline{" + mass + " GeV,}{" + correction_name + "}}").c_str());
	leg.AddEntry(hcentral,"LO Pythia 8","f");
	leg.AddEntry(hcorrection,"NLO Madgraph 5","pl");
	leg.Draw();

	gPad->RedrawAxis();
	can.cd();

	auto botPad = HbbStyle::getRatioBottomPad(0.3);
	botPad->Draw();
	botPad->cd();

	float up_bound = 3, dn_bound = -3;
	auto *botFrame = new TH2D("botFrame","",1,top_frame->GetXaxis()->GetXmin(),top_frame->GetXaxis()->GetXmax(),1,dn_bound,up_bound);
	HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
	botFrame->GetXaxis()->SetTitle(top_frame->GetXaxis()->GetTitle());
	botFrame->GetYaxis()->SetTitle("#frac{LO-NLO}{#sigma(NLO)}");
	botPad->cd();
	auto ratio = static_cast<TH1D*>(hcentral->Clone("ratio"));
	ratio->Sumw2();

	//Minus and Devision
	devideTwoHistograms(ratio,hcentral,hcorrection);
	std::cout<<"Normalisation: m = "<<mass<<" (LO - NLO)/LO*100 = "<<(hcentral->Integral()-hcorrection->Integral())/hcentral->Integral() * 100<<std::endl;
	botFrame->Draw();
	ratio->Draw("Esame");

	TLine *horizLine = new TLine(botFrame->GetXaxis()->GetXmin(),0,botFrame->GetXaxis()->GetXmax(),0);
	horizLine -> SetLineStyle(2);
	horizLine -> Draw();

	can.cd();
	style.drawStandardTitle("out");

	can.SaveAs((output_path + mass + "_" + correction + ".pdf").c_str());

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
		//Assume no correlations
		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb);
		devis 	= difference / ebbnb;
		edevis  = edifference / ebbnb;
		//assume no correlations
		//edevis 	= abs(devis) * sqrt ( pow(edifference/difference,2) + pow(ebbnb/bbnb,2) - 2 * edifference * ebbnb / (difference * bbnb) );
		//Assume the full correlations
		if(devis != devis) { devis = 0; edevis = 0;};
		ratio->SetBinContent(i,devis);
		ratio->SetBinError(i,edevis);
	}
}

std::string GetHardcodedNLOFilePath(const int & mp){
	/*
	 * Function to get path to the NLO file 
	*/
	std::string smass = std::to_string(mp);
	std::string output_n = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/";
	if(mp == 350 || mp == 600 || mp == 900){
		output_n += "MssmHbbSignal_NLO_my_lowM_SUSYGluGluToBBHToBB_M-" + smass + "_cfg_ntuple.root";
	}
	else if (mp == 400 || mp == 1200){
		output_n += "MssmHbbSignal_NLO_my_lowM_rshevche-SUSYGluGluToBBHToBB_M-" + smass + "_cfg_GEN_DIGI80X_RECO80X_MiniAODv2_80X-8ce03fc302488195777eb8e1e985030a.root";
	}
	else if (mp == 1400){
		output_n += "MssmHbbSignal_NLO_my_lowM_chayanit-SUSYGluGluToBBHToBB_M-" + smass + "_extend_cfg_GEN_DIGI80X_RECO80X_MiniAODv2_80X-28028af67189b3de7224b79195bd0e1d.root";
	}
	else throw std::invalid_argument("Mass point " + smass + " has not been produced in NLO. Check the input mass points");
	if (file_exists(output_n)) return output_n;
	else throw std::invalid_argument("FILE: " + output_n + " doesn't exist. Check the spelling,");
}
