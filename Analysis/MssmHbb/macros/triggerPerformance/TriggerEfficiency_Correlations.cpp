/*
 * TriggerEfficiency_Correlations.cpp
 *
 *  Created on: 8 Nov 2017
 *      Author: shevchen
 */

#include "TExec.h"
#include "TStyle.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH1.h"
#include "TEfficiency.h"
#include "TFile.h"

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

struct My2DHisto{
	TH2D * num = nullptr;
	TH2D * den = nullptr;
	TH2D * num_sym = nullptr;
	TH2D * den_sym = nullptr;

	std::string title="";
};

HbbStyle style;
using namespace std;
using namespace mssmhbb;

void CalculateAndPlotCorrelations(My2DHisto & data, My2DHisto& mc, const string& output_path);
TEfficiency * GetSymmetricEfficiency(My2DHisto & data);
TH2D * GetCorrelationCoeff(TEfficiency * eff);
void PlotSymmetric2DEfficiency(TEfficiency * eff, const string& output_path, const bool& data = false);
void Plot2DCorrelations(TH2D * hdata, const string& output_path, const bool& data = false);
void PlotCorrelationSlices(TH2D * hdata, TH2D * hmc, const string& output_path, const bool& data = false);
void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator);

void BayesianErrors(TEfficiency *h);

int main(){
	//Style for the Thesis
	style.set(PRIVATE);
	gStyle->SetOptStat(0000);
	gStyle->SetOptFit(0000);

	//Input TFiles with histograms
	TFile fData( (cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_Data2016_ReReco.root").c_str(),"read");
	TFile fMC( (cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_ReReco_lowM_QCD_Pt.root").c_str(),"read");

	//Extract histograms and they symmetric reflections
	auto h2D_data_numerator 		= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1vspt2");
	auto h2D_data_denumerator 	= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1vspt2");
	auto h2D_data_numerator_sym 		= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt2vspt1");
	auto h2D_data_denumerator_sym 	= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt2vspt1");

	auto h2D_mc_numerator 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1vspt2");
	auto h2D_mc_denumerator 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1vspt2");
	auto h2D_mc_numerator_sym 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt2vspt1");
	auto h2D_mc_denumerator_sym 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt2vspt1");

	//pack them into the My2DHisto container
	My2DHisto data, mc;
	data.num = &h2D_data_numerator; data.den = &h2D_data_denumerator; data.num_sym = &h2D_data_numerator_sym; data.den_sym  = &h2D_data_denumerator_sym;
	mc.num = &h2D_mc_numerator; mc.den = &h2D_mc_denumerator; mc.num_sym = &h2D_mc_numerator_sym; mc.den_sym  = &h2D_mc_denumerator_sym;

	// full title
	data.title = ";p_{T}^{(1)} [GeV];p_{T}^{(2)} [GeV]"; mc.title = ";p_{T}^{(1)} [GeV];p_{T}^{(2)} [GeV]";

	//Output path
	string output_path = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/TriggerPerformance/";
	//Perform Actions ;)
	CalculateAndPlotCorrelations(data,mc,output_path);

}

void CalculateAndPlotCorrelations(My2DHisto & data, My2DHisto& mc, const string& output_path){
	/*
	 * AllInOne method
	 */
	//Get Siametrised 2D efficiency
	auto data_efficiency_symmetric = GetSymmetricEfficiency(data);
	auto mc_efficiency_symmetric = GetSymmetricEfficiency(mc);

	//Plot 2D Symmetric efficiencies
	PlotSymmetric2DEfficiency(data_efficiency_symmetric,output_path,true);
	PlotSymmetric2DEfficiency(mc_efficiency_symmetric,output_path);

//	//Calculate Correlation coefficient and make a 2D histo
	auto data_correlation 	= GetCorrelationCoeff(data_efficiency_symmetric);
	auto mc_correlation 		= GetCorrelationCoeff(mc_efficiency_symmetric);
//
//	//Plot 2D correlations
	Plot2DCorrelations(data_correlation,output_path,true);
	Plot2DCorrelations(mc_correlation,output_path);
//
//	//Plot Correlation Slices
	PlotCorrelationSlices(data_correlation,mc_correlation,output_path,true);
}

TEfficiency * GetSymmetricEfficiency(My2DHisto & data){
	/*
	 * Function to siametrise efficiency
	 */
	TH2D *sym_num = new TH2D(*data.num_sym);
	TH2D *sym_den = new TH2D(*data.den_sym);
	for(int binx = 1; binx <= data.num->GetNbinsX();++binx){
		for(int biny = 1; biny <= data.num->GetNbinsY();++biny){
			//deal with diagonal bins
			if(binx == biny){
				sym_num->SetBinContent(binx, biny,0);
				sym_num->SetBinError(binx, biny,0);
				sym_den->SetBinContent(binx, biny,0);
				sym_den->SetBinError(binx, biny,0);
			}
		}
	}
//	Constract symmetric efficiency
	TEfficiency *eff = new TEfficiency(*data.num,*data.den);
	TEfficiency *eff_sym = new TEfficiency(*sym_num,*sym_den);
	BayesianErrors(eff);
	BayesianErrors(eff_sym);
	eff->Add(*eff_sym);
	eff->SetTitle(data.title.c_str());
	return eff;
}

TH2D * GetCorrelationCoeff(TEfficiency * eff){
	/*
	 * Calculate correlation coefficient bin-by-bin
	 */
	TH2D *correlation = new TH2D(*static_cast<TH2D*>(eff->GetCopyPassedHisto()));
	string xtitle = eff->GetCopyPassedHisto()->GetXaxis()->GetTitle();
	string ytitle = eff->GetCopyPassedHisto()->GetYaxis()->GetTitle();
	correlation->SetTitle( (";" + xtitle + ";" + ytitle).c_str() );
	for(int binx = 1; binx <= eff->GetCopyPassedHisto()->GetNbinsX();++binx){
		for(int biny = 1; biny <= eff->GetCopyPassedHisto()->GetNbinsY();++biny){
//			int xx_bin = binx + (binx ) * (eff->GetCopyPassedHisto()->GetNbinsX() + 2);
//			int yy_bin = biny + (biny ) * (eff->GetCopyPassedHisto()->GetNbinsY() + 2);
//			int global_bin = biny + (binx ) * (eff->GetCopyPassedHisto()->GetNbinsX() + 2);
			int global_bin = eff->GetGlobalBin(binx, biny);
			int xx_bin = eff->GetGlobalBin(binx, binx);
			int yy_bin = eff->GetGlobalBin(biny, biny);

			double xx = eff->GetEfficiency(xx_bin);
			double yy = eff->GetEfficiency(yy_bin);
			double xy = eff->GetEfficiency(global_bin);

			std::cout<<"bins: y = "<<yy_bin<<" x = "<<xx_bin<<" glob. = "<<global_bin<<endl;
			cout<<"x: "<<binx<<" y: "<<biny<<" xx = "<<xx<<" yy = "<<yy<<endl;

			double cxy = 0;
			double ecxy = 0;
			if(xx > 0 && yy>0){
				cxy = xy/sqrt(xx*yy);
				if(xx!=yy && xy != 0){
    				double exy = eff->GetEfficiencyErrorLow(global_bin) / xy;
    				double exx = eff->GetEfficiencyErrorLow(xx_bin) / xx;
    				double eyy = eff->GetEfficiencyErrorLow(yy_bin) / yy;
    		          ecxy = sqrt( exy*exy + 0.25*exx*exx + 0.25*eyy*eyy );
				}
			}
			correlation->SetBinContent(binx,biny,cxy);
			correlation->SetBinError(binx,biny,ecxy);
		}
	}
	return correlation;
}

void PlotSymmetric2DEfficiency(TEfficiency * eff, const string& output_path, const bool& data){
	/*
	 * Plot 2D symmetric efficiency
	 */
	TCanvas can("canv", "histograms", 600, 600);
	can.SetRightMargin(0.15);

	//Setup the frame
	TH2F frame("frame","",1,50.,500.,1,50.,500.);
	style.setFrameStyle(&frame);
	frame.GetXaxis()->SetTitle(eff->GetCopyPassedHisto()->GetXaxis()->GetTitle());
	frame.GetYaxis()->SetTitle(eff->GetCopyPassedHisto()->GetYaxis()->GetTitle());
	frame.GetYaxis()->SetTitleOffset(frame.GetYaxis()->GetTitleOffset()*1.1);
	frame.SetZTitle("#epsilon");
	frame.Draw();

	eff->Draw("COLZ same");
	gPad->RedrawAxis();
	style.drawStandardTitle("right-top","right");

	//Out name
	string out_name = output_path;
	if(data) out_name += "Data_";
	else out_name += "MC_";
	out_name += "TriggerEfficiency_2D_Symmetric_pT" ;

	//Save plots
	can.Print( (out_name + ".pdf").c_str());
	can.Print( (out_name + ".png").c_str());
}

void Plot2DCorrelations(TH2D * hdata, const string& output_path, const bool& data){
	/*
	 * Plot 2D Correlation coefficient
	 */
	TCanvas can("canv", "histograms", 600, 600);
	can.SetRightMargin(0.15);

	//Setup the frame
	TH2F frame("frame","",1,50.,500.,1,50.,500.);
	style.setFrameStyle(&frame);
	frame.SetTitle(hdata->GetTitle());
	frame.SetTitle(";p_{T}^{(1)} [GeV];p_{T}^{(2)} [GeV]");
	frame.SetZTitle("#epsilon");
	frame.GetYaxis()->SetTitleOffset(frame.GetTitleOffset("Y")*1.2);
	frame.Draw();

	hdata->Draw("COLZ same");
	gPad->RedrawAxis();
	style.drawStandardTitle("right-top","right");

	//Out name
	string out_name = output_path;
	if(data) out_name += "Data_";
	else out_name += "MC_";
	out_name += "TriggerEfficiency_2D_Correlations" ;

	//Save plots
	can.Print( (out_name + ".pdf").c_str());
	can.Print( (out_name + ".png").c_str());
}

void PlotCorrelationSlices(TH2D * hdata, TH2D * hmc, const string& output_path, const bool& data){
	/*
	 * Plot slices of correaltions
	 */
	for(int binx = 0; binx < hdata->GetNbinsX();++binx){
		TCanvas can("canv", "histograms", 600, 600);

		//top pad and frame
		auto topPad = HbbStyle::getRatioTopPad(0.7);
		topPad->Draw();
		topPad->cd();

		//Setup the frame
		TH2D frame("frame","",1,50.,500.,1,0.8,1.2);
//		style.setFrameStyle(&frame);
		style.setRatioTopFrame(&frame, can.YtoPixel(can.GetY1()));
		topPad->cd();
		frame.SetTitle(hdata->GetTitle());
		frame.SetTitle(";p_{T}^{(1)} [GeV]; #rho");
		frame.Draw();

		//Project 2D histos
		auto data1D = hdata->ProjectionX( ("Slice_data_binx_" + to_string(binx+1)).c_str(),binx+1,binx+1,"e");
		data1D->SetMarkerStyle(gStyle->GetMarkerStyle());
		data1D->SetMarkerSize(gStyle->GetMarkerSize());
		auto mc1D = hmc->ProjectionX( ("Slice_mc_binx_" + to_string(binx+1)).c_str(),binx+1,binx+1,"e");
		mc1D->SetFillColor(kRed-10);
		mc1D->SetFillStyle(3013);
		mc1D->SetMarkerSize(0);

		TExec er_0("er_0","gStyle->SetErrorX(0)");
		TExec er_1("er_1","gStyle->SetErrorX(0.5)");

		data1D->Draw("Esame");
		er_1.Draw();
		mc1D->Draw("E2 same");
		er_0.Draw();

		//Pt bin
		string binlow = to_string(int(data1D->GetBinLowEdge(binx+1)));
		string binup = to_string(int(data1D->GetBinLowEdge(binx+1) + data1D->GetBinWidth(binx+1)));

		//Setupo TLegend
		TLegend &leg = *HbbStyle::legend("top,right",4);
		leg.Clear();
		leg.SetHeader(("p_{T}^{(2)} #in [" + binlow + ";" + binup + "]").c_str());
		leg.SetY1(leg.GetY1() * 0.88);
		leg.SetY2(leg.GetY2() * 0.93);
		style.setLegendStyle(&leg);
		leg.AddEntry(data1D,"Data, JetHT","pl");
		leg.AddEntry(mc1D,"QCD, #hat{p_{T}}","f");
		leg.Draw();

		gPad->RedrawAxis();
		can.cd();


		//Bottom pad and frame
		auto botPad = HbbStyle::getRatioBottomPad(0.3);
		botPad->Draw();
		botPad->cd();
		TH2D *botFrame = new TH2D("botFrame","",1,frame.GetXaxis()->GetXmin(),frame.GetXaxis()	->GetXmax(),1,-3.9999,3.9999);
		HbbStyle::setRatioBottomFrame(botFrame, can.YtoPixel(can.GetY1()), topPad->YtoPixel(topPad->GetY1()));
		botFrame->GetXaxis()->SetTitle(frame.GetXaxis()->GetTitle());
		botFrame->GetYaxis()->SetTitle("#frac{Data-MC}{#sigma_{MC}}");

		botPad->cd();
		auto ratio = static_cast<TH1D*>(data1D->Clone("ratio"));
		ratio->Sumw2();
		//Minus and Devision
		devideTwoHistograms(ratio,data1D,mc1D);
		botFrame->Draw();
		ratio->Draw("Esame");

		TLine *horizLine = new TLine(frame.GetXaxis()->GetXmin(),0,frame.GetXaxis()->GetXmax(),0);
		horizLine -> SetLineStyle(2);
		horizLine -> Draw();

		can.cd();
		style.drawStandardTitle("out");

		//Out name
		string out_name = output_path;
		if(data) out_name += "Data_";
		else out_name += "MC_";
		out_name += "TriggerCorrelations_Slices_pT_" + binlow + "_" + binup;

		//Save plots
		can.Print( (out_name + ".pdf").c_str());
		can.Print( (out_name + ".png").c_str());
	}
}

void BayesianErrors(TEfficiency *h){
	h->SetStatisticOption(TEfficiency::kBUniform);
	h->SetPosteriorMode(1);
	h->UsesBayesianStat();
	h->SetUseWeightedEvents();
}

void devideTwoHistograms(TH1 * ratio, TH1 * nominator, TH1 * denominator){
	/*
	 * Function to perform:
	 * f = (a - b) / sigma(b)
	 * assume full correlations
	 */
	double difference, edifference, devis, edevis, bbb, ebbb, bbnb, ebbnb;
	for(int i = 1; i < nominator->GetNbinsX(); ++i){
		bbb  = nominator->GetBinContent(i);  ebbb  = nominator->GetBinError(i);
		bbnb = denominator->GetBinContent(i); ebbnb = denominator->GetBinError(i);
		difference  = bbb - bbnb;
		//Assume the full correlations
//		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb - 2 * ebbb * ebbnb );
		//Assume 0 correlation
		edifference = sqrt( ebbb*ebbb + ebbnb * ebbnb);
		devis 	= difference / ebbb;
		//assume the full correlations
		edevis 	= abs(devis) * sqrt ( pow(edifference/difference,2) + pow(ebbb/bbb,2) - 2 * edifference * ebbb / (difference * bbb) );
		//Assume the full correlations
		if(devis != devis) { devis = 0; edevis = 0;};
		ratio->SetBinContent(i,devis);
		ratio->SetBinError(i,edevis);
	}
}
