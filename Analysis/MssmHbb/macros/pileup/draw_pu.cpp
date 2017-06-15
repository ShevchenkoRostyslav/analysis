/*
 * draw_pu.cpp
 *
 *  Created on: 24 May 2017
 *      Author: shevchen
 */

/*
 *
 * Macro to draw PU
 * distribution before
 * and after reweighting
 *
 */

#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

using namespace std;
using namespace analysis;
using namespace mssmhbb;

void Renormalise(TH1&);

int draw_pu(){
	//style
	HbbStyle style;
	style.set(PRELIMINARY);
	string pu_distr_folder = "/src/Analysis/MssmHbb/bin/input_corrections/";
	//ReReco MC
	auto* fMC = new TFile( (cmsswBase + pu_distr_folder + "MC_Autumn16_PU25ns_V1.root").c_str());
	//ReReco Data
	auto *fData_Central  = new TFile( (cmsswBase + pu_distr_folder + "Data_Pileup_2016_ReReco.root").c_str() );				//Central
	auto *fData_2SigmaUp = new TFile( (cmsswBase + pu_distr_folder + "Data_Pileup_2016_ReReco_up2Sigma.root").c_str());   	//Rereco Data 2Sigma Up
	auto *fData_2SigmaDn = new TFile( (cmsswBase + pu_distr_folder + "Data_Pileup_2016_ReReco_down2Sigma.root").c_str());	//Rereco Data 2Sigma Down

	auto* hMC = GetFromTFile<TH1D>(*fMC,"pileup");
	Renormalise(*hMC);
	hMC->SetLineColor(1);

	auto* hData_Central = GetFromTFile<TH1D>(*fData_Central,"pileup");
	hData_Central->SetMarkerStyle(20);
	Renormalise(*hData_Central);

	auto* hData_2SigmaUp = GetFromTFile<TH1D>(*fData_2SigmaUp,"pileup");
	hData_2SigmaUp->SetMarkerStyle(20);
	hData_2SigmaUp->SetMarkerColor(4);
	Renormalise(*hData_2SigmaUp);

	auto* hData_2SigmaDn = GetFromTFile<TH1D>(*fData_2SigmaDn,"pileup");
	hData_2SigmaDn->SetMarkerStyle(20);
	hData_2SigmaDn->SetMarkerColor(2);
	Renormalise(*hData_2SigmaDn);

	auto* leg = style.legend("top,right", 4);
//	auto* leg = new TLegend(0.6,0.6,0.9,0.9);
	leg->AddEntry(hMC,"MC","l");
	leg->AddEntry(hData_Central,"Data, central", "p");
	leg->AddEntry(hData_2SigmaUp,"#sigma(min-bias), 2#sigma up", "p");
	leg->AddEntry(hData_2SigmaDn,"#sigma(min-bias), 2#sigma down", "p");

	auto* canva = new TCanvas("can","can",800,600);
	hMC->SetTitle(";N_{PV};a.u.");
	hMC->Draw("hist");
	hMC->SetMaximum(1.1*hData_2SigmaDn->GetMaximum());
	hMC->SetStats(0);
	hData_Central->Draw("same");
	hData_2SigmaUp->Draw("same");
	hData_2SigmaDn->Draw("same");
	leg->Draw();

	canva->Print("../pictures/PU.pdf");

	return 0;
}

void Renormalise(TH1& h){
	h.Scale(1./h.Integral());
}
