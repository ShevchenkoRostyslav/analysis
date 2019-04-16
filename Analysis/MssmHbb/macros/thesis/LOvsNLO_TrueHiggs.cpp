/*
 * NLO_vs_LO_shapes.cpp
 *
 *  Created on: 20 Nov 2017
 *      Author: shevchen
 */
#include <string>
#include <vector>

#include "TChain.h"
#include "TFileCollection.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "TSystem.h"

using namespace std;

HbbStyle style;

int main(){
	style.setTDRstyle(PRIVATE);
	vector<string> masses = {"350","600","900"};
	vector<string> variables = {"pt[0]"};
//	vector<string> additional_options = {"","","pt[3] > "}
	for(const auto& mass : masses){

		//NLO
		string inputFilelist_NLO = "/afs/desy.de/user/c/chayanit/cms/forRostyslav/21072016/SUSYGluGluToBBHToBB_M-" + mass + "_cfg_ntuple.txt";
		TFileCollection *fNLO = new TFileCollection("fNLO","",inputFilelist_NLO.c_str());
		auto NLO_collection = (TCollection*) fNLO->GetList();

		TChain *cNLO = new TChain("MssmHbb/Events/prunedGenParticles");
//		TChain *cNLO_genJets = new TChain("MssmHbb/Events/slimmedGenJets");
		TChain *cNLO_weight = new TChain("MssmHbb/Events/EventInfo");
		cNLO_weight -> AddFileInfoList(NLO_collection);
		cNLO -> AddFileInfoList(NLO_collection);
//		cNLO_genJets -> AddFileInfoList(NLO_collection);
		cNLO -> AddFriend(cNLO_weight);
//		cNLO -> AddFriend(cNLO_genJets);

		// LO
		string inputFilelist_LO = "/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/test/ForMSSMHbb/2016/Moriond17/SUSYGluGluToBBHToBB_NarrowWidth_M-" + mass + "_TuneCUETP8M1_13TeV-pythia8.txt";
		TFileCollection *fLO = new TFileCollection("fLO","",inputFilelist_LO.c_str());
		auto LO_collection = (TCollection*) fLO->GetList();

		TChain *cLO = new TChain("MssmHbb/Events/prunedGenParticles");
//		TChain *cLO_genJets = new TChain("MssmHbb/Events/slimmedGenJets");
		TChain *cLO_weight = new TChain("MssmHbb/Events/EventInfo");
		cLO_weight -> AddFileInfoList(LO_collection);
		cLO -> AddFileInfoList(LO_collection);
//		cLO_genJets -> AddFileInfoList(NLO_collection);
		cLO -> AddFriend(cLO_weight);
//		cLO -> AddFriend(cNLO_genJets);

		int i = 0;
		for(const auto& variable : variables){
			double min = 0, max =600;
			double ymax = 0.2;
			if(variable.find("eta")!=string::npos) {max = 6.3; min = -6.3; ymax = 0.04;}

			//Make histograms
			TCanvas *can = new TCanvas( ("can_" + mass + "_" + variable).c_str());
			TH2D *topFrame = new TH2D("topFrame","",1,min,max,1,0,ymax);
			auto *topPad = style.getRatioTopPad();
			style.setRatioTopFrame(topFrame,can->YtoPixel(can->GetY1()) );
			topFrame->GetYaxis()->SetTitle("Arbitrary units");
			topPad->Draw();
			topPad->cd();
			auto *leg = style.legend("top,right",4);//new TLegend(0.6,0.6,0.95,0.95);
			leg->SetY1(leg->GetY1()*0.7);
			leg->SetY2(leg->GetY2()*0.95);

			TH1D *hLO 	= new TH1D("hLO","",100,min,max);
			style.applyStandardMCHistoStyle(*hLO);
			hLO->SetFillColor(kRed-10);
			hLO->SetStats(0);
//			topFrame->SetTitle(";Higgs p_{T} [GeV]; 1/#sigma d#sigma / dp_{T} [GeV^{-1}]");
			TH1D *hNLO 	= new TH1D("hNLO","",100,min,max);
			style.applyStandardDataHistoStyle(*hNLO);

			cLO->Draw( (variable + " >> hLO").c_str(),"genWeight && pdg == 36 && status == 62","hist");
			cNLO->Draw((variable + " >> hNLO").c_str(),"genWeight && pdg == 25 && status == 62","Esame");

			leg->SetHeader( ("#font[12]{#splitline{" + mass + " GeV}{NLO vs LO}}").c_str());
			leg->AddEntry(hLO,"LO Pythia 8","f");
			leg->AddEntry(hNLO,"NLO Madgraph 5","p");

			hLO->Scale(1./hLO->GetSumOfWeights());
			hNLO->Scale(1./hNLO->GetSumOfWeights());

			topFrame->Draw();
			hLO->Draw("hist same");
			hNLO->Draw("hist same P");
			leg->Draw();
			gPad->RedrawAxis();

			can->cd();
			auto *botPad = style.getRatioBottomPad();
			botPad->Draw();
			botPad->cd();
			TH2D *frame2 = new TH2D("frame2","",1,min,max,1,0,3);
			HbbStyle::setRatioBottomFrame(frame2->GetXaxis(),frame2->GetYaxis(),can->YtoPixel(can->GetY1()),topPad->YtoPixel(topPad->GetY1()));
			frame2->Draw();
			auto ratio = static_cast<TH1D*>(hNLO->Clone("ratio"));
			ratio->SetLineColor(kBlack);
			if(variable == "pt[0]") frame2->SetTitle( (";Higgs p_{T} [GeV]; #frac{NLO}{LO}"));
			ratio->Divide(hLO);
			ratio->Draw("HISTO same");

			can->cd();
			style.drawStandardTitle("out");

			can->Print(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_8_0_20_patch1/src/Analysis/MssmHbb/macros/pictures/Thesis/NLO_vs_LO_" + variable + "_" + mass + ".pdf").c_str());
			++i;
		}
	}

	return 0;
}


