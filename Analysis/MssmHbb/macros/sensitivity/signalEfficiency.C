#include <memory>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TGraph.h"
#include "../Drawer/HbbStyle.cc"

double xsection(const int & point);
void efficiency(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda, const std::vector<double>& lumis);
void sensitivity(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda);

const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));

HbbStyle style;

void signalEfficiency()
{
   gStyle->SetOptStat(0);
   TH1::SetDefaultSumw2();

   std::map<int,TFile* > inFlowM;
   std::map<int,TFile* > inF2015M;
   std::map<int,TFile* > inFtoComp;
   std::string samples_name_2016 = "";//"ReReco_PU_test_subrange_presc";//
   std::string samples_comp		 = "efficiency_test_nobtag";
   std::string samples_name_2015 = "signalEff_test_SF";//"_signalEff_test";//

   inFlowM[300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[350] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[400] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[500]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[600]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[700] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[900] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[1100] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFlowM[1300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.root").c_str());

   inFtoComp[300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[350] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[400] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[500]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[600]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[700] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[900] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[1100] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inFtoComp[1300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_comp + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.root").c_str());

//   inF2015M[200] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2015 + "_lowM_SUSYGluGluToBBHToBB_M-200_TuneCUETP8M1_13TeV-pythia8.root");
//   inF2015M[250] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2015 + "_lowM_SUSYGluGluToBBHToBB_M-250_TuneCUETP8M1_13TeV-pythia8.root");
//   inF2015M[300] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_SUSYGluGluToBBHToBB_M-300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
//   inF2015M[350] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_SUSYGluGluToBBHToBB_M-350_TuneCUETP8M1_13TeV-pythia8.root").c_str());
//   inF2015M[400] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_SUSYGluGluToBBHToBB_M-400_TuneCUETP8M1_13TeV-pythia8.root").c_str());
//   inF2015M[500]	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_clange-SUSYGluGluToBBHToBB_M-500_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());
//   inF2015M[600]	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_clange-SUSYGluGluToBBHToBB_M-600_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());
//   inF2015M[700] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_chayanit-SUSYGluGluToBBHToBB_M-700_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());
//   inF2015M[900] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_chayanit-SUSYGluGluToBBHToBB_M-900_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());
//   inF2015M[1100] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_chayanit-SUSYGluGluToBBHToBB_M-1100_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());
//   inF2015M[1300] 	= new TFile(("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal" + samples_name_2015 + "_lowM_clange-SUSYGluGluToBBHToBB_M-1300_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root").c_str());

   std::map<int,double> signal_efficiency_8TeV;
   signal_efficiency_8TeV[200] = 2.85e-03;
   signal_efficiency_8TeV[300] = 6.38e-03;
   signal_efficiency_8TeV[350] = 6.32e-03;
   signal_efficiency_8TeV[400] = 6.08e-03;
   signal_efficiency_8TeV[500] = 5.07e-03;
   signal_efficiency_8TeV[600] = 3.85e-03;
   signal_efficiency_8TeV[700] = 2.9e-03;
   signal_efficiency_8TeV[900] = 1.39e-03;

   style.set(PRELIMINARY);
   style.drawStandardTitle();
   gStyle->SetTitleXOffset(1.);

   std::vector<std::map<int,TFile* > > datasets {inFlowM,inFtoComp};
   std::vector<std::string> legends {"2016","no offline b-tag"};
   std::vector<double> lumis {35673.,35673.};

   //Signal Efficiency calculation
   efficiency(datasets,legends,lumis);
   sensitivity(datasets,legends);

   //LowM
   TGraphErrors *Eff = new TGraphErrors(inFlowM.size());
   Eff->SetMarkerStyle(20);
   Eff->SetMarkerSize(1.1);
   Eff->SetTitle("; di-jet M, [GeV]; #epsilon");

   //highM
   TGraphErrors *EffhighM = new TGraphErrors(inF2015M.size());
   EffhighM->SetMarkerStyle(20);
   EffhighM->SetMarkerColor(2);
   EffhighM->SetLineColor(2);
   EffhighM->SetMarkerSize(1.1);
   EffhighM->SetMinimum(0.);
   EffhighM->GetXaxis()->SetRangeUser(0.,1300.);
   EffhighM->SetTitle("; di-jet M, [GeV]; #epsilon");

   // 8 TeV
   TGraphErrors *Eff8TeV = new TGraphErrors(signal_efficiency_8TeV.size());
   Eff8TeV->SetMarkerStyle(20);
   Eff8TeV->SetMarkerColor(3);
   Eff8TeV->SetLineColor(3);
   Eff8TeV->SetMarkerSize(1.1);
   Eff8TeV->SetMaximum(2.2e-02);
   Eff8TeV->SetTitle("; di-jet M, [GeV]; #epsilon");

   TCanvas * c1 = new TCanvas("c1","",700,800);
   TLegend *leg = new TLegend(0.65,0.75,0.92,0.9);
   leg->SetBorderSize(0);

   int i = 0;
   for(const auto & f : inFlowM){

	   double low_border = f.first*0.8;
	   double up_border	 = f.first*1.2;

	   TH1* denum =  (TH1*) f.second->Get("distributions/NumberOfGenEvents_afterMHat");
	   if(denum->GetMean() == 0){
		   denum = (TH1*) f.second->Get("distributions/TotalNumberOfGenEvents");
	   }
	   double nTot = denum->GetMean() * 36620. / denum->GetMean();

	   TH1 *integ = (TH1*) f.second->Get("bbH_Mbb");
	   Eff->SetPoint(i,f.first,(double) integ->Integral()/ (nTot));
	   ++i;
   }
   i=0;

   for(const auto & p : signal_efficiency_8TeV){

	Eff8TeV->SetPoint(i,p.first,p.second);
	++i;
   }
   i = 0;

   for(const auto & f : inF2015M){


	   double low_border = f.first*0.8;
	   double up_border	 = f.first*1.2;

	   TH1* denum =  (TH1*) f.second->Get("distributions/NumberOfGenEvents_afterMHat_rewPU");
	   if(denum->GetMean() == 0){
		   denum = (TH1*) f.second->Get("distributions/TotalNumberOfGenEvents");
	   }
	   TH1 *integ = (TH1*) f.second->Get("bbH_Mbb");
	   double nTot = denum->GetMean() * 2621.181 / denum->GetMean();

	   EffhighM->SetPoint(i,f.first,(double) integ->Integral()/nTot);
	   ++i;
   }
//   Eff8TeV->Draw("AP");
   EffhighM->Draw("APL");
   EffhighM->GetXaxis()->SetRangeUser(0.,1300.);
   EffhighM->Draw("APL");
   Eff->Draw("PL same");
   Eff8TeV->Draw("PL same");
   style.drawStandardTitle();
   c1->Update();
   
   leg->AddEntry(Eff,"Run 2016","p");
   leg->AddEntry(EffhighM,"Run 2015","p");
   leg->AddEntry(Eff8TeV,"Run 2012","p");
   leg->Draw();
   c1->Print("../pictures/SignalEfficiency_16vs15vs12.pdf");


   //Sensitivity

   TH1D	*hBgHighM = (TH1D*) inF2015M[500]->Get("QCD_Mbb");

   //LowM
   TGraphErrors *Sensitivity = new TGraphErrors(inFlowM.size());
   Sensitivity->SetMarkerStyle(20);
   Sensitivity->SetMarkerSize(1.1);
   Sensitivity->SetTitle("Signal Sensitivity = S/#sqrt{Bg}; di-jet M, [GeV]; S/#sqrt{Bg}");
   Sensitivity->GetXaxis()->SetTitleOffset(0.);

   //highM
   TGraphErrors *SensitivityHighM = new TGraphErrors(inF2015M.size());
   SensitivityHighM->SetMarkerStyle(20);
   SensitivityHighM->SetMarkerColor(2);
   SensitivityHighM->SetMarkerSize(1.1);
   SensitivityHighM->SetTitle("Signal Sensitivity = S/#sqrt{Bg}; di-jet M, [GeV]; #epsilon");

   TCanvas * c2 = new TCanvas("c2","",700,800);
   TLegend *leg2 = new TLegend(0.6,0.7,0.85,0.85);
   leg2->SetBorderSize(0);

   i = 0;
   std::cout<<"Low M:"<<std::endl;
   for(const auto & f : inFlowM){
	   ++i;
	   TH1 *s = (TH1*) f.second->Get("bbH_Mbb");
	   TH1D	*hBgLowM = (TH1D*) inFlowM[f.first]->Get("QCD_Mbb");

	   double low_border = f.first*0.8;
	   double up_border	 = f.first*1.2;

	   double signalErr_ = 0;
	   double signal_ = s->IntegralAndError(s->FindBin(low_border),s->FindBin(up_border),signalErr_);

	   double bgErr_ = 0;
	   double bg_ = hBgLowM->IntegralAndError(hBgLowM->FindBin(low_border),hBgLowM->FindBin(up_border),bgErr_);

	   double sens = signal_/std::sqrt(bg_) * xsection(f.first);
	   double sens_err = xsection(f.first) * std::sqrt(1./bg_ * signalErr_ * signalErr_ + signal_*signal_/(4. * bg_*bg_*bg_) * bgErr_* bgErr_);
	   Sensitivity->SetPoint(i,f.first,sens);
	   Sensitivity->SetPointError(i,0,sens_err);

	   std::cout<<"M12 = "<<f.first<<" Eff. = "<<sens<<" Signal = "<<signal_<<" sqrt(Bg). = "<<sqrt(bg_)<<std::endl;
   }
   i = 0;
std::cout<<"High M: "<<std::endl;
/*
   for(const auto & f : inF2015M){
	   ++i;
	   TH1 *s = (TH1*) f.second->Get("bbH_Mbb");

	   double low_border = f.first*0.8;
	   double up_border	 = f.first*1.2;

	   double signalErr_ = 0;
	   double signal_ = s->IntegralAndError(s->FindBin(low_border),s->FindBin(up_border),signalErr_);

	   double bgErr_ = 0;
	   double bg_ = hBgHighM->IntegralAndError(hBgHighM->FindBin(low_border),hBgHighM->FindBin(up_border),bgErr_);

	   double sens = signal_/std::sqrt(bg_) * xsection(f.first);
	   double sens_err = xsection(f.first) * std::sqrt(1./bg_ * signalErr_ * signalErr_ + signal_*signal_/(4. * bg_*bg_*bg_) * bgErr_* bgErr_);
	   SensitivityHighM->SetPoint(i,f.first,sens);
	   SensitivityHighM->SetPointError(i,0,sens_err);

	   std::cout<<"M12 = "<<f.first<<" Eff. = "<<sens<<" Signal = "<<signal_<<" sqrt(Bg). = "<<bg_<<std::endl;
 }
//   Sensitivity->GetYaxis()->SetRangeUser(0,10);
   Sensitivity->SetMinimum(0.);
   Sensitivity->SetTitle("");
   Sensitivity->Draw("AP");
//   SensitivityHighM->Draw("Psame");
   style.drawStandardTitle();

   leg2->AddEntry(Sensitivity,"13 TeV ","p");
//   leg2->AddEntry(SensitivityHighM,"High Mass","p");
   leg2->Draw();


//   c1->SaveAs("../pictures/signalEfficiency.pdf");
   c2->SaveAs("../pictures/signalSensitivity.pdf");
/**/
}

double xsection(const int &point){
	switch (point){
		case 100:
			return 5.23;
		case 120:
			return 41.79;
		case 160:
			return 75.32;
		case 200:
			return 35.42;
		case 250:
			return 15.55;
		case 300:
			return 7.64;
		case 350:
			return 4.1;
		case 400:
			return 2.34;
		case 500:
			return 0.88;
		case 600:
			return 0.377;
		case 700:
			return 0.18;
		case 900:
			return 0.051;
		case 1100:
			return 0.018;
		case 1300:
			return 0.008;
		default:
			std::cerr<<"No such mass point in the list of xsections"<<std::endl;
			exit(-1);
	}
}

void efficiency(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda, const std::vector<double>& lumis){
	std::vector<short int> colours = {1,2,3,4,6,9,12,28,46};
	std::vector<TGraphErrors> Efficiencies;
	TCanvas can("can","",700,800);
	TLegend leg(0.25,0.7,0.4,0.85);
	//Loop  over the datasets:
	int j =0;
	for(const auto& set : files){
		//Create efficiency plot
		TGraphErrors *eff = new TGraphErrors(set.size());
		eff->SetMarkerStyle(20);
		eff->SetMarkerSize(1.1);
		eff->SetMarkerColor(colours.at(j));
		eff->SetLineColor(colours.at(j));
		eff->SetTitle("; di-jet M, [GeV]; #epsilon");
		int i = 0;
		for(const auto& f : set){
			double low_border = f.first*0.8;
			double up_border	 = f.first*1.2;

			TH1* denum =  (TH1*) f.second->Get("distributions/NumberOfGenEvents_afterMHat_rewPU");
			TH1 *integ = (TH1*) f.second->Get("bbH_Mbb");
			double nTot = lumis.at(j);
			double selected = integ->Integral();
			eff->SetPoint(i,f.first,selected/nTot);
			std::cout<<"M12 = "<<f.first<<" Eff. = "<<selected/nTot<<" Tot = "<<nTot<<" Sel. = "<<selected<<std::endl;
			++i;
		}
		i = 0;
		++j;
		Efficiencies.push_back(*eff);
	}

	for(unsigned int i = 0; i != Efficiencies.size(); ++i){
		TGraphErrors& eff = Efficiencies.at(i);
//		eff.SetMaximum(0.021);
		eff.SetMaximum(0.5);
		if(i==0) eff.Draw("APL");
		else eff.Draw("PLsame");

		leg.AddEntry(&eff,(legenda.at(i)).c_str(),"pl");
	}
	style.drawStandardTitle();
	leg.Draw();
	can.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/SignalEfficiency.pdf").c_str());
}

void sensitivity(const  std::vector<std::map<int,TFile* > >& files,const std::vector<std::string>& legenda){
	std::vector<short int> colours = {1,2,3,4,6,9,12,28,46};
	std::vector<TGraphErrors> Efficiencies;
	TCanvas can("can","",700,800);
	TLegend leg(0.65,0.7,0.85,0.85);
	//Loop  over the datasets:
	int j =0;
	for(const auto& set : files){
		//Create efficiency plot
		TGraphErrors *eff = new TGraphErrors(set.size());
		eff->SetMarkerStyle(20);
		eff->SetMarkerSize(1.1);
		eff->SetMarkerColor(colours.at(j));
		eff->SetLineColor(colours.at(j));
		eff->SetTitle("; di-jet M, [GeV]; S/#sqrt{Bg}");
		int i = 0;
		for(const auto& f : set){

			double low_border = f.first*0.8;
			double up_border	 = f.first*1.2;

			TH1 *integ = (TH1*) f.second->Get("bbH_Mbb");
			TH1 *hbg   = (TH1*) f.second->Get("QCD_Mbb");
			double signalErr_ = 0;
			double signal_ = integ->IntegralAndError(integ->FindBin(low_border),integ->FindBin(up_border),signalErr_);

			double bgErr_ = 0;
			double bg_ = hbg->IntegralAndError(hbg->FindBin(low_border),hbg->FindBin(up_border),bgErr_);

			double sens = signal_/std::sqrt(bg_) * xsection(f.first);
			double sens_err = xsection(f.first) * std::sqrt(1./bg_ * signalErr_ * signalErr_ + signal_*signal_/(4. * bg_*bg_*bg_) * bgErr_* bgErr_);
			eff->SetPoint(i,f.first,sens);
			eff->SetPointError(i,0,sens_err);
			++i;
		}
		i = 0;
		++j;
		Efficiencies.push_back(*eff);
	}

	for(unsigned int i = 0; i != Efficiencies.size(); ++i){
		TGraphErrors& eff = Efficiencies.at(i);
		eff.SetMaximum(3.7);
		if(i==0) eff.Draw("APL");
		else eff.Draw("PLsame");

		leg.AddEntry(&eff,(legenda.at(i)).c_str(),"pl");
	}
	style.drawStandardTitle();
	leg.Draw();
	can.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/SignalSensitivity.pdf").c_str());

}
