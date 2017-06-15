#include <memory>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TExec.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/macros/Drawer/RatioPlots.cpp"

RatioPlots style;

template <typename T> struct point{
	point(const T& Nominal,const T& Up, const T& Down) : nominal(Nominal), up(Up), down(Down) {};
	point(const T& Up,const T& Down) : nominal(0), up(Up), down(Down) {};
	point(const point& p) {nominal = p.nominal; up = p.up; down = p.down;}

	T nominal;
	T up;
	T down;
};
double fit(double *x, double *par) {return par[0];}
void SetBottomStyle(TH1 *hRatio);
void sistTemp(const std::map<int,TFile*> & inF,const std::string & regime);
void allInOne(const std::map<int,TFile*> & inF,const std::string & regime);
void selectedTemplates(const std::map<int,TFile*> & inF,const std::string &regime, const std::vector<int> & points,const float& xmin = -1,const float& xmax = -1);
void calcLnN(const std::string& sys, const std::map<int,double > & lnN, const std::map<int,double > & elnN);
double getMean(const double& v1, const double&v2){double f = (v1+v2)/2.;return f;}
double correlatedDivision(const double& v1, const double& ev1,const double& v2, const double& ev2){
	// y = v1/v2
	if(v1 == 0 || v2 == 0) return 0;
	double f = v1/v2 * sqrt( pow(ev1/v1,2) + pow(ev2/v2,2) - 2*(ev1*ev2)/(v2*v1) );
	return f;
}

using namespace std;

void templates()
{
	TH1::SetDefaultSumw2();
   gStyle->SetOptStat(0);

   std::map<int,TFile* > inF;
   std::vector<TH1D*> histo;

   style.set(PRELIMINARY);


   TExec *er_0 = new TExec("er_0","gStyle->SetErrorX(0)");
   TExec *er_1 = new TExec("er_1","gStyle->SetErrorX(0.5)");
   const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));
   std::string samples_name_2016 = "ReReco_35673fb";//"ReReco_PU_test_subrange_presc";//

   inF[300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-300_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[350] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-350_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[400] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-400_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[500]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-500_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[600]		= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-600_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[700] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-700_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[900] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-900_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[1100] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1100_TuneCUETP8M1_13TeV-pythia8.root").c_str());
   inF[1300] 	= new TFile((cmsswBase + "/src/Analysis/MssmHbb/output/MssmHbbSignal_" + samples_name_2016 + "_lowM_SUSYGluGluToBBHToBB_NarrowWidth_M-1300_TuneCUETP8M1_13TeV-pythia8.root").c_str());

//   inF[500]		= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-500_TuneCUETP8M1_13TeV-pythia8.root");
//   inF[600]		= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-600_TuneCUETP8M1_13TeV-pythia8.root");
//   inF[700] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-700_TuneCUETP8M1_13TeV-pythia8.root");
//   inF[900] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-900_TuneCUETP8M1_13TeV-pythia8.root");
//   inF[1100] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-1100_TuneCUETP8M1_13TeV-pythia8.root");
//   inF[1300] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_lowM_SUSYGluGluToBBHToBB_M-1300_TuneCUETP8M1_13TeV-pythia8.root");
/**/
   std::map<int,TFile* > highMF;
   std::vector<TH1D*> histoHighM;
   highMF[300] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_SUSYGluGluToBBHToBB_M-300_TuneCUETP8M1_13TeV-pythia8.root");
   highMF[350] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_SUSYGluGluToBBHToBB_M-350_TuneCUETP8M1_13TeV-pythia8.root");
   highMF[400] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_SUSYGluGluToBBHToBB_M-400_TuneCUETP8M1_13TeV-pythia8.root");
   highMF[500]	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_clange-SUSYGluGluToBBHToBB_M-500_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");
   highMF[600]	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_clange-SUSYGluGluToBBHToBB_M-600_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");
   highMF[700] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_chayanit-SUSYGluGluToBBHToBB_M-700_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");
   highMF[900] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_chayanit-SUSYGluGluToBBHToBB_M-900_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");
   highMF[1100] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_chayanit-SUSYGluGluToBBHToBB_M-1100_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");
   highMF[1300] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_highM_clange-SUSYGluGluToBBHToBB_M-1300_cfg_GEN_DIGI76X_RECO76X_MiniAODv2_76X-17d438ff51ec6b3cada9e499a5a978e0.root");


   //   highMF[500]	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-500_TuneCUETP8M1_13TeV-pythia8.root");
//   highMF[600]	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-600_TuneCUETP8M1_13TeV-pythia8.root");
//   highMF[700] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-700_TuneCUETP8M1_13TeV-pythia8.root");
//   highMF[900] 	= new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-900_TuneCUETP8M1_13TeV-pythia8.root");
//   highMF[1100] = new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-1100_TuneCUETP8M1_13TeV-pythia8.root");
//   highMF[1300] = new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/MssmHbbSignal_1pb_xsection_highM_SUSYGluGluToBBHToBB_M-1300_TuneCUETP8M1_13TeV-pythia8.root");


//   allInOne(inF,"lowM");
//   allInOne(highMF,"highM");

//   sistTemp(inF,"lowM");

   //Wide templates
   std::vector<int> Wide_points = {700,900,1100,1300};
   //Narrow peaks:
   std::vector<int> Narrow_points = {300,350,400};
   // Just some selection of the points
   std::vector<int> Selected_points = {300,500,700,1100,1300};

   selectedTemplates(inF,"Wide peaks low",Wide_points);
   selectedTemplates(inF,"Narrow peaks low",Narrow_points);
   selectedTemplates(inF,"low",Selected_points);


/**/
   
}

void SetBottomStyle(TH1 *hRatio){

	hRatio -> SetNdivisions(505,"YZ");
	hRatio -> SetNdivisions(510,"X");
	hRatio -> SetTickLength(0.04,"YZ");
	hRatio -> SetTickLength(0.07,"X");
	hRatio -> SetTitleSize(0.06*7/3,"XYZ");
	hRatio -> SetTitleFont(42,"XYZ");

	hRatio -> SetLabelFont(42,"XYZ");
	hRatio -> SetLabelSize(0.05*7/3,"XYZ");
	hRatio -> SetLabelOffset(0.01,"XYZ");

	hRatio -> GetYaxis() -> SetTitleOffset(0.6);
	hRatio -> GetXaxis() -> SetTitleOffset(1.);
	hRatio -> GetYaxis() -> SetRangeUser(0.,2.);


}

void selectedTemplates(const std::map<int,TFile*> & inF,const std::string &regime, const std::vector<int> & points,const float& xmin,const float& xmax){
	bool lowM = regime.find("low") != std::string::npos;


	TCanvas *c = new TCanvas(("c"+regime).c_str(),("All mass points " + regime).c_str(),1000,800);
	std::array<int,15> colors = {1,2,4,6,80,95};
	double legXMin = 0.6;
	double legXMax = 0.95;
	double legYMin = 0.4;
	double legYMax = 0.85;
	TLegend *leg = (TLegend*) style.legend("top,right",points.size()+2,0.45);
//	TLegend *leg = new TLegend(0.69,0.4,0.92,0.85);
	leg->SetEntrySeparation(0.01);
	leg->SetLineColor(0);
	float Xmax = xmax, Xmin = xmin;
	float ymax = -1000.;
	if(xmax == -1) Xmax = 1700;

	TH1D *h[15];
	int i=0;
	for(const auto & point : points){
		++i;
		c->cd();
		auto it = inF.find(point);
		h[i-1] = (TH1D*) it->second->Get("templates/bbH_Mbb_VIS");
		h[i-1]->SetLineColor(colors.at(i-1));
		h[i-1]->SetLineWidth(3.5);
		if(ymax <= h[i-1]->GetMaximum()) ymax = 1.1 * h[i-1]->GetMaximum();
		if(i==1){
			h[i-1]->GetYaxis()->SetTitle("a.u.");
			h[i-1]->SetTitle("");
			if(lowM) h[i-1]->SetAxisRange(0.,60.,"y");
			else h[i-1]->SetAxisRange(0.,60.,"y");
			if(xmin == -1) Xmin = h[i-1]->GetXaxis()->GetXmin();
			h[i-1]->SetAxisRange(Xmin,Xmax,"x");
			h[i-1]->Draw("hist");
		}
		else h[i-1]->Draw("hist same");
		leg->AddEntry(h[i-1], (std::to_string(point) + " GeV").c_str(),"l" );
	}
	h[0]->SetAxisRange(0,ymax,"y");
	std::string legHeader;
	if(lowM) legHeader = "Low Mass Scenario";
	else legHeader = "High Mass Scenario";
	leg->SetHeader(legHeader.c_str());
	leg->Draw();
	style.drawStandardTitle();
	std::string cPrinter = "pictures/AllinOne_" + regime + "_Masses.pdf";
	c->Print(cPrinter.c_str());
}

void allInOne(const std::map<int,TFile*> & inF,const std::string & regime){
	double xMin = 0, xMax = 0;
	std::vector<TH1D*> histo;		//vector of syst. histograms.

	TCanvas * c1[15];
	TLegend *leg = new TLegend(0.3,0.3,0.9,0.85);
	leg->SetBorderSize(0);

	int i = 0, j =0;
	   for(const auto & f : inF){
		   ++i;
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_VIS"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVUp"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVDown"));
		   histo.push_back((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVUp"));

		   c1[i-1] = new TCanvas();
		   for(const auto & h : histo){
			   ++j;
			   h->SetTitle((std::to_string(f.first) + "-M mass point").c_str());
			   h->SetMarkerStyle(20);
			   h->SetMarkerColor(j);
			   h->SetLineColor(j);
			   h->SetLineWidth(2);
			   if(j==1){
				   h->Draw("hist");
				   if(i==1) leg->AddEntry(h,h->GetName(),"l");
			   }
			   else {
				   h->Draw("same");
				   if(i==1) leg->AddEntry(h,h->GetName(),"p");
			   }

			   xMin = h->GetXaxis()->GetXmin();
			   xMax = h->GetXaxis()->GetXmax();
		   }
		   j=0;

		   if(i==1) leg->Draw();
		   std::string canva_name = "pictures/AllinOne_" + regime + "_" + std::to_string(f.first) + "_GeV_mass_point.pdf";
		   c1[i-1]->Print(canva_name.c_str());
		   histo.clear();
	   }
}

void sistTemp(const std::map<int,TFile*> & inF,const std::string & regime){
	//One by one
	map<string,map<int,pair<double,double> > > lnN;
	map<string,map<int,pair<double,double> > > elnN;
	map<string,map<int,double> > av_lnN;
	map<string,map<int,double> > av_elnN;
	std::map<std::string,std::pair<TH1D*, TH1D*> > singleH;
	//coordinates for TLegend
	double lxmin = 0.6, lxmax = 0.9, lymin = 0.7, lymax = 0.9;
	//coordinates for pavetext
	double xmin = 0.6, ymin = 0.5, xmax = 0.9, ymax = 0.7;

	int i = 0, j =0, k=0;
	for(const auto & f : inF){
		++i;
		TH1D * first = (TH1D*) f.second->Get("templates/bbH_Mbb_VIS");
		f.second->SetTitle((std::to_string(f.first) + "-M mass point").c_str());
//		first->SetLineWidth(2);
		first->SetLineColor(kRed-9);
		first->SetFillColor(kRed-10);
		first->SetFillStyle(3002);
		first->GetYaxis()->SetTitle("a.u.");

		singleH["SFb"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVUp"));
		singleH["SFl"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVUp"));
		singleH["JER"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVUp"));
		singleH["JES"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVUp"));
		singleH["PU"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVUp"));
		singleH["BTagEff"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVUp"));
		singleH["pT trig Eff"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVUp"));
		for(const auto & p: singleH){
			++j;
			TCanvas can("can","can");

			std::string name2 = "box1_" + std::to_string(i) + "_" + std::to_string(j);
			TPad cPad(name2.c_str(),"",0.,0.1,1,1);
			cPad.SetBottomMargin(0.2);
			cPad.SetLeftMargin(0.16);
			cPad.SetRightMargin(0.05);
			cPad.Draw();
			cPad.cd();
			string leg_position = "top,right";
			double leg_width = 0.47;
			if(f.first > 900) {
				leg_position = "top,left";
				leg_width    = 0.28;
			}

			TLegend legend = (TLegend&)* style.legend(leg_position,3,leg_width);
//			TLegend legend(lxmin,lymin,lxmax,lymax);
			(p.second).second->SetTitle((std::to_string(f.first)+ " GeV mass point, " + p.first+" variation").c_str());

			(p.second).first->SetMarkerColor(kRed);
			(p.second).second->SetMarkerColor(kBlue);
			(p.second).first->SetMarkerStyle(20);
			(p.second).second->SetMarkerStyle(20);
			(p.second).second->GetYaxis()->SetTitle("a.u.");
			(p.second).second->SetMinimum(0.);

			style.SetTopStyle((p.second).second);
			(p.second).second->Draw();
			first->Draw("hist same");
			(p.second).first->Draw("same");
			gPad->RedrawAxis();

			legend.AddEntry(first,"central","fl");
			legend.AddEntry((p.second).first, (" -2#sigma, " + p.first).c_str(),"p");
			legend.AddEntry((p.second).second,("+2#sigma, " + p.first).c_str(),"p");
			legend.Draw();

			can.cd();
			name2 = "box2_" + std::to_string(i) + "_" + std::to_string(j);
			TPad uPad(name2.c_str(),"",0.,0.,1,0.265);
			uPad.SetBottomMargin(0.28);
			uPad.SetTopMargin(0.);
			uPad.SetLeftMargin(cPad.GetLeftMargin());
			uPad.SetRightMargin(cPad.GetRightMargin());
			uPad.Draw();
			uPad.cd();

			TH1D *centrU = (TH1D*) first->Clone("centrU");
			TH1D *centrD = (TH1D*) first->Clone("centrD");
			double cVal, uVal, dVal, ecVal, euVal, edVal;
			cVal = centrU->IntegralAndError(centrU->FindFirstBinAbove(0),centrU->FindLastBinAbove(0),ecVal);
			uVal = (p.second).second->IntegralAndError((p.second).second->FindFirstBinAbove(0),(p.second).second->FindLastBinAbove(0),euVal);
			centrU->Divide((p.second).second);
			double err_u;
			for(auto i = 0; i < centrU->GetNbinsX(); ++i){
				err_u = correlatedDivision(first->GetBinContent(i+1),first->GetBinError(i+1),(p.second).second->GetBinContent(i+1),(p.second).second->GetBinError(i+1));
				centrU->SetBinError(i+1,0);
			}

			double err_d;
			centrD->Divide((p.second).first);
			for(auto i = 0; i < centrD->GetNbinsX(); ++i){
				err_d = correlatedDivision(first->GetBinContent(i+1),first->GetBinError(i+1),(p.second).first->GetBinContent(i+1),(p.second).first->GetBinError(i+1));
				centrD->SetBinError(i+1,0);
			}

			centrD->SetTitle("");
			centrD->GetYaxis()->SetTitle("#frac{central}{down(up)}");
			centrD->SetLineColor(kRed);
			centrD->SetMarkerColor(kRed);
			centrD->SetMinimum(0.5);
			centrD->SetMaximum(1.5);
			centrD->SetMarkerStyle(20);
			style.SetBottomStyle(centrD);
			centrD->GetXaxis()->SetTickLength((p.second).first->GetXaxis()->GetTickLength()*3);
			centrD->Draw("P");

			centrU->SetTitle("");
			centrU->SetLineColor(kBlue);
			centrU->SetMarkerColor(kBlue);
			centrU->SetMarkerStyle(20);
			centrU->Draw("Psame");

			TLine *horizLine = new TLine(centrU->GetXaxis()->GetXmin(),1,centrU->GetXaxis()->GetXmax(),1);
			horizLine -> SetLineStyle(2);
			horizLine -> Draw();

//			can.cd();
//			name2 = "box3_" + std::to_string(i) + "_" + std::to_string(j);
//			TPad dPad(name2.c_str(),"",0.,0.,1,0.2);
//			dPad.SetBottomMargin(0.33);
//			dPad.SetTopMargin(0.);
//			dPad.SetLeftMargin(cPad.GetLeftMargin());
//			dPad.SetRightMargin(cPad.GetRightMargin());
//			dPad.Draw();
//			dPad.cd();
//
//			TH1D *centrD = (TH1D*) first->Clone("centrD");
//			dVal = (p.second).first->IntegralAndError((p.second).first->FindFirstBinAbove(0),(p.second).first->FindLastBinAbove(0),edVal);
//			lnN[p.first][f.first] = make_pair(cVal/dVal,uVal/cVal);
//			elnN[p.first][f.first] = make_pair(correlatedDivision(cVal,ecVal,dVal,edVal),correlatedDivision(uVal,euVal,cVal,ecVal));
//			av_lnN[p.first][f.first] = 1. + std::abs((lnN[p.first][f.first].second + lnN[p.first][f.first].first - 2)/4.);
//			av_elnN[p.first][f.first] = 0.25 * sqrt(elnN[p.first][f.first].second * elnN[p.first][f.first].second + elnN[p.first][f.first].first * elnN[p.first][f.first].first - 2*elnN[p.first][f.first].second*elnN[p.first][f.first].first);
//			//calculate errors of division:
//			double err_d;
//			centrD->Divide((p.second).first);
//			for(auto i = 0; i < centrD->GetNbinsX(); ++i){
//				err_d = correlatedDivision(first->GetBinContent(i+1),first->GetBinError(i+1),(p.second).first->GetBinContent(i+1),(p.second).first->GetBinError(i+1));
//				centrD->SetBinError(i+1,err_d);
//			}
//			centrD->SetTitle("");
//			centrD->SetMarkerStyle(20);
//			SetBottomStyle(centrD);
//			centrD->Draw();
//
//			TLine *horizLine2 = new TLine(centrD->GetXaxis()->GetXmin(),1,centrD->GetXaxis()->GetXmax(),1);
//			horizLine2 -> SetLineStyle(2);
//			horizLine2 -> Draw();
//
//			cPad.cd();
//			if (f.first > 900) {xmin = 0.25; xmax = 0.55;}
//			TPaveText text(xmin,ymin,xmax,ymax,"NDC");
//			text.SetFillColor(0);
//			text.SetShadowColor(0);
//			text.AddText(("lnN +2#sigma: " + std::to_string(lnN[p.first][f.first].second)).c_str());
//			text.AddText(("lnN -2#sigma: " + std::to_string(1. / lnN[p.first][f.first].first)).c_str());
//			text.AddText(("lnN #pm1#sigma: " + std::to_string(av_lnN[p.first][f.first])).c_str());
//			//("lnN +2\sigma: " + std::to_string(lnN[i-1][j-1].first) + " \pm " + std::to_string(elnN[i-1][j-1].first)  + " \splitline lnN -2#sigma: " + std::to_string(lnN[i-1][j-1].second) + " #pm " + std::to_string(elnN[i-1][j-1].second)).c_str());
//			text.Draw();
			std::string name = "pictures/"+std::to_string(f.first)+ "_GeV" + "_mass_point_" +regime + "_" + p.first+"_variation.pdf";
			can.SaveAs(name.c_str());
		}
		j=0;
	}

	//Make plot and fit to lnN of different mass points
	std::vector<string> Syst = {"SFb","SFl","PU"};
	for(const auto& s : Syst){
		calcLnN(s,av_lnN[s],av_elnN[s]);
	}
}


//void sistTemp(const std::map<int,TFile*> & inF,const std::string & regime){
//	//One by one
//	map<string,map<int,pair<double,double> > > lnN;
//	map<string,map<int,pair<double,double> > > elnN;
//	map<string,map<int,double> > av_lnN;
//	map<string,map<int,double> > av_elnN;
//	std::map<std::string,std::pair<TH1D*, TH1D*> > singleH;
//	//coordinates for TLegend
//	double lxmin = 0.6, lxmax = 0.9, lymin = 0.7, lymax = 0.9;
//	//coordinates for pavetext
//	double xmin = 0.6, ymin = 0.5, xmax = 0.9, ymax = 0.7;
//
//	int i = 0, j =0, k=0;
//	for(const auto & f : inF){
//		++i;
//		TH1D * first = (TH1D*) f.second->Get("templates/bbH_Mbb_VIS");
//		f.second->SetTitle((std::to_string(f.first) + "-M mass point").c_str());
////		first->SetLineWidth(2);
//		first->SetLineColor(kRed-9);
//		first->SetFillColor(kRed-10);
//		first->SetFillStyle(3002);
//
//		singleH["SFb"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFb_VIS_13TeVUp"));
////		singleH["SFl"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_SFl_VIS_13TeVUp"));
////		singleH["JER"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JER_VIS_13TeVUp"));
////		singleH["JES"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_JES_VIS_13TeVUp"));
////		singleH["PU"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PU_VIS_13TeVUp"));
////		singleH["BTagEff"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_BTagEff_VIS_13TeVUp"));
////		singleH["pT trig Eff"] = std::make_pair((TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVDown"),(TH1D*) f.second->Get("templates/bbH_Mbb_CMS_PtEff_VIS_13TeVUp"));
//		for(const auto & p: singleH){
//			++j;
//			TCanvas can("can","can");
//
//			std::string name2 = "box1_" + std::to_string(i) + "_" + std::to_string(j);
//			TPad cPad(name2.c_str(),"",0.,0.37,1,1);
//			cPad.SetBottomMargin(0.0);
//			cPad.SetLeftMargin(0.16);
//			cPad.SetRightMargin(0.05);
//			cPad.Draw();
//			cPad.cd();
//			if(f.first > 900) {lxmin = 0.25; lxmax = 0.55;}
//			TLegend legend(lxmin,lymin,lxmax,lymax);
//			(p.second).second->SetTitle((std::to_string(f.first)+ " GeV mass point, " + p.first+" variation").c_str());
//
//			(p.second).first->SetMarkerColor(kRed);
//			(p.second).second->SetMarkerColor(kBlue);
//			(p.second).first->SetMarkerStyle(20);
//			(p.second).second->SetMarkerStyle(20);
//
//			(p.second).second->Draw();
//			first->Draw("hist same");
//			(p.second).first->Draw("same");
//			gPad->RedrawAxis();
//
//			legend.AddEntry(first,"central","fl");
//			legend.AddEntry((p.second).first,("-2#sigma variation " + p.first).c_str(),"p");
//			legend.AddEntry((p.second).second,("+2#sigma variation " + p.first).c_str(),"p");
//			legend.Draw();
//
//			can.cd();
//			name2 = "box2_" + std::to_string(i) + "_" + std::to_string(j);
//			TPad uPad(name2.c_str(),"",0.,0.22,1,0.35);
//			uPad.SetBottomMargin(0.0);
//			uPad.SetTopMargin(0.);
//			uPad.SetLeftMargin(cPad.GetLeftMargin());
//			uPad.SetRightMargin(cPad.GetRightMargin());
//			uPad.Draw();
//			uPad.cd();
//
//			TH1D *centrU = (TH1D*) first->Clone("centrU");
//			double cVal, uVal, dVal, ecVal, euVal, edVal;
//			cVal = centrU->IntegralAndError(centrU->FindFirstBinAbove(0),centrU->FindLastBinAbove(0),ecVal);
//			uVal = (p.second).second->IntegralAndError((p.second).second->FindFirstBinAbove(0),(p.second).second->FindLastBinAbove(0),euVal);
//			centrU->Divide((p.second).second);
//			double err_u;
//			for(auto i = 0; i < centrU->GetNbinsX(); ++i){
//				err_u = correlatedDivision(first->GetBinContent(i+1),first->GetBinError(i+1),(p.second).second->GetBinContent(i+1),(p.second).second->GetBinError(i+1));
//				centrU->SetBinError(i+1,err_u);
//			}
//			centrU->SetTitle("");
//			centrU->SetMarkerStyle(20);
//			SetBottomStyle(centrU);
//
//			centrU->Draw();
//
//			TLine *horizLine = new TLine(centrU->GetXaxis()->GetXmin(),1,centrU->GetXaxis()->GetXmax(),1);
//			horizLine -> SetLineStyle(2);
//			horizLine -> Draw();
//
//			can.cd();
//			name2 = "box3_" + std::to_string(i) + "_" + std::to_string(j);
//			TPad dPad(name2.c_str(),"",0.,0.,1,0.2);
//			dPad.SetBottomMargin(0.33);
//			dPad.SetTopMargin(0.);
//			dPad.SetLeftMargin(cPad.GetLeftMargin());
//			dPad.SetRightMargin(cPad.GetRightMargin());
//			dPad.Draw();
//			dPad.cd();
//
//			TH1D *centrD = (TH1D*) first->Clone("centrD");
//			dVal = (p.second).first->IntegralAndError((p.second).first->FindFirstBinAbove(0),(p.second).first->FindLastBinAbove(0),edVal);
//			lnN[p.first][f.first] = make_pair(cVal/dVal,uVal/cVal);
//			elnN[p.first][f.first] = make_pair(correlatedDivision(cVal,ecVal,dVal,edVal),correlatedDivision(uVal,euVal,cVal,ecVal));
//			av_lnN[p.first][f.first] = 1. + std::abs((lnN[p.first][f.first].second + lnN[p.first][f.first].first - 2)/4.);
//			av_elnN[p.first][f.first] = 0.25 * sqrt(elnN[p.first][f.first].second * elnN[p.first][f.first].second + elnN[p.first][f.first].first * elnN[p.first][f.first].first - 2*elnN[p.first][f.first].second*elnN[p.first][f.first].first);
//			//calculate errors of division:
//			double err_d;
//			centrD->Divide((p.second).first);
//			for(auto i = 0; i < centrD->GetNbinsX(); ++i){
//				err_d = correlatedDivision(first->GetBinContent(i+1),first->GetBinError(i+1),(p.second).first->GetBinContent(i+1),(p.second).first->GetBinError(i+1));
//				centrD->SetBinError(i+1,err_d);
//			}
//			centrD->SetTitle("");
//			centrD->SetMarkerStyle(20);
//			SetBottomStyle(centrD);
//			centrD->Draw();
//
//			TLine *horizLine2 = new TLine(centrD->GetXaxis()->GetXmin(),1,centrD->GetXaxis()->GetXmax(),1);
//			horizLine2 -> SetLineStyle(2);
//			horizLine2 -> Draw();
//
//			cPad.cd();
//			if (f.first > 900) {xmin = 0.25; xmax = 0.55;}
//			TPaveText text(xmin,ymin,xmax,ymax,"NDC");
//			text.SetFillColor(0);
//			text.SetShadowColor(0);
//			text.AddText(("lnN +2#sigma: " + std::to_string(lnN[p.first][f.first].second)).c_str());
//			text.AddText(("lnN -2#sigma: " + std::to_string(1. / lnN[p.first][f.first].first)).c_str());
//			text.AddText(("lnN #pm1#sigma: " + std::to_string(av_lnN[p.first][f.first])).c_str());
//			//("lnN +2\sigma: " + std::to_string(lnN[i-1][j-1].first) + " \pm " + std::to_string(elnN[i-1][j-1].first)  + " \splitline lnN -2#sigma: " + std::to_string(lnN[i-1][j-1].second) + " #pm " + std::to_string(elnN[i-1][j-1].second)).c_str());
//			text.Draw();
//			std::string name = "pictures/"+std::to_string(f.first)+ "_GeV" + "_mass_point_" +regime + "_" + p.first+"_variation.pdf";
//			can.SaveAs(name.c_str());
//		}
//		j=0;
//	}
//
//	//Make plot and fit to lnN of different mass points
//	std::vector<string> Syst = {"SFb","SFl","PU"};
//	for(const auto& s : Syst){
//		calcLnN(s,av_lnN[s],av_elnN[s]);
//	}
//}
/**/
void calcLnN(const std::string& sys, const std::map<int,double > & lnN, const std::map<int,double > & elnN){
	gStyle->SetOptFit(1111);
	TCanvas can("can","can",800,600);
	size_t npoints = lnN.size();
	double x[npoints];
	double y[npoints];
	double ey[npoints];
	int i =0;
	for(const auto& v : lnN) {
		x[i] = v.first;
		y[i] = v.second;
		++i;
	}
	i = 0;
	for(const auto&v: elnN){
		ey[i] = (v.second);
		++i;
	}
	TGraphErrors gr(npoints,x,y,0,ey);
	gr.SetMarkerStyle(20);
	gr.SetMarkerSize(1.2);
	gr.SetTitle(("fit to lnN for " + sys + "; M_{12}, GeV; lnN").c_str());
	gr.Draw("AP");
	gr.Fit("pol0");
	std::string name = "pictures/lnN_" + sys + "_variation.pdf";
	can.SaveAs(name.c_str());

//	can.SaveAs(("template_lnN_" + std::string(h.GetName()) + ".pdf").c_str());
}
