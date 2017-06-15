#include <TH3F.h>
#include <TH3.h>
#include <TGraph.h>
#include <TFile.h>
#include "TSystem.h"
#include <TCanvas.h>
#include <TLegend.h>


#include <string>
#include <vector>

// xsection estimation for MSSM
#include "Analysis/MssmHbb/macros/signal/mssm_xs_tools.h"
#include "Analysis/MssmHbb/macros/signal/mssm_xs_tools.C"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
// Tools
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace std;

const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));

void CompareBrs(const vector<double>&, const vector<double>&, const double& , const vector<string>&, const string&, const string&);
double MssmBr(mssm_xs_tools& my,
		     const double& mA,
		     const double& tanBeta);
double THDMBr(map<string,TH3F*>& f,
				   const double& mA,
				   const double& tanBeta,
				   const double& sinB_A,
				   const string& name);
void DrawGraphs(vector<TGraph>& graphs2HDM, const vector<string> &BRs, const double& sinB_A);
void GetTHDMHistos(TFile& name,map<string,TH3F*>&,const vector<string>&);

string boson = "H";
string mA = "1100";

int main(){
	HbbStyle style;
	style.set(PRIVATE);
//	Comparison should be done
//	Histo with 2HDM results:
	string thdm_production = "production_cosB_A_-1_1_tanB_0p5-100_COMBINATION"; //production_corseBins_cosB_A_-1_1_tanB_1-100 //
	string thdm_type = "type3";
	string thdm_fname = "/nfs/dust/cms/user/shevchen/SusHiScaner/output/" + thdm_production + "/rootFiles/Histograms3D_" + thdm_type + "_mA_mH.root";
	TFile f(thdm_fname.c_str(),"read");
//	Define mass points:
	vector<double> mAs{300,400,500,600,700,900,1100};
//	Define tnaBeta points:
	vector<double> tanBetas{10,20,30,40,50,60};
//	Define sin(beta-alpha) points:
	vector<double> sinB_As;
	for(double i = -0.99; i <= 0.99; i+=0.1) sinB_As.push_back(i);
//	Define types of BRs:
	vector<string> BRs = {"br_" + boson + "ss","br_" + boson + "cc","br_" + boson + "bb","br_" + boson + "tt","br_" + boson + "tautau","br_" + boson + "Zh"};

//	MSSM benchm,ark file
	string mssm_bench = cmsswBase + "/src/Analysis/MssmHbb/macros/signal/mhmodp_mu200_13TeV.root";
	CompareBrs(sinB_As,tanBetas,stoi(mA),BRs,thdm_fname,mssm_bench);
	return 0;
}

void CompareBrs(const vector<double> &sinB_As, const vector<double> &tanBetas, const double&mA, const vector<string>& BRs, const string &thdm_fname,const string &mssm_benchmark){
	//	File for 2HDM scans:
	map<string,TH3F*> thdm_map;
	TFile f(thdm_fname.c_str(),"read");
	GetTHDMHistos(f,thdm_map,BRs);
	// Setup mssm_xs_tools:
	mssm_xs_tools my(mssm_benchmark.c_str(),true,0);
	// Setup mssm_xs_tools:
//	Start mA, tan(beta) and sin(beta-alpha) loops:
	for(const auto& tanBeta : tanBetas){
		// Prepare output TGraphs
		vector<TGraph> graphs2HDM(BRs.size());
		vector<TGraph> graphsMSSM(BRs.size());
		int i = 0;
		for(const auto& BR : BRs){
			TGraph& gr_2hdm = graphs2HDM.at(i);
			TGraph& gr_mssm = graphsMSSM.at(i);
//			Set number of points for TGraph
			gr_2hdm.Set((int)sinB_As.size());
			gr_mssm.Set((int)sinB_As.size());
			++i;
			int j =0;
			for(const auto& sinB_A : sinB_As){
				//Calcualte SigmaBr and fill TGraphs:
//				gr_mssm.SetPoint(j,tanBeta,MssmBr(my,mA,tanBeta));
				gr_2hdm.SetPoint(j,sinB_A,THDMBr(thdm_map,mA,tanBeta,sinB_A,BR));
				++j;
			}
		}
		DrawGraphs(graphs2HDM,BRs,tanBeta);
	}
}

void DrawGraphs(vector<TGraph>& graphs2HDM, const vector<string> &BRs, const double& tanBeta){
	// setup canvas
	TCanvas can("can","can",800,600);
	can.SetGrid();
	double xmin, xmax,ymin,ymax;
	xmin = 0.7;xmax = 0.93;
	ymin = 0.3; ymax = 0.55;
	if(tanBeta <= 20){
		ymin = 0.6; ymax = 0.8;
	}
	if(boson == "H"){
		xmin = 0.74; xmax = 0.94;
		ymin = 0.75; ymax = 0.915;
		if(tanBeta <= 30 ){
			xmin = 0.21; xmax = 0.38;
			ymin = 0.72; ymax = 0.915 ;
		}
	}
	TLegend leg(xmin,ymin,xmax,ymax);
	vector<int> colours = {1,3,2,4,6,9,12,28,46};
	// Loop over the tanBetas
	for (unsigned int i = 0; i != BRs.size(); ++i){
		TGraph& thdm = graphs2HDM.at(i);
		auto br = BRs.at(i);

		//Cosmetics:
		thdm.SetMarkerStyle(20 + i);
		thdm.SetMarkerSize(1.04);
		thdm.SetMarkerColor(colours.at(i));
		thdm.SetLineColor(colours.at(i));
		thdm.SetTitle( ("2HDM(tan#beta = " + to_string_with_precision(tanBeta,2) + ")" + ";cos(#beta-#alpha); BR").c_str());

		leg.AddEntry(&thdm,("2HDM , " + br).c_str(),"pl");

		if(i==0){
			thdm.Draw("APL");
		}
		else thdm.Draw("PLsame");
		if(i==0) thdm.GetYaxis()->SetRangeUser(0.00000001,1);
		if(i==0) thdm.GetXaxis()->SetRangeUser(-1,1.);
		can.SetLogy();
		can.Update();
	}
	leg.Draw();
	can.Print( (cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/2HDM/" + boson + "_thdm_BRs_sinB_A_mA-" + mA + "_tanB-" + to_string_with_precision(tanBeta,2) + ".pdf").c_str() );
}

double THDMBr(map<string,TH3F*>& f,
				   const double& mA,
				   const double& tanBeta,
				   const double& sinB_A,
				   const string& name){
	double BrAbb = f[name]->Interpolate(mA,tanBeta,sinB_A);

	return BrAbb;
}

void GetTHDMHistos(TFile& f, map<string,TH3F*>& histos, const vector<string>& BRs){
	CheckZombie(f);
	for(const auto& br: BRs){
		CheckZombieObjectInTFile(f,br);
		histos[br] = (TH3F*) f.Get(br.c_str());
	}

}

