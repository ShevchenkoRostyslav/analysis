#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TH2.h>
#include <TGraph.h>
#include <TPad.h>
#include <TFile.h>
#include "TSystem.h"

#include "Analysis/MssmHbb/interface/LHCXSGLimits.h"
#include "Analysis/MssmHbb/interface/Limit.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"

using namespace std;
using namespace analysis::mssmhbb;

class Uncertainty {
public:
	std::vector<Limit> central;
	std::vector<Limit> up;
	std::vector<Limit> down;
	std::string name;
};

Uncertainty evaluate(LHCXSGLimits& limits, const std::string& path_to_limits,const std::string& name);
void DrawUncertainty(const Uncertainty&, const std::string& output);
void fillTGraph(TGraph& gr, const std::vector<Limit>& limits);
void divideTGraphs(TGraph& finale, TGraph& up, TGraph& down);
void estimateAverageUnc(TGraph& up, TGraph& down, ofstream & file);

int main(int argc, char **argv) {
	/*
	 * Macro to estimate size of effect of alpha_s+PDF on MSSM xsection
	 * to be used as lnN nuisance for MSSM limits
	 */
	//Set style
	HbbStyle style;
	style.set(PRELIMINARY);
	string path_to_limits = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/datacards/201708/23/unblinded/independent/bias/Hbb.limits";
	string interpretation = "hMSSM_13TeV";//hMSSM_13TeV//mhmodp_mu200_13TeV//lightstopmod_13TeV//lightstau1_13TeV//tauphobic_13TeV//
	string output = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/TheoryUncertainty/" + interpretation + "_";
	string benchmark = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/macros/signal/" + interpretation + ".root";
	LHCXSGLimits limits(mssmhbb::blinded,"both",benchmark);

	//Read output of the combine tool:
	Uncertainty pdfas = evaluate(limits,path_to_limits,"pdfas");
	Uncertainty scale = evaluate(limits,path_to_limits,"scale");

	//Plot
	DrawUncertainty(pdfas,output);
	DrawUncertainty(scale,output);

}

Uncertainty evaluate(LHCXSGLimits& limits, const std::string& path_to_limits,const std::string& name){
	/*
	 * Method to evaluate uncertainty for a given parameter "name"
	 */
	Uncertainty unc;
	limits.ReadCombineLimits(path_to_limits,"");
	unc.central = limits.getLimits();
	limits.ReadCombineLimits(path_to_limits,name,true); //UP
	unc.up = limits.getLimits();
	limits.ReadCombineLimits(path_to_limits,name,false); //DOWN
	unc.down = limits.getLimits();
	unc.name		= name;

	return unc;
}

void DrawUncertainty(const Uncertainty& unc, const std::string& output){
	/*
	 * function to draw graphs dedicated to central, up and down uncertainties
	 */
	TCanvas can("can","can",800,600);
//	gPad->SetLogy();
	TGraph central, up, down;
	fillTGraph(central,unc.central);
	fillTGraph(up,unc.up);
	fillTGraph(down,unc.down);

	central.SetTitle( (unc.name + " variation;M_{A}, [GeV]; tan#beta").c_str() );
	up.SetLineColor(kBlue);
	down.SetLineColor(kRed);

	//Pad1 with limits itself

    TPad pad1("pad1","",0,0.1,1,1);
    pad1.SetBottomMargin(0.2);
    pad1.SetRightMargin(0.05);
    pad1.SetLeftMargin(0.16);
    pad1.Draw();
    pad1.cd();

    double xMin = central.GetXaxis()->GetXmin();
    double xMax = central.GetXaxis()->GetXmax();
    double yMin = 0;
    double yMax = 60;//central.GetYaxis()->GetXmax() + 0.1 * central.GetYaxis()->GetXmax();
	TH2F frame("frame","",2,xMin,xMax,2,yMin,yMax);
	frame.GetXaxis()->SetTitle(central.GetXaxis()->GetTitle());
	frame.GetYaxis()->SetTitle(central.GetYaxis()->GetTitle());
	frame.GetXaxis()->SetNdivisions(505);
	frame.GetYaxis()->SetNdivisions(206);
	frame.GetYaxis()->SetTitleOffset(1.3);
	frame.GetXaxis()->SetTitleOffset(999);
	frame.GetXaxis()->SetLabelOffset(999);
	frame.GetYaxis()->SetTitleSize(0.048);
	frame.SetStats(0);
	frame.Draw();

	central.Draw("lsame");
	up.Draw("lsame");
	down.Draw("lsame");
//	central.GetYaxis()->SetRangeUser(0,60);
	central.Draw("lsame");

	TLegend leg(0.6,0.3,0.9,0.5);
	leg.SetBorderSize(0);
	leg.AddEntry(&central,"Central","l");
	leg.AddEntry(&up,("+1 #sigma, " + unc.name).c_str(),"l");
	leg.AddEntry(&down,("-1 #sigma, " + unc.name).c_str(),"l");
	leg.Draw();

	//Pad2 with ratios

    can.cd();
    TPad pad2("pad2","",0,0.0,1,0.265);
    pad2.SetTopMargin(1);
    pad2.SetBottomMargin(0.33);
    pad2.SetLeftMargin(pad1.GetLeftMargin());
    pad2.SetRightMargin(pad1.GetRightMargin());
    pad2.SetGridy();
    pad2.Draw();
    pad2.cd();

    TH2F frame2("frame2","",2,xMin,xMax,2,0.9,1.1);
    frame2.SetTitle("");
    frame2.GetXaxis()->SetTitle(frame.GetXaxis()->GetTitle());
    frame2.GetXaxis()->SetTitleSize(0.15);
    frame2.GetXaxis()->SetTitleOffset(1.06);
    frame2.GetXaxis()->SetLabelSize(0.215);
    frame2.GetXaxis()->SetLabelOffset(0.010);
    frame2.GetXaxis()->SetNdivisions(505);
    frame2.GetXaxis()->SetTickLength(frame.GetXaxis()->GetTickLength()*3);
    frame2.SetYTitle("#frac{up(down)}{central}");
    frame2.GetYaxis()->CenterTitle(kTRUE);
    frame2.GetYaxis()->SetTitleSize(0.14);
    frame2.GetYaxis()->SetTitleOffset(0.4);
    frame2.GetYaxis()->SetNdivisions(206);
    frame2.GetYaxis()->SetLabelSize(0.115);
    frame2.GetYaxis()->SetLabelOffset(0.011);
    frame2.Draw();

    TGraph rat_c_up,rat_d_c;
    divideTGraphs(rat_c_up,up,central);
    divideTGraphs(rat_d_c,down,central);


    rat_c_up.SetLineColor(up.GetLineColor());
    rat_d_c.SetLineColor(down.GetLineColor());
    rat_c_up.Draw("lsame");
    rat_d_c.Draw("lsame");
    cout<<"Average for "<<unc.name<<endl;
    //output .txt file with uncertainties
    ofstream file;
    file.open( output + unc.name + ".txt");
    file<<"Average for "<<unc.name<<"\n";
    estimateAverageUnc(rat_c_up,rat_d_c,file);
    file.close();

	can.Print( (output + unc.name + ".pdf").c_str() );
}

void fillTGraph(TGraph& gr, const std::vector<Limit>& limits){
	/*
	 * Function to fill TGraph from vector of Limits
	 */
	int i = 0;
	for(const auto& l : limits){
		gr.SetPoint(i,l.getX(),l.getMedian());
		++i;
	}
}

void divideTGraphs(TGraph& finale, TGraph& gr1, TGraph& gr2){
	/*
	 * Function to divide two TGraphs with the same binning and x-range
	 */
	for(int i = 0;i<gr1.GetN();++i){
		finale.SetPoint(i,gr1.GetX()[i],gr1.GetY()[i]/gr2.GetY()[i]);
	}

}

void estimateAverageUnc(TGraph& up, TGraph& down, ofstream & file){
	/*
	 * Method to estimate the avarage lnN for each mass point
	 */
	for(int i = 0; i < up.GetN(); ++i){
		cout<<"x = "<<up.GetX()[i]<<" lnN: "<<1. + (1. - up.GetY()[i] + down.GetY()[i] - 1.) / 2.<<endl;
		file<<"x = "<<up.GetX()[i]<<" lnN: "<<1. + (1. - up.GetY()[i] + down.GetY()[i] - 1.) / 2.<<"\n";
	}
}
