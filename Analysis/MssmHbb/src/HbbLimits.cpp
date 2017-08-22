/*
 * HbbLimits.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/HbbLimits.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

namespace analysis {
namespace mssmhbb {

HbbLimits::HbbLimits() :
		blindData_(true),
		TEST_(false),
		xMin_(200),
		xMax_(1300),
		yMin_(0.1),
		yMax_(30)
		{
};
HbbLimits::HbbLimits(const bool& blindData, const double& xMin, const double& xMax, const double& yMin, const double& yMax, const bool& test) :
		blindData_(blindData),
		TEST_(test),
		xMin_(xMin),
		xMax_(xMax),
		yMin_(yMin),
		yMax_(yMax) {
};

HbbLimits::~HbbLimits() {
	// TODO Auto-generated destructor stub
}

void HbbLimits::ReadCombineLimits(const std::string& file_name){
	/*
	 * Method to read Limits calculated with Asymptotic method in the combine tool
	 * and written to the file.
	 */
//	input stream:
	std::ifstream inputList(file_name);
	if(is_empty(inputList)) std::logic_error("ERROR in HbbLimits::ReadCombineLimits. WRONG INPUT FILE NAME");
//	TFile name
	std::string tfileName;
	std::vector<Limit> GxBr_limits{};
//	Loop over the lines in the input file:
	while (inputList >> tfileName) {
		Limit limit = ReadCombineLimit(tfileName,blindData_);
		GxBr_limits.push_back(limit);
	}
	setLimits(GxBr_limits);
}

void HbbLimits::LimitPlotter(
		TLegend& leg,
		const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){
	/*
	 * Method to avoid null_vec creation by user
	 */
	LimitsToCompare differ_limits;
	LimitPlotter(differ_limits,leg,output,Lumi,xtitle,ytitle,logY);
}

void HbbLimits::LimitPlotter(
		const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){
	/*
	 * Method to avoid null_vec creation by user
	 */
	LimitsToCompare differ_limits;
	TLegend leg(0.62,0.6,0.9,0.85);
	LimitPlotter(differ_limits,leg,output,Lumi,xtitle,ytitle,logY);
}

void HbbLimits::LimitPlotter(const LimitsToCompare& differ_limits,
		TLegend& leg,
		const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){


	if(limits_.size() == 0) {
		std::cerr<<"Error in HbbLimits::LimitPlotter: No limits with this name. Please check spelling";
		exit(-1);
	}
	HbbStyle style;

	gStyle->SetPadTopMargin   (0.08);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);

	//Write value of limits in file
//	CheckOutputDir(output);
//	Write(output);

	int nPointsX = limits_.size();
	std::vector<double> X;
	std::vector<double> obs;
	std::vector<double> median;
	std::vector<double> minus1;
	std::vector<double> minus2;
	std::vector<double> plus1;
	std::vector<double> plus2;
	std::vector<double> zero;

	std::cout<<"SIZE: "<<nPointsX<<std::endl;
	for(const auto& l : limits_){
		std::cout<<"x = "<<l.getX()<<" obs = "<<l.getObserved()<<" exp = "<<l.getMedian()<<std::endl;
		X.push_back(l.getX());
		obs.push_back(l.getObserved());
		median.push_back(l.getMedian());
		minus2.push_back(l.getMedian() - l.getMinus2G());
		minus1.push_back(l.getMedian() - l.getMinus1G());
		plus2.push_back(l.getPlus2G() - l.getMedian());
		plus1.push_back(l.getPlus1G() - l.getMedian());
		zero.push_back(0);

	}

	TGraph * differ_exp = nullptr;
	bool compare_limits = false;
	if(differ_limits.limits.size() != 0) compare_limits = true;
	std::vector<double> X_2;
	std::vector<double> median_2;
	if(compare_limits){
		for(const auto& l: differ_limits.limits){
			X_2.push_back(l.getX());
			median_2.push_back(l.getMedian());
		}
		differ_exp = new TGraph(differ_limits.limits.size(),X_2.data(),median_2.data());
		differ_exp->SetLineWidth(3);
		differ_exp->SetLineColor(1);
		differ_exp->SetLineStyle(3);
	}

	TGraph * obsG = new TGraph(nPointsX, X.data(), obs.data());
	obsG->SetLineWidth(3);
	obsG->SetLineColor(1);
	obsG->SetLineWidth(2);
	obsG->SetMarkerColor(1);
	obsG->SetMarkerStyle(20);
	obsG->SetMarkerSize(1.4);

	TGraph * expG = new TGraph(nPointsX, X.data(),median.data());
	expG->SetLineWidth(3);
	expG->SetLineColor(2);
	expG->SetLineStyle(2);

	TGraphAsymmErrors * innerBand = new TGraphAsymmErrors(nPointsX, X.data(), median.data(), zero.data(), zero.data(), minus1.data(), plus1.data());
	innerBand->SetFillColor(kGreen);
	innerBand->SetLineColor(kGreen);

	TGraphAsymmErrors * outerBand = new TGraphAsymmErrors(nPointsX, X.data(), median.data(), zero.data(), zero.data(), minus2.data(), plus2.data());
	outerBand->SetFillColor(kYellow);
	outerBand->SetLineColor(kYellow);

	TH2F * frame = nullptr;

	frame = new TH2F("frame","",2,xMin_,xMax_,2,yMin_,yMax_);
	frame->GetXaxis()->SetTitle(xtitle.c_str());
	frame->GetYaxis()->SetTitle(ytitle.c_str());
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetNdivisions(206);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetXaxis()->SetTitleOffset(1.05);
	frame->GetYaxis()->SetTitleSize(0.048);
	frame->SetStats(0);

	TCanvas *canv = new TCanvas("canv", "histograms", 600, 600);

	frame->Draw();

	leg.SetFillColor(0);
	leg.SetTextSize(0.035);
	leg.SetLineWidth(2);
	leg.SetBorderSize(1);

	outerBand->Draw("3same");
	innerBand->Draw("3same");
	expG->Draw("lsame");
	if(compare_limits) differ_exp->Draw("lsame");

	if(compare_limits) {
		leg.AddEntry(differ_exp,differ_limits.legend.c_str(),"l");
	}
	AddPlottingObjects(*frame,leg,*obsG,*expG,*innerBand,*outerBand);

	TPad * pad = (TPad*)canv->GetPad(0);
	pad->RedrawAxis();
	leg.Draw();
	style.drawStandardTitle();

	if(logY) canv->SetLogy();
	canv->Update();
	canv->Print( (output+".pdf").c_str() ,"Portrait pdf");
	canv->Print( (output+".png").c_str()) ;
}

void HbbLimits::AddPlottingObjects(TH2F &frame, TLegend &leg, TGraph& obs, TGraph& exp, TGraphAsymmErrors& inner_band, TGraphAsymmErrors& outer_band){
	/*
	 * Virtual function that is actually plotting objects on the frame
	 */
	if (!blindData_){
		obs.Draw("lpsame");
		leg.AddEntry(&obs,"Observed","lp");
	}

	leg.AddEntry(&exp,"Expected","l");
	leg.AddEntry(&inner_band,"#pm1#sigma Expected","f");
	leg.AddEntry(&outer_band,"#pm2#sigma Expected","f");

}

const Limit HbbLimits::ReadCombineLimit(const std::string& tfile_name, const bool& blindData){
	/*
	 * Method to get information about limit from root file
	 */
	TFile file(tfile_name.c_str(),"read");
//	Check whether file is Ok.
	CheckZombie(file);
//	Get TTree from this file
	TTree& tree = *((TTree*)file.Get("limit"));
	double LIMIT;
	tree.SetBranchAddress("limit",&LIMIT);
	double MH;
	tree.SetBranchAddress("mh",&MH);
//	Get iinformation from TTree:
    Limit limit;

    tree.GetEntry(0);
    limit.setX(double(MH));
    limit.setMinus2G(double(LIMIT));

    tree.GetEntry(1);
    limit.setMinus1G(double(LIMIT));

    tree.GetEntry(2);
    limit.setMedian(double(LIMIT));

    tree.GetEntry(3);
    limit.setPlus1G(double(LIMIT));

    tree.GetEntry(4);
    limit.setPlus2G(double(LIMIT));

    tree.GetEntry(5);
    limit.setObserved(double(LIMIT));
//    if(blindData) limit.setObserved(limit.getExpected());

    return limit;
}

void HbbLimits::Write(const std::string& name){
	/*
	 * Method to write limits to a .txt file
	 */
	std::fstream fs;
	fs.open( (name + ".txt").c_str(), std::fstream::out );
	fs<<"X     -2s     -1s  median     +1s     +2s     obs \n";
	for(const auto& limit: limits_){
		char strOut[400];
		sprintf(strOut,"%4f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f",
				(limit.getX()),limit.getMinus2G(),limit.getMinus1G(),limit.getExpected(),limit.getPlus1G(),limit.getPlus2G(),limit.getObserved());
		fs<<strOut;
		fs<<"\n";
	}
	fs.close();
	std::cout<<"File: "<<name + ".txt"<<" has been written."<<std::endl;
}

} /* namespace MssmHbb */
} /* namespace Analysis */
