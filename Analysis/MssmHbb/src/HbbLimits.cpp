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

bool is_empty(std::ifstream&);

HbbLimits::HbbLimits() : blindData_(true), boson_("both"), TEST_(false) {
//	HbbStyle* style = new HbbStyle();
	style_.set(PRIVATE);
//	style->SetStyle();
//	gStyle->SetOptFit(0000);
//	gStyle->SetErrorX(0.5);
};
HbbLimits::HbbLimits(const bool& blindData, const bool& test) : blindData_(blindData), TEST_(test) {
	style_.set(PRIVATE);
//	HbbStyle* style = new HbbStyle();
//	style->set(PRIVATE);
//	gStyle->SetOptFit(0000);
//	gStyle->SetErrorX(0.5);
};

HbbLimits::HbbLimits(const bool& blindData, const std::string& boson, const bool& test) : blindData_(blindData), boson_(boson), TEST_(test) {
	/*
	 * Constructor with 2 arguments.
	 * Higgs boson name has to be:
	 * A or H or both
	 */
	CheckHiggsBoson();
	style_.set(PRIVATE);
//	HbbStyle* style = new HbbStyle();
//	style->set(PRIVATE);
//	gStyle->SetOptFit(0000);
//	gStyle->SetErrorX(0.5);
};

HbbLimits::~HbbLimits() {
	// TODO Auto-generated destructor stub
}

std::vector<Limit> HbbLimits::ReadCombineLimits(const std::string& file_name){
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
	return GxBr_limits;
}

const std::vector<Limit> HbbLimits::GetMSSMLimits(const std::vector<Limit>& GxBR_limits, const std::string& benchmark_path, const std::string& uncert, const bool& UP, const std::string& benchmark_ref_path, const double& tanBref){
	/*
	 * Method to get MSSM interpretation of GxBR limits
	 */

	std::vector<Limit> limits;
	for(unsigned int i = 0; i != GxBR_limits.size(); ++i){
		Limit gxbr_limit = GxBR_limits.at(i);
		Limit tan_b_limit;
//		Set Mass point
		tan_b_limit.setX(gxbr_limit.getX());
//		Set limits
		tan_b_limit.setMinus2G(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getMinus2G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setMinus1G(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getMinus1G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setMedian(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getMedian(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setPlus1G(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getPlus1G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setPlus2G(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getPlus2G(),uncert,UP,benchmark_ref_path,tanBref)));
		tan_b_limit.setObserved(double(MSSMTanBeta(benchmark_path,tan_b_limit.getX(),gxbr_limit.getObserved(),uncert,UP,benchmark_ref_path,tanBref)));

		limits.push_back(tan_b_limit);
	}
	return limits;
}

TH3D HbbLimits::Get2HDM_GxBR_3D(const std::string& benchmark_path){
	/*
	 * Method to get 2HDM 3D GxBR
	 */
	//Output histo
	TH3D *hOut = nullptr;
	//map with histos:
	std::map<std::string,TH3D*> histos;
	//map with projections
	std::map<std::string,TH2*> proj;
	//Open TFile with 2HDM scans:
	TFile f(benchmark_path.c_str(),"READ");
	CheckZombie(f);
	SetTHDMHistos(f,histos);

	if(boson_ != "both"){
		histos["GxBR_" + boson_] = (TH3D*) histos["xs_bb" + boson_]->Clone();
		histos["GxBR_" + boson_]->Multiply(histos["br_" + boson_ + "bb"]);
		hOut = histos["GxBR_" + boson_];
	}
	else {
		histos["GxBR_A"] = (TH3D*) histos["xs_bbA"]->Clone();
		histos["GxBR_A"]->Multiply(histos["br_Abb"]);
		histos["GxBR_H"] = (TH3D*) histos["xs_bbH"]->Clone();
		histos["GxBR_H"]->Multiply(histos["br_Hbb"]);

		//Add HxBR + AxBR
		histos["GxBR_HandA"] = (TH3D*) histos["GxBR_A"]->Clone("GxBR_HandA");
		histos["GxBR_HandA"]->Add(histos["GxBR_H"],1.);
		hOut = histos["GxBR_HandA"];
	}

	if(TEST_){
		std::cout<<" HbbLimits::Get2HDM_GxBR_3D - DONE"<<std::endl;
	}

	gStyle->SetPadLeftMargin(0.08);
	gStyle->SetPadRightMargin(0.12);

	return *hOut;

}

TH2D HbbLimits::Get2HDM_GxBR_2D(const TH3& h_in, const double& var, const std::string& axis){
	/*
	 * Method to get 2D 2HDM GxBR for particular value of variable VAR, that correspondse to the axis - axis:
	 * X - mA
	 * Y - tan(b)
	 * Z - cos(b-a)
	 */
	//Copy input h
	TH3D h = *(TH3D*) h_in.Clone();
	//Out h;
	TH2D *proj = nullptr;
	//Var Range has to be setuped
	int bin = -10;
	if(axis == "x" || axis == "X") {
		bin = h.GetXaxis()->FindBin(var);
		h.GetXaxis()->SetRange(bin,bin);
		//Project
		proj = (TH2D*)  h.Project3D("yz");
	}
	if(axis == "y" || axis == "Y") {
		bin = h.GetYaxis()->FindBin(var);
		h.GetYaxis()->SetRange(bin,bin);
		//Project
		proj = (TH2D*)  h.Project3D("zx");
	}
	if(axis == "z" || axis == "Z") {
		bin = h.GetZaxis()->FindBin(var);
		h.GetZaxis()->SetRange(bin,bin);
		//Project
		proj = (TH2D*)  h.Project3D("yx");
	}

	if(TEST_){
		std::cout<<" HbbLimits::Get2HDM_GxBR_2D - DONE"<<std::endl;
	}

	//Special Pad margins to plot 2D mu scans
	gStyle->SetPadTopMargin   (0.09);
	gStyle->SetPadLeftMargin(0.14);
	gStyle->SetPadRightMargin(0.15);

	return *proj;
}

HbbLimits::THDMScan HbbLimits::Get2HDMmuScan(const TH2& GxBR_2hdm, const Limit& GxBR_95CL){
	/*
	 * Method to make a Signal strength scan of 2HDM model for particular value of VAR
	 */
	//Clone input histo
	TH2D exp = *(TH2D*) GxBR_2hdm.Clone();
	TH2D obs = *(TH2D*) GxBR_2hdm.Clone();
	double GxBR_exp  = GxBR_95CL.getMedian();
	double GxBR_obs  = GxBR_95CL.getMedian();

	//Expected
	exp.Scale(1./GxBR_exp);
	//Observed
	obs.Scale(1./GxBR_obs);

	exp.SetTitle("#mu = (#sigma x BR)_{theory} / (#sigma x BR)_{95%C.L.}");
	exp.SetLabelSize(0.050,"XYZ");
//	exp.SetLabelSize(0.050,"Y");
	exp.SetTitleSize(0.055,"XYZ");
//	exp.SetTitleSize(0.055,"Y");

	THDMScan out;
	out.expected = (TH2D) exp;
	out.observed = (TH2D) obs;

	return out;
}

std::vector<Limit> HbbLimits::Get2HDM_Limits(const TH2& h_in, const Limit& GxBR_95CL, const double& min, const double& max){
	/*
	 * Method to calculate 2HDM Limits
	 * and plot them
	 */

	if(TEST_){
		std::cout<<"START HbbLimits::Get2HDM_Limits"<<std::endl;
	}

	//Clone input h
	TH2D GxBR_2hdm = *(TH2D*) h_in.Clone();
	GxBR_2hdm.GetXaxis()->SetRangeUser(min,max);
	double Limit_95CL[6];
	Limit_95CL[0] = GxBR_95CL.getMinus2G();
	Limit_95CL[1] = GxBR_95CL.getMinus1G();
	Limit_95CL[2] = GxBR_95CL.getMedian();
	Limit_95CL[3] = GxBR_95CL.getPlus1G();
	Limit_95CL[4] = GxBR_95CL.getPlus2G();
	Limit_95CL[5] = GxBR_95CL.getObserved();

	TCanvas canva("canva","",800,600);
	GxBR_2hdm.DrawCopy("colz");
	GxBR_2hdm.SetContour(6,Limit_95CL);
	GxBR_2hdm.Draw("cont Z LIST");
	canva.Update();

	// Get Contours
	TObjArray *conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
	gROOT->GetListOfSpecials()->Print();

	if (conts == nullptr){
		std::logic_error("ERROR in HbbLimits::Get2HDM_Limits: *** No Contours Were Extracted!");
	}

	TGraph *innerBand_up 	= nullptr;
	TGraph *innerBand_down = nullptr;
	TGraph *outerBand_up = nullptr;
	TGraph *outerBand_down = nullptr;
	TGraph *expected = nullptr;
	TGraph *observed = nullptr;

	outerBand_down	= (TGraph*) ((TList*)conts->At(0))->First();
	innerBand_down	= (TGraph*) ((TList*)conts->At(1))->First();
	expected		= (TGraph*) ((TList*)conts->At(2))->First();
	innerBand_up 	= (TGraph*) ((TList*)conts->At(3))->First();
	outerBand_up 	= (TGraph*) ((TList*)conts->At(4))->First();
	observed		= (TGraph*) ((TList*)conts->At(5))->First();

	if(!outerBand_down || !outerBand_up ){
		std::logic_error("ERROR in HbbLimits::Get2HDM_Limits: ***No 2G up or Down contours***");
	}
	if(!innerBand_down || !innerBand_up){
		std::logic_error("ERROR in HbbLimits::Get2HDM_Limits: ***No 1G up or Down contours***");
	}
	//Define future binning:
	int np = 201;
//	double max = GxBR_2hdm.GetXaxis()->GetXmax();
//	double min = GxBR_2hdm.GetXaxis()->GetXmin();

	double step = (max - min) / (np-1);
	//limit map:
	std::vector<Limit> limits_2hdm;
	//get values of limits:
	double in_up[202], in_down[202], exp[202], x[202], out_up[202], out_down[202], obs[202];
	for(int i = 0;i < np;++i) {
		x[i] = min + i*step;
		exp[i] = expected->Eval(x[i]);
		obs[i] = observed->Eval(x[i]);
		in_up[i] = innerBand_up->Eval(x[i]);
		in_down[i] = innerBand_down->Eval(x[i]);
		out_up[i] = outerBand_up->Eval(x[i]);
		out_down[i] = outerBand_down->Eval(x[i]);

		Limit l;
		l.setX(x[i]);
		l.setMedian(exp[i]);
		l.setObserved(obs[i]);
		l.setMinus1G(in_down[i]);
		l.setMinus2G(out_down[i]);
		l.setPlus1G(in_up[i]);
		l.setPlus2G(out_up[i]);
		limits_2hdm.push_back(l);
	}

	//Clean memory
//	delete conts;
//	delete innerBand_up;
//	delete innerBand_down;
//	delete outerBand_up;
//	delete outerBand_down;
//	delete expected;
//	delete observed;
//	delete []Limit_95CL;
//	delete []in_up;
//	delete []in_down;
//	delete []exp;
//	delete []obs;
//	delete []x;
//	delete []zero;
//	delete []out_up;
//	delete []out_down;

	if(TEST_){
		std::cout<<" HbbLimits::Get2HDM_Limits - DONE"<<std::endl;
	}

	return limits_2hdm;
}

const std::vector<Limit> HbbLimits::Get2HDM_1D_Limits(const TH2& h_in,const std::vector<Limit>& GxBR_limits){
	/*
	 * Method to get 2HDM 1D interpretation of GxBR limits
	 */
	//Clone input h
	TH2D GxBR_2hdm = *(TH2D*) h_in.Clone();

	std::vector<Limit> limits;
	for(unsigned int i = 0; i != GxBR_limits.size(); ++i){
		Limit gxbr_limit = GxBR_limits.at(i);
		Limit tan_b_limit;
//		Set Mass point
		tan_b_limit.setX(gxbr_limit.getX());
//		Set limits
		tan_b_limit.setMinus2G(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getMinus2G())));
		tan_b_limit.setMinus1G(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getMinus1G())));
		tan_b_limit.setMedian(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getMedian())));
		tan_b_limit.setPlus1G(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getPlus1G())));
		tan_b_limit.setPlus2G(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getPlus2G())));
		tan_b_limit.setObserved(double(THDMTanBeta(GxBR_2hdm,tan_b_limit.getX(),gxbr_limit.getObserved())));

		limits.push_back(tan_b_limit);
	}
	return limits;
}

void HbbLimits::LimitPlotter(const std::vector<Limit>& limits,
		const std::vector<Limit>& differ_limits,
		TLegend& leg,
		const std::string& output,
		const float& yMin,
		const float& yMax,
		const float& xMin,
		const float& xMax,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){


	if(limits.size() == 0) {
		std::cerr<<"Error in HbbLimits::LimitPlotter: No limits with this name. Please check spelling";
		exit(-1);
	}

	gStyle->SetPadTopMargin   (0.08);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);

	//Write value of limits in file
//	CheckOutputDir(output);
	Write(limits,output);

	int nPointsX = limits.size();
	std::vector<double> X;
	std::vector<double> obs;
	std::vector<double> median;
	std::vector<double> minus1;
	std::vector<double> minus2;
	std::vector<double> plus1;
	std::vector<double> plus2;
	std::vector<double> zero;

	std::cout<<"SIZE: "<<nPointsX<<std::endl;
	for(const auto& l : limits){
		std::cout<<"x = "<<l.getX()<<std::endl;
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
	if(differ_limits.size() != 0) compare_limits = true;
	std::vector<double> X_2;
	std::vector<double> median_2;
	if(compare_limits){
		for(const auto& l: differ_limits){
			X_2.push_back(l.getX());
			median_2.push_back(l.getMedian());
		}
		differ_exp = new TGraph(differ_limits.size(),X_2.data(),median_2.data());
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

	frame = new TH2F("frame","",2,xMin,xMax,2,yMin,yMax);
	frame->GetXaxis()->SetTitle(xtitle.c_str());
	frame->GetYaxis()->SetTitle(ytitle.c_str());
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetNdivisions(206);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetXaxis()->SetTitleOffset(1.05);
	frame->GetYaxis()->SetTitleSize(0.048);
	frame->SetStats(0);

//	TGraph atlas_zhll;
//	bool thdm_type3 = false;
//	if(output.find("type3") != std::string::npos) {
//		atlas_zhll = GetAtlasZhll_flipped();
//		thdm_type3 = true;
//	}


	TCanvas *canv = new TCanvas("canv", "histograms", 600, 600);

	frame->Draw();

	outerBand->Draw("3same");
	innerBand->Draw("3same");
	expG->Draw("lsame");
	if(compare_limits) differ_exp->Draw("lsame");
	if (!blindData_)
		obsG->Draw("lpsame");

//	if(thdm_type3) {
//		atlas_zhll.Draw("PLsame");
//	}

	std::string entry_comp = "Standard";//"Prev. version";
	if(compare_limits) {
		if(Lumi.find("24.6") != std::string::npos) entry_comp = "Expected 7 + 8 TeV";
		else if (Lumi.find("19.7") != std::string::npos) entry_comp = "Expected 13 TeV";
		else if (Lumi.find("2.62") != std::string::npos) entry_comp = "Expected 2015";
		leg.AddEntry(differ_exp,entry_comp.c_str(),"l");
	}
	if (!blindData_)
		leg.AddEntry(obsG,"Observed","lp");
	leg.AddEntry(expG,"Expected","l");
	leg.AddEntry(innerBand,"#pm1#sigma Expected","f");
	leg.AddEntry(outerBand,"#pm2#sigma Expected","f");

	TPad * pad = (TPad*)canv->GetPad(0);
//	Luminosity lum;
//	lum.writeExtraText = true;
//	lum.lumi_13TeV = Lumi;
//	lum.extraText = "Private Work";
//	int offs = 33;
//	if(ytitle.find("tan") != std::string::npos) offs = 11;
//	lum.CMS_lumi(pad,4,offs);

	pad->RedrawAxis();

	leg.Draw();
	style_.drawStandardTitle();

	if(logY) canv->SetLogy();
	canv->Update();
	canv->Print( (output+".pdf").c_str() ,"Portrait pdf");
	canv->Print( (output+".png").c_str()) ;
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

void HbbLimits::Write(const std::vector<Limit>& limits, const std::string& name){
	/*
	 * Method to write limits to a .txt file
	 */
	std::fstream fs;
	fs.open( (name + ".txt").c_str(), std::fstream::out );
	fs<<"X     -2s     -1s  median     +1s     +2s     obs \n";
	for(const auto& limit: limits){
		char strOut[400];
		sprintf(strOut,"%4f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %6.2f",
				(limit.getX()),limit.getMinus2G(),limit.getMinus1G(),limit.getExpected(),limit.getPlus1G(),limit.getPlus2G(),limit.getObserved());
		fs<<strOut;
		fs<<"\n";
	}
	fs.close();
	std::cout<<"File: "<<name + ".txt"<<" has been written."<<std::endl;
}

double HbbLimits::MSSMTanBeta(const std::string& benchmark_path, double mA, double xsection, const std::string& uncert, const bool& UP, const std::string& benchmark_ref_path, const double& tanBref){
	/*
	 * Method to translate GxBR value to tanBeta for MSSM
	 */
	mssm_xs_tools my(benchmark_path.c_str(),true,0);
	//for combination - check reference path:
	bool combination = false;
	mssm_xs_tools ref(benchmark_ref_path.c_str(),true,0);
	if(benchmark_ref_path != ""){
		combination = true;
	}
	double minimalDifference = 1e+10;
	bool rangeExceeded = true;
	double tanBetaTarget = -1;
	double xsecTarget    = -1;
	double maxTanBeta = 60;
	double minTanBeta = 1;

	double sigmaBBA = -100,sigmaBBH = -100;
	double tanBeta, BrAbb, BrHbb, totXSec, difference;
	int ibmax = 10*int(maxTanBeta-minTanBeta);
	for (int ib=0; ib<ibmax; ++ib) {

		tanBeta = minTanBeta + 0.1*double(ib);

		//Check systyematics:
		if(uncert == ""){
		    sigmaBBA = my.bbHSantander_A(mA,tanBeta);
		    sigmaBBH = my.bbHSantander_H(mA,tanBeta);
		}
		else if (uncert == "pdfas"){
		    sigmaBBA = my.bbHSantander_A_pdfas(mA,tanBeta,UP);
		    sigmaBBH = my.bbHSantander_H_pdfas(mA,tanBeta,UP);
		}
		else if (uncert == "scale"){
		    sigmaBBA = my.bbHSantander_A_scale(mA,tanBeta,UP);
		    sigmaBBH = my.bbHSantander_H_scale(mA,tanBeta,UP);
		}

	    BrAbb = my.br_Abb(mA,tanBeta);
	    BrHbb = my.br_Hbb(mA,tanBeta);

	    totXSec = sigmaBBA*BrAbb + sigmaBBH*BrHbb;
	    // for the combination of 7 + 8 + 13 TeV crosssection should be divided by the xsection at reference tanB
	    if(combination){
	    	totXSec /= (ref.bbH5F_A(mA, tanBref) * ref.br_Abb(mA, tanBref) + ref.bbH5F_H(mA, tanBref) *  ref.br_Hbb(mA, tanBref));
	    }
	    difference = TMath::Abs(totXSec-xsection);

	    if (difference<minimalDifference) {
	      minimalDifference = difference;
	      tanBetaTarget = tanBeta;
	      xsecTarget = totXSec;
	    }

	    if (totXSec>xsection) {
	      rangeExceeded = false;
	      break;
	    }

	}

	if (rangeExceeded)
		tanBetaTarget = tanBetaTarget*TMath::Sqrt(xsection/xsecTarget);

	return tanBetaTarget;
}

double HbbLimits::THDMTanBeta(const TH2& h_in, double mA, double xsection){
	/*
	 * Method to translate GxBR limits to tanBeta limits
	 */
	//Clone histo
	TH2D GxBR_2hdm = *(TH2D*) h_in.Clone();
	double minimalDifference = 1e+10;
	bool rangeExceeded = true;
	double tanBetaTarget = -1;
	double xsecTarget    = -1;
	double maxTanBeta = 100;
	double minTanBeta = 1;
	int ibmax = 10*int(maxTanBeta-minTanBeta);
	for (int ib=0; ib<ibmax; ++ib) {

		double tanBeta = minTanBeta + 0.1*double(ib);

//		std::cout<<""
	    double totXSec = GxBR_2hdm.Interpolate(mA,tanBeta);
	    double difference = TMath::Abs(totXSec-xsection);

	    if (difference<minimalDifference) {
	      minimalDifference = difference;
	      tanBetaTarget = tanBeta;
	      xsecTarget = totXSec;
	    }

	    if (totXSec>xsection) {
	      rangeExceeded = false;
	      break;
	    }

	}

	if (rangeExceeded)
		tanBetaTarget = tanBetaTarget*TMath::Sqrt(xsection/xsecTarget);

	return tanBetaTarget;
}

void HbbLimits::CheckHiggsBoson(){
	if(boson_ != "A" && boson_ != "H" && boson_ != "both"){
		std::invalid_argument("ERROR: in HbbLimits::CheckHiggsBoson() - wrong name of Higgs boson was specified. Available configs: A/H/both");
	}
}

void HbbLimits::SetTHDMHistos(TFile& f,std::map<std::string,TH3D*>& histos){
	/*
	 * Method to set map of 2HDM histograms according to boson name
	 */
	if(boson_ != "both"){
		CheckZombieObjectInTFile(f,"xs_bb" + boson_);
		histos["xs_bb" + boson_] = (TH3D*) f.Get(("xs_bb" + boson_).c_str());

		CheckZombieObjectInTFile(f,"br_" + boson_ + "bb");
		histos["br_" + boson_ + "bb"] = (TH3D*) f.Get(("br_" + boson_ + "bb").c_str());
	}
	else{
		CheckZombieObjectInTFile(f,"xs_bbA");
		histos["xs_bbA"] = (TH3D*) f.Get("xs_bbA");
		CheckZombieObjectInTFile(f,"br_Abb");
		histos["br_Abb"] = (TH3D*) f.Get("br_Abb");

		CheckZombieObjectInTFile(f,"xs_bbH");
		histos["xs_bbH"] = (TH3D*) f.Get("xs_bbH");
		CheckZombieObjectInTFile(f,"br_Hbb");
		histos["br_Hbb"] = (TH3D*) f.Get("br_Hbb");
	}
}

bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

TGraph HbbLimits::GetAtlasZhll_flipped(){

	class point {
	public:
		point(const double& v1, const double& v2) : x(v1), y(v2) {};
		double x;
		double y;
	};

	std::vector<point> lower_left= {
				point(-0.02,0.5),
				point(-0.021,0.6),
				point(-0.022,0.7),
				point(-0.023,0.8),
				point(-0.024,0.9),
				point(-0.025,1.),
				point(-0.026,1.3),
				point(-0.05,2.),
				point(-0.054,2.5),
				point(-0.1,3.),
				point(-0.15,3.5),
				point(-0.21,4),
				point(-0.3,4.5),
				point(-0.45,5),
				point(-0.64,5.5),
				point(-0.95,6)
		};

		std::vector<double> cos_lower_left, tan_lower_left;
		for(const auto& val : lower_left) {
			cos_lower_left.push_back(val.x);
			tan_lower_left.push_back(val.y);
		}

		std::vector<point> upper_left = {
				point(-0.95,8.5),
				point(-0.75,9.),
				point(-0.7,9.5),
				point(-0.675,10.),
				point(-0.65,11),
				point(-0.6,12),
				point(-0.55,12.5),
				point(-0.5,13.5),
				point(-0.45,14),
				point(-0.4,18),
				point(-0.375,25.),
				point(-0.36,40.),
				point(-0.35,50.)
		};

		std::vector<double> cos_upper_left, tan_upper_left;
		for(const auto& val : upper_left) {
			cos_upper_left.push_back(val.x);
			tan_upper_left.push_back(val.y);
		}

		std::vector<point> upper_right = {
				point(0.35,50),
				point(0.36,40.),
				point(0.375,25.),
				point(0.4,20),
				point(0.45,14),
				point(0.5,13.5),
				point(0.55,12.5),
				point(0.6,12),
				point(0.65,11),
				point(0.675,10.),
				point(0.7,9.7),
				point(0.75,9.),
				point(0.95,8.5),
		};

		std::vector<double> cos_upper_right, tan_upper_right;
		for(const auto& val : upper_right) {
			cos_upper_right.push_back(val.x);
			tan_upper_right.push_back(val.y);
		}

		std::vector<point> lower_right= {
				point(0.95,6),
				point(0.725,5.5),
				point(0.575,5),
				point(0.5,4.5),
				point(0.425,4),
				point(0.4,3.7),
				point(0.35,3.5),
				point(0.3,3.5),
				point(0.275,3.5),
				point(0.25,3.4),
				point(0.15,3.2),
				point(0.12,3.),
				point(0.1,2.7),
				point(0.07,2.5),
				point(0.05,2.),
				point(0.026,1.3),
				point(0.025,1.),
				point(0.024,0.9),
				point(0.023,0.8),
				point(0.022,0.7),
				point(0.021,0.6),
				point(0.02,0.5),
		};

		std::vector<double> cos_lower_right, tan_lower_right;
		for(const auto& val : lower_right) {
			cos_lower_right.push_back(val.x);
			tan_lower_right.push_back(val.y);
		}

//		TGraph *fl_lower_left = new TGraph(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());
//		TGraph *fl_upper_left = new TGraph(cos_upper_left.size(),cos_upper_left.data(),tan_upper_left.data());
//		TGraph *fl_upper_right = new TGraph(cos_upper_right.size(),cos_upper_right.data(),tan_upper_right.data());
//		TGraph *fl_lower_right = new TGraph(cos_lower_right.size(),cos_lower_right.data(),tan_lower_right.data());

		//One graph for all of them:
		int total_size = cos_upper_left.size() + cos_lower_left.size() + cos_upper_right.size() + cos_lower_right.size();
		std::vector<double> cos_all, tan_all;
		cos_all.reserve(total_size);
		cos_all.insert(cos_all.end(),cos_lower_left.begin(),cos_lower_left.end());
		cos_all.insert(cos_all.end(),cos_upper_left.begin(),cos_upper_left.end());
		cos_all.insert(cos_all.end(),cos_upper_right.begin(),cos_upper_right.end());
		cos_all.insert(cos_all.end(),cos_lower_right.begin(),cos_lower_right.end());

		tan_all.reserve(total_size);
		tan_all.insert(tan_all.end(),tan_lower_left.begin(),tan_lower_left.end());
		tan_all.insert(tan_all.end(),tan_upper_left.begin(),tan_upper_left.end());
		tan_all.insert(tan_all.end(),tan_upper_right.begin(),tan_upper_right.end());
		tan_all.insert(tan_all.end(),tan_lower_right.begin(),tan_lower_right.end());

		TGraph all(total_size,cos_all.data(),tan_all.data());

		return all;

}

} /* namespace MssmHbb */
} /* namespace Analysis */
