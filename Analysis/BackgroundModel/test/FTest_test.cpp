/*
 * FTest_test.cpp
 *
 *  Created on: 12 Apr 2017
 *      Author: shevchen
 *
 *      Unit test of FTest class
 */

#include "Analysis/BackgroundModel/interface/FTest.h"

#include "TFile.h"
#include "TTree.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

using namespace analysis::backgroundmodel;
using namespace std;

struct sub_range{
	sub_range(const double& xmin, const double& xmax, const int& nbins, const string& pdf) :
		xmin(xmin), xmax(xmax), nbins(nbins), pdf(pdf) {};
	double xmin; double xmax;
	int nbins;
	std::string pdf;
};

int main(){
	//Input data file
	TFile f("data/2016DataRereco_05_01_2017/TripleBTagReverseSelectionBtoH2016_prescale_13TeV_G4.root","READ");
	TTree& t = *GetFromTFile<TTree>(f,"MssmHbb_13TeV");
	//PDF family to be studied
	string pdf = "exppolynom";//"supernovosibirsk";//"supernovoeffprod";//"superdijetlinearprod";//"supernovoeffprod";//"superdijeteffprod";//"berneffprod";//
	//Sub-range
//	sub_range sr(200,650,45,pdf);
	sub_range sr(350,1190,42,pdf);
//	sub_range sr(500,1700,48,pdf);
	RooRealVar x("mbb","mbb",sr.xmin,sr.xmax,"GeV");
	x.setBins(sr.nbins);
	RooDataSet data("data","data",RooArgSet(x),RooFit::Import(t));
	data.Print();


	FTest test(data,x,pdf);
	string dir_path = "test/FTest/" + test.DefineSubrange() + "/" + pdf;
	CheckOutputDir(dir_path);
	test.PerformFTest(9,dir_path);

	return 0;
}
