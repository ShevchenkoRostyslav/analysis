#include "TEfficiency.h"
#include "TPaveStats.h"
#include "TVirtualFitter.h"
#include "TBinomialEfficiencyFitter.h"
#include "TGraphAsymmErrors.h"
#include "TClass.h"
#include "TKey.h"
#include "TApplication.h"
#include "TROOT.h"

#include "TMath.h"

#include "TH2.h"
#include "TH2F.h"
#include "TF1.h"

#include "Analysis/MssmHbb/interface/HbbStyleClass.h"
#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

double sigmoid_sigmoid(double *x, double *par);
double func_erf(double *x, double *par);
double test_fcn(double *x, double  *par);
double SF_function(double *x, double *par);

HbbStyle style;
using namespace std;
using namespace mssmhbb;

class MyTEfficiency : public TEfficiency{
public:

	MyTEfficiency(TH2 & num, TH2 & den, const int& bin_number, const bool& bayesian_errs = true);
	void SetNiceStyle(const bool& data);
	void SetFitFunction(TF1 &);
	TF1 *  GetFitFunction() {return f_;};

private:
	void SetBayesianErrors();
	TF1 *f_ = nullptr;

};

MyTEfficiency::MyTEfficiency(TH2 & num, TH2 & den, const int& bin_number, const bool& bayesian_errs) :
	TEfficiency(*num.ProjectionX( ("Slice_num_" + to_string(bin_number)).c_str(),bin_number + 1, bin_number + 1,"E"),
			*den.ProjectionX( ("Slice_den_" + to_string(bin_number)).c_str(),bin_number + 1, bin_number + 1,"E"))
		{
	this->SetName( (string(this->GetName()) + "_" + to_string(bin_number)).c_str());
	if(bayesian_errs) SetBayesianErrors();
}

void MyTEfficiency::SetBayesianErrors(){
	this->SetStatisticOption(TEfficiency::kBUniform);
	this->SetPosteriorMode(1);
	this->UsesBayesianStat();
	this->SetUseWeightedEvents();
}

void MyTEfficiency::SetNiceStyle(const bool& data){
	if(data){
		this->SetMarkerStyle(gStyle->GetMarkerStyle());
		this->SetMarkerSize(gStyle->GetMarkerSize());
	}
	else{
		this->SetMarkerStyle(20);
		this->SetMarkerColor(kMagenta);
	}
}

void MyTEfficiency::SetFitFunction(TF1 &f){
//	f_ = static_cast<TF1*>(f.Clone(f.GetName()));
	f_=&f;
}

void CalculateAndPlotEfficienciesAndSFs(TH2D & data_num, TH2D& data_den, TH2D & mc_num, TH2D & mc_den,const std::vector<std::string> var2_bins, const std::string& var_name = "|#eta^{(1)}|", const string& output_path = "");
pair<MyTEfficiency,MyTEfficiency> CalculateAndPlotEfficiency(TH2D & data_num, TH2D& data_den, TH2D & mc_num, TH2D & mc_den, const int& index, const std::string& var2_bins, const std::string& var_name, const string& output_path);
TF1 * CalculateAndPlotSF(MyTEfficiency *data, MyTEfficiency *mc, const int& index, const string& var_name, const string&var2_bin, const string& output_path);
void DrawTriggerEfficiency(MyTEfficiency *data, MyTEfficiency *mc,const string& var_name, const string&var2_bin, const string& output_path);
void DrawEfficiencyScaleFactor(TF1 *SF,const string& var_name, const string&var2_bin, const string& output_path);
TF1 * GetFitFunction(TEfficiency &h, const string& name, const bool& data);
TF1 * CalculateEfficiencyScaleFactor(MyTEfficiency *data, MyTEfficiency *mc, const int& index);
void PlotEfficiencyVar2Dependence( vector<MyTEfficiency> & data_v, const string& var_name, const std::vector<std::string>& var2_bins, const string& output_path, const bool& data);


int main(){
	style.set(PRIVATE);
	gStyle->SetOptStat(0000);
	gStyle->SetOptFit(0000);
	// Set stat options
	gStyle->SetStatY(0.5);
	// Set y-position (fraction of pad size)
	gStyle->SetStatX(0.8);
	// Set x-position (fraction of pad size)
	//Data and MC TFiles
	TFile fData( (cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_Data2016_ReReco.root").c_str(),"read");
	TFile fMC( (cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_ReReco_lowM_QCD_Pt.root").c_str(),"read");
	// TH2D with numerator and denumerator
	auto h2D_data_numerator 		= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1eta1");
	h2D_data_numerator.GetXaxis()->SetTitle("p_{T}^{(1)} [GeV]");
	auto h2D_data_denumerator 	= *GetFromTFile<TH2D>(fData	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1eta1");
	auto h2D_mc_numerator 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1eta1");
	auto h2D_mc_denumerator 		= *GetFromTFile<TH2D>(fMC	,"TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1eta1");
	//Output path
	string output_path = cmsswBase + "/src/Analysis/MssmHbb/macros/pictures/TriggerPerformance/";
	// Eta bins stored in the histogram
	vector<string> etas = {"[0;0.6]","[0.6;1.7]","[1.7;2.2]"};
	// Calculation and plotting
	CalculateAndPlotEfficienciesAndSFs(h2D_data_numerator,h2D_data_denumerator,h2D_mc_numerator,h2D_mc_denumerator,etas,"|#eta^{(1)}|",output_path);
}

void CalculateAndPlotEfficienciesAndSFs(TH2D & data_num, TH2D& data_den, TH2D & mc_num, TH2D & mc_den,const std::vector<std::string> var2_bins, const std::string& var_name, const string& output_path){
	/*
	 * Regions of the 2 variable
	 */
	//Check output path
	vector<MyTEfficiency> data, mc;
	vector< TF1 > SFs;
	CheckOutputDir(output_path);
	for(const auto& var2 : var2_bins){
		auto index = &var2 - &var2_bins[0];
		auto data_mc_pair = CalculateAndPlotEfficiency(data_num,data_den,mc_num,mc_den,index,var2,var_name,output_path);
		auto SF = CalculateAndPlotSF(&data_mc_pair.first,&data_mc_pair.second,index,var_name,var2,output_path);

		data.push_back(data_mc_pair.first);
		mc.push_back(data_mc_pair.second);
		SFs.push_back(*SF);
	}
	PlotEfficiencyVar2Dependence(data,var_name,var2_bins,output_path,true);
	PlotEfficiencyVar2Dependence(mc,var_name,var2_bins,output_path,false);
}

void PlotEfficiencyVar2Dependence( vector<MyTEfficiency> & data_v, const string& var_name, const std::vector<std::string>& var2_bins, const string& output_path, const bool& data){
	/*
	 * Plot dependence of the trigger efficiency on Var2
	 */
	TCanvas can("canv", "histograms", 600, 600);

	//Setup the frame
	TH2F frame("frame","",1,50.,500.,1,0.,1.2);
	style.setFrameStyle(&frame);
	frame.SetYTitle("#epsilon");
	frame.Draw();

	//Setupo TLegend
	TLegend &leg = *HbbStyle::legend("top,right",data_v.size() + 2);
	leg.Clear();
	leg.SetY1(leg.GetY1() * 0.75);
	leg.SetY2(leg.GetY2() * 0.8);
	style.setLegendStyle(&leg);

	for(int i = 0; i < int(data_v.size()); ++i){
		auto data = &data_v.at(i);
		auto var2_bin	  = var2_bins.at(i);
		frame.SetXTitle(data->TEfficiency::GetCopyPassedHisto()->GetXaxis()->GetTitle());

		auto f = data->GetFitFunction();
		leg.AddEntry(f,(var_name + " #in  " + var2_bin).c_str(),"l");
		f->SetLineStyle(i+1);
		f-> Draw("same");
	}
	leg.Draw("same");
	style.drawStandardTitle("out");
	string out_name = output_path + "TriggerEfficiency_" + var_name + "_dependence";
	if(data) out_name += "_data";
	else out_name += "_mc";

	//Save plots
	can.Print( (out_name + ".pdf").c_str());
	can.Print( (out_name + ".png").c_str());
}

TF1 * CalculateAndPlotSF(MyTEfficiency *data, MyTEfficiency *mc, const int& index, const string& var_name, const string&var2_bin, const string& output_path){
	/*
	 * Calcualte and plot SFs
	 */
	auto SF = CalculateEfficiencyScaleFactor(data,mc,index);
	DrawEfficiencyScaleFactor(SF,var_name,var2_bin,output_path);

	return SF;
}

std::pair<MyTEfficiency,MyTEfficiency> CalculateAndPlotEfficiency(TH2D & data_num, TH2D& data_den, TH2D & mc_num, TH2D & mc_den, const int& index, const std::string& var2_bin, const std::string& var_name, const string& output_path){
	/*
	 * AllInOne method to calcualte and draw Kinematic Trigger Efficiencies
	 */
	MyTEfficiency data(data_num,data_den,index,true);
	data.SetNiceStyle(true);
	MyTEfficiency mc(mc_num,mc_den,index,true);
	mc.SetNiceStyle(false);

	//Trigger Efficiencies
	DrawTriggerEfficiency(&data,&mc,var_name,var2_bin,output_path);

	return {data,mc};
}

TF1 * CalculateEfficiencyScaleFactor(MyTEfficiency *data, MyTEfficiency *mc, const int& index){
	/*
	 * Function to define the moethod of SFs calculation
	 */
	TF1 *sf = nullptr;
	auto fData = data->GetFitFunction();
	if(mc) {
		auto fMC = mc->GetFitFunction();
		sf = new TF1( ("SF" + to_string(index)).c_str(), SF_function,fData->GetXaxis()->GetXmin(),fData->GetXaxis()->GetXmax(),fMC->GetNpar() + fData->GetNpar());
		sf->SetLineColor(gStyle->GetLineColor());
		sf->SetLineWidth(fData->GetLineWidth());
		string xaxis_title = data->TEfficiency::GetCopyPassedHisto()->GetXaxis()->GetTitle();
		sf->SetTitle( (";" + xaxis_title + ";").c_str());
		for(int i = 0; i< fData->GetNpar(); ++i) sf->SetParameter(i, fData->GetParameter(i));
		for(int i = 0; i< fMC->GetNpar(); ++i) sf->SetParameter(i+fData->GetNpar(), fMC->GetParameter(i));
	}
	return sf;
}

void DrawEfficiencyScaleFactor(TF1 *SF,const string& var_name, const string&var2_bin, const string& output_path){
	/*
	 * Draw Trigger Efficiency Data/MC SFs
	 */
	TCanvas can("canv", "histograms", 600, 600);

	//Setup the frame
	TH2F frame("frame","",1,50.,500.,1,0.5,1.5);
	style.setFrameStyle(&frame);
	frame.SetXTitle(SF->GetXaxis()->GetTitle());
	frame.SetYTitle("Scale Factor");
	frame.Draw();

	//Setupo TLegend
	TLegend &leg = *HbbStyle::legend("top,right",2);
	leg.Clear();
	leg.SetHeader((var_name + " #in  " + var2_bin).c_str());
//	leg.SetY1(leg.GetY1() * 0.75);
//	leg.SetY2(leg.GetY2() * 0.8);
	style.setLegendStyle(&leg);
	leg.Draw();

	SF->Draw("same");
	style.drawStandardTitle("out");

	//Save plots
	string out_name = output_path + "TriggerSF_" + var_name + "_" + var2_bin;
	can.Print( (out_name + ".pdf").c_str());
	can.Print( (out_name + ".png").c_str());
}

void DrawTriggerEfficiency(MyTEfficiency *data, MyTEfficiency *mc,const string& var_name, const string&var2_bin, const string& output_path){
	/*
	 * Draw TriggerEfficiency plots
	 *
	 * XAxis is Titled equaly to the Data histo;
	 */
	TCanvas can("canv", "histograms", 600, 600);

	//Setup the frame
	TH2F frame("frame","",1,50.,500.,1,0.,1.2);
	style.setFrameStyle(&frame);
	frame.SetXTitle(data->TEfficiency::GetCopyPassedHisto()->GetXaxis()->GetTitle());
	frame.SetYTitle("#epsilon");
	frame.Draw();

	//Setupo TLegend
	TLegend &leg = *HbbStyle::legend("top,right",4);
	leg.Clear();
	leg.SetHeader((var_name + " #in  " + var2_bin).c_str());
	leg.SetY1(leg.GetY1() * 0.75);
	leg.SetY2(leg.GetY2() * 0.8);
	style.setLegendStyle(&leg);

	data->Draw("Esame");
	leg.AddEntry(data,"Data, JetHT","pl");
	//Fits to the Trigger Efficiency
	string fit_name = "fit";
	auto fitData = GetFitFunction(*data,("data_" + fit_name),true);
	data->SetFitFunction(*fitData);
	TF1 *fitMC = nullptr;
	fitData->Draw("same");
	if(mc){
		leg.AddEntry(mc,"QCD, #hat{p_{T}}","pl");
		mc->Draw("E same");
		fitMC = GetFitFunction(*mc,("mc_" + fit_name),false);
		mc->SetFitFunction(*fitMC);
		fitMC->Draw("same");
	}
	leg.Draw("same");
	style.drawStandardTitle("out");
	string out_name = output_path + "TriggerEfficiency_" + var_name + "_" + var2_bin;

	//Save plots
	can.Print( (out_name + ".pdf").c_str());
	can.Print( (out_name + ".png").c_str());
}

TF1 * GetFitFunction(TEfficiency &h, const string& name, const bool& data){
	/*
	 * Function creates and return TF1 object
	 *
	 * properties are taken from the histogram
	 */
	TF1 *f = nullptr;
	if(data){
		f = new TF1(name.c_str(),sigmoid_sigmoid,100,500,4);
		f->SetParameters(0.99,100.,0.99,100.);
	}
	else {
		f = new TF1(name.c_str(),test_fcn,100,500,7);
		f->SetParameters(4.01543e+01,-108.4,0.1182,108,0.99,-66.06,84.45);
	}
	f->SetLineColor(h.GetMarkerColor());
	h.Fit(f,"R EX0");

	return f;
}

double sigmoid_sigmoid(double *x, double *par){
	double sigmoid = 1./(1.+exp( -par[0] * (x[0] - par[1]) ));
	double sigmoid2 = 1./ (1.+exp( -par[2] * (x[0] - par[3]) ));
	return sigmoid * sigmoid2;
}

double func_erf(double *x, double *par){
	return 1 / 2. * (1 + par[0]*TMath::Erf((x[0]-par[1])/(sqrt(2)*par[2])));
}

double test_fcn(double *x, double  *par){
	return sigmoid_sigmoid(x,par)  * func_erf(x,par+4);
}

double SF_function(double *x,double *par){
	return sigmoid_sigmoid(x,par)/test_fcn(x,par+4);
}
