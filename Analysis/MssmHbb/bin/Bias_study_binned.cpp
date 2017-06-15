#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <memory>
#include <ostream>
#include <iostream>
#include <string>
#include <exception>

//Root includes
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TH1.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TIterator.h"

//RooFit includes
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooDataHist.h"
#include "RooHist.h"
#include "RooNumIntConfig.h"
#include "RooMinimizer.h"
#include "RooMinuit.h"
#include "Math/Math.h"
#include "Math/MinimizerOptions.h"
//include boost
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/Tools/interface/RooFitUtils.h"
#include "Analysis/BackgroundModel/interface/RooFitQuality.h"
#include "Analysis/MssmHbb/src/namespace_mssmhbb.cpp"

namespace fs = boost::filesystem;
using namespace boost::program_options;
using namespace analysis;
using namespace analysis::backgroundmodel;

using namespace std;
using namespace RooFit;

void Bias_study(const int & mass_point, const int & n_obs_, int & signal_strength, const int & toy_max, const std::string & pdf, const std::string & generator, const std::string & output, const bool& bg_only = false, const bool & test = false,const bool & store_plots = true, const bool& fix_turnon = false);
void make_Plot(RooPlot &xframe1, RooPlot &xframe2, const int& i, const std::string & output);

const boost::program_options::variables_map ParseUserInput(const int& argc, const char* argv[]);
const std::string SetupOutputDir(const variables_map& user_input);

double AssignNSignalGen(const int& mass);
double GetEffectiveNSignal(const double& n_gen, const int& signal_strength);

double GetGausRandomisation(TRandom3& r, const double& mean, const double&sigma, const double& max, const double& min);
void RandomiseVarsInWorkspace(RooWorkspace &w, const std::vector<RooRealVar>& init_vars, const int& nsigma);
void FixTurnOn(RooWorkspace &w, const std::string& pdf_type);
void UnFixTurnOn(RooWorkspace &w, const std::string& pdf_type);
class Printer{

public:
	Printer(const bool& test) : test_(test){}; //constructor with value of "test"
	void PrintInputMessage(const double& mass,
						const double& mu,
						const std::string& pdf,
						const std::string& generator){
		/*
		 * Information about the user input
		 */
		if(test_){
			std::cout<<" Starting Bias test with M = "<<mass<<" and signal strength = "<<mu<<std::endl;
			std::cout<<" Selected PDF is: "<<pdf<<std::endl;
			std::cout<<" Selected Generator is: "<<generator<<std::endl;
		}
	};

	void PrintPdfInformation(const std::string& pdf, const int& npars, const std::string& outfile_path) const{
		/*
		 * Information about selected PDF
		 */
		if(test_){
			std::cout<<" PDF: "<<pdf<<" with npar = "<<npars<<" is from the file: "<<outfile_path<<std::endl;
		}
	};

	void PrintNparsInfo(const int& npars_pdf, const int& npars_generator) const{
		/*
		 * Information about the NDFs of generator and pdf
		 */
		std::cout<<"NDFs for PDF      : "<<npars_pdf<<std::endl;
		std::cout<<"NDFs for Generator: "<<npars_generator<<std::endl;
	};

	void PrintRooFitResults(RooFitResult & r) const{
		/*
		 * Information about pre- and post-fit pars
		 */
		std::cout << "\nconstant parameters:" << std::endl;
  		r.constPars().Print("v");
  		std::cout << "\nfloating parameters (init):" << std::endl;
  		r.floatParsInit().Print("v");
  		std::cout << "\nfloating parameters (final):" << std::endl;
  		r.floatParsFinal().Print("v");
	}


private:
	bool test_;
};


std::string GetPdfPath(const std::string & pdf, const Printer& p, int & npars);
std::string GetPdfPath(const std::string & pdf, const Printer& p);
std::vector<RooRealVar> GetRooRealVars(RooWorkspace& w, const std::vector<std::string>& names);
RooRealVar* AdjustMbbVariable(RooWorkspace& w, const double& mass, const std::string& name = "mbb");
void AddOutBranches(TTree &outtree,
		int&    minuit_status,
		double& edm,
		double& minNLL,
		double& nbg_prefit,
		double& nbg_postfit,
		double& nbg_err_postfit,
		double& ns_prefit,
		double& ns_postfit,
		double& ns_err_postfit,
		double& ns_gen,
		int&    r,
		double& chi_2,
		double& bias,
		double& dN_post_pre,
		vector<string>& sg_par_names,
		double *sg_pars,
		double *esg_pars,
		vector<string>& bg_par_names,
		double *bg_pars,
		double *ebg_pars,
		double& ntot
);
void SetPoissonBins(TH1 & h);
vector<string> GetRooRealVarNames(RooWorkspace & w);
void FillArrayOfPdfPars(const vector<string>& names, RooWorkspace & w, double *pars, double *epars);
void SetVarsInWorkspace(RooWorkspace &w, const std::vector<RooRealVar>& vars);

void ChangeParameterRange(RooWorkspace &w, const std::string& pdf_name);

int main(int argc, const char * argv[])
{

	//Parse User input
	variables_map user_input = ParseUserInput(argc,argv);
	//Setup output file
	string output_file = SetupOutputDir(user_input);

	int mass_ 		= user_input["mass_point"].as<int>();
	int n_obs_		= user_input["n_obs"].as<int>();
	int signal_strength_ 	= user_input["signal_strength"].as<int>();
	int toy_max_		= user_input["toy_max"].as<int>();
	string pdf_ 	= user_input["pdf"].as<std::string>();
	string generator_  = user_input["generator"].as<std::string>();
	string signal_	= std::to_string(mass_);
	string output_	= user_input["output_file"].as<std::string>();
	bool bg_only	= user_input["bg_only"].as<bool>();
	bool test_		= user_input["test"].as<bool>();
	bool plots_     = user_input["plots"].as<bool>();
	auto fix_turnon_= user_input["fix_turnon"].as<bool>();

	Bias_study(mass_,n_obs_,signal_strength_,toy_max_,pdf_,generator_,output_file,bg_only,test_,plots_,fix_turnon_);

	return  0;
}

void Bias_study(const int & mass,
				const int & n_obs,
				int & signal_strength,
				const int & toy_max,
				const std::string & pdf,
				const std::string & generator,
				const std::string & output,
				const bool& bg_only,
				const bool & test,
				const bool & store_plots,
				const bool & fix_turnon){
	/*
	 * Function to perform a bias study.
	 */
//	RooAbsPdf::defaultIntegratorConfig()->setEpsRel(1e-10) ;
//	RooAbsPdf::defaultIntegratorConfig()->setEpsAbs(1e-10) ;
//	ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-05);
//	ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);

	std::cout << std::fixed;
//	TH1::SetDefaultSumw2();
	Printer p(test);
	//Input message
	p.PrintInputMessage(mass,signal_strength,pdf,generator);
	//Get number of generated signal
	double n_signal_generated = AssignNSignalGen(mass);
	//Translate it to the total number of signal according to specified signal strength
	double ns_prefit = GetEffectiveNSignal(n_signal_generated,signal_strength);
	string signal_point		= to_string(mass);
	int npars_pdf = 0, npars_generator = 0;
	//Get path to the workspaces:
	string p_to_pdf 		= GetPdfPath(pdf,p,npars_pdf);
	string p_to_generator 	= GetPdfPath(generator,p,npars_generator);
	string p_to_signal 		= GetPdfPath(signal_point,p);
	p.PrintNparsInfo(npars_pdf,npars_generator);
	//Get workspaces from the root files:
	RooWorkspace &wPDF = *GetRooWorkspace(p_to_pdf);
	RooWorkspace &wGenerator = *GetRooWorkspace(p_to_generator);
	RooWorkspace &wSignal	= *GetRooWorkspace(p_to_signal);
	//Adjust x-axis (mbb)
	RooRealVar &mbb = *AdjustMbbVariable(wSignal,mass,"mbb");
	//Output TFile and TTree
	TFile outputFile((output+ "/" + output + ".root").c_str(), "RECREATE");
	TTree outtree("Toys","Toys");

	double nbg_prefit = n_obs, nbg_postfit = 0 ,ns_postfit=0.,ns_err_postfit=0.,nbg_err_postfit=0.,chi_2=0.,bias=0., dN_post_pre = 0, ntot = 0;
	double edm = -100, minNLL = 100;
	int minuit_status = -100;
	double signal_pars[10], esignal_pars[10], bg_pars[10], ebg_pars[10];
	vector<string> bg_par_names = GetRooRealVarNames(wPDF);
	vector<string> sg_par_names = GetRooRealVarNames(wSignal);
	//Make parameters range wider
	ChangeParameterRange(wPDF,pdf);
	//Get initial values of the parameters - needed to reset value at each step of the Toy production
	vector<RooRealVar> bg_roovars = GetRooRealVars(wPDF,bg_par_names);
	vector<RooRealVar> sg_roovars = GetRooRealVars(wSignal,sg_par_names);
	//Output branches
	AddOutBranches(outtree,minuit_status,edm,minNLL,nbg_prefit,nbg_postfit,nbg_err_postfit,ns_prefit,ns_postfit,ns_err_postfit,n_signal_generated,signal_strength,chi_2,bias,dN_post_pre,sg_par_names,signal_pars,esignal_pars,bg_par_names,bg_pars,ebg_pars,ntot);
	//pdfs to be used
	auto &pdf_bkg_fit 	= *GetFromRooWorkspace<RooAbsPdf>(wPDF,"background");
	auto &pdf_bkg_alt 	= *GetFromRooWorkspace<RooAbsPdf>(wGenerator,"background");
	auto &pdf_sgn 		= *GetFromRooWorkspace<RooAbsPdf>(wSignal,"signal");

	//S+Bg pdf

	//random generator for Poisson
	for(int i =1; i<=toy_max;++i){
		/*
		 * Toy is produced by alternative function
		 * and fitted with backgorund pdf
		 */

		//Bias loop
		cout<<"\n\nToy number "<<i<<"\n\n"<<endl;
		string toy_number = to_string(i);

		//Reset Vars
		SetVarsInWorkspace(wPDF, bg_roovars);
		if(fix_turnon) FixTurnOn(wPDF,pdf);
		SetVarsInWorkspace(wSignal, sg_roovars);

		TH1D &data_toy = *(TH1D*) pdf_bkg_alt.createHistogram( ("data_toy" + toy_number).c_str(), mbb);
		data_toy.Scale(n_obs/data_toy.Integral());
		TH1D &data_toy_sgn = *(TH1D*) pdf_sgn.createHistogram( ("data_toy_sgn" + toy_number).c_str(), mbb);
		data_toy_sgn.Scale(ns_prefit/data_toy_sgn.Integral());


//		SetPoissonBins(data_toy);
		nbg_prefit = data_toy.Integral();
		RooDataHist data_bg_toy( ("data_bg_toy"+toy_number).c_str(),"data_bg_toy", mbb, RooFit::Import(data_toy));
		data_toy.Add(&data_toy_sgn); //Sumup
		SetPoissonBins(data_toy);
		RooDataHist data_toy_sb( ("data_toy_sb_"+toy_number).c_str(),"data_toy_sb", mbb, RooFit::Import(data_toy));

		RooFitResult *fitResult;//, *bgFitResult;

		//Set vars
		RooRealVar roo_nsignal  ("roo_nsignal","roo_nsignal",ns_prefit); roo_nsignal.setConstant(false);
		RooRealVar roo_ntot     ("roo_ntot","roo_ntot",data_toy.Integral()); roo_ntot.setConstant(true);
		RooFormulaVar roo_nbkg ("roo_nbkg","roo_nbkg","roo_ntot - roo_nsignal",RooArgList(roo_ntot,roo_nsignal));

		RooAddPdf  pdf_SpBg             ("pdf_SpBg","pdf_SpBg",RooArgList(pdf_sgn,pdf_bkg_fit),RooArgList(roo_nsignal,roo_nbkg));

		if(!bg_only){
			//try Bg only fit before the S+Bg one:
//			if(fix_turnon) UnFixTurnOn(wPDF,pdf);
//			bgFitResult = pdf_bkg_fit.fitTo(data_bg_toy, RooFit::Save(), RooFit::PrintLevel(0), RooFit::Verbose(0), RooFit::Minimizer("Minuit2", "migrad"), RooFit::Hesse(1), RooFit::Offset(1));
//			p.PrintRooFitResults(*bgFitResult);
//			roo_ntot.Print(); roo_nsignal.Print(); roo_nbkg.Print();
//			if(fix_turnon) FixTurnOn(wPDF,pdf);
//			fitResult = pdf_SpBg.fitTo(data_toy_sb, RooFit::Save(), RooFit::PrintLevel(3), RooFit::Verbose(1), RooFit::Minimizer("Minuit2", "migrad"), RooFit::Hesse(1) , RooFit::Offset(1) );

			/*
			//RooMinimizer fit
			//bg fit
//			if(fix_turnon) UnFixTurnOn(wPDF,pdf);
			RooAbsReal &nll_bg = *pdf_bkg_fit.createNLL(data_bg_toy);
			RooMinimizer m_bg(nll_bg);
			m_bg.setMinimizerType("Minuit2");
			m_bg.setOffsetting(true);
			m_bg.setPrintLevel(3);
//			m_bg.setWarnLevel(3);
			m_bg.setVerbose(kTRUE);
			m_bg.setEps(1e-09);
			m_bg.optimizeConst(kTRUE);
			m_bg.minimize("Minuit2","Migrad") ;
			m_bg.hesse();
			bgFitResult = m_bg.save();
			p.PrintRooFitResults(*bgFitResult);
//			if(fix_turnon) FixTurnOn(wPDF,pdf);
			 */

			//s+bg fit
			RooAbsReal &nll = *pdf_SpBg.createNLL(data_toy_sb);
			RooMinimizer m(nll);
			m.setMinimizerType("Minuit2");
			m.setOffsetting(true);
			m.setPrintLevel(3);
//			m.setWarnLevel(3);
			m.setVerbose(kTRUE);
			m.setEps(1e-09);
			m.optimizeConst(kTRUE);
			m.minimize("Minuit2","Migrad") ;
			m.hesse();
			fitResult = m.save();
//			p.PrintRooFitResults(*fitResult);
		}
		else {
//			fitResult = pdf_bkg_fit.fitTo(data_toy_sb, RooFit::Save(), RooFit::PrintLevel(1), RooFit::Verbose(1), RooFit::Minimizer("Minuit2", "migrad"), RooFit::Hesse(1), RooFit::Offset(1));
			RooAbsReal &nll = *pdf_bkg_fit.createNLL(data_toy_sb);
			RooMinimizer m(nll);
			m.setMinimizerType("Minuit2");
			m.setOffsetting(true);
			m.setPrintLevel(3);
//			m.setWarnLevel(3);
			m.setVerbose(kTRUE);
			m.setEps(1e-09);
			m.optimizeConst(kTRUE);
			m.minimize("Minuit2","Migrad") ;
			m.hesse();
			fitResult = m.save();
		}
	  	p.PrintRooFitResults(*fitResult);
	  	minuit_status = fitResult->status();

		cout<<"AFter: nBg = "<<roo_nbkg.getValV()<<" nsg = "<<roo_nsignal.getValV()<<" n_tot = "<<roo_nbkg.getValV() + roo_nsignal.getValV()<<endl;

		RooPlot &xframe1 = *mbb.frame(Title("Bias Test")) ;
	  	data_toy_sb.plotOn(&xframe1, RooFit::Name("data_curve"), RooFit::MarkerSize(0.1));
	  	if(!bg_only){
	  		pdf_SpBg.plotOn(&xframe1, RooFit::Name("fit_curve"));
	  	}
	  	else pdf_bkg_fit.plotOn(&xframe1, RooFit::Name("fit_curve"));

//	  	chi_2 = xframe1.chiSquare("fit_curve", "data_curve", npars_pdf + 1);	//floating parameters from sig+bkg fit including n_b,n_s
	  	RooFitQuality fit_quality;
	  	Chi2Ndf chi2ndf = fit_quality.chiSquare(xframe1,"fit_curve", "data_curve",npars_pdf + 1);
	  	chi_2 = chi2ndf.chi2 / chi2ndf.ndf;
	  	std::cout<<"Chi2 / ndf = "<<chi2ndf.chi2<<" / "<<chi2ndf.ndf<<" = "<<chi_2<<std::endl;
	  	std::cout << "Chi^2/df from the fit = " << chi_2 << std::endl;

		ns_postfit = roo_nsignal.getValV();
		nbg_postfit = roo_nbkg.getValV();
		//nbg_err_postfit = roo_nbkg.getError();
		ns_err_postfit = roo_nsignal.getError();
		bias   = (ns_postfit - ns_prefit)/ns_err_postfit;
		dN_post_pre = ns_postfit - ns_prefit;
		edm = fitResult->edm();
		minNLL = fitResult->minNll();
		ntot =  roo_ntot.getValV();

		std::cout << "ns_prefit      = " << ns_prefit << std::endl; 	// number of signal before the fit
		std::cout << "nbg_prefit     = " << roo_ntot.getValV() - ns_prefit<< std::endl;
		std::cout <<" ntot_pretfit   = " << roo_ntot.getValV()<<endl;
	  	std::cout << "ns_postfit     = " << ns_postfit << std::endl; 	// number of signal after the fit
	  	std::cout << "nbg_postfit    = " << nbg_postfit<< std::endl;
	  	std::cout <<" ntot_postfit   = "<<roo_nbkg.getValV() + roo_nsignal.getValV()<<endl;
	  	std::cout <<" ntot_fromHisto = "<<data_toy.Integral()<<std::endl;

	  	std::cout << "err_s_postfit  = " << ns_err_postfit << std::endl; // error of n_s from the fit

	  	std::cout << "Bias   = " << bias << std::endl;

	  	if(i <= 5 && store_plots){
	  		RooHist &hpull = *xframe1.pullHist();
	  		hpull.SetMarkerSize(0.1);
	  		RooPlot &xframe2 = *mbb.frame();
	  		xframe2.addPlotable(&hpull,"P");
	  		xframe2.SetMinimum(-5.);
	  		xframe2.SetMaximum(+5.);

	  		pdf_SpBg.plotOn(&xframe1, RooFit::Name("bkg"), Components("background"), LineStyle(kDashed), LineColor(kRed), LineWidth(4));
	  		pdf_SpBg.plotOn(&xframe1, RooFit::Name("sgn"), Components("signal"), LineStyle(kDashed), LineColor(kGreen+2), LineWidth(4));

	  		make_Plot(xframe1,xframe2,i,output);
	  	}

		FillArrayOfPdfPars(sg_par_names,wSignal,signal_pars,esignal_pars);
		FillArrayOfPdfPars(bg_par_names,wPDF,bg_pars,ebg_pars);

		outtree.Fill();

	}


	outputFile.cd();
	outtree.Write();
	outputFile.Close();
	cout << "End of code :)" << endl;

}

void make_Plot(RooPlot &xframe1, RooPlot &xframe2, const int& i, const std::string & output){
	std::string ntoy = std::to_string(i);

	TCanvas canvas;
        canvas.SetCanvasSize(500,500);
        canvas.cd();

	TPad* pad1;
	pad1 = new TPad("pad1","",0,0.1,1,1);
	pad1->SetBottomMargin(0.2);
        pad1->SetRightMargin(0.05);
        pad1->SetLeftMargin(0.16);
        pad1->Draw();
        pad1->cd();
        xframe1.GetXaxis()->SetTitleOffset(999); //Effectively turn off x axis title on main plot
        xframe1.GetXaxis()->SetLabelOffset(999); //Effectively turn off x axis label on main plot
        xframe1.GetYaxis()->SetTitleSize(0.038);
        xframe1.GetYaxis()->SetTitleOffset(1.6);
        xframe1.GetYaxis()->SetLabelSize(0.033);
	xframe1.Draw();
        TLegend *leg1 = new TLegend(0.67,0.63,0.93,0.89);
        leg1->SetFillColor(kWhite);
        leg1->SetLineColor(kWhite);
        leg1->AddEntry("data_curve","Toy Data", "P");
        leg1->AddEntry("fit_curve","Signal + background","L");
        leg1->AddEntry("bkg","Background only", "L");
        leg1->AddEntry("sgn","Signal only", "L");
        leg1->Draw("same");
	
	canvas.cd();
        TPad *pad2 = new TPad("pad2","",0,0.0,1,0.265);
        pad2->SetTopMargin(1);
        pad2->SetBottomMargin(0.33);
        pad2->SetLeftMargin(pad1->GetLeftMargin());
        pad2->SetRightMargin(pad1->GetRightMargin());
        pad2->SetGridy();
        pad2->Draw();
        pad2->cd();
        xframe2.SetTitle("");
        xframe2.GetXaxis()->SetTitleSize(0.15);
        xframe2.GetXaxis()->SetTitleOffset(0.9);
        xframe2.GetXaxis()->SetLabelSize(0.115);
        xframe2.GetXaxis()->SetLabelOffset(0.010);
        xframe2.SetYTitle("Pulls");
        xframe2.GetYaxis()->CenterTitle(kTRUE);
        xframe2.GetYaxis()->SetTitleSize(0.14);
        xframe2.GetYaxis()->SetTitleOffset(0.4);
        xframe2.GetYaxis()->SetNdivisions(3,5,0);
        xframe2.GetYaxis()->SetLabelSize(0.115);
        xframe2.GetYaxis()->SetLabelOffset(0.011);
        xframe2.Draw();

	canvas.SaveAs((output+std::string("/")+std::string("Toy_")+ntoy+std::string("_linear.pdf")).c_str());
	pad1->SetLogy();
        xframe1.GetYaxis()->SetRangeUser(0.1, xframe1.GetMaximum()*100);
        canvas.Modified();
        canvas.Update();
        canvas.SaveAs((output+std::string("/")+std::string("Toy_")+ntoy+std::string("_log.pdf")).c_str());
}

const boost::program_options::variables_map ParseUserInput(const int& argc, const char* argv[]){
	// firts arg is mass point (int ), second - signal strength, third - number of toys
	// list of the pdfs to fit toy S+B MC
	// list of the pdfs to generate toy S+B MC

	options_description generalOptions("General options");
    	generalOptions.add_options()
				("help,h", "Produce help message.")
				("verbose,v", value<int>()->default_value(0), "More verbose output.");

        // Declare a group of options that will be allowed only on command line
		options_description cmdLineOptions("CMD Optional arguments");
		cmdLineOptions.add_options()
			("mass_point,m", value<int>()->required(),"Select mass point from the list of available: \n"
			" \t300;\n"
			" \t350;\n"
			" \t400;\n"
			" \t500;\n"
			" \t600;\n"
			" \t700;\n"
			" \t900;\n"
			" \t1100;\n"
			" \t1300;\n")
			("n_obs,n",value<int>()->required(),"Number of observable in each subrange: \n"
			" \t259399 -- SR1;\n"
			" \t105053 -- SR2;\n"
			" \t26760  -- SR3;\n")
			("signal_strength,r",value<int>()->required(),"Signal strength to be injected to toys fit")
			("toy_max,t",value<int>()->required(),"Number of toys")
			("pdf,f",value<std::string>()->required(),"Select nominal PDF to fit toys data: \n"
			" \textnovoeffprod_G4_R1;\n"
			" \tsupernovoeffprod2_G4_R1;\n"
			" \tsupernovoeffprod3_G4_R1;\n"
			" \tberneffprod9par_G4_R1;\n"
			" \tberneffprod8par_G4_R1;\n"
			" \tdijetv3logprod_G4_R1;\n"
			" \tsuperdijeteffprod2_G4_R1;\n"
			" \tnovoeffprod_G7_R1;\n"
			" \tberneffprod7par_G7_R1;\n"
			" \tnovosibirsk_G4_R2;\n"
			" \tsupernovosibirsk1_G4_R2;\n"
			" \tsupernovosibirsk2_G4_R2;\n"
			" \tdijetexp_G4_R2;\n"
			" \tnovosibirsk_G7_R2;\n"
			" \tdijetexp_G7_R2;\n"
			" \tnovosibirsk_G4_R3;\n"
			" \tsupernovosibirsk1_G4_R3;\n"
			" \tsupernovosibirsk2_G4_R3;\n"
			" \tdijetexp_G4_R3;\n"
		 	" \tnovosibirsk_G7_R3;\n"
			" \tdijetexp_G7_R3;\n")
			("generator,g",value<std::string>()->required(),"Select alternative PDF to produce toys data: \n"
			" \textnovoeffprod_G4_R1;\n"
			" \tsupernovoeffprod2_G4_R1;\n"
			" \tsupernovoeffprod3_G4_R1;\n"
			" \tberneffprod9par_G4_R1;\n"
			" \tberneffprod8par_G4_R1;\n"
			" \tdijetv3logprod_G4_R1;\n"
			" \tsuperdijeteffprod2_G4_R1;\n"
			" \tnovoeffprod_G7_R1;\n"
			" \tberneffprod7par_G7_R1;\n"
			" \tnovosibirsk_G4_R2;\n"
			" \tdijetexp_G4_R2;\n"
			" \tnovosibirsk_G7_R2;\n"
			" \tdijetexp_G7_R2;\n"
			" \tnovosibirsk_G4_R3;\n"
			" \tdijetexp_G4_R3;\n"
			" \tnovosibirsk_G7_R3;\n"
			" \tdijetexp_G7_R3;\n")
			("output_file,o", value<std::string>()->default_value(""), "Output file name, if not specified will be created automatically")
			("plots,p", bool_switch()->default_value(false),"if specified - no plots will be stored")
			("bg_only",bool_switch()->default_value(false),"if specified - background only fit would be performed")//value<bool>()->default_value(0),"if ");
			("fix_turnon", bool_switch()->default_value(false),"if specified - Turn-on parameters will be fixed")
		;


        // Hidden options, will be allowed both on command line and
        // in config file, but will not be shown to the user.
        options_description hidden("Hidden options");
        hidden.add_options()
            ("test", value<bool>()->default_value(1),"if test specify 1. Default 0.")
            ;

        options_description all_options;
        all_options.add(generalOptions).add(cmdLineOptions).add(hidden);

        options_description visible("Allowed options");
        visible.add(generalOptions).add(cmdLineOptions);


   	variables_map output_vm;
    	store(command_line_parser(argc, argv).options(generalOptions).allow_unregistered().run(), output_vm);
    	notify(output_vm);

    	if (output_vm.count("help")) {
	    std::cout << visible << std::endl;
	    exit(0);
	}

    	store(parse_command_line(argc, argv,all_options), output_vm);
    	try {
    		notify(output_vm);
    	}
    	catch (const required_option& ex) {
    		std::cerr << ex.what() << std::endl;
        	exit(-1);
    	}

    	return output_vm;
}

const std::string SetupOutputDir(const variables_map& user_input){
	/*
	 * Function to setup output firectory path.
	 *
	 * Setup is done according to the user input AND
	 * defined template
	 */
	const string user_path 	= user_input["output_file"].as<string>();
	auto fix_turnon 		= user_input["fix_turnon"].as<bool>();
	auto bg_only			= user_input["bg_only"].as<bool>();
	std::string generator_  = user_input["generator"].as<std::string>();
	std::string pdf_ 		= user_input["pdf"].as<std::string>();
	int mass_ 				= user_input["mass_point"].as<int>();
	int signal_strength_ 	= user_input["signal_strength"].as<int>();


	std::string oDir 		= user_path;
		//Pre-defined template is used
	oDir = "BiasTest_toy"+generator_+"_fit"+pdf_+"_"+std::to_string(mass_)+"GeV_r"+std::to_string(signal_strength_);
	if(bg_only) oDir += "_BGonly";
	if(fix_turnon) oDir += "_fixedTurnOn";
	if(user_path != "") oDir = user_path + "_" + oDir;
	CheckOutputDir(oDir);
	std::cout << "Writing results to " << oDir << std::endl;
	return oDir;
}

double AssignNSignalGen(const int& mass){
	/*
	 * Function to assign and get number of generated signal
	 */
	switch(mass){
	case 300:
		return 1187;
	case 350:
		return 1681;
	case 400:
		return 968;
	case 500:
		return 888;
	case 600:
		return 521;
	case 700:
		return 299;
	case 900:
		return 141;
	case 1100:
		return 82;
	case 1300:
		return 59;
	default:
		throw invalid_argument("Invalid mass: " + to_string(mass) + " has been provided.");
		return 0;
	}
}
double GetEffectiveNSignal(const double& n_gen, const int& signal_strength){
	/*
	 * Function to get effective number of signal events.
	 *
	 * N_effective = signal_srength * N_generated
	 */
	return signal_strength*n_gen;
}

RooRealVar* AdjustMbbVariable(RooWorkspace& w, const double& mass, const std::string& mbb_name){
	/*
	 * Function to adjust Mbb x-variable
	 */
	vector<int> sr1 = {300,350,400,500};
	vector<int> sr2 = {600,700,900};
	vector<int> sr3 = {1100,1300};
	int nbins = 500;
	auto *mbb = GetFromRooWorkspace<RooRealVar>(w,mbb_name);
	if( find(sr3.begin(),sr3.end(),mass) != sr3.end()) {
		mbb->setRange(500,1700);
	}
	else if( find(sr2.begin(),sr2.end(),mass) != sr2.end() ){
		mbb->setRange(350,1190);
	}
	else if ( find(sr1.begin(),sr1.end(),mass) != sr1.end() ) {
		mbb->setRange(200,650);
	}
	else {
		throw invalid_argument("AdjustMbbVariable:: No sub-range for mass: " + to_string(mass));
	}
	mbb->setBins(nbins);

	TIterator &s_iter = *w.componentIterator();
	auto *s_var = s_iter.Next();
		while(s_var){
			string class_name = s_var->ClassName();
			string name = s_var->GetName();
			if(class_name == "RooRealVar"){
				auto &var = *GetFromRooWorkspace<RooRealVar>(w,name);
				if(name != "mbb"){
					var.setConstant();
					}
			}
			s_var = s_iter.Next();
		}

	return mbb;
}

std::string GetPdfPath(const std::string & pdf, const Printer& p){
	/*
	 * Function to get path for particular pdf
	 */
	int npars;
	return GetPdfPath(pdf,p,npars);
}

std::string GetPdfPath(const std::string & pdf, const Printer& p, int & npar){
	/*
	 * Function to get path for particular pdf and number of pars
	 */
	std::string outname;
//	int npar = -1;
	if(findStrings(pdf,"extnovoeffprod_G4_R1")){
		npar = 6;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/extnovoeffprod_200to650_10GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if (findStrings(pdf,"extnovoeffprod_R1")){
		npar = 6;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/fullBBnB/extnovoeffprod_200to650_10GeV__45bins/workspace/FitContainer_workspace.root";
	}
	else if (findStrings(pdf,"supernovoeffprod2_G4_R1")){
		npar = 7;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr1/supernovoeffprod/pdf2/workspace/FitContainer.root";
	}
	else if (findStrings(pdf,"supernovoeffprod3_G4_R1")){
		npar = 8;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr1/supernovoeffprod/pdf3/workspace/FitContainer.root";
	}
	else if(findStrings(pdf,"berneffprod9par_G4_R1")) {
		npar = 11;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/Prescale_v5/berneff9_200_to_650_10GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if(findStrings(pdf,"berneffprod8par_G4_R1")){
		npar = 8;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/Prescale_v5/berneff8_200_to_650_10GeV/workspace/FitContainer_workspace.root";
	}
	else if(findStrings(pdf,"berneffprod8par_R1")){
			npar = 8;
			outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/fullBBnB/berneffprod,8_200to650_10GeV__45bins/workspace/FitContainer_workspace.root";
		}
	else if(findStrings(pdf,"dijetv3logprod_G4_R1")) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/dijetv3logprod_200to650_10GeV_G4/workspace/FitContainer_workspace_corr.root";
	}
	else if(findStrings(pdf,"superdijeteffprod2_G4_R1")){
		npar = 7;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/Prescale_v5/superdijeteffprod2_200_to_650_10GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if(findStrings(pdf,"superdijeteffprod2_R1")){
			npar = 7;
			outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/fullBBnB/superdijeteffprod,2_200to650_10GeV__45bins/workspace/FitContainer_workspace.root";
		}
	else if(findStrings(pdf,"berneffprod7par_G7_R1")) {
		npar = 9;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/berneffprod7par_200to650_10GeV_G7/workspace/FitContainer_workspace.root";
        }
	else if(findStrings(pdf,"novosibirsk_G4_R2")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_350to1190_20GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if (findStrings(pdf,"supernovosibirsk1_G4_R2") ){
		npar = 4;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr2/supernovosibirsk/pdf1/workspace/FitContainer.root";
	}
	else if (findStrings(pdf,"supernovosibirsk2_G4_R2") ){
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr2/supernovosibirsk/pdf2/workspace/FitContainer.root";
	}
	else if(findStrings(pdf,"dijetexp_G4_R2")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/dijetv2_350to1190_20GeV_G4/workspace/FitContainer_workspace.root";
        }
	else if(findStrings(pdf,"novosibirsk_G7_R2")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_350to1190_20GeV_G7/workspace/FitContainer_workspace.root";
        }
	else if(findStrings(pdf,"dijetexp_G7_R2")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/dijetv2_350to1190_20GeV_G7/workspace/FitContainer_workspace.root";
        }

	else if(findStrings(pdf,"novosibirsk_G4_R3")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_500to1700_25GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if (findStrings(pdf,"supernovosibirsk1_G4_R3") ){
		npar = 4;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr3/supernovosibirsk/pdf1/workspace/FitContainer.root";
	}
	else if (findStrings(pdf,"supernovosibirsk2_G4_R3") ){
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/FTest/sr3/supernovosibirsk/pdf2/workspace/FitContainer.root";
	}
	else if(findStrings(pdf,"dijetexp_G4_R3")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/dijetv2_500to1700_25GeV_G4/workspace/FitContainer_workspace.root";
	}
	else if(findStrings(pdf,"novosibirsk_G7_R3")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/novosibirsk_500to1700_25GeV_G7/workspace/FitContainer_workspace.root";
        }
	else if(findStrings(pdf,"dijetexp_G7_R3")) {
		npar = 3;
		outname = mssmhbb::cmsswBase + "/src/Analysis/BackgroundModel/test/dijetv2_500to1700_25GeV_G7/workspace/FitContainer_workspace.root";
        }
	else if(pdf == "300" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-300/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "350" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-350/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "400" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-400/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "500" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-500/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "600" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-600/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "700") {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-700/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "900" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-900/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "1100" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1100/workspace/FitContainer_workspace.root";
	}
	else if(pdf == "1300" ) {
		npar = 5;
		outname = mssmhbb::cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-1300/workspace/FitContainer_workspace.root";
	}
	else{
		throw invalid_argument("Invalid pdf: " + pdf + " has been provided.");
	}
	p.PrintPdfInformation(pdf,npar,outname);

	return outname;
}

void AddOutBranches(TTree &outtree,
					int&    minuit_status,
					double& edm,
					double& minNLL,
					double& nbg_prefit,
					double& nbg_postfit,
					double& nbg_err_postfit,
					double& ns_prefit,
					double& ns_postfit,
					double& ns_err_postfit,
					double& ns_gen,
					int&    r,
					double& chi_2,
					double& bias,
					double& dN_post_pre,
					vector<string>& sg_par_names,
					double *sg_pars,
					double *esg_pars,
					vector<string>& bg_par_names,
					double *bg_pars,
					double *ebg_pars,
					double& ntot
					){
	/*
	 * Function to add special branches to the TTree
	 */
	outtree.Branch("minuit_status",     &minuit_status,     "minuit_status/I");
	outtree.Branch("edm",				&edm,				"edm/D");
	outtree.Branch("minNLL",			&minNLL,			"minNLL/D");
	outtree.Branch("nbg_prefit",		&nbg_prefit,		"nbg_prefit/D");
	outtree.Branch("ntot",				&ntot,				"ntot/D");
	outtree.Branch("nbg_postfit", 		&nbg_postfit, 		"nbg_postfit/D");
	outtree.Branch("nbg_err_postfit",	&nbg_err_postfit,	"nbg_err_postfit/D");
	outtree.Branch("signal_strength",	&r,					"signal_strength/I");
	outtree.Branch("ns_prefit",			&ns_prefit,			"ns_prefit/D");
	outtree.Branch("ns_postfit", 		&ns_postfit, 		"ns_postfit/D");
	outtree.Branch("ns_err_postfit",	&ns_err_postfit,	"ns_err_postfit/D");
	outtree.Branch("ns_gen",			&ns_gen,	"ns_gen/D");
	outtree.Branch("chi_2",				&chi_2,		"chi_2/D");
	outtree.Branch("bias",				&bias,		"bias/D");
	outtree.Branch("dN_post_pre",		&dN_post_pre,		"dN_post_pre/D");
	//branches with bg and signal parameters TODO: introduce lambda-expression
	auto lambda_n = [](TTree & tree, auto vec, auto pars) { for(auto & v : vec) {auto index = &v - &vec[0]; tree.Branch(v.c_str(), &pars[index] , (v + "/D").c_str());}};
	auto lambda_e = [](TTree & tree, auto vec, auto pars) { for(auto & v : vec) {auto index = &v - &vec[0]; tree.Branch((v + "_err").c_str(), &pars[index] , (v + "_err/D").c_str());}};
	lambda_n(outtree,sg_par_names,sg_pars);
	lambda_e(outtree,sg_par_names,esg_pars);
	lambda_n(outtree,bg_par_names,bg_pars);
	lambda_e(outtree,bg_par_names,ebg_pars);
}

void SetPoissonBins(TH1 & h){
	/*
	 * Function to set Poisson distr. bin content
	 */
	TRandom3 r(0);
	for(int j = 1; j <= h.GetNbinsX(); ++j)
	{
		int n = r.Poisson(h.GetBinContent(j));
		h.SetBinContent(j,n);
		h.SetBinError(j,TMath::Sqrt(n));
//		h.SetBinError(j,TMath::Sqrt(h.GetBinContent(j)));
	}
}

//void GetRooRealVarMeans(RooWorkspace & w, vector<pair<string, double > >& v, vector<pair<string, double > >& ev){
vector<string> GetRooRealVarNames(RooWorkspace & w){
	/*
	 * Function to get a pair of NAME : DOUBLE & from RooRealVar
	 *
	 * To be written in TTree
	 */
	vector<string> par_names;
	TIterator &s_iter = *w.componentIterator();
	auto *s_var = s_iter.Next();
		while(s_var){
			string class_name = s_var->ClassName();
			string name = s_var->GetName();
			if(class_name == "RooRealVar"){
				if(name != "mbb" && name != "weight"){
					par_names.push_back(name);
				}
			}
			s_var = s_iter.Next();
		}
		return par_names;
}

void FillArrayOfPdfPars(const vector<string>& names, RooWorkspace & w, double *pars, double *epars){
	/*
	 * Function to fill the arrays of pdf parameters
	 */
	for(const auto& v : names){
		auto index = &v - &names[0];
		pars[index]  = w.var(v.c_str())->getValV();
		epars[index] = w.var(v.c_str())->getError();
	}
}

std::vector<RooRealVar> GetRooRealVars(RooWorkspace& w, const std::vector<std::string>& names){
	/*
	 * Function to get copies of RooRealVars from the
	 * Workspace
	 */
	vector<RooRealVar> v;
	for(const auto& n : names){
		RooRealVar var(*w.var(n.c_str()),n.c_str());
		v.push_back(var);
	}
	return v;
}

void SetVarsInWorkspace(RooWorkspace &w, const std::vector<RooRealVar>& vars){
	/*
	 * Function to set vars in the workspace according to vars
	 */
	std::cout<<"\nWorkspace content: "<<w.GetName()<<" VARS"<<std::endl;
	w.Print("v");
	for(const auto& v : vars){
		RooRealVar &init = *w.var(v.GetName());
		init.setVal(v.getValV());
		init.setError(v.getError());
		init.setMax(v.getMax());
		init.setMin(v.getMin());
		init.Print();
	}
}

void ChangeParameterRange(RooWorkspace &w, const std::string& pdf_name){
	/*
	 * Function to change parameter range
	 * hard coded
	 */
	if(pdf_name == "extnovoeffprod_G4_R1"){
		w.var("tail")->setRange(-10, 10);
	}
	else if (pdf_name == "novosibirsk_G4_R2"){
		w.var("tail1")->setRange(-10, 10);
	}
	else if (pdf_name == "berneffprod9par_G4_R1"){
//		for(const string& s : {"0","1","2","3","4","5","6","7","8"}){
//			w.var( ("bernstein_coefficient_0" + s).c_str())->setRange(-20, 20);
//		}
	}
	else if (pdf_name == "berneffprod8par_G4_R1"){
	}
	else if (pdf_name == "dijetv3logprod_G4_R1"){
	}
	else if (pdf_name == "superdijeteffprod2_G4_R1"){
		w.var("parb_dijet")->setRange(-1,1);
	}
}

void FixTurnOn(RooWorkspace &w, const std::string& pdf_name){
	/*
	 * Function to fix turn-on curve
	 */
	if(pdf_name == "extnovoeffprod_G4_R1" || pdf_name == "extnovoeffprod_R1" || pdf_name == "supernovoeffprod2_G4_R1" || pdf_name == "supernovoeffprod3_G4_R1"){
		w.var("turnon_novoeff")->setConstant();
		w.var("slope_novoeff")->setConstant();
	}
	else if (pdf_name == "novosibirsk_G4_R2" || pdf_name == "supernovosibirsk1_G4_R2" || pdf_name == "supernovosibirsk2_G4_R2"){
	}
	else if (pdf_name == "novosibirsk_G4_R3" || pdf_name == "supernovosibirsk1_G4_R3" || pdf_name == "supernovosibirsk2_G4_R3"){
	}
	else if (pdf_name == "berneffprod9par_G4_R1" || pdf_name == "berneffprod9par_R1"){
		w.var("slope_novoeff")->setConstant();
		w.var("turnon_novoeff")->setConstant();
	}
	else if (pdf_name == "berneffprod8par_G4_R1" || pdf_name == "berneffprod8par_R1"){
		w.var("slope_novoeff")->setConstant();
		w.var("turnon_novoeff")->setConstant();
	}
	else if (pdf_name == "dijetv3logprod_G4_R1"){
	}
	else if (pdf_name == "superdijeteffprod2_G4_R1" || pdf_name == "superdijeteffprod2_R1"){
		w.var("slope_novoeff")->setConstant();
		w.var("turnon_novoeff")->setConstant();
	}
	else throw invalid_argument("Wrong pdf name: " + pdf_name + " at FixTurnOn::FixTurnOn");
}

void UnFixTurnOn(RooWorkspace &w, const std::string& pdf_name){
	/*
	 * Function to un-fix turn-on curve
	 */
	if(pdf_name == "extnovoeffprod_G4_R1" || pdf_name == "extnovoeffprod_R1"){
		w.var("turnon_novoeff")->setConstant(0);
		w.var("slope_novoeff")->setConstant(0);
	}
	else if (pdf_name == "novosibirsk_G4_R2"){
	}
	else if (pdf_name == "berneffprod9par_G4_R1" || pdf_name == "berneffprod9par_R1"){
		w.var("slope_novoeff")->setConstant(0);
		w.var("turnon_novoeff")->setConstant(0);
	}
	else if (pdf_name == "berneffprod8par_G4_R1" || pdf_name == "berneffprod8par_R1"){
		w.var("slope_novoeff")->setConstant(0);
		w.var("turnon_novoeff")->setConstant(0);
	}
	else if (pdf_name == "dijetv3logprod_G4_R1"){
	}
	else if (pdf_name == "superdijeteffprod2_G4_R1" || pdf_name == "superdijeteffprod2_R1"){
		w.var("slope_novoeff")->setConstant(0);
		w.var("turnon_novoeff")->setConstant(0);
	}
	else throw invalid_argument("Wrong pdf name: " + pdf_name + " at FixTurnOn::FixTurnOn");
}

double GetGausRandomisation(TRandom3& r, const double& mean, const double&sigma, const double& max, const double& min){
	/*
	 * Gaus randomisation within the limits
	 */
	double res = r.Gaus(mean,sigma);
	if(res >= max || res <= min) res = GetGausRandomisation(r,mean,sigma,max,min);
	return res;
}

void RandomiseVarsInWorkspace(RooWorkspace &w, const std::vector<RooRealVar>& init_vars, const int& nsigma){
	/*
	 * Function to randomise vars in workspace
	 */
	TRandom3 r(0);
	for(const auto& v : init_vars){
		auto init_val = v.getValV();
		auto init_err = v.getError();
		auto up_range = v.getMax();
		auto low_range= v.getMin();

		auto &var = *GetFromRooWorkspace<RooRealVar>(w,v.GetName());
		var.setVal(GetGausRandomisation(r,init_val,init_err * nsigma, up_range, low_range));
		cout<<"Randomise var: "<<var.GetName()<<" ";
		var.Print();
	}
}

