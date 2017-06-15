#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/macros/Drawer/RatioPlots.cpp"
#include "Analysis/MssmHbb/src/utilLib.cpp"
#include "TEfficiency.h"
#include "TPaveStats.h"
#include "TVirtualFitter.h"
#include "TBinomialEfficiencyFitter.h"
#include "TGraphAsymmErrors.h"
#include "TClass.h"
#include "TKey.h"
#include "TApplication.h"
#include "TROOT.h"

double errorPropagation(const double & num, const double &e_num, const double &den, const  double &e_den);

TEfficiency * GetEfficiency(TFile * f, const std::string & name_num, const std::string & name_den);
void BayesianErrors(TEfficiency *h);
double finaleFunction(double *x, double *par);
double sigmoid(double *x, double *par);
double sigmoid_sigmoid(double *x, double *par);
double func_erf(double *x, double *par);
double erfc(double *x, double*par);

double ratio_erf(double *x, double *par);
double ratio_sigmoid(double *x, double *par);
double ratio_sigmoid2(double *x, double *par);

double sgm_sgm_erf(double *x, double *par);
double test_fcn(double *x, double  *par);

double ratio_finale(double *x,double *par);


int fit_correction(){
	gStyle->SetOptStat(0000);
	gStyle->SetOptFit(0000);
	// Set stat options
	gStyle->SetStatY(0.5);
	// Set y-position (fraction of pad size)
	gStyle->SetStatX(0.8);
	// Set x-position (fraction of pad size)

//	ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
//	ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-06);

	TH1::SetDefaultSumw2();
	const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));

	TFile *fDataLowM = new TFile( ( cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_Data2016_ReReco.root").c_str());
	CheckZombie(*fDataLowM);
//	TFile *fDataHighM = new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/TriggerEff_highM_Run2015D-16Dec2015-v1.root");

   //BG MC:
   TFile *BgMC = new TFile( ( cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_ReReco_lowM_QCD_Pt.root").c_str());
   CheckZombie(*BgMC);
//   TFile *BgMC_high = new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/TriggerEff_highM_QCD.root");

   //Setup style
   TH1::SetDefaultSumw2();

   //Setup Ratio
   RatioPlots ratio(PRELIMINARY);
   ratio.SetRatioRange(0.5,1.5);
   ratio.SetRatioTitle("Data/MC");
//   ratio.set(PRIVATE);
   TH1 *nullH = nullptr;

   /****************************LOW M*************************/
  TPaveText *eta_def[3];
//  TLegend *leg_low_m = (TLegend*) ratio.legend("bottom,right",3,0.6);
  TLegend *leg_low_m = new TLegend(0.65,0.3,0.9,0.45);
  TH2D *h_lowM_pt_eta_num 		= (TH2D*)	BgMC->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1eta1");
  TH2D *h_lowM_pt_eta_num_data 	= (TH2D*)	fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1eta1");
  TH2D *h_lowM_pt_eta_den 		= (TH2D*)	BgMC->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1eta1");
  TH2D *h_lowM_pt_eta_den_data 	= (TH2D*)	fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1eta1");

  TH1D *h_lowM_pt_num_slice[3], *h_lowM_pt_num_slice_data[3], *h_lowM_pt_den_slice[3], *h_lowM_pt_den_slice_data[3];
  TEfficiency *eff_lowM_pt_slice[3], *eff_lowM_pt_slice_data[3];
  TF1 *lowm_fit_sigmoid2[3], *lowm_fit_finale[3];
  TF1 *lowm_fit_sigmoid2_mc[3], *lowm_fit_finale_mc[3];
  TGraphAsymmErrors *ratio_lowM_pt_slice[3], *gr_lowM_eff[3], *gr_lowM_eff_data[3], *gr_lowM_data_vs_fit[3], *gr_lowM_mc_vs_fit[3];
  TCanvas *can_lowM_pt_slice[3];
  TPaveStats *st[3],*st_data[3];
  TPad *pad_lowM_fit_top[3], *pad_lowM_fit_bot[3];
  std::string eta_bins[3]={"0-0.6","0.6-1.7","1.7-2.2"};

  double *num_x;
  double xMin = 50;
  double xMax = 500;
  double yMin = 0;//total->GetMinimum();
  double yMax = 1.2;
  for(int i = 0; i<3;++i){
	  can_lowM_pt_slice[i] = new TCanvas(("can_lowM_pt_slice_"+std::to_string(i)).c_str(),"canv",1000,800);
	  pad_lowM_fit_top[i] = new TPad(("pad_lowM_fit_top"+std::to_string(i+1)).c_str(),"pad1",0,0.1,1,1);
	  pad_lowM_fit_top[i] -> SetBottomMargin(0.21);
	  pad_lowM_fit_top[i] -> SetRightMargin(0.05);
	  pad_lowM_fit_top[i] -> SetLeftMargin(0.16);
	  pad_lowM_fit_top[i] -> Draw();
	  pad_lowM_fit_top[i] -> cd();

	  //Frame to adjust drawing style!
	  TH2F *frame = new TH2F("frame","",2,xMin,xMax,2,yMin,yMax);
	  frame->SetTitle(";p^{(1)}_{T}, [GeV]; #epsilon");
	  ratio.SetTopStyle(frame);
	  frame->Draw();

	  h_lowM_pt_num_slice[i] = (TH1D*)  h_lowM_pt_eta_num->ProjectionX(("Slice_num_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  h_lowM_pt_den_slice[i] = (TH1D*)  h_lowM_pt_eta_den->ProjectionX(("Slice_den_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  eff_lowM_pt_slice[i] = new TEfficiency(*h_lowM_pt_num_slice[i],*h_lowM_pt_den_slice[i]);
	  eff_lowM_pt_slice[i]->SetMarkerStyle(21);
	  eff_lowM_pt_slice[i]->SetMarkerColor(kMagenta);
	  BayesianErrors(eff_lowM_pt_slice[i]);

	  h_lowM_pt_num_slice_data[i] = (TH1D*)  h_lowM_pt_eta_num_data->ProjectionX(("Slice_num_data_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  h_lowM_pt_den_slice_data[i] = (TH1D*)  h_lowM_pt_eta_den_data->ProjectionX(("Slice_den_data_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  eff_lowM_pt_slice_data[i] = new TEfficiency(*h_lowM_pt_num_slice_data[i],*h_lowM_pt_den_slice_data[i]);
	  BayesianErrors(eff_lowM_pt_slice_data[i]);
	  eff_lowM_pt_slice_data[i]->SetMarkerColor(1);
	  eff_lowM_pt_slice_data[i]->SetMarkerStyle(20);
	  eff_lowM_pt_slice[i]->SetTitle(";p^{(1)}_{T}, [GeV]; #epsilon");
	  eta_def[i] = new TPaveText(0.7,0.6,0.9,0.68,"NDC");
	  eta_def[i]->SetFillColor(0);
	  eta_def[i]->SetBorderSize(0);
	  eta_def[i]->AddText(("|#eta| = "+eta_bins[i]).c_str());
	  eff_lowM_pt_slice[i]->Draw("Psame");
	  eff_lowM_pt_slice_data[i]->Draw("P same");

	  gPad->Update();
	  gr_lowM_eff_data[i] = (TGraphAsymmErrors*) eff_lowM_pt_slice_data[i]->GetPaintedGraph();
	  gr_lowM_eff[i] = (TGraphAsymmErrors*) eff_lowM_pt_slice[i]->GetPaintedGraph();
	  gr_lowM_eff[i]->GetYaxis()->SetRangeUser(0.,1.1);
	  gr_lowM_eff[i]->GetXaxis()->SetTitleOffset(1.1);
	  gr_lowM_eff[i]->Draw("Psame");
	  gr_lowM_eff_data[i]->Draw("P same");
	  ratio.drawStandardTitle();
	  eta_def[i]->Draw();

	  lowm_fit_sigmoid2[i] = new TF1(("fit_"+std::to_string(i+1)).c_str(),sigmoid_sigmoid,100,500,4);
	  lowm_fit_sigmoid2[i]->SetParameters(0.99,100.,0.99,100.);
	  lowm_fit_sigmoid2[i]->SetLineColor(1);
	  gr_lowM_eff_data[i]->Fit(("fit_"+std::to_string(i+1)).c_str(),"R EX0");
	  can_lowM_pt_slice[i]->Modified();
	  can_lowM_pt_slice[i]->Update();

	  //Uncomment if stat boxes are needed!!!!!

	  st_data[i] = (TPaveStats*)(gr_lowM_eff_data[i]->GetListOfFunctions()->FindObject("stats"));
	  if(st_data[i]){
		  st_data[i]->SetX1NDC(0.6);
	  	  st_data[i]->SetX2NDC(0.9);
	  	  st_data[i]->SetY1NDC(0.2);
	  	  st_data[i]->SetY2NDC(0.45);
	  }
	  pad_lowM_fit_top[i]->Modified();
	  pad_lowM_fit_top[i]->Update();

	  lowm_fit_sigmoid2_mc[i] = new TF1(("fit_mc_"+std::to_string(i+1)).c_str(),test_fcn,100,500,7);
	  lowm_fit_sigmoid2_mc[i]->SetParameters(4.01543e+01,-108.4,0.1182,108,0.99,-66.06,84.45);
	  lowm_fit_sigmoid2_mc[i]->SetLineColor(2);
	  gr_lowM_eff[i]->Fit(lowm_fit_sigmoid2_mc[i],"+R EX0");
	  pad_lowM_fit_top[i]->Modified();
	  pad_lowM_fit_top[i]->Update();
	  st[i] = (TPaveStats*)(gr_lowM_eff[i]->GetListOfFunctions()->FindObject("stats"));
	  if(st[i]){
		  st[i]->SetX1NDC(0.6);
		  st[i]->SetX2NDC(0.9);
		  st[i]->SetY1NDC(0.45);
		  st[i]->SetY2NDC(0.7);
		  pad_lowM_fit_top[i]->Modified();
		  pad_lowM_fit_top[i]->Update();
	  }


	  if(i == 0 ){
		  leg_low_m->AddEntry(eff_lowM_pt_slice_data[i],"JetHT, 2016","pl");
		  leg_low_m->AddEntry(eff_lowM_pt_slice[i],"QCD, #hat{p_{T}}","pl");
	  }
	  leg_low_m->Draw();
	  ///////////////Pad2////////////////////
	  can_lowM_pt_slice[i]->cd();
	  pad_lowM_fit_bot[i] = new TPad("pad2","pad2",0,0.0,1,0.265);
	  pad_lowM_fit_bot[i] -> SetTopMargin(0.0);
	  pad_lowM_fit_bot[i] -> SetLeftMargin( pad_lowM_fit_top[i]->GetLeftMargin());
	  pad_lowM_fit_bot[i] -> SetRightMargin( pad_lowM_fit_top[i]->GetRightMargin());
	  pad_lowM_fit_bot[i] -> SetBottomMargin(0.28);
	  pad_lowM_fit_bot[i] -> Draw();
	  pad_lowM_fit_bot[i] -> cd();

	  TH2F *frame2 = new TH2F("frame2","",2,xMin,xMax,2,0.5,1.5);
	  frame2->SetTitle("");
	  frame2->GetXaxis()->SetTitle(frame->GetXaxis()->GetTitle());
	  frame2->GetXaxis()->SetTickLength(frame->GetXaxis()->GetTickLength()*3);
	  ratio.SetBottomStyle(frame2);
	  frame2->SetMinimum(0.5);
	  frame2->SetMaximum(1.5);
	  frame2->Draw();

	  gr_lowM_data_vs_fit[i] = (TGraphAsymmErrors*) gr_lowM_eff_data[i]->Clone(("gr_lowM_data_vs_fit_"+std::to_string(i+1)).c_str());
	  gr_lowM_mc_vs_fit[i] = (TGraphAsymmErrors*) gr_lowM_eff[i]->Clone(("gr_lowM_mc_vs_fit_"+std::to_string(i+1)).c_str());

	  num_x = gr_lowM_eff_data[i]->GetX();
	  double *num_y = gr_lowM_eff_data[i]->GetY();
	  for( int j = 1; j <= gr_lowM_eff_data[i]->GetN(); ++j){
		double den_y = lowm_fit_sigmoid2[i]->Eval(num_x[j-1]);
		double val = num_y[j-1]/den_y;
		gr_lowM_data_vs_fit[i]->SetPoint(j-1,num_x[j-1],val);
		double el = gr_lowM_eff_data[i]->GetErrorYlow(j-1)/den_y;
		gr_lowM_data_vs_fit[i]->SetPointEYlow(j-1,el);
		double eu = gr_lowM_eff_data[i]->GetErrorYhigh(j-1)/den_y;
		gr_lowM_data_vs_fit[i]->SetPointEYhigh(j-1,eu);
	  }

	  num_y = gr_lowM_eff[i]->GetY();
	  for( int j = 1; j <= gr_lowM_eff_data[i]->GetN(); ++j){
		//MC
		double den_y = lowm_fit_sigmoid2_mc[i]->Eval(num_x[j-1]);
		double val = num_y[j-1]/den_y;
		gr_lowM_mc_vs_fit[i]->SetPoint(j-1,num_x[j-1],val);
		double el = gr_lowM_eff[i]->GetErrorYlow(j-1)/den_y;
		gr_lowM_mc_vs_fit[i]->SetPointEYlow(j-1,el);
		double eu = gr_lowM_eff[i]->GetErrorYhigh(j-1)/den_y;
		gr_lowM_mc_vs_fit[i]->SetPointEYhigh(j-1,eu);

//		std::cout<<"WTF "<<gr_lowM_eff_data[i]->GetN()<<std::endl;
	  }
	  gr_lowM_data_vs_fit[i]->GetListOfFunctions()->RemoveAll();
	  gr_lowM_mc_vs_fit[i]->GetListOfFunctions()->RemoveAll();
	  ratio.SetBottomStyle(gr_lowM_data_vs_fit[i]);
	  gr_lowM_data_vs_fit[i] -> Draw("Psame");
	  gr_lowM_data_vs_fit[i] -> GetYaxis() -> SetTitle("Val./Fit");
	  gr_lowM_data_vs_fit[i] -> GetXaxis() -> SetTitle(gr_lowM_eff[i]->GetXaxis()->GetTitle());
	  gr_lowM_data_vs_fit[i] -> GetXaxis() -> SetTitleOffset(gr_lowM_eff[i]->GetXaxis()->GetTitleOffset());
	  gr_lowM_mc_vs_fit[i] -> Draw("P same");

	  can_lowM_pt_slice[i]->cd();
	  can_lowM_pt_slice[i]->Modified();
	  can_lowM_pt_slice[i]->Update();
	  can_lowM_pt_slice[i]->Print(("../pictures/TriggerPerformance/FittedTriggerEff_lowM_eta_bin_"+std::to_string(i+1) + "_PF60_PF100_pt1_eta1.pdf").c_str());
  }

  //Draw All 3 eta bins fit at 1 plot
  TCanvas *can3in1_fits_lowM_data = new TCanvas("can3in1_fits_lowM_data","can3in1_fits in DATA",1000,800);
  ratio.drawStandardTitle();
  TH1D *ax_plot = new TH1D("ax_plot","",10,60.,500.);
  ax_plot->GetYaxis()->SetRangeUser(0.,1.1);
  ax_plot->GetXaxis()->SetRangeUser(60.,500.);
  ax_plot->GetXaxis()->SetTitleOffset(1.);
  ax_plot->SetTitle(";p^{(1)}_{T}, [GeV]; #epsilon");
  ax_plot->Draw();
  TLegend *leg3in1_leg_data = new TLegend(0.4,0.4,0.85,0.8);
  lowm_fit_sigmoid2[0]->SetLineColor(1);
  lowm_fit_sigmoid2[0]->Draw("same");
  lowm_fit_sigmoid2[1]->SetLineColor(2);
  lowm_fit_sigmoid2[1]->SetLineStyle(kDashed);
  lowm_fit_sigmoid2[1]->Draw("same");
  lowm_fit_sigmoid2[2]->SetLineColor(6);
  lowm_fit_sigmoid2[2]->SetLineStyle(9);
  lowm_fit_sigmoid2[2]->Draw("same");

  leg3in1_leg_data->SetHeader("JetHT lowM, ReReco-2016");
  leg3in1_leg_data->AddEntry(lowm_fit_sigmoid2[0],"0<|#eta|<0.6","l");
  leg3in1_leg_data->AddEntry(lowm_fit_sigmoid2[1],"0.6<|#eta|<1.7","l");
  leg3in1_leg_data->AddEntry(lowm_fit_sigmoid2[2],"1.7<|#eta|<2.2","l");
  leg3in1_leg_data->Draw();
  ratio.drawStandardTitle();
  can3in1_fits_lowM_data->Print("../pictures/TriggerPerformance/FittedTriggerEff_lowM_Data_AllEta_fits_PF60_PF100_pt1_eta1.pdf");

  TCanvas *can3in1_fits_mc = new TCanvas("can3in1_fits_mc","can3in1_fits in MC",1000,800);
  ax_plot->Draw();
  TLegend *leg3in1_leg_mc = new TLegend(0.4,0.4,0.85,0.8);
  lowm_fit_sigmoid2_mc[0]->SetLineColor(1);
  lowm_fit_sigmoid2_mc[0]->Draw("same");
  lowm_fit_sigmoid2_mc[1]->SetLineColor(2);
  lowm_fit_sigmoid2_mc[1]->SetLineStyle(kDashed);
  lowm_fit_sigmoid2_mc[1]->Draw("same");
  lowm_fit_sigmoid2_mc[2]->SetLineColor(6);
  lowm_fit_sigmoid2_mc[2]->SetLineStyle(9);
  lowm_fit_sigmoid2_mc[2]->Draw("same");

  leg3in1_leg_mc->SetHeader("QCD, #hat{p_{T}} MC             ");
  leg3in1_leg_mc->AddEntry(lowm_fit_sigmoid2_mc[0],"0<|#eta|<0.6","l");
  leg3in1_leg_mc->AddEntry(lowm_fit_sigmoid2_mc[1],"0.6<|#eta|<1.7","l");
  leg3in1_leg_mc->AddEntry(lowm_fit_sigmoid2_mc[2],"1.7<|#eta|<2.2","l");
  leg3in1_leg_mc->Draw("same");
  ratio.drawStandardTitle();
  can3in1_fits_mc->Print("../pictures/TriggerPerformance/FittedTriggerEff_lowM_MC_AllEta_fits_PF60_PF100_pt1_eta1.pdf");

  //Finale
  TF1 *f1_ratio_lowM_erf[3], *f1_ratio_lowM_sigmoid[3], *f1_ratio_lowM_sigmoid2[3];
  TCanvas *can_lowM_correction_funct[3];
  for(int i=0;i<3;++i){
  can_lowM_correction_funct[i] = new TCanvas();
  	ratio_lowM_pt_slice[i] = (TGraphAsymmErrors*) ratio.DrawRatio((TGraphAsymmErrors*) gr_lowM_eff_data[i]->Clone(),(TGraphAsymmErrors*) gr_lowM_eff[i]->Clone(), can_lowM_correction_funct[i],nullptr);
  	can_lowM_correction_funct[i]->cd();
  	ratio_lowM_pt_slice[i]->GetXaxis()->SetTitleOffset(1.1);
  	ratio_lowM_pt_slice[i]->GetXaxis()->SetTitle("p^{(1)}_{T}, [GeV]");
  	ratio_lowM_pt_slice[i]->Draw("AP");
  	eta_def[i]->Draw();
  	ratio.drawStandardTitle();

  	f1_ratio_lowM_sigmoid2[i]= new TF1(("ratio_sigmoid2_"+std::to_string(i+1)).c_str(),ratio_finale,100,500,11);
  	f1_ratio_lowM_sigmoid2[i]->SetParameters(lowm_fit_sigmoid2[i]->GetParameter(0),
  										lowm_fit_sigmoid2[i]->GetParameter(1),
										lowm_fit_sigmoid2[i]->GetParameter(2),
										lowm_fit_sigmoid2[i]->GetParameter(3),
										lowm_fit_sigmoid2_mc[i]->GetParameter(0),
										lowm_fit_sigmoid2_mc[i]->GetParameter(1),
										lowm_fit_sigmoid2_mc[i]->GetParameter(2),
										lowm_fit_sigmoid2_mc[i]->GetParameter(3),
										lowm_fit_sigmoid2_mc[i]->GetParameter(4),
										lowm_fit_sigmoid2_mc[i]->GetParameter(5),
										lowm_fit_sigmoid2_mc[i]->GetParameter(6));
  	f1_ratio_lowM_sigmoid2[i]->SetNpx(500);
  	f1_ratio_lowM_sigmoid2[i]->SetLineColor(kBlack);
  	ratio_lowM_pt_slice[i]->GetListOfFunctions()->RemoveAll();
  	ratio_lowM_pt_slice[i]->Draw("P same");
  	f1_ratio_lowM_sigmoid2[i]->Draw("same");
  	can_lowM_correction_funct[i]->Modified();
  	can_lowM_correction_funct[i]->Update();
  	std::cout<<"VALUE at 100: "<<f1_ratio_lowM_sigmoid2[i]->Eval(100)<<std::endl;
  	can_lowM_correction_funct[i]->Print(("../pictures/TriggerPerformance/CorrectionFunction_PF60_PF100_pt1_eta1_eta_bin"+std::to_string(i+1) + ".pdf").c_str());
  }


  /************************************************************/

  /**************************high M****************************/
  /*
  TLegend *leg_highM = new TLegend(0.3,0.5,0.5,0.7);
  TH2D *h_highM_pt_eta_num 		= (TH2D*)	BgMC_high->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF80_PF160_pt1eta1");
  TH2D *h_highM_pt_eta_num_data 	= (TH2D*)	fDataHighM->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF80_PF160_pt1eta1");
  TH2D *h_highM_pt_eta_den 		= (TH2D*)	BgMC_high->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF80_PF160_pt1eta1");
  TH2D *h_highM_pt_eta_den_data 	= (TH2D*)	fDataHighM->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF80_PF160_pt1eta1");

  TH1D *h_highM_pt_num_slice[3], *h_highM_pt_num_slice_data[3], *h_highM_pt_den_slice[3], *h_highM_pt_den_slice_data[3];
  TEfficiency *eff_highM_pt_slice[3], *eff_highM_pt_slice_data[3];
  TF1 *highM_fit_sigmoid2[3], *highM_fit_finale[3];
  TF1 *highM_fit_sigmoid2_mc[3], *highM_fit_finale_mc[3];
  TGraphAsymmErrors *ratio_highM_pt_slice[3], *gr_highM_eff[3], *gr_highM_eff_data[3], *gr_highM_data_vs_fit[3], *gr_highM_mc_vs_fit[3];
  TCanvas *can_highM_pt_slice[3];
  TPad *pad_highM_fit_top[3], *pad_highM_fit_bot[3];

  for(int i = 0; i<3;++i){
	  can_highM_pt_slice[i] = new TCanvas(("can_highM_pt_slice_"+std::to_string(i)).c_str(),"canv",1000,800);
	  pad_highM_fit_top[i] = new TPad(("pad_highM_fit_top"+std::to_string(i+1)).c_str(),"pad1",0,0.3,1,1);
	  pad_highM_fit_top[i] -> SetBottomMargin(0.01);
	  pad_highM_fit_top[i] -> Draw();
	  pad_highM_fit_top[i] -> cd();

	  h_highM_pt_num_slice[i] = (TH1D*)  h_highM_pt_eta_num->ProjectionX(("Slice_num_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  h_highM_pt_den_slice[i] = (TH1D*)  h_highM_pt_eta_den->ProjectionX(("Slice_den_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  eff_highM_pt_slice[i] = new TEfficiency(*h_highM_pt_num_slice[i],*h_highM_pt_den_slice[i]);
	  eff_highM_pt_slice[i]->SetMarkerStyle(21);
	  eff_highM_pt_slice[i]->SetMarkerColor(kMagenta);
	  BayesianErrors(eff_highM_pt_slice[i]);

	  h_highM_pt_num_slice_data[i] = (TH1D*)  h_highM_pt_eta_num_data->ProjectionX(("Slice_num_data_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  h_highM_pt_den_slice_data[i] = (TH1D*)  h_highM_pt_eta_den_data->ProjectionX(("Slice_den_data_"+std::to_string(i)).c_str(),i+1,i+1,"E");
	  eff_highM_pt_slice_data[i] = new TEfficiency(*h_highM_pt_num_slice_data[i],*h_highM_pt_den_slice_data[i]);
	  BayesianErrors(eff_highM_pt_slice_data[i]);
	  eff_highM_pt_slice_data[i]->SetMarkerColor(1);
	  eff_highM_pt_slice_data[i]->SetMarkerStyle(20);
	  eff_highM_pt_slice[i]->SetTitle(";p^{(1)}_{T}, [GeV]; #epsilon");
	  eta_def[i] = new TPaveText(0.7,0.7,0.9,0.8,"NDC");
	  eta_def[i]->SetFillColor(0);
	  eta_def[i]->SetBorderSize(0);
	  eta_def[i]->AddText(("|#eta| = "+eta_bins[i]).c_str());
	  eff_highM_pt_slice[i]->Draw("AP");
	  eff_highM_pt_slice_data[i]->Draw("P same");

	  gPad->Update();
	  gr_highM_eff_data[i] = (TGraphAsymmErrors*) eff_highM_pt_slice_data[i]->GetPaintedGraph();
	  gr_highM_eff[i] = (TGraphAsymmErrors*) eff_highM_pt_slice[i]->GetPaintedGraph();
	  gr_highM_eff[i]->GetYaxis()->SetRangeUser(0.,1.1);
	  gr_highM_eff[i]->GetXaxis()->SetTitleOffset(1.1);
	  gr_highM_eff[i]->Draw("AP");
	  gr_highM_eff_data[i]->Draw("P same");
	  ratio.drawStandardTitle();
	  eta_def[i]->Draw();

	  highM_fit_sigmoid2[i] = new TF1(("fit_highM_"+std::to_string(i+1)).c_str(),sigmoid_sigmoid,160,500,4);
//	  highM_fit_sigmoid2[i] = new TF1(("fit_"+std::to_string(i+1)).c_str(),test_fcn,100,500,7);
	  highM_fit_sigmoid2[i]->SetParameters(0.99,150.,0.99,150.);
	  highM_fit_sigmoid2[i]->SetLineColor(1);
	  gr_highM_eff_data[i]->Fit(highM_fit_sigmoid2[i],"R EX0 N");
	  highM_fit_sigmoid2[i]->SetParameters(highM_fit_sigmoid2[i]->GetParameters());
	  gr_highM_eff_data[i]->Fit(highM_fit_sigmoid2[i],"R EX0");
	  can_highM_pt_slice[i]->Modified();
	  can_highM_pt_slice[i]->Update();
	  st_data[i] = (TPaveStats*)(gr_highM_eff_data[i]->GetListOfFunctions()->FindObject("stats"));
	  if(st_data[i]){
		  st_data[i]->SetX1NDC(0.6);
	  	  st_data[i]->SetX2NDC(0.9);
	  	  st_data[i]->SetY1NDC(0.2);
	  	  st_data[i]->SetY2NDC(0.45);
	  }
	  pad_highM_fit_top[i]->Modified();
	  pad_highM_fit_top[i]->Update();

//
//	  highM_fit_sigmoid2_mc[i] = new TF1(("fit_mc_"+std::to_string(i+1)).c_str(),sigmoid_sigmoid,100,500,4);
	  highM_fit_sigmoid2_mc[i] = new TF1(("fit_highM_mc_"+std::to_string(i+1)).c_str(),test_fcn,160,500,7);
	  highM_fit_sigmoid2_mc[i]->SetParameters(4.01543e+01,-100.4,0.1182,100,0.99,-66.06,84.45);
	  highM_fit_sigmoid2_mc[i]->SetLineColor(2);
	  gr_highM_eff[i]->Fit(highM_fit_sigmoid2_mc[i],"+R EX0 N");
	  highM_fit_sigmoid2_mc[i]->SetParameters(highM_fit_sigmoid2_mc[i]->GetParameters());
	  gr_highM_eff[i]->Fit(highM_fit_sigmoid2_mc[i],"R EX0");
	  pad_highM_fit_top[i]->Modified();
	  pad_highM_fit_top[i]->Update();
	  st[i] = (TPaveStats*)(gr_highM_eff[i]->GetListOfFunctions()->FindObject("stats"));
	  if(st[i]){
		  st[i]->SetX1NDC(0.6);
		  st[i]->SetX2NDC(0.9);
		  st[i]->SetY1NDC(0.45);
		  st[i]->SetY2NDC(0.7);
		  pad_highM_fit_top[i]->Modified();
		  pad_highM_fit_top[i]->Update();
	  }
	  if(i==0){
		  leg_highM->AddEntry(eff_highM_pt_slice_data[i],"JetHT, 2015D","pl");
		  leg_highM->AddEntry(eff_highM_pt_slice[i],"QCD, #hat{p_{T}}","pl");
	  }
	  leg_highM->Draw();
	  ///////////////Pad2////////////////////
	  can_highM_pt_slice[i]->cd();
	  pad_highM_fit_bot[i] = new TPad("pad2","pad2",0,0.0,1,0.3);
	  pad_highM_fit_bot[i] -> SetTopMargin(0.0);
	  pad_highM_fit_bot[i] -> SetBottomMargin(0.35);
	  pad_highM_fit_bot[i] -> Draw();
	  pad_highM_fit_bot[i] -> cd();

	  gr_highM_data_vs_fit[i] = (TGraphAsymmErrors*) gr_highM_eff_data[i]->Clone(("gr_highM_data_vs_fit_"+std::to_string(i+1)).c_str());
	  gr_highM_mc_vs_fit[i] = (TGraphAsymmErrors*) gr_highM_eff[i]->Clone(("gr_highM_mc_vs_fit_"+std::to_string(i+1)).c_str());

	  num_x = gr_highM_eff_data[i]->GetX();
	  double *num_y = gr_highM_eff_data[i]->GetY();
	  for( int j = 1; j <= gr_highM_eff_data[i]->GetN(); ++j){
		double den_y = highM_fit_sigmoid2[i]->Eval(num_x[j-1]);
		double val = num_y[j-1]/den_y;
		gr_highM_data_vs_fit[i]->SetPoint(j-1,num_x[j-1],val);
		double el = gr_highM_eff_data[i]->GetErrorYlow(j-1)/den_y;
		gr_highM_data_vs_fit[i]->SetPointEYlow(j-1,el);
		double eu = gr_highM_eff_data[i]->GetErrorYhigh(j-1)/den_y;
		gr_highM_data_vs_fit[i]->SetPointEYhigh(j-1,eu);
	  }

	  num_y = gr_highM_eff[i]->GetY();
	  for( int j = 1; j <= gr_highM_eff_data[i]->GetN(); ++j){
		//MC
		double den_y = highM_fit_sigmoid2_mc[i]->Eval(num_x[j-1]);
		double val = num_y[j-1]/den_y;
		gr_highM_mc_vs_fit[i]->SetPoint(j-1,num_x[j-1],val);
		double el = gr_highM_eff[i]->GetErrorYlow(j-1)/den_y;
		gr_highM_mc_vs_fit[i]->SetPointEYlow(j-1,el);
		double eu = gr_highM_eff[i]->GetErrorYhigh(j-1)/den_y;
		gr_highM_mc_vs_fit[i]->SetPointEYhigh(j-1,eu);

//		std::cout<<"WTF "<<gr_highM_eff_data[i]->GetN()<<std::endl;
	  }
	  gr_highM_data_vs_fit[i]->GetListOfFunctions()->RemoveAll();
	  gr_highM_mc_vs_fit[i]->GetListOfFunctions()->RemoveAll();
	  ratio.SetBottomStyle(gr_highM_data_vs_fit[i]);
	  gr_highM_data_vs_fit[i] -> Draw("AP");
	  gr_highM_data_vs_fit[i] -> GetYaxis() -> SetTitle("Val./Fit");
	  gr_highM_data_vs_fit[i] -> GetXaxis() -> SetTitle(gr_highM_eff[i]->GetXaxis()->GetTitle());
	  gr_highM_data_vs_fit[i] -> GetXaxis() -> SetTitleOffset(gr_highM_eff[i]->GetXaxis()->GetTitleOffset());
	  gr_highM_mc_vs_fit[i] -> Draw("P same");

	  can_highM_pt_slice[i]->cd();
	  can_highM_pt_slice[i]->Modified();
	  can_highM_pt_slice[i]->Update();
	  can_highM_pt_slice[i]->Print(("../pictures/TriggerPerformance/FittedTriggerEff_highM_eta_"+std::to_string(i+1) + "_PF80_PF160_pt1_eta1.pdf").c_str());
  }

  //Draw All 3 eta bins fit at 1 plot
  TCanvas *can3in1_highM_fits_highM_data = new TCanvas("can3in1_highM_fits_highM_data","can3in1_highM_fits in DATA",1000,800);
  TH1D *ax_plot_high = new TH1D("ax_plot_high","",10,60.,500.);
  ax_plot_high->GetYaxis()->SetRangeUser(0.,1.1);
  ax_plot_high->GetXaxis()->SetRangeUser(100.,500.);
  ax_plot_high->GetXaxis()->SetTitleOffset(1.);
  ax_plot_high->SetTitle(";p^{(1)}_{T}, [GeV]; #epsilon");
  ax_plot_high->Draw();
  TLegend *leg3in1_highM_leg_data = new TLegend(0.5,0.5,0.8,0.8);
  highM_fit_sigmoid2[0]->SetLineColor(1);
  highM_fit_sigmoid2[0]->Draw("same");
  highM_fit_sigmoid2[1]->SetLineColor(2);
  highM_fit_sigmoid2[1]->SetLineStyle(kDashed);
  highM_fit_sigmoid2[1]->Draw("same");
  highM_fit_sigmoid2[2]->SetLineColor(6);
  highM_fit_sigmoid2[2]->SetLineStyle(9);
  highM_fit_sigmoid2[2]->Draw("same");

  leg3in1_highM_leg_data->SetHeader("JetHT highM, 2015D");
  leg3in1_highM_leg_data->AddEntry(highM_fit_sigmoid2[0],"0<|#eta|<0.6","l");
  leg3in1_highM_leg_data->AddEntry(highM_fit_sigmoid2[1],"0.6<|#eta|<1.7","l");
  leg3in1_highM_leg_data->AddEntry(highM_fit_sigmoid2[2],"1.7<|#eta|<2.2","l");
  leg3in1_highM_leg_data->Draw();
  ratio.drawStandardTitle();
  can3in1_highM_fits_highM_data->Print("../pictures/TriggerPerformance/FittedTriggerEff_highM_Data_AllEta_fits_PF80_PF160_pt1_eta1.pdf");

  TCanvas *can3in1_highM_fits_mc = new TCanvas("can3in1_highM_fits_mc","can3in1_highM_fits in MC",1000,800);
  ax_plot_high->Draw();
  TLegend *leg3in1_highM_leg_mc = new TLegend(0.5,0.5,0.7,0.7);
  highM_fit_sigmoid2_mc[0]->SetLineColor(1);
  highM_fit_sigmoid2_mc[0]->Draw("same");
  highM_fit_sigmoid2_mc[1]->SetLineColor(2);
  highM_fit_sigmoid2_mc[1]->SetLineStyle(kDashed);
  highM_fit_sigmoid2_mc[1]->Draw("same");
  highM_fit_sigmoid2_mc[2]->SetLineColor(6);
  highM_fit_sigmoid2_mc[2]->SetLineStyle(9);
  highM_fit_sigmoid2_mc[2]->Draw("same");

  leg3in1_highM_leg_mc->SetHeader("QCD, #hat{p_{T}} MC");
  leg3in1_highM_leg_mc->AddEntry(highM_fit_sigmoid2_mc[0],"0<|#eta|<0.6","l");
  leg3in1_highM_leg_mc->AddEntry(highM_fit_sigmoid2_mc[1],"0.6<|#eta|<1.7","l");
  leg3in1_highM_leg_mc->AddEntry(highM_fit_sigmoid2_mc[2],"1.7<|#eta|<2.2","l");
  leg3in1_highM_leg_mc->Draw("same");
  ratio.drawStandardTitle();
  can3in1_highM_fits_mc->Print("../pictures/TriggerPerformance/FittedTriggerEff_highM_MC_AllEta_fits_PF80_PF160_pt1_eta1.pdf");

  //Finale
  TF1 *f1_ratio_highM_erf[3], *f1_ratio_highM_sigmoid[3], *f1_ratio_highM_sigmoid2[3];
  ratio.drawStandardTitle();

  TCanvas *can_highM_correction_funct[3];
  for(int i=0;i<3;++i){
  can_highM_correction_funct[i] = new TCanvas();
  	ratio_highM_pt_slice[i] = (TGraphAsymmErrors*) ratio.DrawRatio((TGraphAsymmErrors*) gr_highM_eff_data[i]->Clone(),(TGraphAsymmErrors*) gr_highM_eff[i]->Clone(), can_highM_correction_funct[i],nullptr);
  	can_highM_correction_funct[i]->cd();
  	ratio_highM_pt_slice[i]->GetXaxis()->SetTitleOffset(1.1);
  	ratio_highM_pt_slice[i]->GetXaxis()->SetTitle("p^{(1)}_{T}, [GeV]");
  	ratio_highM_pt_slice[i]->Draw("AP");
  	eta_def[i]->Draw();
  	ratio.drawStandardTitle();

  	f1_ratio_highM_sigmoid2[i]= new TF1(("ratio_sigmoid2_"+std::to_string(i+1)).c_str(),ratio_finale,160,500,11);
  	f1_ratio_highM_sigmoid2[i]->SetParameters(highM_fit_sigmoid2[i]->GetParameter(0),
  										highM_fit_sigmoid2[i]->GetParameter(1),
										highM_fit_sigmoid2[i]->GetParameter(2),
										highM_fit_sigmoid2[i]->GetParameter(3),
										highM_fit_sigmoid2_mc[i]->GetParameter(0),
										highM_fit_sigmoid2_mc[i]->GetParameter(1),
										highM_fit_sigmoid2_mc[i]->GetParameter(2),
										highM_fit_sigmoid2_mc[i]->GetParameter(3),
										highM_fit_sigmoid2_mc[i]->GetParameter(4),
										highM_fit_sigmoid2_mc[i]->GetParameter(5),
										highM_fit_sigmoid2_mc[i]->GetParameter(6));
  	f1_ratio_highM_sigmoid2[i]->SetNpx(500);
  	f1_ratio_highM_sigmoid2[i]->SetLineColor(kBlack);
  	ratio_highM_pt_slice[i]->GetListOfFunctions()->RemoveAll();
  	ratio_highM_pt_slice[i]->Draw("P same");
  	f1_ratio_highM_sigmoid2[i]->Draw("same");
  	can_highM_correction_funct[i]->Modified();
  	can_highM_correction_funct[i]->Update();
  	std::cout<<"VALUE at 100: "<<f1_ratio_highM_sigmoid2[i]->Eval(160)<<std::endl;
  	can_highM_correction_funct[i]->Print(("../pictures/TriggerPerformance/CorrectionFunction_PF80_PF160_pt1_eta1_eta_bin"+std::to_string(i+1) + ".pdf").c_str());
  }

  /************************************************************/

	return 0;


}

void BayesianErrors(TEfficiency *h){
	h->SetStatisticOption(TEfficiency::kBUniform);
	h->SetPosteriorMode(1);
	h->UsesBayesianStat();
	h->SetUseWeightedEvents();
}

TEfficiency * GetEfficiency(TFile * f, const std::string & name_num, const std::string & name_den){
	if(!f){
		std::cerr<<"empty file were provided. Interupt."<<std::endl;
		exit(0);
	}

	//check 1D or 2D
	bool twoD = false;
	if(name_num.find("2D") != std::string::npos) twoD = true;

	TEfficiency *out;
	TH1D *num, *den;
	TH2D *num2D, *den2D;
	//1D case:
	if(!twoD){
		num = (TH1D*) f->Get(name_num.c_str());
		den = (TH1D*) f->Get(name_den.c_str());
		out = new TEfficiency(*num,*den);
		out->SetTitle(num->GetTitle());
	}
	else {
		num2D = (TH2D*) f->Get(name_num.c_str());
		den2D = (TH2D*) f->Get(name_den.c_str());
		out = new TEfficiency(*num2D,*den2D);
		out->SetTitle(num2D->GetTitle());
	}
	BayesianErrors(out);

	return out;
}

double finaleFunction(double *x, double *par)
{

	double sigmoid = 1./(1.+exp( -par[0] * (x[0] - par[1]) ));
	double erf = 1 / 2. * (1 + par[0]*TMath::Erf((x[0]-par[1])/(sqrt(2)*par[2])));
	double erfc = 1+ TMath::Erfc((par[5]-x[0])/(sqrt(2)*par[6]));
	double result =  0.5 * (1. + TMath::Erf((par[2]*x[0]-par[1])/(sqrt(2)*par[0])));
	double sigmoid2 = 1./ (1.+exp( -par[2] * (x[0] - par[3]) ));
	return sigmoid * sigmoid2;

}

double sigmoid(double *x, double *par){
	return 1./(1.+exp( -par[0] * (x[0] - par[1]) ));
}
double sigmoid_sigmoid(double *x, double *par){
	double sigmoid = 1./(1.+exp( -par[0] * (x[0] - par[1]) ));
	double sigmoid2 = 1./ (1.+exp( -par[2] * (x[0] - par[3]) ));
	return sigmoid * sigmoid2;
}
double func_erf(double *x, double *par){
	return 1 / 2. * (1 + par[0]*TMath::Erf((x[0]-par[1])/(sqrt(2)*par[2])));
}
double erfc(double *x, double*par){
	return 1 + TMath::Erfc((par[1]-x[0])/(sqrt(2)*par[0]));
}

double ratio_erf(double *x, double *par){
	return func_erf(x,par)/func_erf(x,par+3);
}

double ratio_sigmoid(double *x, double *par){
	return sigmoid(x,par)/sigmoid(x,par+2);
}

double ratio_sigmoid2(double *x, double *par){
	return sigmoid_sigmoid(x,par)/sigmoid_sigmoid(x,par+4);
}

double sgm_sgm_erf(double *x, double *par){
	return 1./(1.+exp( -par[0] * (x[0] - par[1]) ) + par[2] * (x[0]*x[0] + par[3]) + par[4]*(x[0]*x[0]*x[0] + par[5])) * par[6];
}

double test_fcn(double *x, double  *par){
	return sigmoid_sigmoid(x,par)  * func_erf(x,par+4);
}

double errorPropagation(const double & num, const double &e_num, const double &den, const  double &e_den){
	double err = 1./(den*den) * std::sqrt((den - num)*(den-num)*e_num*e_num + num*num*(e_den*e_den - e_num*e_num) );
	return err;
}

double ratio_finale(double *x,double *par){
	return sigmoid_sigmoid(x,par)/test_fcn(x,par+4);
}
