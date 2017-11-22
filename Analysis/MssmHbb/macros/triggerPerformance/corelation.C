#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/macros/Drawer/RatioPlots.cpp"
#include "Analysis/MssmHbb/src/utilLib.cpp"
#include "TEfficiency.h"

TEfficiency * GetEfficiency(TFile * f, const std::string & name_num, const std::string & name_den);
void BayesianErrors(TEfficiency *h);
void drawRatio(TH1D *, TF1 *,TCanvas *,std::string );
double finaleFunction(double *x, double *par);
double dEtaFitFunction(double *x, double *par);

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

int corelation(){

	TH1::SetDefaultSumw2();

	const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));
	TFile *fDataLowM = new TFile( ( cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_Data2016_ReReco.root").c_str());
	CheckZombie(*fDataLowM);
//	TFile *fDataHighM = new TFile("/afs/desy.de/user/s/shevchen/cms/cmssw-analysis/CMSSW_7_6_3_patch2/src/Analysis/MssmHbb/output/TriggerEff_highM_Run2015D-16Dec2015-v1.root");

   //BG MC:
   TFile *BgMC = new TFile( ( cmsswBase + "/src/Analysis/MssmHbb/output/TriggerEff_ReReco_lowM_QCD_Pt.root").c_str());
   CheckZombie(*BgMC);


   std::string error_mode = "";

   //Setup style
   HbbStyle style;
   style.set(PRIVATE);
   TH1::SetDefaultSumw2();

   //Setup Ratio
   RatioPlots ratio(PRIVATE);
   ratio.SetRatioRange(0.5,1.5);
   ratio.SetRatioTitle("Data/MC");
   TH1 *nullH = nullptr;

	//***************2D Efficiency in Data*****************//
	TCanvas *twoDplotLowM = new TCanvas("twoDplotLowM","Two dimensional Efficiency");
	twoDplotLowM->SetRightMargin(0.13);
    TH2D *TwoDEff_Num = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1vspt2");
    TH2D *TwoDEff_Den = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1vspt2");
	TEfficiency *TwoDEff_data = new TEfficiency(*TwoDEff_Num,*TwoDEff_Den);
	TH2D *TwoDEff_NumSym = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt2vspt1");
    TH2D *TwoDEff_DenSym = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt2vspt1");
    for(int binx = 1; binx <= TwoDEff_NumSym->GetNbinsX(); binx++)
    {
    	for(int biny = 1; biny <= TwoDEff_NumSym->GetNbinsY(); biny++)
    	{
    		if(binx == biny){
    			TwoDEff_NumSym->SetBinContent(binx,biny,0);
    			TwoDEff_NumSym->SetBinError(binx,biny,0);
    			TwoDEff_DenSym->SetBinContent(binx,biny,0);
    			TwoDEff_DenSym->SetBinError(binx,biny,0);
    		}
    	}
    }
    TEfficiency *TwoDEff_data_sym = new TEfficiency(*TwoDEff_NumSym,*TwoDEff_DenSym);
    BayesianErrors(TwoDEff_data);
    BayesianErrors(TwoDEff_data_sym);
    TwoDEff_data->Add(*TwoDEff_data_sym);
    TwoDEff_data->SetTitle(";sub-leading jet p_{T}, [GeV];leading jet p_{T}, [GeV]; #epsilon");
    TwoDEff_data->Draw("COLZ");
    gPad->Update();
    TwoDEff_data->GetPaintedHistogram()->GetXaxis()->SetTitleOffset(1.1);
    TwoDEff_data->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.,1.);
    TwoDEff_data->GetPaintedHistogram()->GetZaxis()->SetTitleOffset(0.8);
    style.drawStandardTitle();
	twoDplotLowM->Print("../pictures/TriggerPerformance/2DEff_Data_PF60_PF100_LowM_COLZ.pdf");
    //***************2D Efficiency in QCD*****************//
	TCanvas *twoDplotLowM_qcd = new TCanvas("twoDplotLowM_qcd","Two dimensional Efficiency");
	twoDplotLowM_qcd->SetRightMargin(0.13);
    TH2D *TwoDEff_QCD_Num = (TH2D*) BgMC->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt1vspt2");
    TH2D *TwoDEff_QCD_Den = (TH2D*) BgMC->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt1vspt2");
	TEfficiency *TwoDEff_qcd = new TEfficiency(*TwoDEff_QCD_Num,*TwoDEff_QCD_Den);
	TH2D *TwoDEff_QCD_NumSym = (TH2D*) BgMC->Get("TriggerEfficiencies/KinTrigEff_Num_2D_PF60_PF100_pt2vspt1");
	TH2D *TwoDEff_QCD_DenSym = (TH2D*) BgMC->Get("TriggerEfficiencies/KinTrigEff_Den_2D_PF60_PF100_pt2vspt1");
	    for(int binx = 1; binx <= TwoDEff_QCD_NumSym->GetNbinsX(); binx++)
	    {
	    	for(int biny = 1; biny <= TwoDEff_QCD_NumSym->GetNbinsY(); biny++)
	    	{
	    		if(binx == biny){
	    			TwoDEff_QCD_NumSym->SetBinContent(binx,biny,0);
	    			TwoDEff_QCD_NumSym->SetBinError(binx,biny,0);
	    			TwoDEff_QCD_DenSym->SetBinContent(binx,biny,0);
	    			TwoDEff_QCD_DenSym->SetBinError(binx,biny,0);
	    		}
	    	}
	    }
	    TEfficiency *TwoDEff_qcd_sym = new TEfficiency(*TwoDEff_QCD_NumSym,*TwoDEff_QCD_DenSym);
	    BayesianErrors(TwoDEff_qcd);
	    BayesianErrors(TwoDEff_qcd_sym);
	    TwoDEff_qcd->Add(*TwoDEff_qcd_sym);
	    TwoDEff_qcd->SetTitle(";sub-leading jet p_{T}, [GeV];leading jet p_{T}, [GeV]; #epsilon");
	    TwoDEff_qcd->Draw("COLZ");
	    gPad->Update();
	    TwoDEff_qcd->GetPaintedHistogram()->GetXaxis()->SetTitleOffset(1.1);
	    TwoDEff_qcd->GetPaintedHistogram()->GetZaxis()->SetTitleOffset(0.8);
    	TwoDEff_qcd->GetPaintedHistogram()->GetZaxis()->SetRangeUser(0.,1.);
	    style.drawStandardTitle();
		twoDplotLowM_qcd->Print("../pictures/TriggerPerformance/2DEff_QCD_PF60_PF100_LowM_COLZ.pdf");

		//******************Corelation****************//
		TCanvas *can_correl = new TCanvas("can_correl","Correlation in data");
		can_correl->SetRightMargin(0.13);
		int nBinsPt1Pt2 = TwoDEff_QCD_NumSym->GetNbinsY();
		TH1D *h_pt_den_slice[nBinsPt1Pt2], *h_pt_num_slice[nBinsPt1Pt2], *h_pt_den_slice_sym[nBinsPt1Pt2], *h_pt_num_slice_sym[nBinsPt1Pt2];
		TEfficiency *eff_pt_slice[nBinsPt1Pt2], *eff_pt_slice_sym[nBinsPt1Pt2];
		TH1D *h_qcd_pt_den_slice[nBinsPt1Pt2], *h_qcd_pt_num_slice[nBinsPt1Pt2], *h_qcd_pt_den_slice_sym[nBinsPt1Pt2], *h_qcd_pt_num_slice_sym[nBinsPt1Pt2];
		TEfficiency *eff_qcd_pt_slice[nBinsPt1Pt2], *eff_qcd_pt_slice_sym[nBinsPt1Pt2];
		//data for ratio
		TH2D *Correlation_PF60_LowM = new TH2D(*TwoDEff_Num);
	    for(int binx = 1; binx <= TwoDEff_Num->GetNbinsX(); binx++)
	    {
	    	//Data
    		h_pt_num_slice[binx-1] 	= (TH1D*) TwoDEff_Num->ProjectionX(("Slice_num_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		h_pt_den_slice[binx-1] 	= (TH1D*) TwoDEff_Den->ProjectionX(("Slice_den_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		eff_pt_slice[binx-1] 	= new TEfficiency(*h_pt_num_slice[binx-1],*h_pt_den_slice[binx-1]);
    		BayesianErrors(eff_pt_slice[binx-1]);
    		h_pt_num_slice_sym[binx-1] 	= (TH1D*) TwoDEff_NumSym->ProjectionX(("Slice_num_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		h_pt_den_slice_sym[binx-1] 	= (TH1D*) TwoDEff_DenSym->ProjectionX(("Slice_den_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		eff_pt_slice_sym[binx-1] 	= new TEfficiency(*h_pt_num_slice_sym[binx-1],*h_pt_den_slice_sym[binx-1]);
    		BayesianErrors(eff_pt_slice_sym[binx-1]);
    		eff_pt_slice[binx-1]->Add(*eff_pt_slice_sym[binx-1]);

    		//MC
    		h_qcd_pt_num_slice[binx-1] 	= (TH1D*) TwoDEff_QCD_Num->ProjectionX(("Slice_qcdnum_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		h_qcd_pt_den_slice[binx-1] 	= (TH1D*) TwoDEff_QCD_Den->ProjectionX(("Slice_qcdden_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		eff_qcd_pt_slice[binx-1] 	= new TEfficiency(*h_qcd_pt_num_slice[binx-1],*h_qcd_pt_den_slice[binx-1]);
    		BayesianErrors(eff_qcd_pt_slice[binx-1]);
    		h_qcd_pt_num_slice_sym[binx-1] 	= (TH1D*) TwoDEff_QCD_NumSym->ProjectionX(("Slice_qcdnum_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		h_qcd_pt_den_slice_sym[binx-1] 	= (TH1D*) TwoDEff_QCD_DenSym->ProjectionX(("Slice_qcdden_"+to_string_with_precision<int>(binx,1)).c_str(),binx,binx,"E");
    		eff_qcd_pt_slice_sym[binx-1] 	= new TEfficiency(*h_qcd_pt_num_slice_sym[binx-1],*h_qcd_pt_den_slice_sym[binx-1]);
    		BayesianErrors(eff_qcd_pt_slice_sym[binx-1]);
    		eff_qcd_pt_slice[binx-1]->Add(*eff_qcd_pt_slice_sym[binx-1]);
	    }

	    for(int binx = 0; binx < TwoDEff_Num->GetNbinsX(); binx++)
	    {
	    	for(int biny = 0; biny < TwoDEff_Num->GetNbinsY(); biny++)
	    	{

	    		double xx = eff_pt_slice[binx]->GetEfficiency(binx+1);
	    		double yy = eff_pt_slice[biny]->GetEfficiency(biny+1);
	    		double xy = eff_pt_slice[binx]->GetEfficiency(biny+1);

	    		std::cout<<"bins: y = "<<biny+1<<" x = "<<binx+1<<" glob. = "<<binx+1 + biny*(TwoDEff_Num->GetNbinsY() + 2)<<std::endl;
	    		cout<<"x: "<<binx+1<<" y: "<<biny+1<<" xx = "<<xx<<" yy = "<<yy<<endl;

	    		double cxy = 0;
	    		double ecxy = 0;
	    		if(xx > 0 && yy>0){
	    			cxy = xy/sqrt(xx*yy);
	    			if(xx!=yy && xy != 0){
	    				double exy = eff_pt_slice[binx]->GetEfficiencyErrorLow(biny+1) / xy;
	    				double exx = eff_pt_slice[binx]->GetEfficiencyErrorLow(binx+1) / xx;
	    				double eyy = eff_pt_slice[biny]->GetEfficiencyErrorLow(biny+1) / yy;
	    		          ecxy = sqrt( exy*exy + 0.25*exx*exx + 0.25*eyy*eyy );
	    			}
	    		}
	    		Correlation_PF60_LowM->SetBinContent(binx+1,biny+1,cxy);
	    		Correlation_PF60_LowM->SetBinError(binx+1,biny+1,ecxy);
	    	}
	    }
	    Correlation_PF60_LowM->SetTitle(";sub-leading jet p_{T}, [GeV];leading jet p_{T}, [GeV]; #epsilon");
	    Correlation_PF60_LowM->GetXaxis()->SetTitleOffset(1.1);
	    Correlation_PF60_LowM->GetYaxis()->SetTitleOffset(1.5);
	    Correlation_PF60_LowM->GetZaxis()->SetTitleOffset(0.8);
	    Correlation_PF60_LowM->Draw("COLZ");
	    style.drawStandardTitle();
	    can_correl->Print("../pictures/TriggerPerformance/Correlation_PF60_PF100_LowM_COLZ.pdf");

	    TCanvas *can_correl_qcd = new TCanvas("can_correl_qcd","Correlation in QCD");
	    can_correl_qcd->SetRightMargin(0.13);
	    TH2D *Correlation_qcd_PF60_LowM = new TH2D(*TwoDEff_QCD_Num);
	    for(int binx = 0; binx < TwoDEff_Num->GetNbinsX(); binx++)
	    {
	    	for(int biny = 0; biny < TwoDEff_Num->GetNbinsY(); biny++)
	    	{

	    		double xx = eff_qcd_pt_slice[binx]->GetEfficiency(binx+1);
	    		double yy = eff_qcd_pt_slice[biny]->GetEfficiency(biny+1);
	    		double xy = eff_qcd_pt_slice[binx]->GetEfficiency(biny+1);

	    		double cxy = 0;
	    		double ecxy = 0;
	    		if(xx > 0 && yy>0){
	    			cxy = xy/sqrt(xx*yy);
	    			if(xx!=yy && xy != 0){
	    				double exy = eff_qcd_pt_slice[binx]->GetEfficiencyErrorLow(biny+1) / xy;
	    				double exx = eff_qcd_pt_slice[binx]->GetEfficiencyErrorLow(binx+1) / xx;
	    				double eyy = eff_qcd_pt_slice[biny]->GetEfficiencyErrorLow(biny+1) / yy;
	    		          ecxy = sqrt( exy*exy + 0.25*exx*exx + 0.25*eyy*eyy );
	    			}
	    		}
	    		std::cout<<xx<<" "<<yy<<" "<<xy<<" val: "<<cxy<<" eval: "<<ecxy<<" up: "<<eff_qcd_pt_slice[binx-1]->GetEfficiencyErrorUp(biny)<<" low: "<<eff_qcd_pt_slice[binx-1]->GetEfficiencyErrorLow(biny)<<std::endl;
	    		Correlation_qcd_PF60_LowM->SetBinContent(binx+1,biny+1,cxy);
	    		Correlation_qcd_PF60_LowM->SetBinError(binx+1,biny+1,ecxy);
	    	}
	    }
	    Correlation_qcd_PF60_LowM->SetTitle(";sub-leading jet p_{T}, [GeV];leading jet p_{T}, [GeV]; #epsilon");
	    Correlation_qcd_PF60_LowM->GetXaxis()->SetTitleOffset(1.1);
	    Correlation_qcd_PF60_LowM->GetYaxis()->SetTitleOffset(1.5);
	    Correlation_qcd_PF60_LowM->GetZaxis()->SetTitleOffset(0.8);
	    Correlation_qcd_PF60_LowM->Draw("COLZ");
	    style.drawStandardTitle();
	    can_correl_qcd->Print("../pictures/TriggerPerformance/Correlation_QCD_PF60_PF100_LowM_COLZ.pdf");
	    
	    /****************Projection of Correlations*******************/
	    	    //Slices:
	    TCanvas *canl, *canb;
	    TLegend *slice_leg = nullptr;
	    TPaveText *pt_def[TwoDEff_QCD_Num->GetNbinsY()];
	    TH1D ** slix_corel_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];
	    TH1D ** slix_qcd_corel_lowM = new TH1D*[TwoDEff_QCD_Num->GetNbinsY()];
	    TH1D ** slix_eff_corel_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];
	    std::string name;
	    for(int biny = 1; biny < TwoDEff_Num->GetNbinsY(); biny++){
	    	canb = new TCanvas("canb");
	    	pt_def[biny] = new TPaveText(0.5,0.4,0.7,0.6,"NDC");
	    	slix_qcd_corel_lowM[biny]	= Correlation_qcd_PF60_LowM->ProjectionX(("Slice_qcd_Corel_LowM_ybin_" + to_string_with_precision<int>(biny+1,1) ).c_str(),biny+1,biny+1,"e");
	    	slix_qcd_corel_lowM[biny]->SetLineColor(2);
	    	slix_qcd_corel_lowM[biny]->Draw("E");
	    	slix_corel_lowM[biny]	= Correlation_PF60_LowM->ProjectionX(("Slice_Corel_LowM_ybin_" + to_string_with_precision<int>(biny+1,1) ).c_str(),biny+1,biny+1,"e");
	    	slix_corel_lowM[biny]->Draw("E same");

	    	std::cout<<"\n"<<std::endl;
	    	for(int i=1;i<TwoDEff_Num->GetNbinsY();++i){
	    		std::cout<<"MC: "<<slix_qcd_corel_lowM[biny]->GetBinContent(i)<<" Data: "<<slix_corel_lowM[biny]->GetBinContent(i)<<std::endl;
	    	}
	    	slix_eff_corel_lowM[biny] = (TH1D*) ratio.DrawRatio(slix_corel_lowM[biny],slix_qcd_corel_lowM[biny],canb,slice_leg,nullH);
	    	name = "../pictures/TriggerPerformance/slices/" + std::string(slix_qcd_corel_lowM[biny]->GetName()) + ".pdf";

	    	pt_def[biny]->AddText(("p_{T} = "+ to_string_with_precision<double>(slix_qcd_corel_lowM[biny]->GetBinLowEdge(biny+1),3) + " - " + to_string_with_precision<double>(slix_qcd_corel_lowM[biny]->GetBinLowEdge(biny+1) + slix_qcd_corel_lowM[biny]->GetBinWidth(biny+1),3)).c_str());
	    	pt_def[biny]->SetTextSize(0.04);
	    	style.drawStandardTitle();
	    	canb->cd();
			pt_def[biny]->SetFillColor(0);
			pt_def[biny]->SetBorderSize(0);
			pt_def[biny]->Draw();
	    	canb->Print(name.c_str());
	    }
	    

/*
	    //MC
		TH2D *Correlation_QCD_PF80_LowM = new TH2D(*TwoDEff_QCD_Num);
	    for(int binx = 1; binx <= TwoDEff_QCD_Num->GetNbinsX(); binx++)
	    {
	    	for(int biny = 1; biny <= TwoDEff_QCD_Num->GetNbinsY(); biny++)
	    	{
	    		double xx = TwoDEff_QCD_Num->GetBinContent(binx,binx);
	    		double yy = TwoDEff_QCD_Num->GetBinContent(biny,biny);
	    		double xy = TwoDEff_QCD_Num->GetBinContent(binx,biny);
	    		double cxy = 0;
	    		if(xx > 0 && yy>0) std::cout<<"x: "<<binx<<" y: "<<biny<<" val: "<<xy/sqrt(xx*yy)<<std::endl;
	    		double ecxy = 0;
	    		if(xx > 0 && yy>0){
	    			cxy = xy/sqrt(xx*yy);
	    			if(xx!=yy && xy != 0){
	    		          double exy = TwoDEff_QCD_Num->GetBinError(binx,biny) / xy;
	    		          double exx = TwoDEff_QCD_Num->GetBinError(binx,binx) / xx;
	    		          double eyy = TwoDEff_QCD_Num->GetBinError(biny,biny) / yy;
	    		          ecxy = sqrt( exy*exy + 0.25*exx*exx + 0.25*eyy*eyy );
	    			}
	    		}
	    		Correlation_QCD_PF80_LowM->SetBinContent(binx,biny,cxy);
	    		Correlation_QCD_PF80_LowM->SetBinError(binx,biny,ecxy);
	    	}
	    }
	    Correlation_QCD_PF80_LowM->SetTitle("Correlation coefficient");
	    TCanvas *canCorelLowM = new TCanvas("canCorelLowM","Correlation for LowM");
	    canCorelLowM->SetLeftMargin(0.08);
	    canCorelLowM->SetRightMargin(0.1);
	    Correlation_QCD_PF80_LowM->Draw("COLZ");
	    canCorelLowM->Print("../pictures/TriggerPerformance/2DEff_QCD_Correlation_PF80_PF100_LowM_COLZ.pdf");
	    Correlation_QCD_PF80_LowM->Draw("BOX");
	    canCorelLowM->Print("../pictures/TriggerPerformance/2DEff_QCD_Correlation_PF80_PF100_LowM_BOX.pdf");
	    Correlation_QCD_PF80_LowM->Draw("CONT1");
	    canCorelLowM->Print("../pictures/TriggerPerformance/2DEff_QCD_Correlation_PF80_PF100_LowM_CONT1.pdf");
	    Correlation_QCD_PF80_LowM->Draw("CONT4");
	    canCorelLowM->Print("../pictures/TriggerPerformance/2DEff_QCD_Correlation_PF80_PF100_LowM_CONT4.pdf");
	    Correlation_QCD_PF80_LowM->Draw("E");
	    canCorelLowM->Print("../pictures/TriggerPerformance/2DEff_QCD_Correlation_PF80_PF100_LowM_E.pdf");

	    //Slices:
	    TCanvas *canl, *canb;
	    TLegend *slice_leg = nullptr;
	    TH1D ** slix_eff_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];
	    TH1D ** slix_corel_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];

	    TH1D ** slix_qcd_eff_lowM = new TH1D*[TwoDEff_QCD_Num->GetNbinsY()];
	    TH1D ** slix_qcd_corel_lowM = new TH1D*[TwoDEff_QCD_Num->GetNbinsY()];

	    TH1D ** slix_eff_ratio_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];
	    TH1D ** slix_eff_corel_lowM = new TH1D*[TwoDEff_Num->GetNbinsY()];
	    std::string name;

	    for(int biny = 0; biny < TwoDEff_Num->GetNbinsY(); biny++){
	    	canl = new TCanvas("canl");
	    	slix_qcd_eff_lowM[biny] 	= TwoDEff_QCD_Num->ProjectionX(("Slice_qcd_Eff_LowM_ybin_" + std::to_string(biny) ).c_str(),biny+1,biny+1,"e");
	    	slix_qcd_eff_lowM[biny]->SetLineColor(2);
	    	slix_qcd_eff_lowM[biny]->Draw("E");
	    	slix_eff_lowM[biny] 	= TwoDEff_Num->ProjectionX(("Slice_Eff_LowM_ybin_" + std::to_string(biny) ).c_str(),biny+1,biny+1,"e");
	    	slix_eff_lowM[biny]->Draw("E same");
	    	slix_eff_ratio_lowM[biny] = (TH1D*) ratio.DrawRatio(slix_eff_lowM[biny],slix_qcd_eff_lowM[biny],nullH,slice_leg,canl);
	    	name = "../pictures/TriggerPerformance/slices/" + std::string(slix_qcd_eff_lowM[biny]->GetName()) + ".pdf";
	    	canl->Print(name.c_str());
	    	canb = new TCanvas("canb");
	    	slix_qcd_corel_lowM[biny]	= Correlation_QCD_PF80_LowM->ProjectionX(("Slice_qcd_Corel_LowM_ybin_" + std::to_string(biny) ).c_str(),biny+1,biny+1,"e");
	    	slix_qcd_corel_lowM[biny]->SetLineColor(2);
	    	slix_qcd_corel_lowM[biny]->Draw("E");
	    	slix_corel_lowM[biny]	= Correlation_PF80_LowM->ProjectionX(("Slice_Corel_LowM_ybin_" + std::to_string(biny) ).c_str(),biny+1,biny+1,"e");
	    	slix_corel_lowM[biny]->Draw("E");
	    	slix_eff_corel_lowM[biny] = (TH1D*) ratio.DrawRatio(slix_corel_lowM[biny],slix_qcd_corel_lowM[biny],nullH,slice_leg,canb);
	    	name = "../pictures/TriggerPerformance/slices/" + std::string(slix_qcd_corel_lowM[biny]->GetName()) + ".pdf";
	    	canb->Print(name.c_str());
	    }
	    
/*

	TCanvas *oneDFitLowM = new TCanvas("oneDFitLowM","Factor Low M",1000,800);
	TH1D *PF100Eff = (TH1D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_1D_PF80_PF100_pt1");
	PF100Eff->SetMarkerStyle(20);
	TH1D *PF100Eff_mc = (TH1D*) fMCLowM->Get("TriggerEfficiencies/KinTrigEff_1D_PF80_PF100_pt1");
	auto ratioDataMC_lowM_PF80 = (TH1D*) ratio.DrawRatio(PF100Eff,PF100Eff_mc,nullH,leg,oneDFitLowM);

	/*
	oneDFitLowM->SaveAs("../pictures/TriggerPerformance_MC_1DEff_PF80_LowM.pdf");

	/*
	TCanvas *twoDplotLowM = new TCanvas("twoDplotLowM","Two dimensional Efficiency");
	TH2D *TwoDEff_Num = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_2D_PF80_PF100_pt1vspt2");
    TH2D *TwoDEff_NumSym = (TH2D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_2D_PF80_PF100_pt2vspt1");
    for(int binx = 1; binx <= TwoDEff_NumSym->GetNbinsX(); binx++)
    {
    	for(int biny = 1; biny <= TwoDEff_NumSym->GetNbinsX(); biny++)
    	{
    		if(binx == biny){
    			TwoDEff_NumSym->SetBinContent(binx,biny,0);
    			TwoDEff_NumSym->SetBinError(binx,biny,0);
    		}
    	}
    }
    TwoDEff_Num->Add(TwoDEff_NumSym);

	TwoDEff_Num ->SetMarkerStyle(20);
	TwoDEff_Num ->SetTitle("2D Efficiency;sub-leading jet p_{T};leading jet p_{T}; #epsilon");
    gPad->SetRightMargin(0.15);
    TwoDEff_Num->GetZaxis()->SetTitleOffset(0.8);
//    gStyle->SetPadLeftMargin(0.85);
    TwoDEff_Num->GetXaxis()->SetTitleOffset(1.1);
	TwoDEff_Num->Draw("E");

	TCanvas *canSave_low;
	TwoDEff_Num->Draw("COLZ");
	TwoDEff_Num->SaveAs("../pictures/TriggerPerformance_2DEff_PF80_PF100_LowM.pdf");
//	TwoDEff_Num->Write();

*/
	auto dEtaCan = new TCanvas("dEtaCan","dEta Trig Eff",1000,800);
	style.standardTitle(PRIVATE);

	TH1D *dEtaEff_Num = (TH1D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Num_1D_PF60_PF100_dEta");
	TH1D *dEtaEff_Den = (TH1D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_Den_1D_PF60_PF100_dEta");
        TEfficiency *dEtaEff = new TEfficiency(*dEtaEff_Num,*dEtaEff_Den);
//	TH1D *dEtaEff = (TH1D*) fDataLowM->Get("TriggerEfficiencies/KinTrigEff_1D_PF60_PF100_dEta");
	dEtaEff->Draw("AP");
	dEtaEff->SetMarkerStyle(20);

//	TF1 *fitEta = new TF1("fitEta",finaleFunction,0,2.5,4);
//	fitEta-> SetParameters(-1.15698e+02,1.59862e+00,9.93479e-02,-4.49790e+01);
//	auto ratiodEta = (TH1D*) ratio.DrawRatio(dEtaEff,fitEta,"Fit:",dEtaCan);
//	ratiodEta->Draw();
	dEtaCan->SaveAs("../pictures/TriggerPerformance/dEtaEff_lowM.pdf");
/*


	//............................	HIGH M ....................................

	TCanvas *oneDFitHighM = new TCanvas("oneDFitHighM","Factor Low M",1000,800);
	TH1D *PF160Eff = (TH1D*) fDataHighM->Get("TriggerEfficiencies/KinTrigEff_1D_PF80_PF160_pt1");
	PF160Eff->SetMarkerStyle(20);
	PF160Eff->Draw("E1");
	TF1 *fit1DHighM = new TF1("fit1DHighM",finaleFunction,160,500,4);
	fit1DHighM->SetParameters(1.16893e-01,1.05584e+02,1.16893e-01,1.05584e+02);
	auto hHighMRatio = (TH1D*) ratio.DrawRatio(PF160Eff,fit1DHighM,"Sigmoid*Sigmoid",oneDFitHighM);

	oneDFitHighM->SaveAs("../pictures/TriggerPerformance_1DEff_PF80_HighM.pdf");

	TCanvas *twoDplotHighM = new TCanvas("twoDplotHighM","Two dimensional Efficiency");

	TH2D *TwoDEffPF160_Num = (TH2D*) fDataHighM->Get("TriggerEfficiencies/KinTrigEff_2D_PF80_PF160_pt1vspt2");
    TH2D *TwoDEffPF160_NumSym = (TH2D*) fDataHighM->Get("TriggerEfficiencies/KinTrigEff_2D_PF80_PF160_pt2vspt1");
    for(int binx = 1; binx <= TwoDEffPF160_NumSym->GetNbinsX(); binx++)
    {
    	for(int biny = 1; biny <= TwoDEffPF160_NumSym->GetNbinsX(); biny++)
    	{
    		if(binx == biny){
    			TwoDEffPF160_NumSym->SetBinContent(binx,biny,0);
    			TwoDEffPF160_NumSym->SetBinError(binx,biny,0);
    		}
    	}
    }
    TwoDEffPF160_Num->Add(TwoDEffPF160_NumSym);

	TwoDEffPF160_Num ->SetMarkerStyle(20);
	TwoDEffPF160_Num ->SetTitle("2D Efficiency;sub-leading jet p_{T};leading jet p_{T}; #epsilon");
    gPad->SetRightMargin(0.15);
    TwoDEffPF160_Num->GetZaxis()->SetTitleOffset(0.8);
//    gStyle->SetPadLeftMargin(0.85);
    TwoDEffPF160_Num->GetXaxis()->SetTitleOffset(1.1);
	TwoDEffPF160_Num->Draw("E");

	TCanvas *canSave_high;
	TwoDEffPF160_Num->Draw("COLZ");
	TwoDEffPF160_Num->SaveAs("../pictures/TriggerPerformance_2DEff_PF80_PF160_HighM.pdf");
//	TwoDEffPF160_Num->Write();

	/*
	TFile *outputHM = new TFile("TwoDPtHighMassEfficiency.root","recreate");

	TCanvas *can_high_M = new TCanvas("can_high_M","Two D Pt trig Eff",1000,800);

	const float Bins_H_mass[] = {120.,140.,150,160,170,180,190,200,220,240,260,300,350,400,500};
	size = sizeof(Bins_H_mass)/sizeof(float) - 1;

	TH2F *TwoDEff_H_mass_Num = new TH2F("TwoDEff_H_mass_Num","2D Efficiency",size,Bins_H_mass,size,Bins_H_mass);
	TwoDEff_H_mass_Num ->SetMarkerStyle(20);
	TH2F *TwoDEff_H_mass_Denum = new TH2F("TwoDEff_H_mass_Denum","2D Efficiency",size,Bins_H_mass,size,Bins_H_mass);
	TwoDEff_H_mass_Num ->SetTitle("2D Efficiency;sub-leading jet p_{T};leading jet p_{T}; #epsilon");
	data ->Draw("LeadPt[0]:LeadPt[1]>>TwoDEff_H_mass_Num","LeadPt[0] > 140 && LeadPt[1] > 140 && LeadMatch160[0] == 1 && LeadMatch160[1] == 1 && doubleJetTopology == 1 && dPhiFS>2.7");
	data ->Draw("LeadPt[0]:LeadPt[1]>>TwoDEff_H_mass_Denum","LeadPt[0] > 140 && LeadPt[1] > 140 && PFJet80 == 1 && doubleJetTopology == 1 && dPhiFS>2.7");
	TwoDEff_H_mass_Num ->Divide(TwoDEff_H_mass_Num,TwoDEff_H_mass_Denum,1,1,error_mode.c_str());
	TwoDEff_H_mass_Num ->Draw("E");

    TH2F *TwoDEff_H_mass_NumSym = new TH2F("TwoDEff_H_mass_NumSym","2D Efficiency with True method",size,Bins_H_mass,size,Bins_H_mass);
    TwoDEff_H_mass_NumSym ->SetMarkerStyle(20);
    TH2F *TwoDEff_H_massTrueMC_DenumSym = new TH2F("TwoDEff_H_massTrueMC_DenumSym","2D Efficiency with True method",size,Bins_H_mass,size,Bins_H_mass);
    TwoDEff_H_mass_NumSym ->SetTitle("2D Efficiency with True method;sub-leading jet p_{T};leading jet p_{T}");
    data ->Draw("LeadPt[1]:LeadPt[0]>>TwoDEff_H_mass_NumSym","LeadPt[0] > 140 && LeadPt[1] > 140 && LeadMatch160[0] == 1 && LeadMatch160[1] == 1 && doubleJetTopology == 1 && dPhiFS>2.7");
    data ->Draw("LeadPt[1]:LeadPt[0]>>TwoDEff_H_massTrueMC_DenumSym","LeadPt[0] > 140 && LeadPt[1] > 140 && PFJet80 == 1 && doubleJetTopology == 1 && dPhiFS>2.7");
    TwoDEff_H_mass_NumSym ->Divide(TwoDEff_H_mass_NumSym,TwoDEff_H_massTrueMC_DenumSym,1,1,error_mode.c_str());

    for(int binx = 1; binx <= size; binx++)
    {
    	for(int biny = 1; biny <= size; biny++)
    	{
    		if(binx == biny){
    			std::cout<<"bin = "<<binx<<std::endl;
    			TwoDEff_H_mass_NumSym->SetBinContent(binx,biny,0);
				TwoDEff_H_mass_NumSym->SetBinError(binx,biny,0);
    		}
    	}
    }

    TwoDEff_H_mass_Num->SetTitle("Work in progress, 2.2 fb^{-1}, #sqrt{s} = 13 TeV; p^{(1)}_{T}, (GeV); p^{(2)}_{T}, (GeV)");
    TwoDEff_H_mass_Num->Add(TwoDEff_H_mass_NumSym);
    gPad->SetRightMargin(0.15);
    TwoDEff_H_mass_Num->GetZaxis()->SetTitleOffset(0.8);
//    gStyle->SetPadLeftMargin(0.85);
    TwoDEff_H_mass_Num->GetXaxis()->SetTitleOffset(1.1);
	TwoDEff_H_mass_Num->Draw("E");

	TwoDEff_H_mass_Num->Draw("COLZ");
	TwoDEff_H_mass_Num->SaveAs("pictures/HighMTwoDPtTriggerEff.pdf");
//	TwoDEff_H_mass_Num->Write("TwoDEff_Num");


	TCanvas *oneDFitHighM = new TCanvas("oneDFitHighM","Factor High M",1000,800);
	TH1D *PF160Eff = new TH1D("PF160Eff",";p_{T}, [GeV]; #epsilon",100,100.,500.);
	PF160Eff->SetMarkerStyle(20);
	TH1D *PF160EffDenum = new TH1D("PF160EffDenum","",100,100.,500.);
	data->Draw("LeadPt[0]>>PF160Eff","LeadMatch160[0] == 1");
	data->Draw("LeadPt[0]>>PF160EffDenum","PFJet80 == 1");
	PF160Eff->Divide(PF160Eff,PF160EffDenum,1,1,error_mode.c_str());
	PF160Eff->Draw("E");

	fit->SetRange(160,500);
	fit->SetParameters(1.16893e-01,1.05584e+02,1.16893e-01,1.05584e+02);
//	PF160Eff->Fit(fit);
//	PF160Eff->Draw("E");
//	drawRatio(PF160Eff,fit,oneDFitLowM,"bla");
	auto hHighMRatio = (TH1D*) ratio.DrawRatio(PF160Eff,fit,"Sigmoid*Sigmoid",oneDFitHighM);
	oneDFitLowM->SaveAs("pictures/OneDEfficiencyHighM.pdf");

//	outputHM->Close();
//	output->Close();


/**/
	return 0;


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

double dEtaFitFunction(double *x, double *par)
{
	double sigmoid = 1./(1.+exp( - par[0] * (x[0] - par[1]) ));
	return sigmoid;
}

void drawRatio(TH1D *histo,TF1 *fit, TCanvas *canva, std::string fitName)
{

   //TCanvas *canva = new TCanvas("canva","",1000,800);
   TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
   pad1 -> SetBottomMargin(0.05);
   pad1 -> Draw();
   pad1 -> cd();

//   histo -> SetTitle("HLT_PFJet80 100 GeV Pt treshold;; #epsilon");
//   histo -> SetMarkerColor(1);
//   histo -> SetMarkerStyle(20);
//   histo -> Draw("ep");

   TH1D * hRatio = (TH1D*)histo->Clone("hRatio");
   histo -> Fit(fit,"R+");
   char name[200];
   sprintf(name,"}{#chi^{2}/ndf = %.2f}",fit->GetChisquare()/fit->GetNDF());
   fitName = "#splitline{ " + fitName + name;
   TLegend *fitLegend = new TLegend(0.4,0.45,0.9,0.65);
   fitLegend -> SetTextSize(0.04);
   fitLegend -> AddEntry(fit,fitName.c_str(),"l");
   histo -> Draw();
   fitLegend->Draw();
   canva -> cd();

   TPad *pad2 = new TPad("pad2","pad2",0,0.0,1,0.3);
   pad2 -> SetTopMargin(0);
   pad2 -> SetBottomMargin(0.2);
   pad2 -> Draw();
   pad2 -> cd();

   hRatio -> Sumw2();

   hRatio -> Divide(fit);
   hRatio -> Draw();

   TLine *horizLine = new TLine(120,1,500,1);
   horizLine -> SetLineStyle(2);
   horizLine -> Draw();


//   cc2->cd();
//   cc2->Update();
   hRatio -> SetTitle("");
   hRatio -> GetYaxis() -> SetRangeUser(0.6,1.4);
   hRatio -> GetYaxis() -> SetTitle("Ratio histo/fit ");
   hRatio -> GetYaxis() -> SetNdivisions(505);
   hRatio -> GetYaxis() -> SetTitleSize(20);
   hRatio -> GetYaxis() -> SetTitleFont(43);
   hRatio -> GetYaxis() -> SetTitleOffset(1.55);
   hRatio -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
   hRatio -> GetYaxis() -> SetLabelSize(20);

   // X axis ratio plot settings
   hRatio -> GetXaxis() -> SetTitle(histo->GetXaxis() -> GetTitle());
   hRatio -> GetXaxis() -> SetTitleSize(20);
   hRatio -> GetXaxis() -> SetTitleFont(43);
   hRatio -> GetXaxis() -> SetTitleOffset(3.);
   hRatio -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
   hRatio -> GetXaxis() -> SetLabelSize(20);

//   std::cout<<" Chi2/ndf: = "<<fit1D -> GetChisquare()/fit1D -> GetNDF()<<" sigmoid: "<<fit2D -> GetChisquare() / fit2D -> GetNDF()<<" Finale: "<<fit3D -> GetChisquare() / fit3D -> GetNDF()<<std::endl;

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
	}
	else {
		TH2D *num2D = (TH2D*) f->Get(name_num.c_str());
		TH2D *den2D = (TH2D*) f->Get(name_den.c_str());
		out = new TEfficiency(*num2D,*den2D);
	}
	BayesianErrors(out);

	return out;
}

