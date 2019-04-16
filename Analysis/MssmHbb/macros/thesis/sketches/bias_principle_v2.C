#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"

using namespace std;

void bias_principle_v2(){

	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.13);
	gStyle->SetPadLeftMargin(0.17);
	gStyle->SetPadRightMargin(0.02);

	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasColor(kWhite);
	gStyle->SetCanvasDefH(600); //Height of canvas
	gStyle->SetCanvasDefW(600); //Width of canvas
	gStyle->SetCanvasDefX(0);   //POsition on screen
	gStyle->SetCanvasDefY(0);

	// For the Pad:
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadColor(kWhite);
	gStyle->SetPadGridX(false);
	gStyle->SetPadGridY(false);
	gStyle->SetGridColor(0);
	gStyle->SetGridStyle(3);
	gStyle->SetGridWidth(1);

	// For the frame:
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(1);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameFillStyle(0);
	gStyle->SetFrameLineColor(1);
	gStyle->SetFrameLineStyle(1);
	gStyle->SetFrameLineWidth(1);

	gStyle->SetOptTitle(0);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleColor(1);
	gStyle->SetTitleTextColor(1);
	gStyle->SetTitleFillColor(10);
	gStyle->SetTitleFontSize(0.05);
	gStyle->SetTitleColor(1, "XYZ");
	gStyle->SetTitleFont(42, "XYZ");
	gStyle->SetTitleSize(0.06, "XYZ");
	gStyle->SetTitleXOffset(0.9);
	gStyle->SetTitleYOffset(1.25);

	TCanvas *can = new TCanvas();
	TH2F *frame = new TH2F("frame",";x [a.u.]; Entries",1,55,80,1,0,16000);
	frame->SetStats(0);
	frame->Draw();

	/*
	 * Generate the data according to the Bg-only distribution
	 */
	TF1 *datagen_bg = new TF1("datagen_bg","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
	datagen_bg->SetParameters(1,50,10);
	// generate the data according to this function
	TH1D *databg = new TH1D("databg","",200,0,100);
        databg->SetMarkerStyle(20);
        int nevents_bg = 960000;
        for(int i = 1; i <nevents_bg; ++i){
                double n = datagen_bg->GetRandom();
                databg->Fill(n);
        }
	databg->Draw("E same");
	can->Print("bias_scatch_bg.pdf");
	frame->Draw();
	

	// Data from signal only
	TF1 *datagen_sg = new TF1("datagen_sg","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
	datagen_sg->SetParameters(1,68,1.5);
	TH1D *datasg = new TH1D("datasg","",200,0,100);
	datasg->SetMarkerStyle(20);
	datasg->SetMarkerColor(kRed);
	for(int i = 1; i <4000; ++i){
                double n = datagen_sg->GetRandom();
                datasg->Fill(n);
        }
        datasg->Draw("E same");

	// S + Bg
	TH1D *data = static_cast<TH1D*>(databg->Clone("data"));
	data->Add(datasg);
	data->Draw("E same");
	can->Print("bias_scatch_sgpbg.pdf");
        frame->Draw();
	data->Draw("E same");

	TF1 *bg1 = new TF1("bg1","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
        bg1->SetLineColor(kBlue);
        bg1->SetLineStyle(2);

	TF1 *sg = new TF1("sg","[0]*exp(-0.5*((x-[1])/[2])**2) ",0,100);
        sg->SetLineColor(kRed);
        sg->SetFillStyle(3004);
        sg->SetFillColor(kRed);
        sg->SetLineStyle(kDashed);

	TF1 *spbg1 = new TF1("spbg1","[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)",00,100);
        spbg1->SetParameters(30000,50,10,0.4,68,1.5);

	//databg->Draw("E same");
	
        data->Fit(spbg1,"N");
        bg1->SetParameters(spbg1->GetParameter(0),spbg1->GetParameter(1),spbg1->GetParameter(2));
        sg->SetParameters(spbg1->GetParameter(3),spbg1->GetParameter(4),spbg1->GetParameter(5));
        spbg1->Draw("same");
        bg1->Draw("same");
        sg->Draw("same");

        //number of signal
        TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
        leg->SetBorderSize(0);
        leg->AddEntry(spbg1,"S+Bg","l");
        leg->AddEntry(bg1,"Bg","l");
        leg->AddEntry(sg,"Signal","f");
        leg->Draw("same");

	can->Print("bias_sketch1.pdf");
	/*

	TF1 *datagen = new TF1("datagen","[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)",0,100);
	datagen->SetParameters(1,50,10,0.04,68,1.5);

	TH1D *data = new TH1D("data","",100,0,100);
	data->SetMarkerStyle(20);
	int nevents = 100000;
	for(int i = 1; i <nevents; ++i){
		double n = datagen->GetRandom();
		data->Fill(n);
	}

	TF1 *sg = new TF1("sg","[0]*exp(-0.5*((x-[1])/[2])**2) ",0,100);
	sg->SetLineColor(kRed);
	sg->SetFillStyle(3004);
	sg->SetFillColor(kRed);
	sg->SetLineStyle(kDashed);
	TF1 *sg2 = static_cast<TF1*>(sg->Clone("sg2"));
	TF1 *bg1 = new TF1("bg1","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
	bg1->SetLineColor(kBlue);
	bg1->SetLineStyle(2);
	TF1 *spbg1 = new TF1("spbg1","[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)",00,100);
	spbg1->SetParameters(3000,50,10,0.04,68,1.5);

	TF1 *bg2 = new TF1("bg2","[0] +[1]*x + [2]*x*x",60,80);
	bg2->SetLineColor(kBlue);
	bg2->SetLineStyle(2);
	TF1 *spbg2 = new TF1("spbg2","[0] +[1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)",55,80);
	spbg2->SetParameters(3.97136e+04,-1.00176e+03,6.32775e+00,0.04,68,1.5);

	data->Draw("E same");
	data->Fit(spbg1,"N");
	bg1->SetParameters(spbg1->GetParameter(3),spbg1->GetParameter(4),spbg1->GetParameter(5));
	sg->SetParameters(spbg1->GetParameter(0),spbg1->GetParameter(1),spbg1->GetParameter(2));
	spbg1->Draw("same");
	bg1->Draw("same");
	sg->Draw("same");

	//number of signal
	TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
	leg->SetBorderSize(0);
	leg->AddEntry(spbg1,"S+Bg","l");
	leg->AddEntry(bg1,"Bg","l");
	leg->AddEntry(sg,"Signal","f");
	leg->Draw("same");

	TCanvas *can2 = new TCanvas();
	frame->Draw();
	data->Draw("E same");
	data->Fit(spbg2,"Nr");
	bg2->SetParameters(spbg2->GetParameter(0),spbg2->GetParameter(1),spbg2->GetParameter(2));
	sg2->SetParameters(spbg2->GetParameter(3),spbg2->GetParameter(4),spbg2->GetParameter(5));
	spbg2->Draw("same");
	bg2->Draw("same");
	sg2->Draw("same");

	leg->Clear();
	leg->AddEntry(spbg2,"S+Bg","l");
	leg->AddEntry(bg2,"Bg","l");
	leg->AddEntry(sg2,"Signal","f");
	leg->Draw("same");

	cout<<"signal with bg1: = "<<sg->Integral(60,80)<<endl;
	cout<<"signal with bg2: = "<<sg2->Integral(60,80)<<endl;

	can->Print("../../pictures/Thesis/bias_sketch1.pdf");
	can2->Print("../../pictures/Thesis/bias_sketch2.pdf");
/*
	TF1 *sg = new TF1("sg","[0]*exp(-0.5*((x-[1])/[2])**2) ",0,100);
	sg->SetParameters(100,68,1.5);
	sg->SetLineStyle(kDashed);
	TF1 *spbg1 = new TF1("spbg1","[0]*exp(-0.5*((x-[1])/[2])**2) + [3]*exp(-0.5*((x-[4])/[5])**2)",0,100);
	spbg1->SetParameters(1,50,10,0.04,68,1.5);
	TF1 *bg1  = new TF1("bg1","[0]*exp(-0.5*((x-[1])/[2])**2)",0,100);
	bg1->SetParameters(1,50,10);
	bg1->SetLineStyle(2);
	bg1->SetLineColor(kBlue);
	TF1 *bg2 = new TF1("bg2","[0] +[1]*x + [2]*x*x",60,80);
	bg2->SetLineStyle(2);
	TF1 *spbg2 = new TF1("spbg1","[0] +[1]*x + [2]*x*x + [3]*exp(-0.5*((x-[4])/[5])**2)",60,80);
	spbg2->SetParameters(3.97136e+04,-1.00176e+03,6.32775e+00,0.04,68,1.5);
	spbg2->SetLineColor(kRed);
	bg2->SetLineColor(6);

	TH1D *data = new TH1D("data","",100,0,100);
	int nevents = 100000;
	for(int i = 1; i <nevents; ++i){
		double n = spbg1->GetRandom();
		data->Fill(n);
	}


	data->SetMarkerStyle(20);
	data->Draw("E same");
	data->Fit(spbg1);
//	data->Fit(spbg2,"r");
//	data->Fit(bg1);
//	data->Fit(bg2,"r");
//	bg1->Draw("same");
//	bg2->Draw("same");
//	sg->Draw("same");

	TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
	leg->SetBorderSize(0);
	leg->AddEntry(bg1,"bg: f_1(x)","l");
	leg->AddEntry(bg2,"bg: f_2(x)","l");
	leg->AddEntry(sg,"signal","f");
	leg->Draw("same");


//	can->Print("../../pictures/Thesis/p_values.pdf");

	/*
	 * Plot a chi^2 distribution
	 * to illustrate 95% CL
	 */
	/*
	TF1 *f_chi2 = new TF1("f_chi2_1","TMath::Prob(x,1)",4,10);
	f_chi2->SetLineColor(1);
	f_chi2->GetXaxis()->SetTitle("#tilde{q}_{#mu}");
	f_chi2->GetYaxis()->SetTitle("#chi^{2}_{1}(#tilde{q}_{#mu}|#mu,#tilde{#theta}^{obs}_{#mu})");
	f_chi2->Draw();

	TF1 *f_chi2_5 = new TF1("f_chi2_5","TMath::Prob(x,1)",5.1,10);
	f_chi2_5->SetLineColor(1);
	f_chi2_5->SetFillStyle(3002);
	f_chi2_5->SetFillColor(kRed-4);
	f_chi2_5->Draw("same");

	TLine *l_5 = new TLine(f_chi2_5->GetXmin(),0,f_chi2_5->GetXmin(),f_chi2_5->Eval(f_chi2_5->GetXmin())+0.004);
	l_5->SetLineStyle(2);
	l_5->Draw("same");

	TLatex *txt_5 = new TLatex();
	txt_5->DrawLatex(f_chi2_5->GetXmin() + 0.2, l_5->GetY2() - 0.002,
	                     ("#tilde{q}_{#mu3}^{obs}"));

	TF1 *f_chi2_2 = new TF1("f_chi2_10","TMath::Prob(x,1)",6,10);
	f_chi2_2->SetLineColor(1);
	f_chi2_2->SetFillStyle(3002);
	f_chi2_2->SetFillColor(kBlue-4);
	f_chi2_2->Draw("same");

	TLine *l_2 = new TLine(f_chi2_2->GetXmin(),0,f_chi2_2->GetXmin(),f_chi2_2->Eval(f_chi2_2->GetXmin())+0.004);
	l_2->SetLineStyle(2);
	l_2->Draw("same");

	TLatex *txt_2 = new TLatex();
	txt_2->DrawLatex(f_chi2_2->GetXmin() + 0.2, l_2->GetY2() - 0.002,
	                     ("#tilde{q}_{#mu2}^{obs}"));

	TF1 *f_chi2_0p5 = new TF1("f_chi2_0p5","TMath::Prob(x,1)",8.3,10);
	f_chi2_0p5->SetLineColor(1);
	f_chi2_0p5->SetFillStyle(3002);
	f_chi2_0p5->SetFillColor(kGreen-4);
	f_chi2_0p5->Draw("same");

	TLine *l_0p5 = new TLine(f_chi2_0p5->GetXmin(),0,f_chi2_0p5->GetXmin(),f_chi2_0p5->Eval(f_chi2_0p5->GetXmin())+0.004);
	l_0p5->SetLineStyle(2);
	l_0p5->Draw("same");

	TLatex *txt_0p5 = new TLatex();
	txt_0p5->DrawLatex(f_chi2_0p5->GetXmin() + 0.2, l_0p5->GetY2() - 0.002,
	                     ("#tilde{q}_{#mu1}^{obs}"));

	TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
	leg->SetBorderSize(0);
	leg->AddEntry(f_chi2_0p5,"p_{#mu1} = 0.5% | #tilde{q}_{#mu1}^{obs}","f");
	leg->AddEntry(f_chi2_2,"p_{#mu2} = 2% | #tilde{q}_{#mu2}^{obs}","f");
	leg->AddEntry(f_chi2_5,"p_{#mu3} = 5% | #tilde{q}_{#mu3}^{obs}","f");
	leg->Draw("same");


	can->Print("../../pictures/Thesis/p_values.pdf");

	/*

	//Full gaus function
	TF1 *full_gaus = new TF1("full_gaus","TMath::Gaus(x, 5, 0.8, 1)",2,8);
	full_gaus->SetTitle("");
	full_gaus->SetLineColor(1);
	full_gaus->GetXaxis()->SetTitle("#tilde{#theta}");
	full_gaus->GetYaxis()->SetTitle("f(#tilde{#theta},#theta)");
	full_gaus->SetMaximum(0.7);
	full_gaus->Draw();

	//Right area to be hashed
	TF1 *right_gaus_part = new TF1("right_gaus_part","TMath::Gaus(x, 5, 0.8, 1)",6.2,8);
	right_gaus_part->SetLineColor(1);
	right_gaus_part->SetLineStyle(2);
	right_gaus_part->SetFillStyle(3004);
	right_gaus_part->SetFillColor(kRed-4);
	right_gaus_part->Draw("same");

	TLine *l_right = new TLine(right_gaus_part->GetXmin(),0,right_gaus_part->GetXmin(),0.4);
	l_right->SetLineStyle(2);
	l_right->Draw("same");

	TLatex *txt_falpha = new TLatex();
//	txt_falpha->SetTextColor(kBlue+2);
	txt_falpha->SetTextSize(txt_falpha->GetTextSize()*1.2);
	txt_falpha->DrawLatex(right_gaus_part->GetXmin()*1.02, 0.4,
	                     ("f_{#alpha}(#theta)"));

	//Left area to be hashed
	TF1 *left_gaus_part = new TF1("left_gaus_part","TMath::Gaus(x, 5, 0.8, 1)",2,4.2);
	left_gaus_part->SetLineColor(1);
	left_gaus_part->SetLineStyle(2);
	left_gaus_part->SetFillStyle(3005);
	left_gaus_part->SetFillColor(kBlue-4);
	left_gaus_part->Draw("same");

	TLine *l_left = new TLine(left_gaus_part->GetXmax(),0,left_gaus_part->GetXmax(),0.4);
	l_left->SetLineStyle(2);
	l_left->Draw("same");

	TLatex *txt_fbeta = new TLatex();
//	txt_fbeta->SetTextColor(kBlue+2);
	txt_fbeta->SetTextSize(txt_fbeta->GetTextSize()*1.2);
	txt_fbeta->DrawLatex(left_gaus_part->GetXmax()*0.8, 0.4,
	                     ("f_{#beta}(#theta)"));

	TLegend *leg = new TLegend(0.6,0.7,0.96,0.93);
	leg->SetBorderSize(0);
	leg->AddEntry(left_gaus_part,"#beta = P(#tilde{#theta} #leq f_{#beta}(#theta))","f");
	leg->AddEntry(right_gaus_part,"#alpha = P(#tilde{#theta} #geq f_{#alpha}(#theta))","f");
	leg->Draw("same");

	can->Print("../../pictures/Thesis/gaus_with_porbs.pdf");
	/**/

}
