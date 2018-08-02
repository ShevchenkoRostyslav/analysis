/*
 * Weights.cpp
 *
 *  Created on: Feb 25, 2016
 *      Author: shevchen
 */

#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include "TMath.h"

#include "cmath"

#include "Analysis/MssmHbb/interface/Weights.h"

double sigmoid(const double &x, double *par){
	return 1./(1.+exp( -par[0] * (x - par[1]) ));
}

double func_erf(const double &x, double *par){
	return 1 / 2. * (1 + par[0]*TMath::Erf((x-par[1])/(sqrt(2)*par[2])));
}

double Weights::FactorizationPtWeight(const double & pt1, const double & pt2){

	if(pt1 <= 0 || pt2<=0){
		std::cerr<<"Error: in Weights::FactorizationPtWeight wrong pt (<0)"<<std::endl;
		exit(0);
	}


	double weight = 0;

	/*
	// NUMBERS FOR 74X
	if(lowM_){
		weight = 1./(1.+exp( -8.87968e-02 * (pt1 - 1.03879e+02) )) * 1./ (1.+exp( -1.88620e-01 * (pt1 - 1.02891e+02) ));
		weight *= 1./(1.+exp( -8.87968e-02 * (pt2 - 1.03879e+02) )) * 1./ (1.+exp( -1.88620e-01 * (pt2 - 1.02891e+02) ));
	}
	else {
		weight = 1./(1.+exp( -1.60938e-01 * (pt1 - 1.64206e+02) )) * 1./ (1.+exp( -1.53883e-01 * (pt1 + 3.50936e+01) ));
		weight *= 1./(1.+exp( -1.60938e-01 * (pt2 - 1.64206e+02) )) * 1./ (1.+exp( -1.53883e-01 * (pt2 + 3.50936e+01) ));
	}
	*/
	if(lowM_){
		weight = 1./(1.+exp( -1.13097e-01 * (pt1 - 1.05075e+02) )) * 1./ (1.+exp( -1.13097e-01 * (pt1 - 1.05075e+02) ));
		weight *= 1./(1.+exp( -1.13097e-01 * (pt2 - 1.05075e+02) )) * 1./ (1.+exp( -1.13097e-01 * (pt2 - 1.05075e+02) ));
	}
	else {
		weight = 1./(1.+exp( -1.46385e-01 * (pt1 - 1.66919e+02) )) * 1./ (1.+exp( -5.99105e-02 * (pt1 + 1.29321e+02) ));
		weight *= 1./(1.+exp( -1.46385e-01 * (pt2 - 1.66919e+02) )) * 1./ (1.+exp( -5.99105e-02 * (pt2 + 1.29321e+02) ));
	}

	return weight;
}

double Weights::triggerCorrectionFunction(const double & pt, const double & eta){

	if(pt <= 0 || std::abs(eta)>100){
		std::cerr<<"Error: in Weights::triggerCorrectionFunction wrong pt (<0) or |eta|>100"<<std::endl;
		exit(0);
	}

	double par_data[4];
	double par_mc[7];
	double data;
	double mc;
	//Low M case
	if(lowM_){
		if(std::abs(eta) >= 0 && std::abs(eta) <=0.6){
			par_data[0] = 4.86604e-02; par_data[1] = 6.25882e+01; par_data[2] = 1.40447e-01; par_data[3] = 1.03946e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 4.01543e+01; par_mc[1] = -1.08400e+02; par_mc[2] = 1.43523e-01; par_mc[3] = 1.03724e+02; par_mc[4] = 9.99452e-01; par_mc[5] = 5.93944e+01; par_mc[6] = 4.11896e+01;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else if(std::abs(eta) >0.6 && std::abs(eta) <=1.7){
			par_data[0] = 5.87019e-02; par_data[1] = 8.25334e+01; par_data[2] = 1.49105e-01; par_data[3] = 1.04546e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 4.01543e+01; par_mc[1] = -1.08400e+02; par_mc[2] = 1.31936e-01; par_mc[3] = 1.04404e+02; par_mc[4] = 9.99052e-01; par_mc[5] = 5.71751e+01; par_mc[6] = 4.32183e+01;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else if(std::abs(eta) > 1.7 && std::abs(eta) < 2.2){
			par_data[0] = 7.67114e-02; par_data[1] = 1.00578e+02; par_data[2] = 1.84631e-01; par_data[3] = 1.07048e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 6.15856e+01; par_mc[1] = -1.31705e+02; par_mc[2] = 1.29581e-01; par_mc[3] = 1.06427e+02; par_mc[4] = 9.99821e-01; par_mc[5] = -3.76594e+00; par_mc[6] = 6.44400e+01;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else{
			std::cerr<<"Error: WRONG ETA in Weights::triggerCorrectionFunction"<<std::endl;
			exit(-1);
		}
	}
	//High Mass
	else{
		if(std::abs(eta) >= 0 && std::abs(eta) <=0.6){
			par_data[0] = 0.008951; par_data[1] = -427.9; par_data[2] = 0.1461; par_data[3] = 158.3;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 99.51; par_mc[1] = -385.; par_mc[2] = 0.1891; par_mc[3] = 159.1; par_mc[4] = 0.9998; par_mc[5] = -115.6; par_mc[6] = 119.1;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else if(std::abs(eta) >0.6 && std::abs(eta) <=1.7){
			par_data[0] = 0.01579; par_data[1] = -117.1; par_data[2] = 0.1415; par_data[3] = 161.5;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 122.2; par_mc[1] = -324.3; par_mc[2] = 0.1563; par_mc[3] = 161.7; par_mc[4] = 0.9998; par_mc[5] = -140.8; par_mc[6] = 128.6;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else if(std::abs(eta) > 1.7 && std::abs(eta) < 2.2){
			par_data[0] = 0.2366; par_data[1] = 86.47; par_data[2] = 0.1485; par_data[3] = 170.2;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
			par_mc[0] = 40.48; par_mc[1] = -103.7; par_mc[2] = 0.03285; par_mc[3] = 43.11; par_mc[4] = 0.9988; par_mc[5] = 167.4; par_mc[6] = 8.899;
			mc   = sigmoid(pt,par_mc)*sigmoid(pt,par_mc+2)*func_erf(pt,par_mc+4);
		}
		else{
			std::cerr<<"Error: WRONG ETA in Weights::triggerCorrectionFunction"<<std::endl;
			exit(-1);
		}
	}

	return data/mc;
}

double Weights::PtTriggerEfficiency(const double &pt, const double &eta){
	if(pt <= 0 || std::abs(eta)>100){
		std::cerr<<"Error: in Weights::triggerCorrectionFunction wrong pt (<0) or |eta|>100"<<std::endl;
		exit(0);
	}

	double par_data[4];
	double data;
	//Low M case
	if(lowM_){
		if(std::abs(eta) >= 0 && std::abs(eta) <=0.6){
			par_data[0] = 4.86604e-02; par_data[1] = 6.25882e+01; par_data[2] = 1.40447e-01; par_data[3] = 1.03946e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else if(std::abs(eta) >0.6 && std::abs(eta) <=1.7){
			par_data[0] = 5.87019e-02; par_data[1] = 8.25334e+01; par_data[2] = 1.49105e-01; par_data[3] = 1.04546e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else if(std::abs(eta) > 1.7 && std::abs(eta) < 2.2){
			par_data[0] = 7.67114e-02; par_data[1] = 1.00578e+02; par_data[2] = 1.84631e-01; par_data[3] = 1.07048e+02;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else{
			std::cerr<<"Error: WRONG ETA in Weights::triggerCorrectionFunction"<<std::endl;
			exit(-1);
		}
	}
	//High Mass
	else{
		if(std::abs(eta) >= 0 && std::abs(eta) <=0.6){
			par_data[0] = 0.008951; par_data[1] = -427.9; par_data[2] = 0.1461; par_data[3] = 158.3;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else if(std::abs(eta) >0.6 && std::abs(eta) <=1.7){
			par_data[0] = 0.01579; par_data[1] = -117.1; par_data[2] = 0.1415; par_data[3] = 161.5;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else if(std::abs(eta) > 1.7 && std::abs(eta) < 2.2){
			par_data[0] = 0.2366; par_data[1] = 86.47; par_data[2] = 0.1485; par_data[3] = 170.2;
			data = sigmoid(pt,par_data)*sigmoid(pt,par_data+2);
		}
		else{
			std::cerr<<"Error: WRONG ETA in Weights::triggerCorrectionFunction"<<std::endl;
			exit(-1);
		}
	}

	return data;
}

std::pair<double,double> Weights::BTagTriggerEfficiency(const double &pt, const double& lumi_data){
	if(pt <= 0){
		std::cerr<<"Error: in Weights::BTagTriggerEfficiency wrong pt (<0)"<<std::endl;
		exit(0);
	}

	//  CSVM GtoH   nPar=2
	double SF_GtoH = 0;  //Value
	double vSF_GtoH = 0; //Variance
	double eSF_GtoH = 0; //an Error
	//  CSVM BtoF   nPar=3
	double SF_BtoF = 0;  //Value
	double vSF_BtoF = 0; //Variance
	double eSF_BtoF = 0; //an Error
	// Total:
	double SF = 0;
	double eSF = 0;
	//Lumis
	double Lumi_GtoH = 15956;//16382;
	double Lumi_BtoF = 19717;//20233;
	double Lumi_tot = Lumi_GtoH + Lumi_BtoF;
//	if(Lumi_tot != lumi_data){
//		std::cerr<<"Error: in Weights::BTagTriggerEfficiency wrong LUMI!!!!!"<<std::endl;
//		exit(0);
//	}

	double pt_valid = pt;
	if(pt <= 100) pt_valid = 100;
	else if (pt >= 900) pt_valid = 900;

	if(lowM_){
		SF_GtoH  = (0.989572) + (3.723e-05) * std::pow(pt_valid,1);
		vSF_GtoH = (0.000175237) + (-1.05947e-06) * std::pow(pt_valid,1) + (2.20153e-09) * std::pow(pt_valid,2);
		eSF_GtoH = sqrt( vSF_GtoH );

		SF_BtoF  = (0.990768) + (-0.000544311) * std::pow(pt_valid,1) + (3.93919e-07) * std::pow(pt_valid,2);
		vSF_BtoF = (0.000362084) + (-4.15513e-06) * std::pow(pt_valid,1) + (1.79752e-08) * std::pow(pt_valid,2) + (-3.08284e-11) * std::pow(pt_valid,3)  + (1.8781e-14) * std::pow(pt_valid,4);
		eSF_BtoF = sqrt( vSF_BtoF );

		SF = SF_GtoH * Lumi_GtoH / (Lumi_tot) + SF_BtoF * Lumi_BtoF / (Lumi_tot);
		eSF = (eSF_GtoH * Lumi_GtoH / (Lumi_tot) + eSF_BtoF * Lumi_BtoF / (Lumi_tot) ); //
	}

	return std::make_pair(SF,eSF);
}

double Weights::TwoDPtWeight(TH2 *histo, const double & pt1, const double & pt2){

	if(histo==0){
		std::cerr<<"Error: in Weights::TwoDPtWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	if(pt1 <= 0 || pt2<=0){
		std::cerr<<"Error: in Weights::TwoDPtWeight wrong pt (<0)"<<std::endl;
		exit(0);
	}

	double weight = 0;

	if (pt1 > 500) {
		if (pt2 > 500){
			weight = histo->Interpolate(499,499);
		}
		else weight = histo->Interpolate(pt2,499);
	}
	else weight = histo->Interpolate(pt2,pt1);

	return weight;

}

double Weights::dEtaWeight(const double & dEta){
	double weight = 0;

	if(lowM_){
		weight = 1./(1.+exp( 1.41132e+02 * (dEta - 1.59991e+00) )) * 1./ (1.+exp( -2.06352e-02 * (dEta + 2.27738e+02) ));
	}
	else weight = 1;

	return weight;
}

double Weights::BTagWeight(TH1 * low_eta, TH1 * middle_eta, TH1 * high_eta, const double & pt, const double &eta){

	if(low_eta == 0 || middle_eta == 0 || high_eta == 0){
		std::cerr<<"Error: in Weights::BTagWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	if(pt <= 0){
		std::cerr<<"Error: in Weights::BTagWeight wrong pt (<0)"<<std::endl;
		exit(0);
	}

	double weight = 0;

	if(pt <1000) {
		if(eta < 0.9) weight = low_eta -> Interpolate(pt);
		else if (eta > 0.9 &&eta < 1.4) weight = middle_eta -> Interpolate(pt);
		else weight = high_eta -> Interpolate(pt);

	}
	else {
		if(eta < 0.9) weight = low_eta -> Interpolate(995.);
		else if (eta > 0.9 &&eta < 1.4) weight = middle_eta -> Interpolate(995.);
		else weight = high_eta -> Interpolate(995.);
	}

	return weight;
}

double Weights::BTagWeight(TH2 *pt_eta, const double &pt, const double & eta){

	if(pt_eta == 0){
		std::cerr<<"Error: in Weights::BTagWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	double weight = 0;

	if(pt <1000) {
		weight = pt_eta->Interpolate(pt,eta);
	}
	else {
		weight = pt_eta->Interpolate(995.,eta);
	}

	return weight;
}

double Weights::PileUpWeight( TH1& data, TH1& mc, const double & NTruePileUpMC){
//	if(data==0 || mc ==0){
//		std::cerr<<"Error: in Weights::PileUpWeight empty weight histogram"<<std::endl;
//		exit(1);
//	}
	double weight = 0;

	double data_integral = data.Integral();
	double mc_integral   = mc.Integral();

	double N_data = data.GetBinContent(data.FindBin(NTruePileUpMC)) / data_integral;
	double N_mc = mc.GetBinContent(mc.FindBin(NTruePileUpMC)) / mc_integral;
	weight = N_data/N_mc;

//	std::cout<<" True = "<<NTruePileUpMC<<std::endl;
//	std::cout<<" Data: "<<hData->FindBin(NTruePileUpMC)<<" "<<hData->GetBinContent(hData->FindBin(NTruePileUpMC))<<std::endl;
//	std::cout<<" MC:   "<<hMc->FindBin(NTruePileUpMC)<<" "<<hMc->GetBinContent(hMc->FindBin(NTruePileUpMC))<<std::endl;
//	std::cout<<" weight = "<<weight<<" "<<N_data<<" "<<N_mc<<std::endl;

	return weight;
}

double Weights::PileUpWeight(TH1 *data, const double & N_mc, const double & NTruePileUpMC){
	if(data==0){
		std::cerr<<"Error: in Weights::PileUpWeight empty weight histogram"<<std::endl;
		exit(1);
	}
	double weight = 0;
	double N_data = data->GetBinContent(data->FindBin(NTruePileUpMC));
	weight = N_data/N_mc;

	std::cout<<"PileUp = "<<N_data<<" Nmc = "<<N_mc<<" "<<NTruePileUpMC<<std::endl;

	return weight;
}

double Weights::LOEfficiencyScaleToNLO(const int& mass){
	/*
	 * Scale factors for the signal efficiency according to the LO to NLO
	 */
	double sf = 1;
	if(mass == 300) sf = 1;
	else if (mass == 350) sf = 0.98;
	else if (mass == 400) sf = 0.97;
	else if (mass == 500) sf = 0.95;
	else if (mass == 600) sf = 0.93;
	else if (mass == 700) sf = 0.91;
	else if (mass == 900) sf = 0.9;
	else if (mass == 1100) sf = 0.9;
	else if (mass == 1300) sf = 0.92;
	else sf = 1;
	return sf;
}

/*
double Weights::PileUpWeight(TH1 *data, TH1 *mc, const double & NTruePileUpMC){

	if(data==0 || mc ==0){
		std::cerr<<"Error: in Weights::PileUpWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	NVertex_mc_ = NTruePileUpMC;

	double weight = 0;

	PileUp pu;
	pu.setDataHisto(data);
	pu.setMCHisto(mc);
	weight = pu.getPUweight(NVertex_mc_);

	return weight;
}

double Weights::PileUpWeight(TH1 *data, TH1 *mc){

	if(data==0 || mc ==0){
		std::cerr<<"Error: in Weights::PileUpWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	double weight = 0;

	PileUp pu;
	pu.setDataHisto(data);
	pu.setMCHisto(mc);
	weight = pu.getPUweight(NVertex_mc_);

	return weight;
}
*/

double Weights::HtWeight(TH1 * histo, const double &Ht){

	if(histo==0){
		std::cerr<<"Error: in Weights::HtWeight empty weight histogram"<<std::endl;
		exit(1);
	}

	if(Ht == 0){
		std::cerr<<"Error: in Weights::HtWeight Ht = 0"<<std::endl;
		exit(0);
	}

	double weight = 0;

	if(Ht < 2600){
		weight = histo->Interpolate(Ht);
	}
	else {
		weight = histo->Interpolate(2550);
	}

	return weight;
}

double Weights::dRWeight(const double& dR){
	/*
	 * Experimental weight
	 */
	double weight = 1.33441e+00 + (-8.19299e-01)*dR + 4.53292e-01 * dR * dR + (-7.19293e-02) * dR * dR * dR;
	return weight;
}

double Weights::M12Weight(const double& m12){
	double weight = 1;
	if(m12 < 400) weight= 1./(1.+exp( -0.0540249 * (m12 - 188.762) ));

	return weight;
}
double Weights::Ht_bbxWeight(TH1& h, const double& Ht){
	double weight = 0;
	if(Ht < 2600){
		weight = h.Interpolate(Ht);
	}
	else {
		weight = h.Interpolate(2550);
	}
	return weight;
}
//.............Functions...................

