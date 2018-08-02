/*
 * Weights.h
 *
 *  Created on: Feb 25, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_SRC_WEIGHTS_H_
#define MSSMHBB_SRC_WEIGHTS_H_

#include <iostream>
#include <utility>

class Weights {

public:
	Weights(const bool & selection): lowM_(selection) {};
	virtual ~Weights() { 
      //std::cout<<"I'm at ~Weights"<<std::endl;
   };

	// Factorized Pt trigger efficiency weight in 74X
	double FactorizationPtWeight(const double & pt1, const double &pt2);

	//New pT Correction function for 2015 trigger efficiency
	double triggerCorrectionFunction(const double & pt, const double & eta);

	//Systematics for pT trigger efficiency:
	double PtTriggerEfficiency(const double &pt, const double &eta);

	//Systematics for pT trigger efficiency:
	std::pair<double,double> BTagTriggerEfficiency(const double &pt, const double& lumi_data);

	// 2D Pt trigger Efficiency weight
	double TwoDPtWeight(TH2 *histo, const double &pt1, const double &pt2);

	// dEta trigger efficiency weight
	double dEtaWeight(const double & dEta);

	// Online BTag trigger efficiency
	double BTagWeight(TH1 * low_eta, TH1 * middle_eta, TH1 * high_eta, const double & pt, const double & eta);
	double BTagWeight(TH2 * pt_eta, const double & pt, const double & eta);

	// PileUp weight
	double PileUpWeight(TH1& data, TH1& mc, const double & NTruePileUpMC);
	double PileUpWeight(TH1 * data, const double & N_mc, const double &NTruePileUpMC);
//	double PileUpWeight(TH1 * data, TH1 * mc);																	// Could be used after the first call

	// Ht weight
	double HtWeight(TH1 * histo, const double & Ht);

	// Scale signal efficiency of the LO to the NLO
	double LOEfficiencyScaleToNLO(const int& mass);

	// Luminosity normalisation
	double LumiWeight(const double & data_lumi, const double & mc_lumi);

	//Experimental weights
	double dRWeight(const double& dR);
	double M12Weight(const double& m12);
	double Ht_bbxWeight(TH1& h, const double& Ht);

private:
	const bool lowM_;
	double NVertex_mc_;

};

inline double Weights::LumiWeight(const double & data_lumi, const double & mc_lumi){
	if(mc_lumi ==0 ){
		std::cerr<<"Error in Weights::LumiWeight - 0 MCLuminosity were provided. Interrupt."<<std::endl;
		exit(0);
	}
	double weight = data_lumi/mc_lumi;
	return weight;
}

#endif /* MSSMHBB_SRC_WEIGHTS_H_ */
