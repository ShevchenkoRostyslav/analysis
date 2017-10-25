/*
 * THDMLimits.h
 *
 *  Created on: 20 авг. 2017 г.
 *      Author: rostyslav
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_THDMLIMITS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_THDMLIMITS_H_

#include "Analysis/MssmHbb/interface/LimitsInterpretation.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

//ROOT includes
#include "TCanvas.h"
#include "TH1.h"
#include "TH3.h"
#include "TROOT.h"

namespace analysis {
namespace mssmhbb {

class THDMLimits : public LimitsInterpretation{
public:

	THDMLimits(const bool& blindData, const std::string& boson, const double& xMin = 200, const double& xMax = 900, const double& yMin = 0, const double& yMax = 60, const bool& test = false) :
			LimitsInterpretation(blindData,boson,xMin,xMax,yMin,yMax,test) {};

	struct THDMScan{
		TH2D expected;
		TH2D observed;
	};
	//Overwrite ReadCombineLimits method in order to save GxBrLimits
	void ReadCombineLimits(const std::string& file_name);
	
	//	Method to get 2HDM 3D GxBR
		TH3D Get2HDM_GxBR_3D(const std::string& benchmark_path);
	//	Method to get 2HDM 2D GxBR for particular value of VAR
		TH2D Get2HDM_GxBR_2D(const TH3& GxBR, const double& var, const std::string& axis = "X");
		TH2D Get2HDM_GxBR_2D(const std::string& benchmark_path, const double& var, const std::string& axis = "X");
	//	Method to plot Brazil 2HDM limits
		const std::vector<Limit> Get2HDM_1D_Limits(const TH2& GxBR_2hdm);
		const std::vector<Limit> Get2HDM_1D_Limits(const std::string& benchmark_path, const double& var, const std::string& axis = "X");

	//	Method to get 2D limits of G_95%CL / G_pred for particular mass point
		THDMScan Get2HDMmuScan(const TH2& GxBR_2hdm, const Limit& GxBR_95CL);
	//	Method to calculate 2HDM limits
		std::vector<Limit> Get2HDM_Limits(const TH2& GxBR_2hdm, const Limit& GxBR_95CL, const double& xmin = -1, const double& xmax = 1);

		// Method to receive tanBeta value for 1D Siga x BR in 2HDMinterpretation
		double THDMTanBeta(const TH2& GxBR_2hdm, double mA, double xsection);
		
		//
		void AddPlottingObjects(TH2F &frame, TLegend &leg, TGraph& obs, TGraph& exp, TGraphAsymmErrors& inner_band, TGraphAsymmErrors& outer_band, TCanvas &can);

protected:
		void SetTHDMHistos(TFile& file,std::map<std::string,TH3D*>&);
		std::string var_axis_ {"z"};
		float var_point_ {300};
};

} /* namespace mssmhbb */
} /* namespace analysis */

#endif /* ANALYSIS_MSSMHBB_INTERFACE_THDMLIMITS_H_ */
