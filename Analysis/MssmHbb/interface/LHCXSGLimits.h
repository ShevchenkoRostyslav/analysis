/*
 * LHCXSGLimits.h
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_LHCXSGLIMITS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_LHCXSGLIMITS_H_

#include "Analysis/MssmHbb/interface/LimitsInterpretation.h"

namespace analysis {
namespace mssmhbb {

/*
 * @author shevchen
 *
 * 18 Aug 2017
 */
 
class LHCXSGLimits : public LimitsInterpretation{
public:
	LHCXSGLimits(const bool& blindData, const std::string& boson, const std::string& benchmark_path, const double& xMin = 300, const double& xMax = 900, const double& yMin = 0, const double& yMax = 60, const bool& test = false);

	//	Method to get tanBeta limits according to theoretical Br and Sigma
	virtual void ReadCombineLimits(const std::string& file_name, const std::string& uncert = "", const bool& UP = false, const std::string& benchmark_ref_path = "", const double& tanBref = -1);
	//	Method to receive tanBeta value from Sigma x BR in MSSM interpretation
	double TanBeta(double mA, double xsection, const std::string& uncert = "", const bool& UP = false, const std::string& benchmark_ref_path = "", const double& tanBref = -1);

	//Additional plotting objects, i.e. box with model
	virtual void AddPlottingObjects(TH2F &frame, TLegend &leg, TGraph& obs, TGraph& exp, TGraphAsymmErrors& inner_band, TGraphAsymmErrors& outer_band);
	// Plotter method with adjusted x and y axis according to MSSM interpretation
//	void LimitPlotter(const HbbLimits::LimitsToCompare& differ_limits,
//			TLegend& leg,
//			const std::string& output = "",
//			const std::string& Lumi = "35.7 fb^{-1}",
//			const std::string& xtitle = "m_{A/H} [GeV]",
//			const std::string& ytitle = "tan(#beta)",
//			const bool& logY = false);

protected:
	mssm_xs_tools xs_tool_;
};

} /* namespace mssmhbb */
} /* namespace analysis */

#endif /* ANALYSIS_MSSMHBB_INTERFACE_LHCXSGLIMITS_H_ */
