/*
 * LimitsInterpretation.h
 *
 *  Created on: Dec 13, 2016
 *      Author: shevchen
 */
#ifndef MSSMHBB_INTERFACE_LIMITSINTERPRITATION_H_
#define MSSMHBB_INTERFACE_LIMITSINTERPRITATION_H_

#include <string>
#include <memory>
#include <iostream>
#include <map>

#include "TH3.h"
#include <TH2.h>
#include <TH1.h>
#include "TLegend.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"

// Root includes
#include <TTree.h>
#include <TMath.h>
#include <TH1.h>
#include <TH3F.h>

// cpp includes
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

//my includes
#include "Analysis/MssmHbb/interface/HbbLimits.h"
#include "Analysis/MssmHbb/interface/LHCXSGScenarious.h"
// MSSM tools
#include "Analysis/MssmHbb/interface/mssm_xs_tools.h"
//#include "Analysis/MssmHbb/macros/signal/mssm_xs_tools.C"

typedef std::unique_ptr<analysis::mssmhbb::Scenario> pScenario;

namespace analysis {
namespace mssmhbb {

class LimitsInterpretation : public HbbLimits {
	/*
	 * General class for the limits interpretation
	 */
public:
	LimitsInterpretation(const bool& blindData, const std::string& boson, const bool& test = false);
	LimitsInterpretation(const bool& blindData, const std::string& boson, const double& xMin = 200, const double& xMax = 900, const double& yMin = 0, const double& yMax = 60, const bool& test = false);

//	Method to set higgs boson to be used: A/H/both
	void SetHiggsBoson(const std::string& boson);
	virtual void setScenario(AvailableScenarios scenario) {scenario_ = std::move(Scenario::Create(scenario)) /* scenario_ = std::make_unique<Scenario>(Scenario::Create(scenario)) */; xMax_ = scenario_->getXMax();}
	void compareWithPrevious(const std::string& pas) {compareWithPrevious_ = pas;}
	void drawExcludedRegion(const int& deltaM) {drawExcludedRegion_ = deltaM;}

	//getters
	std::vector<Limit> getLimits() const;
	std::vector<Limit> getGxBrLimits() const;

	// Get (tanb,mA) parameter space incompatible with mh or mH = 125+/-1 GeV
	// IncompatibleAreas points have BinContent = 1, compatible points BinContent = 0
	std::vector<TGraphAsymmErrors*> getH125IncompatibleAreas(mssm_xs_tools& xs_banchmark, const double & delta_m = 1) const;

	//Function to define whether for this particular mA and tanBeta
	// phase scpae is compatible with SM Higgs of mH = 125 +- dm
	bool isCompatibleWithH125(const double& mA, const double& tb, mssm_xs_tools& xs_banchmark, const double & delta_m = 1) const;
	bool isCompatibleWithH125(const double& mA, const double& tb, const std::string& banchmark, const double & delta_m = 1) const;

	void LimitPlotter(const HbbLimits::LimitsToCompare& differ_limits,
			const std::string& output = "",
			const std::string& Lumi = "35.7 fb^{-1}",
			const std::string& xtitle = "m_{A/H} [GeV]",
			const std::string& ytitle = "tan(#beta)",
			const bool& logY = true);

	void LimitPlotter(const std::string& output = "",
				const std::string& Lumi = "35.7 fb^{-1}",
				const std::string& xtitle = "m_{A/H} [GeV]",
				const std::string& ytitle = "tan(#beta)",
				const bool& logY = true);

	void LimitPlotter(
				TLegend leg,
				const std::string& output = "",
				const std::string& Lumi = "35.7 fb^{-1}",
				const std::string& xtitle = "m_{A/H} [GeV]",
				const std::string& ytitle = "95%C.L. limit on #sigma x BR [pb]",
				const bool& logY = true);

protected:

	std::vector<Limit> GxBr_limits_;
	std::string boson_;
	pScenario scenario_;
	int drawExcludedRegion_{1};
	std::string compareWithPrevious_{""};
	void CheckHiggsBoson();

};

inline std::vector<Limit> LimitsInterpretation::getLimits() const {
	return limits_;
}

inline std::vector<Limit> LimitsInterpretation::getGxBrLimits() const {
	return GxBr_limits_;
}

inline void LimitsInterpretation::SetHiggsBoson(const std::string& boson) {
	this->boson_ = boson;
}

} /* namespace MssmHbb */
} /* namespace Analysis */

#endif /* MSSMHBB_INTERFACE_LIMITSINTERPRITATION_H_ */
