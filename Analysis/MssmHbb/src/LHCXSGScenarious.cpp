/*
 * LHCXSGScenarious.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/LHCXSGScenarious.h"

namespace analysis {
namespace mssmhbb {

Scenario* Scenario::Create(AvailableScenarios scenario){
	/*
	* Realisation of a factory pattern
	*/
	if(scenario == MHMODP_200) return new mhmodp_200();
	else if(scenario == LIGHT_STOP) return new light_stop();
	else if(scenario == LIGHT_STAU) return new light_stau();
	else if(scenario == HMSSM) return new hMSSM();
	else if(scenario == TAU_PHOBIC) return new tau_phobic();
	else return NULL;
}

// LHCXSGScenario* LHCXSGScenario::Create(AvailableLHCXSGScenarios scenario){
// 	/*
// 	* Realisation of a factory pattern
// 	*/
// 	if(scenario == mhmodp_m200) return new mhmodp_m200();
// 	else if(scenario == light_stop) return new light_stop();
// 	else if(scenario == light_stau) return new light_stau();
// 	else if(scenario == hMSSM) return new hMSSM();
// 	else if(scenario == tau_phobic) return new tau_phobic();
// 	else return NULL;
// }

TGraph * mhmodp_200::getPreviousResults() const{
	std::array<double,8> x = {{100,140,160,200,300,350,400,500}};
	std::array<double,8> y = {{19.8,19.6,18.0,19.2,24.8,30.0,36.5,53.5}};
	TGraph *gr = new TGraph(8,&x[0],&y[0]);
	return gr;
}

TPaveText * mhmodp_200::getPreviousResultsLabel() const{
	/*
	* Return a TPaveText box with a sign for a previous results
	* when available
	*/
	TPaveText *text = new TPaveText(300,24.8,400,36.5);
	text->AddText("7+8 TeV expected");
	text->SetTextAngle(30);
	return text;
}

} /* namespace mssmhbb */
} /* namespace analysis */
