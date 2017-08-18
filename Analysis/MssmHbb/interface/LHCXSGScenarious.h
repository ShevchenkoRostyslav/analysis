/*
 * LHCXSGScenarious.h
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_

#include <string>
#include "TPaveText.h"
#include "TGraph.h"

namespace analysis {
namespace mssmhbb {

/*
 * @author shevchen
 *
 * 18 Aug 2017
 */
 
//  enum Available2HDMScenarious{typeI,typeII,flipped,lepton_specific};
 enum AvailableScenarios{MHMODP_200,LIGHT_STOP,LIGHT_STAU,HMSSM,TAU_PHOBIC,TYPEI,TYPEII,FLIPPED,LEPTON_SPECIFIC};
//  enum AvailableScenarious{AvailableLHCXSGScenarios,Available2HDMScenarious};
 
class Scenario{
	/*
	* Abstract interface for a general scenarious
	*/
public:
// 	virtual Scenario(AvailableScenarious scenario) {return Create(scenario);}
	virtual std::string getLabel() const = 0;
	static Scenario* Create(AvailableScenarios scenario);
}; 

class LHCXSGScenario : public Scenario {
	/*
	* Abstract interface class to specific LHCXSG scenarios
	*/
public:
	virtual std::string getLabel() const = 0;
	virtual TGraph* getPreviousResults() const = 0;
	virtual TPaveText* getPreviousResultsLabel() const = 0;
// 	static LHCXSGScenario* Create(AvailableLHCXSGScenarios scenario);
};

//LHCXSGScenario implementations

class mhmodp_200 : public LHCXSGScenario{
public:
	std::string getLabel() const {return "m_{h}^{mod+} scenario,  #mu = +200";}
	TGraph* getPreviousResults() const;
	TPaveText* getPreviousResultsLabel() const;
};

class light_stop : public LHCXSGScenario{
public:
	std::string getLabel() const {return "Light-t scenario,  #mu = +200";}
	TGraph* getPreviousResults() const;
	TPaveText* getPreviousResultsLabel() const;
};

class light_stau : public LHCXSGScenario{
public:
	std::string getLabel() const {return "Light-#tau scenario,  #mu = +200";}
	TGraph* getPreviousResults() const;
	TPaveText* getPreviousResultsLabel() const;
};

class hMSSM : public LHCXSGScenario{
public:
	std::string getLabel() const {return "hMSSM scenario,  #mu = +200";}
	TGraph* getPreviousResults() const;
	TPaveText* getPreviousResultsLabel() const;
};

class tau_phobic : public LHCXSGScenario{
public:
	std::string getLabel() const {return "#tau-phobic scenario";}
	TGraph* getPreviousResults() const;
	TPaveText* getPreviousResultsLabel() const;
};

} /* namespace mssmhbb */
} /* namespace analysis */

#endif /* ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_ */
