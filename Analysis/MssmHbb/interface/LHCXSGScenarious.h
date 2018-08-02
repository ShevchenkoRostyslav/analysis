/*
 * LHCXSGScenarious.h
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_

#include <memory>

#include <vector>
#include <string>
#include "TText.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include <boost/algorithm/string/join.hpp>

#include "Analysis/MssmHbb/interface/utilLib.h"

namespace analysis {
namespace mssmhbb {

/*
 * @author shevchen
 *
 * 18 Aug 2017
 */
 
 enum AvailableScenarios{MHMODP_200,LIGHT_STOP,LIGHT_STAU,HMSSM,TAU_PHOBIC,TYPE1,TYPE2,FLIPPED,LEPTON_SPECIFIC};
 std::string AvailableScenariosToString(AvailableScenarios scenario);
 AvailableScenarios AvailableScenariosFromString(const std::string& scenario_string);
 
class Scenario{
	/*
	* Abstract interface for a general scenarious
	*/
public:
	 Scenario() = default;
	virtual std::string getLabel() const = 0;
	virtual std::vector<TGraph> getPreviousResults(const std::string& PAS) const = 0;
	virtual std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const = 0;
	virtual TText getPreviousResultsLabel(const std::string& PAS) const = 0;
	virtual double getXMax() const = 0;
	virtual bool previousExists() const = 0;
	void checkResultsToCompareWith(const std::string& pas);
	
	static std::unique_ptr<Scenario> Create(AvailableScenarios scenario);
	
protected:
	/*
	 * Available results are encoded in terms of PAS/CONFNOTES/PAPERS/ARXIVs etc
	 */
	std::vector<std::string> availableResults_ = {""};

}; 

class LHCXSGScenario : public Scenario {
	/*
	* Abstract interface class to specific LHCXSG scenarios
	*/
public:
	LHCXSGScenario() = default;
	virtual std::string getLabel() const = 0;
	virtual std::vector<TGraph> getPreviousResults(const std::string& PAS) const = 0;
	virtual TText getPreviousResultsLabel(const std::string& PAS) const = 0;
	virtual std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const = 0;
	virtual double getXMax() const = 0;
	virtual bool previousExists() const = 0;
};

//LHCXSGScenario implementations

class mhmodp_200 : public LHCXSGScenario{
public:
	mhmodp_200();
	std::string getLabel() const {return "m_{h}^{mod+} scenario,  #mu = +200 GeV";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const;
	TText getPreviousResultsLabel(const std::string& PAS) const;
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class light_stop : public LHCXSGScenario{
public:
	light_stop();
	std::string getLabel() const {return "Light-#tilde{t} scenario";}	
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const;
	TText getPreviousResultsLabel(const std::string& PAS) const;
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class light_stau : public LHCXSGScenario{
public:
	light_stau();
	std::string getLabel() const {return "Light-#tilde{#tau} scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const;
	TText getPreviousResultsLabel(const std::string& PAS) const;
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class hMSSM : public LHCXSGScenario{
public:
	hMSSM();
	std::string getLabel() const {return "hMSSM scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const;
	TText getPreviousResultsLabel(const std::string& PAS) const;
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return false;}
};

class tau_phobic : public LHCXSGScenario{
public:
	tau_phobic();
	std::string getLabel() const {return "#tau-phobic scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const;
	TText getPreviousResultsLabel(const std::string& PAS) const;
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class type1 : public Scenario{
public:
	type1();
	std::string getLabel() const {return "2HDM type-I scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	TText getPreviousResultsLabel(const std::string& PAS) const {TText tx; return tx;}
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class lepton_specific : public Scenario{
public:
	lepton_specific();
	std::string getLabel() const {return "2HDM Lepton Specific scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	TText getPreviousResultsLabel(const std::string& PAS) const {TText tx; return tx;}
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class type2 : public Scenario{
public:
	type2();
	std::string getLabel() const {return "2HDM type-II scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	TText getPreviousResultsLabel(const std::string& PAS) const {TText tx; return tx;}
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class flipped : public Scenario{
public:
	flipped();
	std::string getLabel() const {return "2HDM flipped scenario";}
	std::vector<TGraph> getPreviousResults(const std::string& PAS) const {std::vector<TGraph> gr; return gr;}
	TText getPreviousResultsLabel(const std::string& PAS) const {TText tx; return tx;}
	std::vector<TGraph> getPreviousResults(const std::string& var,const std::string& PAS) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

} /* namespace mssmhbb */
} /* namespace analysis */

#endif /* ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_ */
