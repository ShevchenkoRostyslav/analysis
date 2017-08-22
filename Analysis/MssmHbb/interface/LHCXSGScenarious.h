/*
 * LHCXSGScenarious.h
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#ifndef ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_
#define ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_

#include <memory>

#include <string>
#include "TText.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"

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
	virtual std::string getLabel() const = 0;
	virtual TGraph getPreviousResults() const = 0;
	virtual std::vector<TGraph> getAtlasResults(const std::string& var) const = 0;
	virtual TText getPreviousResultsLabel() const = 0;
	virtual double getXMax() const = 0;
	virtual bool previousExists() const = 0;
	
	static std::unique_ptr<Scenario> Create(AvailableScenarios scenario);
	
}; 

class LHCXSGScenario : public Scenario {
	/*
	* Abstract interface class to specific LHCXSG scenarios
	*/
public:
	virtual std::string getLabel() const = 0;
	virtual TGraph getPreviousResults() const = 0;
	virtual TText getPreviousResultsLabel() const = 0;
	virtual std::vector<TGraph> getAtlasResults(const std::string& var) const = 0;
	virtual double getXMax() const = 0;
	virtual bool previousExists() const = 0;
};

//LHCXSGScenario implementations

class mhmodp_200 : public LHCXSGScenario{
public:
	std::string getLabel() const {return "m_{h}^{mod+} scenario,  #mu = +200";}
	TGraph getPreviousResults() const;
	TText getPreviousResultsLabel() const;
	std::vector<TGraph> getAtlasResults(const std::string& var) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class light_stop : public LHCXSGScenario{
public:
	std::string getLabel() const {return "Light-#tilde{t} scenario";}	
	TGraph getPreviousResults() const;
	TText getPreviousResultsLabel() const;
	std::vector<TGraph> getAtlasResults(const std::string& var) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class light_stau : public LHCXSGScenario{
public:
	std::string getLabel() const {return "Light-#tilde{#tau} scenario";}
	TGraph getPreviousResults() const;
	TText getPreviousResultsLabel() const;
	std::vector<TGraph> getAtlasResults(const std::string& var) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return true;}
};

class hMSSM : public LHCXSGScenario{
public:
	std::string getLabel() const {return "hMSSM scenario";}
	TGraph getPreviousResults() const;
	TText getPreviousResultsLabel() const;
	std::vector<TGraph> getAtlasResults(const std::string& var) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return false;}
};

class tau_phobic : public LHCXSGScenario{
public:
	std::string getLabel() const {return "#tau-phobic scenario";}
	TGraph getPreviousResults() const;
	TText getPreviousResultsLabel() const;
	std::vector<TGraph> getAtlasResults(const std::string& var) const {std::vector<TGraph> gr; return gr;}
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class type1 : public Scenario{
public:
	std::string getLabel() const {return "2HDM Type-I scenario";}
	TGraph getPreviousResults() const {TGraph gr; return gr;}
	TText getPreviousResultsLabel() const {TText tx; return tx;}
	std::vector<TGraph> getAtlasResults(const std::string& var) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class type2 : public Scenario{
public:
	std::string getLabel() const {return "2HDM Type-II scenario";}
	TGraph getPreviousResults() const {TGraph gr; return gr;}
	TText getPreviousResultsLabel() const {TText tx; return tx;}
	std::vector<TGraph> getAtlasResults(const std::string& var) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

class flipped : public Scenario{
public:
	std::string getLabel() const {return "2HDM Flipped scenario";}
	TGraph getPreviousResults() const {TGraph gr; return gr;}
	TText getPreviousResultsLabel() const {TText tx; return tx;}
	std::vector<TGraph> getAtlasResults(const std::string& var) const;
	double getXMax() const {return 900;}
	bool previousExists() const {return false;} //exists but too small mH/A
};

} /* namespace mssmhbb */
} /* namespace analysis */

#endif /* ANALYSIS_MSSMHBB_INTERFACE_LHCXSGSCENARIOUS_H_ */
