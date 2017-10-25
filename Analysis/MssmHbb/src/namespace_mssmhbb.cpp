#include <algorithm>
#include <string>
#include <map>
#include <vector>

#include "Analysis/MssmHbb/interface/namespace_mssmhbb.h"
/*
 * Namespace that contains all basic variables
 * used in the analysis flow
 */

namespace mssmhbb{
	const std::string getModelForTheMass(const int& mass){
		auto findResult = std::find_if(std::begin(signal_models),std::end(signal_models),[&](const std::pair<std::string,std::vector<int> >&pair){
			return (std::find(pair.second.begin(),pair.second.end(),mass) != std::end(pair.second));
		});
		if(findResult != std::end(signal_models)) return findResult->first;
		else throw std::logic_error("No model for mass " + std::to_string(mass));
	}

};
