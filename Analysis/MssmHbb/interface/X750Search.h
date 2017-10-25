/*
 * doubleB.h
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_INTERFACE_X750SEARCH_H_
#define MSSMHBB_INTERFACE_X750SEARCH_H_

#include <iostream>
#include <memory>
#include <string>

#include "Analysis/MssmHbb/interface/TriggerEfficiency.h"
#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"


namespace analysis{
	namespace selection{
		class X750Search : public analysis::selection::TriggerEfficiency {
//			friend class MssmHbbSignal;
//			friend class ReverseB;
		public:
			X750Search(const std::string & inputFilelist,
							const double & dataLumi,
							const bool & lowM = true,
							const bool & test = true);
			virtual ~X750Search();

			virtual void applySelection();

		protected:

		};
	}
}

#endif /* MSSMHBB_INTERFACE_X750SEARCH_H_ */
