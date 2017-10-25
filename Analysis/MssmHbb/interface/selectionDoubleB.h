/*
 * doubleB.h
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_INTERFACE_SELECTIONDOUBLEB_H_
#define MSSMHBB_INTERFACE_SELECTIONDOUBLEB_H_

#include <iostream>
#include <memory>
#include <string>

#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"
#include "Analysis/MssmHbb/interface/CutFlow.h"
//#include "Analysis/MssmHbb/interface/MssmHbbSignal.h"


namespace analysis{
	namespace selection{
		class selectionDoubleB : public analysis::selection::JetAnalysisBase {
//			friend class MssmHbbSignal;
//			friend class ReverseB;
		public:
			selectionDoubleB(const std::string & inputFilelist,
							const double & dataLumi,
							const bool & lowM = true,
							const bool & test = true);
			virtual ~selectionDoubleB();

		protected:

			//Overwrite Leading jet selection from JetAnalysisBase class
			virtual const bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);
			//Overwrite Methods to work with histograms:
			virtual void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > &offlineJets, const double & weight);
			//Overwrite assignWeight method:
			const double assignWeight();

		};
	}
}

#endif /* MSSMHBB_INTERFACE_SELECTIONDOUBLEB_H_ */
