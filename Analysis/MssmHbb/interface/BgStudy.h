/*
 * BgStudy.h
 *
 *  Created on: Apr 23, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_INTERFACE_BGSTUDY_H_
#define MSSMHBB_INTERFACE_BGSTUDY_H_

#include <iostream>
#include <memory>
#include <string>
#include <map>

#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"
#include "Analysis/MssmHbb/interface/MssmHbbSignal.h"

namespace analysis{
	namespace selection{
		class BgStudy : public analysis::selection::MssmHbbSignal {
			public:
				BgStudy(const std::string & inputFilelist,
								const double & dataLumi,
								const bool & lowM = true,
								const bool & test = true);
				virtual ~BgStudy();

				//Overwrite Leading jet selection from JetAnalysisBase class
				const bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);
				//Overwrite Methods to work with histograms:
				void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > &offlineJets, const double & weight);
				//Overwrite assignWeight method:
				const double assignWeight();


			protected:

		        double pt3_  = 40;
		        double eta3_ = 2.2;
		        double btag3_;

			private:

		};
	}
}

#endif /* MSSMHBB_INTERFACE_BGSTUDY_H_ */
