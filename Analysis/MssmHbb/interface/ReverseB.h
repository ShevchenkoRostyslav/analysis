/*
 * ReverseB.h
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_INTERFACE_REVERSEB_H_
#define MSSMHBB_INTERFACE_REVERSEB_H_

#include <iostream>
#include <memory>
#include <string>

#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"
#include "Analysis/MssmHbb/interface/selectionDoubleB.h"

namespace analysis{
	namespace selection{
	class ReverseB : public analysis::selection::selectionDoubleB {
	public:
		ReverseB(const std::string & inputFilelist,
				const bool & lowM = true,
					const bool & test = true);
		virtual ~ReverseB();

		//Overwrite Leading jet selection from JetAnalysisBase class
		const bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);
		//Overwrite Methods to work with histograms:
		void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > &offlineJets,const double & weight);
		//Overwrite assignWeight method:
		const double assignWeight();

	protected:

        double pt3_;
        double eta3_ = 2.2;
        double btag3_;

	private:

	};
	}
}

#endif /* MSSMHBB_INTERFACE_REVERSEB_H_ */
