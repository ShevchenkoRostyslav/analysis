/*
 * doubleB.h
 *
 *  Created on: Apr 24, 2016
 *      Author: shevchen
 */

#ifndef MSSMHBB_INTERFACE_DATAMCCOMPARISON_H_
#define MSSMHBB_INTERFACE_DATAMCCOMPARISON_H_

#include <iostream>
#include <memory>
#include <string>

#include "TLorentzVector.h"

#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"
#include "Analysis/MssmHbb/interface/selectionDoubleB.h"
#include "Analysis/Tools/interface/Jet.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
#include "Analysis/MssmHbb/interface/CutFlow.h"


namespace analysis{
	namespace selection{
		class DataMcComparison : public analysis::selection::selectionDoubleB {
//			friend class MssmHbbSignal;
//			friend class ReverseB;
		public:
			DataMcComparison(const std::string & inputFilelist,
							const double & dataLumi,
							const bool & lowM = true,
							const int & njets = 2,
							const bool & test = true);
			virtual ~DataMcComparison();

			const virtual bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);
			//Overwrite Leading jet selection from JetAnalysisBase class
			virtual void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > &offlineJets, const double & weight);
			//Overwrite assignWeight method:
			const double assignWeight();
			//Overwrite runAnalysis method
			void runAnalysis(const std::string &json, const std::string &output, const int &size);
			//Overwrite writeHistograms method
			void writeHistograms();
			void SetupStandardOutputFile(const std::string & outputFileName);
		protected:
			//Overwrite modifyJetCollection method
			std::shared_ptr<tools::Collection<tools::Jet> > modifyJetCollection(tools::Jet & jet,
									std::shared_ptr<tools::Collection<tools::Jet> > & initialJets
									);

			void histoToFill(const std::string & syst, const std::string & variation, const double & weight);

			tools::Jet jet1_{};
			tools::Jet jet2_{};
			tools::Jet jet3_{};
			TLorentzVector diJetObject_{};

		private:

			double jet_size__{};


		};
	}
}

#endif /* MSSMHBB_INTERFACE_DATAMCCOMPARISON_H_ */
