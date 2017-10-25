/*
 * JetAnalysisBase.h
 *
 *  Created on: 17 апр. 2016 г.
 *      Author: rostyslav
 */

#ifndef ANALYSIS_MSSMHBB_SRC_TRIGGEREFFICIENCY_H_
#define ANALYSIS_MSSMHBB_SRC_TRIGGEREFFICIENCY_H_

#include <iostream>
#include <memory>
#include <string>
#include <map>

#include "Analysis/MssmHbb/interface/JetAnalysisBase.h"
#include "Analysis/MssmHbb/interface/selectionDoubleB.h"

namespace analysis{
	namespace selection{
		class TriggerEfficiency : public analysis::selection::selectionDoubleB {
			public:
			TriggerEfficiency(const std::string & inputFilelist,
								const double & dataLumi,
								const bool & lowM = true,
								const bool & test = true);
				virtual ~TriggerEfficiency();

				//Overwrite Leading jet selection from JetAnalysisBase class
				const bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);
				//Overwrite Methods to work with histograms:
				void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > &offlineJets, const double & weight);
				//Overwrite assignWeight method:
				const double assignWeight();
				//Overwrite writeHistograms method
				void writeHistograms();
				//Overwrite runner
				void runAnalysis(const std::string &json, const std::string &output = "", const int &size = 100);
				//Calculate Signal Efficiency
				void signalEfficiency();
				//Overwrite naming for the output file!
				void SetupStandardOutputFile(const std::string & outputFileName);

				//Trigger matching
				//PF100
		        bool matchToPF60(const analysis::tools::Jet & jet);
		        bool matchToPF100_PF60(const analysis::tools::Jet & jet);
		        bool matchToPF80(const analysis::tools::Jet & jet);
		        bool matchToPF100_PF80(const analysis::tools::Jet & jet);
		        bool matchToPF100dEta1p6_PF60(const analysis::tools::Jet & jet1, const analysis::tools::Jet & jet2);
		        bool matchToPF100L3_PF60(const analysis::tools::Jet & jet);
		        bool matchToPF100dEta1p6_PF80(const analysis::tools::Jet & jet1, const analysis::tools::Jet & jet2);
		        bool matchToPF100L3_PF80(const analysis::tools::Jet & jet);

		        bool matchToPF160_PF60(const analysis::tools::Jet & jet);
		        bool matchToPF160_PF80(const analysis::tools::Jet & jet);

			protected:

			private:

		};
	}
}



#endif /* ANALYSIS_MSSMHBB_SRC_TRIGGEREFFICIENCY_H_ */
