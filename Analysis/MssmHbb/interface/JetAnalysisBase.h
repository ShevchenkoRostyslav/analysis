/*
 * JetAnalysisBase.h
 *
 *  Created on: 17 апр. 2016 г.
 *      Author: rostyslav
 */

#ifndef ANALYSIS_MSSMHBB_SRC_JETANALYSISBASE_H_
#define ANALYSIS_MSSMHBB_SRC_JETANALYSISBASE_H_

#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"

// system include files
#include <memory>
#include <vector>
#include <string>
#include <type_traits>
#include <stdlib.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/case_conv.hpp>
#include "TSystem.h"
// PDF uncertainties
#include "LHAPDF/LHAPDF.h"
//
// user include files

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/Tools/interface/Collection.h"
#include "Analysis/Tools/interface/Jet.h"

#include "Analysis/MssmHbb/interface/Weights.h"
#include "Analysis/MssmHbb/interface/Histograms.h"
#include "Analysis/MssmHbb/interface/BTagScaleFactor.h"
#include "Analysis/MssmHbb/interface/utilLib.h"

#include "Analysis/MssmHbb/interface/CutFlow.h"


//TODO: Add template arguments to unique_ptr depending on the TH type
typedef std::unique_ptr<TFile> pTFile;

namespace analysis{
	namespace selection{
	class JetAnalysisBase : public analysis::tools::Analysis {
	public:

		//Default constructor
		JetAnalysisBase(const std::string & inputFilelist, const double & dataLumi, const bool & lowM = true, const bool & test = true);
		JetAnalysisBase(const std::string & inputFilelist, const bool & lowM = true);
		virtual ~JetAnalysisBase();

        //Add trigger Object Trees
        void addTriggerObjects(const std::vector<std::string> & triggerObjectName = {"hltL1sL1DoubleJetC100",
                   																	"hltDoubleJetsC100",
                   																	"hltDoublePFJetsC100",
                   																	"hltDoubleBTagCSV0p9",
                   																	"hltDoublePFJetsC100MaxDeta1p6"},
     		   	   	   	   	  const std::string & path = "MssmHbb/Events/selectedPatTrigger/");

        void loadCorrections();

        //Method that make an event loop
        virtual void applySelection();
        //setup all corrections, ttrees etc
        virtual void setupAnalysis(const std::string & json);
        //run everything
        virtual void runAnalysis(const std::string &json, const std::string &output = "", const int &size = 100);
        // virtual method that create histogrmas
        virtual void makeHistograms(const int &size = 100);
        //Create standart output file name
        virtual void SetupStandardOutputFile(const std::string & outputFileName = "");
        void createOutputFile(const std::string &name );

        const virtual bool OnlineSelection(const analysis::tools::Jet &fLeadOfflineJet,const analysis::tools::Jet &sLeadOfflineJet);

        //write histograms and close the output file.
        virtual void writeHistograms();

        //Return std::vector that contains names of all trigger objects
        std::vector<std::string> getTriggerObjectNames();
        TFile * getOutputFile();
        const bool & getLowM();
        const bool & isSignalMC();
        int returnMassPoint() const;

	protected:

        double dataLumi_{};
        const bool lowM_{};
        std::string triggerLogicName_{};
        int nJets_{};
        Histograms histo_{};
        const bool TEST{};
        std::map<std::string, double > weight_{};
        std::vector<std::string> triggerObjectName_{};
        bool signalMC_ = false;
        pTFile outputFile_{};
        double JESshift_{};
        double JERshift_{};
        std::string baseOutputName_{};
        CutFlow cuts_;

        //generated jetcollection:
        std::shared_ptr<tools::Collection<tools::Jet > > genJets_;

        //Protected methods
        virtual void fillHistograms(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets, const double & weight);
        //virtual void fillHistograms(const tools::Jet &LeadingJet, const tools::Jet &sub_LeadingJet);
        //method to assign total weight
        virtual const double assignWeight();
        //Modify jet Collection
        virtual std::shared_ptr<tools::Collection<tools::Jet> > modifyJetCollection(tools::Jet & jet,
        																	  std::shared_ptr<tools::Collection<tools::Jet> > & initialJets
																			  );
//        const virtual bool leadingJetSelection(const int & iJet, const tools::Jet & Jet);
        const virtual bool leadingJetSelection(const std::shared_ptr<tools::Collection<tools::Jet> > & offlineJets);

//        Method to combine BTag SFs from different jets:
        virtual void combineBTagSFs(const std::array<BTagScaleFactor::ScaleFactor,3>& sf);
        virtual void combineBTagSFs(const std::array< std::vector<BTagScaleFactor::ScaleFactor> , 3 >& sf, const std::vector<double>& lumis, const double& lumi_tot);

				// PDF uncertainties
				const void pdfUncertainties(const double& genScale, const analysis::tools::PDF& pdf);

        void Ht(const double & Ht);
        const double & Ht();

		//Selection constants
        double pt1_ = 100;
        double pt2_ = 100;
        double pt3_ =  40;
        double eta1_ = 2.2;
        double eta2_ = 2.2;
        double eta3_ = 2.2;
        double btag1_= 0.935;
        double btag2_= 0.935;
        double btag3_= 0.46;
        double dEta_ = 1.55;
        double dR_ = 1.;
        int btagOP1_ = 2; //OPerating points for BTag. 0 - Loose, 1 - Mid, 2 - Tight
        int btagOP2_ = 2;
        int btagOP3_ = 0;
        double mHat_ = 0;

        std::unique_ptr<Weights> pWeight_;
        BTagScaleFactor BTagLib_;
        std::map<std::string,TH1D *>  hCorrections1D_;
        std::map<std::string,TH2D *>  hCorrections2D_;
        //Signal xsections
        std::map<int,double> xsection_;
	const std::string defaultInputDir_ = std::string(gSystem->Getenv("CMSSW_BASE"))+"/src/Analysis/MssmHbb/";
	private:

        std::map<std::string,pTFile> fCorrections_;

        double Ht__{};

        //mHat calculation
        const double mHat();
        //Setup xsections
        void setupXSections();
				// PDF uncertainties
				std::vector<LHAPDF::PDF*> pdfs;



	};

	inline std::vector<std::string> JetAnalysisBase::getTriggerObjectNames() {return triggerObjectName_;}
	inline TFile * JetAnalysisBase::getOutputFile(){ return outputFile_.get();}
	inline const bool & JetAnalysisBase::getLowM(){ return lowM_;}
	inline const bool & JetAnalysisBase::isSignalMC(){ return signalMC_;}
	inline const double & JetAnalysisBase::Ht(){ return Ht__;}
	inline void JetAnalysisBase::Ht(const double & Ht){ Ht__ = Ht;}

	}
}



#endif /* ANALYSIS_MSSMHBB_SRC_JETANALYSISBASE_H_ */
