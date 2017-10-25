#ifndef Analysis_MssmHbb_MssmHbb_h
#define Analysis_MssmHbb_MssmHbb_h 1

// -*- C++ -*-
//
// Package:    Analysis/MssmHbb
// Class:      Analysis
// 
/**\class Analysis MssmHbb.cc Analysis/MssmHbb/src/MssmHbb.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
// Modified by: Rostyslav Shevchenko
//

// system include files
#include <memory>
#include <vector>
#include <string>
#include <boost/filesystem.hpp>
// 
// user include files

#include "Analysis/Tools/interface/Analysis.h"
#include "Analysis/MssmHbb/interface/BasicTree.h"
#include "Analysis/MssmHbb/interface/utilLib.h"
//Macros that provide information about the current CMSSW version
#define CMSSW std::string(PROJECT_VERSION)

//
// class declaration
//

namespace analysis {
   namespace selection {

      class MssmHbb : public analysis::tools::Analysis, public BasicTree {
         public:
            MssmHbb(const std::string & inputFilelist, const std::string & evtinfo = "MssmHbb/Events/EventInfo");
           ~MssmHbb();
           
           //Create standart output file name
           void SetupStandardOutputFile(const std::string & outputFileName);
           //Set Output Branches
           void setBranches();
           //Clean variables
           void cleanVariables();
           //BTag SF calculation
           void calculateBTagSF(const BTagCalibrationReader & reader, const BTagCalibrationReader &reader_up, const BTagCalibrationReader &reader_down);
           //Compute Flavour composition
           void calculateFlavourComposition();
           //Trigger Selection
           bool OnlineSelection(const analysis::tools::Jet &fLeadOfflineJet,const analysis::tools::Jet &sLeadOfflineJet);

           //Add trigger Object Trees
           void addTriggerObjects(const std::vector<std::string> & triggerObjectName = {"hltL1sL1DoubleJetC100","hltDoubleJetsC100","hltDoublePFJetsC100",
											"hltDoubleBTagCSV0p9","hltDoublePFJetsC100MaxDeta1p6"},
        		   	   	   	   	  				const std::string & path = "MssmHbb/Events/selectedPatTrigger/");

           //define selection type for MC
           void setLowMSelection(const bool & type);

           //Add Vertex variables to analysis
           void addVertexInfo(analysis::tools::Vertex & vert);

           // MC weights
           void calculateWeights(TH1F *btag_low_eta,TH1F * btag_mid_eta, TH1F * btag_high_eta, TH2F *pt,const double & dataLumi);
           void setBTagWeight(const double & weight);
           void setHtWeight(const double & weight);
           void setPtFactorizationWeight(const double &weight);
           void setPt2DWeight(const double & weight);
           void setdEtaWeight(const double &weight);
           void setPileUpWeight(const double & central, const double &up, const double & down);

           //returns
           void ShowMCInformation(); // Show MC information
           std::vector<std::string> getTriggerObjectNames();
           std::string getOutPutFileName();
           const char * getTriggerLogicName();
           int getTriggerResult();
           const bool & getLowMSelection();
           bool isSignalMc() const;
           double getMHat()     const;
           bool mHatCut(const double & value) const;

           //Default cuts
           const double & Pt1Cut();
           const double & Pt2Cut();
           const double & Pt3Cut();
           const double & dEtaCut();
           const double & BTag1Cut();
           const double & BTag2Cut();
           const double & BTag3Cut();

           void setSignalMc(bool signalMc);
           //void setRun2016(bool Run2016);
           int returnMassPoint() const;
                               
            // ----------member data ---------------------------
         protected:
           std::vector<std::string> triggerObjectName_;
           std::string outPutName_;
           std::string triggerLogicName_;
           bool lowMSelection_;

           //Setup cuts according to the trigger
           void SetupConstants();
               
         private:

           double maxBJetPt__ = 670.; // For BTagCalibrationStandAlone
           double pt1_cut__;
           double pt2_cut__;
           double pt3_cut__;
           double deta_cut__;
           double btag1_cut__;
           double btag2_cut__;
           double btag3_cut__;
           bool signalMC_;
           int Mpoint__;

           double factorizationPtWeight1D(const double &);
           double dEtaWeight(const double &);
           double BTagWeight(TH1F*region1,TH1F* region2,TH1F* region3, const double &pt, const double &eta);
           double twoDPtWeight(TH2F *, const double &, const double &);




      };

      //Sets
      inline void MssmHbb::setHtWeight(const double &weight) {WeightHt_ = weight;}
      inline void MssmHbb::setBTagWeight(const double & weight){ BTagWeight_ = weight;}
      inline void MssmHbb::setPtFactorizationWeight(const double &weight){ FactorizationPtWeight_ = weight; }
      inline void MssmHbb::setPt2DWeight(const double & weight) {TwoDPtWeight_ = weight;}
      inline void MssmHbb::setdEtaWeight(const double &weight){ dEtaWeight_ = weight;}
      inline void MssmHbb::setPileUpWeight(const double &central, const double &up, const double &down){
    	  WeightPileUp_[0] = central;
    	  WeightPileUp_[1] = up;
    	  WeightPileUp_[2] = down;
      }

      inline void MssmHbb::setLowMSelection(const bool & type){lowMSelection_ = type; this->SetupConstants();}
      inline void MssmHbb::setSignalMc(bool signalMc) {signalMC_ = signalMc;    }

      //Gets
      inline std::vector<std::string> MssmHbb::getTriggerObjectNames() {return triggerObjectName_;}
      inline const char * MssmHbb::getTriggerLogicName(){ return triggerLogicName_.c_str();	}
      inline std::string MssmHbb::getOutPutFileName(){ return outPutName_; }
      inline int MssmHbb::getTriggerResult(){ return this->triggerResult(triggerLogicName_.c_str());}
      inline bool MssmHbb::isSignalMc() const { return signalMC_;}
      inline double MssmHbb::getMHat() const {return mHat_;}
      inline bool MssmHbb::mHatCut(const double & value) const {
          if (!isSignalMc()) return false;
          if (mHat_ < value * returnMassPoint()) return false;
          else return true;
      }

      //Cuts
      inline const double & MssmHbb::Pt1Cut(){ return pt1_cut__;}
      inline const double & MssmHbb::Pt2Cut(){ return pt2_cut__;}
      inline const double & MssmHbb::Pt3Cut(){ return pt3_cut__;}

      inline const double & MssmHbb::dEtaCut(){ return deta_cut__;}

      inline const double & MssmHbb::BTag1Cut(){ return btag1_cut__;}
      inline const double & MssmHbb::BTag2Cut(){ return btag2_cut__;}
      inline const double & MssmHbb::BTag3Cut(){ return btag3_cut__;}

      inline const bool & MssmHbb::getLowMSelection(){ return lowMSelection_;}




   }
}

#endif  // Analysis_MssmHbb_MssmHbb_h
