#ifndef Analysis_Objects_Btagging_h
#define Analysis_Objects_Btagging_h 1

// -*- C++ -*-
//
// Package:    Analysis/Objects
// Class:      Analysis
// 
/**\class Analysis Btagging.cc Analysis/Objects/src/Btagging.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Roberval Walsh Bastos Rangel
//         Created:  Mon, 20 Oct 2014 14:24:08 GMT
//
//

// system include files
#include <memory>
#include <vector>
#include <string>
// 
// user include files
#include "TH1.h" 
#include "TH2.h" 

#include "Analysis/Objects/interface/Object.h"

//
// class declaration
//

namespace analysis {
   namespace objects {

      class Btagging : public analysis::objects::Object {
         public:
            Btagging(const std::string & inputFilelist, const std::string & evtinfo = "MssmHbb/Events/EventInfo");
           ~Btagging();
           
            void efficiencies();
            void efficiencies(const std::string &);
            void jetsCollection(const std::string &);
            void jetsCollection(const std::string &, const std::string &);
            void efficiencyHistograms();
            void workingPoint(const float &);
            void ptBinning(const int & , const float *);
            void etaBinning(const  int & , const float * );
            void ptMin(const float & ptmin = 20.);
            void etaMax(const float & etamax = 2.5);
            void flavourDefinition(const std::string & flavdef = "Hadron");
           
         
            // ----------member data ---------------------------
         protected:
               
         private:
            std::string jets_;
            float wp_;
            float ptmin_;
            float etamax_;
            std::string flavdef_;
            std::map<std::string,TH2F *> h2d_eff_;
            int nptbins_;
            const float * ptbins_;
            int netabins_;
            const float * etabins_;
            

      };
      
      inline void  Btagging::workingPoint(const float & wp) { wp_ = wp; }
      inline void  Btagging::ptBinning(const int & nbins,  const float * bins) { nptbins_ = nbins; ptbins_ = bins; }
      inline void  Btagging::etaBinning(const int & nbins,  const float * bins) { netabins_ = nbins; etabins_ = bins; }
      inline void  Btagging::ptMin(const float & ptmin) { ptmin_ = ptmin; }
      inline void  Btagging::etaMax(const float & etamax) { etamax_ = etamax; }
      inline void  Btagging::flavourDefinition(const std::string & flavdef) { flavdef_ = flavdef; }
   }
}

#endif  // Analysis_Objects_Btagging_h
