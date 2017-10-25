#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/MssmHbb/interface/MssmHbb.h"

#include "Analysis/Tools/interface/Jets.h"
#include "Analysis/Tools/interface/Muons.h"
#include "Analysis/Tools/interface/METs.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;
using namespace analysis::selection;


// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   std::string inputList = "rootFileList.txt";
   MssmHbb analysis(inputList);
   
   analysis.AddPhysicsObject("Jets","MssmHbb/Events/slimmedJetsPuppi");
   analysis.AddPhysicsObject("Muons","MssmHbb/Events/slimmedMuons");
   analysis.AddPhysicsObject("METs","MssmHbb/Events/slimmedMETsPuppi"); 
   
   Jets  jets(analysis.PhysicsObject("Jets"));
   Muons muons(analysis.PhysicsObject("Muons"));
   METs  mets(analysis.PhysicsObject("METs"));
   
   // Analysis of events
   for ( int i = 0 ; i < analysis.NumberOfEvents() ; ++i )
   {
      analysis.Event(i);
      for ( int i = 0; i < jets.size(); ++i )
      {
         Jet jet = jets.at(i);
         std::cout << jet.pt() << "  " << jet.eta() << std::endl;
      }
      
   }
   
   // cross sections
   analysis.CrossSections("MssmHbb/Metadata/CrossSections");
   analysis.ListCrossSections();
   
   //generator filter
   analysis.GeneratorFilter("MssmHbb/Metadata/GeneratorFilter");
   analysis.ListGeneratorFilter();
   
      
//    
}

