#include <string>
#include <iostream>
#include <vector>

#include "TFile.h" 
#include "TFileCollection.h"
#include "TChain.h"
#include "TH1.h" 

#include "Analysis/Tools/interface/Analysis.h"

using namespace std;
using namespace analysis;
using namespace analysis::tools;


// =============================================================================================   
int main(int argc, char * argv[])
{
   TH1::SetDefaultSumw2();  // proper treatment of errors when scaling histograms
   
   // Input files list
   std::string inputList = "rootFileList.txt";
   Analysis analysis(inputList);
   
   // Analysis of events
   std::cout << "This analysis has " << analysis.size() << " events" << std::endl;
      std::cout << "Printing information from first 5 entries... \n\n";
            
   for ( int i = 0 ; i < analysis.size() ; ++i )
   {
      analysis.event(i);
      
      if ( i < 5 ) 
      {
         std::cout << "++++++    ENTRY  " << i;
         std::cout << ", Run = " << analysis.run();
         std::cout << ", Event = " << analysis.event();
         std::cout << ", LumiSection = " << analysis.lumiSection();
         std::cout << ", nPileup = " << analysis.nPileup();
         std::cout << ", nTruePileup = " << analysis.nTruePileup();
         std::cout << std::endl;
         std::cout << "                  ";
         std::cout << "  Gen Weight = " << analysis.genWeight();
         std::cout << ", Gen Scale = " << analysis.genScale();
         std::cout << ", PDF id 1 = " << analysis.pdf().id.first;
         std::cout << ", PDF id 2 = " << analysis.pdf().id.second;
         std::cout << ", PDF x 1  = " << analysis.pdf().x.first;
         std::cout << ", PDF x 2  = " << analysis.pdf().x.second;
         std::cout << "\n" << std::endl;
         
      }            
      
   }
   
//    
}

