/*
 * LimitsInterpretation.cpp
 *
 *  Created on: Dec 13, 2016
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/LimitsInterpretation.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"

namespace analysis {
namespace mssmhbb {

bool is_empty(std::ifstream&);

LimitsInterpretation::LimitsInterpretation(const bool& blindData, const std::string& boson, const bool& test) :
		HbbLimits(blindData,200,900,0,60,test),
		boson_(boson){
	/*
	 * Constructor with 2 arguments.
	 * Higgs boson name has to be:
	 * A or H or both
	 */
	CheckHiggsBoson();
};

LimitsInterpretation::LimitsInterpretation(const bool& blindData, const std::string& boson, const double& xMin, const double& xMax, const double& yMin, const double& yMax, const bool& test) :
		HbbLimits(blindData,xMin,xMax,yMin,yMax,test),
		boson_(boson) {
	/*
	 * Constructor with 2 arguments.
	 * Higgs boson name has to be:
	 * A or H or both
	 */
	CheckHiggsBoson();
};

void LimitsInterpretation::CheckHiggsBoson(){
	if(boson_ != "A" && boson_ != "H" && boson_ != "both"){
		std::invalid_argument("ERROR: in LimitsInterpretation::CheckHiggsBoson() - wrong name of Higgs boson was specified. Available configs: A/H/both");
	}
}

void LimitsInterpretation::LimitPlotter(TLegend legenda,
		const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){
	/*
	 * Overwriting HbbLimits plotter method with new defaults
	 */
	HbbLimits::LimitPlotter(legenda,output,Lumi,xtitle,ytitle,logY);
}

void LimitsInterpretation::LimitPlotter(const HbbLimits::LimitsToCompare& differ_limits,
		const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){
	/*
	 * Overwriting HbbLimits plotter method with new defaults
	 */
	TLegend legenda(0.65,0.17,0.92,0.44);
	HbbLimits::LimitPlotter(differ_limits,legenda,output,Lumi,xtitle,ytitle,logY);
}

void LimitsInterpretation::LimitPlotter(const std::string& output,
		const std::string& Lumi,
		const std::string& xtitle,
		const std::string& ytitle,
		const bool& logY){
	/*
	 * Overwriting HbbLimits plotter method with new defaults
	 */
	LimitsToCompare differ_limits;
	LimitsInterpretation::LimitPlotter(differ_limits,output,Lumi,xtitle,ytitle,logY);
}

bool LimitsInterpretation::isCompatibleWithH125(const double& mA, const double& tb, mssm_xs_tools& xs_banchmark, const double & delta_m) const{
	/*
	 * Function to define whether for this particular mA and tanBeta
	 * phase space is compatible with SM Higgs of mH = 125 +- delta
	 */
	double xMinss = 125 - delta_m;
	double xMax_ss = 125 + delta_m;
	double mh = xs_banchmark.mh(mA, tb);
	bool mhInRange = mh > xMinss && mh < xMax_ss;
	double mH = xs_banchmark.mH(mA, tb);
	bool mHInRange = mH > xMinss && mH < xMax_ss;
	return mhInRange || mHInRange;

}

bool LimitsInterpretation::isCompatibleWithH125(const double& mA, const double& tb, const std::string& banchmark, const double & delta_m) const{
	mssm_xs_tools xs_banchmark(banchmark.c_str(),true,0);
	return isCompatibleWithH125(mA,tb,xs_banchmark,delta_m);
}

std::vector<TGraphAsymmErrors*> LimitsInterpretation::getH125IncompatibleAreas(mssm_xs_tools& xs_banchmark, const double & delta_m) const{
	/*
	 * set of a TGraphs with tanBeta correspond. to mH=125 +- delta GeV
	 */
	  const double stepTB = 1.;
	  const double tanBMax_ = yMax_;
	  const double tanBMin_ = yMin_;
	  const size_t nStepsTB = static_cast<size_t>(tanBMax_-tanBMin_)/stepTB;

	  std::vector<Limit> borderL;   // from low tb
	  std::vector<Limit> borderH;   // from high tb
	  double mA = xMin_;
	  while( mA < xMax_ ) {         // scan mA values

	    // 1) scan from low tb
	    if( !isCompatibleWithH125(mA,tanBMin_,xs_banchmark,delta_m) ) { // assume that incompatible area starts at edge
	      if( borderL.size() == 0 ) {
	        borderL.push_back(Limit(mA,tanBMin_));
	      }
	      size_t iTB = 0;
	      for(; iTB < nStepsTB; ++iTB) {
	        const double tb = tanBMin_ + iTB*stepTB;
	        if( isCompatibleWithH125(mA,tb,xs_banchmark,delta_m) ) {
	          borderL.push_back( Limit( mA, tb-stepTB ) );
	          break;
	        }
	      }
	      if( iTB == nStepsTB ) {
	        borderL.push_back( Limit( mA, tanBMax_ ) );
	      }
	    }

	    // 2) scan from high tb
	    if( !isCompatibleWithH125(mA,tanBMax_,xs_banchmark,delta_m) ) { // assume that incompatible area starts at edge
	      if( borderH.size() == 0 ) {
	        borderH.push_back(Limit(mA,tanBMax_));
	      }
	      size_t iTB = 0;
	      for(; iTB < nStepsTB; ++iTB) {
	        const double tb = tanBMax_ - iTB*stepTB;
	        if( isCompatibleWithH125(mA,tb,xs_banchmark,delta_m) ) {
	        	borderH.push_back( Limit( mA, tb+stepTB ) );
	        	break;
	        }
	      }
	      if( iTB == nStepsTB ) {
	        borderH.push_back( Limit( mA, tanBMin_ ) );
	      }
	    }

	    if(      mA < 150 ) mA += 5.;
	    else if( mA < 300 ) mA += 5.;
	    else if( mA < 700 ) mA += 15.;
	    else                mA += 50.;

	  } // end of mA scan
	  if( borderL.size() > 0 ) {
	    borderL.push_back( Limit( xMax_, borderL.back().getMedian() ) );
	    borderL.push_back( Limit( xMax_, tanBMin_ ) );
	  }
	  if( borderH.size() > 0 ) {
	    borderH.push_back( Limit( xMax_, borderH.back().getMedian() ) );
	    borderH.push_back( Limit( xMax_, tanBMax_ ) );
	  }

	  // create TGraphs
	  std::vector<TGraphAsymmErrors*> graphs;
	  for(int direction = 0; direction < 2; ++direction) {
	    std::vector<Limit>* border = &borderL;
	    if( direction == 1 ) {
	      if( borderH.size() > 0 ) border = &borderH;
	      else break;
	    }
	    std::vector<double> xs;
	    std::vector<double> ys;
	    std::vector<double> xel;
	    std::vector<double> xeh;
	    std::vector<double> yel;
	    std::vector<double> yeh;
	    for(size_t i = 0; i < border->size(); ++i) {
	      xs.push_back( border->at(i).getX() );
	      ys.push_back( border->at(i).getMedian() );
	      xel.push_back( i>0 ? 0.5*(border->at(i).getX()-border->at(i-1).getX()) : 0. );
	      xeh.push_back( i<border->size()-1 ? 0.5*(border->at(i+1).getX()-border->at(i).getX()) : 0. );
	      yel.push_back( 2. );
	      yeh.push_back( 0. );
	    }
	    graphs.push_back( new TGraphAsymmErrors(xs.size(),
	                                            &(xs.front()),&(ys.front()),
	                                            &(xel.front()),&(xeh.front()),
	                                            &(yel.front()),&(yeh.front())) );
	  }

	  return graphs;

}

} /* namespace MssmHbb */
} /* namespace Analysis */
