/*
 * RooFitQuality.cpp
 *
 *  Created on: 20 Apr 2017
 *      Author: shevchen
 */

#include "Analysis/BackgroundModel/interface/RooFitQuality.h"

namespace analysis {
namespace backgroundmodel {

void RooFitQuality::PrintParametersInfo(const RooFitResult& res){
	/*
	 * Method to print information about pre-/post-fit parameters.
	 */
	  std::cout << "\nconstant parameters:" << std::endl;
	  res.constPars().Print("v");
	  std::cout << "\nfloating parameters (init):" << std::endl;
	  res.floatParsInit().Print("v");
	  std::cout << "\nfloating parameters (final):" << std::endl;
	  res.floatParsFinal().Print("v");
}

Chi2Ndf RooFitQuality::chiSquare(RooPlot& frame,const std::string& curvename, const std::string& histname,const int& nFitParam, const double& blind_lowEdge, const double& blind_highEdge, const double& chi2_lowEdge, const double& chi2_highEdge){
	/*
	 * Method to get a chi^2
	 */
	auto &curve = *GetFromRooPlot<RooCurve>(frame, curvename);
	auto  &hist  = *GetFromRooPlot<RooHist> (frame, histname);

	//Curve range
	auto curve_xstart = CurveRange(curve).first;
	auto curve_xstop = CurveRange(curve).second;

	//Calculate chi2
	double x = 0,y = 0, chi2 = 0;
	int nbins = 0;
	for(int i = 0; i < hist.GetN(); ++i){
		//point-by-point loop

		//Retrieve histogram content
		hist.GetPoint(i,x,y);
		//Check whether point is in the range of a curve
		if(x < curve_xstart || x > curve_xstop) continue;
		//Check whether point is in a blinded region
		if (x > blind_lowEdge && x < blind_highEdge) continue;
		if (x < chi2_lowEdge || x > chi2_highEdge) continue;
//		std::cout<<"WTF: x = "<<x<<" "<<chi2_lowEdge<<" "<<chi2_highEdge<<std::endl;

		auto eyl = hist.GetEYlow()[i] ;
		auto eyh = hist.GetEYhigh()[i] ;
		auto exl = hist.GetEXlow()[i] ;
		auto exh = hist.GetEXhigh()[i] ;
		auto curve_yavg = curve.average(x-exl,x+exh);

		if(y!=0){
			double pull = (y > curve_yavg) ? ((y - curve_yavg)/eyl) : ((y - curve_yavg)/eyh);
			chi2 += pull*pull;
			++nbins;
		}
	}
	return Chi2Ndf(chi2,nbins-nFitParam);
}


std::pair<double,double> RooFitQuality::CurveRange(const RooCurve &curve){
	/*
	 * Method to retrieve curve range
	 */
	double curve_xstart, curve_xstop;
	int npoints = curve.GetN();
	curve_xstart = curve.GetX()[0];
	curve_xstop  = curve.GetX()[npoints-1];
	return std::make_pair(curve_xstart,curve_xstop);
}

} /* namespace backgroundmodel */
} /* namespace analysis */
