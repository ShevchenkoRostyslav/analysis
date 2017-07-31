#include <iostream>
#include <cmath>
#include <math.h>
#include "Riostream.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "Analysis/BackgroundModel/interface/RooDoubleBukin.h"

double EvaluateSingleBukin(const RooRealProxy &x ,const RooRealProxy& lambda,const RooRealProxy& sigma,const RooRealProxy& xg);

using namespace analysis::backgroundmodel;

ClassImp(RooDoubleBukin)


RooDoubleBukin::RooDoubleBukin(const char *name, const char *title,
        RooAbsReal& x,
        RooAbsReal& xg1,
		RooAbsReal& xg2,
		RooAbsReal& ksi,
		RooAbsReal& sigma1,
		RooAbsReal& lambda1,
		RooAbsReal& sigma2,
		RooAbsReal& lambda2
		) :
  RooAbsPdf(name, title),
  x_("x","x",this,x),
  xg1_("xg1","xg1",this,xg1),
  xg2_("xg2","xg2",this,xg2),
  ksi_("ksi","ksi",this,ksi),
  sigma1_("sigma1","sigma1",this,sigma1),
  sigma2_("sigma2","sigma2",this,sigma2),
  lambda1_("lambda1","lambda1",this,lambda1),
  lambda2_("lambda2","lambda2",this,lambda2)
{
}


RooDoubleBukin::RooDoubleBukin(const RooDoubleBukin& other, const char* name) :
  RooAbsPdf(other, name),
  x_("x",this,other.x_),
  xg1_("xg1",this,other.xg1_),
  xg2_("xg2",this,other.xg2_),
  ksi_("ksi",this,other.ksi_),
  sigma1_("sigma1",this,other.sigma1_),
  sigma2_("sigma2",this,other.sigma2_),
  lambda1_("lambda1",this,other.lambda1_),
  lambda2_("lambda2",this,other.lambda2_)
{
}


TObject* RooDoubleBukin::clone(const char* newname) const {
  return new RooDoubleBukin(*this, newname);
}


double RooDoubleBukin::evaluate() const {
	double fit_result = 0, bukin1 = 0, bukin2 = 0;
//	erfx = - (x_ - xm_) * lambda1_ / (sigma1_ * abs(lambda1_) *sqrt(2)) + sigma1_ / (abs(lambda1_) * sqrt(2));
//	expx = - (x_ - xm_)/lambda1_ + sigma1_*sigma1_ / (2 * lambda1_*lambda1_);
//	bukin1 = 0.5 * lambda1_ * (1 - TMath::Erf(erfx)) * exp(expx);

	bukin1 = EvaluateSingleBukin(x_,lambda1_,sigma1_,xg1_);
	bukin2 = EvaluateSingleBukin(x_,lambda2_,sigma2_,xg2_);
	fit_result = cos(ksi_)*cos(ksi_) * bukin1 + sin(ksi_)*sin(ksi_)*bukin2;
//	fit_result = bukin1;


//	double fit_result = 0;
//	fit_result = EvaluateSingleBukin(x_,lambda1_,sigma1_,xm_);

	return fit_result;
}

double EvaluateSingleBukin(const RooRealProxy &x ,const RooRealProxy& lambda,const RooRealProxy& sigma,const RooRealProxy& xg){
	double fit_result = 0, erfx = 0, expx = 0;
	erfx = - (x - xg) * lambda / (sigma * abs(lambda) *sqrt(2)) + sigma / (abs(lambda) * sqrt(2));
	expx = - (x - xg)/lambda + sigma*sigma / (2 * lambda*lambda);
	fit_result = 0.5 * lambda * (1 - TMath::Erf(erfx)) * exp(expx);

	return fit_result;
}


