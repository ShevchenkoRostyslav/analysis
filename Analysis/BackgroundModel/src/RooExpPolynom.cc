#include <iostream>
#include <cmath>
#include <math.h>
#include "Riostream.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "Analysis/BackgroundModel/interface/RooExpPolynom.h"


using namespace analysis::backgroundmodel;

ClassImp(RooExpPolynom)


RooExpPolynom::RooExpPolynom(const char *name, const char *title,
		RooAbsReal& x,
		RooAbsReal& para,
		const RooArgList& coefList) :
  RooAbsPdf(name, title),
  x_("x", "x", this, x),
  para_("para", "para", this, para),
  coefList_("coefList","List of coefficients",this){
	RooFIter coefIter = coefList.fwdIterator();
	RooAbsArg *coef;
	while((coef = static_cast<RooAbsArg*>(coefIter.next()))) {
		if (!static_cast<RooAbsReal*>(coef)) {
			throw std::invalid_argument("RooExpPolynom::" + std::string(GetName()) + ", coeficient " + std::string(coef->GetName()) + " is not of type RooAbsReal");
		}
		coefList_.add(*coef);
	}
}


RooExpPolynom::RooExpPolynom(const RooExpPolynom& other, const char* name) :
  RooAbsPdf(other, name),
  x_("x", this, other.x_),
  para_("para", this, other.para_),
  coefList_("coefList",this,other.coefList_) {
}


TObject* RooExpPolynom::clone(const char* newname) const {
  return new RooExpPolynom(*this, newname);
}


double RooExpPolynom::evaluate() const {
	double fit_result = 0, arg = 0;
	double std = x_ / 13000.;
	RooFIter it = coefList_.fwdIterator();
	RooAbsReal* c;
	int i = 0;
//	arg = para_*TMath::Log(std);
	while ((c = (RooAbsReal*) it.next())) {
		arg += c->getVal() * TMath::Power(std,i);
		++i;
	}
//	int i = 0;
//	while ((c = (RooAbsReal*) it.next())) {
//		arg += c->getVal() * TMath::Power(x_,i);
//		++i;
//	}
	fit_result = TMath::Exp(-arg);
	fit_result *= (1+para_ * std);

	return fit_result;
}
