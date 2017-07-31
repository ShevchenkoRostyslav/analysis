#include <iostream>
#include <cmath>
#include <math.h>
#include "Riostream.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "TMath.h"
#include "Analysis/BackgroundModel/interface/RooBukinExp.h"


using namespace analysis::backgroundmodel;

ClassImp(RooBukinExp)


RooBukinExp::RooBukinExp(const char *name, const char *title,
							RooAbsReal& x,
		                    RooAbsReal& Xp,
		                    RooAbsReal& sigp,
		                    RooAbsReal& xi,
		                    RooAbsReal& rho1,
							RooAbsReal& rho2,
							RooAbsReal& Xt,
							RooAbsReal& tau) :
  RooAbsPdf(name, title),
  x_("x","x",this,x),
  Xp_("Xp","Xp",this,Xp),
  sigp_("sigp","sigp",this,sigp),
  xi_("xi","xi",this,xi),
  rho1_("rho1","rho1",this,rho1),
  rho2_("rho2","rho2",this,rho2),
  Xt_("Xt","Xt",this,Xt),
  tau_("tau","tau",this,tau){
}


RooBukinExp::RooBukinExp(const RooBukinExp& other, const char* name) :
  RooAbsPdf(other, name),
  x_("x",this,other.x_),
  Xp_("Xp",this,other.Xp_),
  sigp_("sigp",this,other.sigp_),
  xi_("xi",this,other.xi_),
  rho1_("rho1",this,other.rho1_),
  rho2_("rho2",this,other.rho2_),
  Xt_("Xt",this,other.Xt_),
  tau_("tau",this,other.tau_) {
}


TObject* RooBukinExp::clone(const char* newname) const {
  return new RooBukinExp(*this, newname);
}


double RooBukinExp::evaluate() const {
	double r1=0,r2=0,r3=0,r4=0,r5=0,hp=0, r2_tail = 0;
	double x1 = 0,x2 = 0;
	double fit_result = 0;
	double consts = 2*sqrt(2*log(2.));
	hp=sigp_*consts;
	r3=log(2.);
	r4=sqrt(TMath::Power(xi_,2)+1);
	r1=xi_/r4;

	if(TMath::Abs(xi_) > exp(-6.)){
		r5=xi_/log(r4+xi_);
	}
	else
		r5=1;
	x1 = Xp_ + (hp / 2) * (r1-1);
	x2 = Xp_ + (hp / 2) * (r1+1);
	//--- Left Side
	if(x_ < x1){
		r2		= rho1_*TMath::Power((x_-x1)/(Xp_-x1),2)-r3 + 4 * r3 * (x_-x1)/hp * r5 * r4/TMath::Power((r4-xi_),2);
//		r2_tail = rho1_*TMath::Power((Xt_-x1)/(Xp_-x1),2)-r3 + 4 * r3 * (Xt_-x1)/hp * r5 * r4/TMath::Power((r4-xi_),2);
	}


	//--- Center
	else if(x_ < x2) {
		if(TMath::Abs(xi_) > exp(-6.)) {
			r2		= log(1 + 4 * xi_ * r4 * (x_-Xp_)/hp)/log(1+2*xi_*(xi_-r4));
			r2=-r3*(TMath::Power(r2,2));
		}
		else{
			r2		=-4*r3*TMath::Power(((x_-Xp_)/hp),2);
		}
	}

	//--- Right Side before the Expo tail

	else {
		r2		= rho2_*TMath::Power((x_-x2)/(Xp_-x2),2)-r3 - 4 * r3 * (x_-x2)/hp * r5 * r4/TMath::Power((r4+xi_),2);
	}

	if(TMath::Abs(r2) > 100){
		fit_result = 0;
	}
	else{
		//---- Normalize the result
		fit_result = exp(r2);
	}
//	if(Xt_ > x2)
	r2_tail = rho2_*TMath::Power((Xt_-x2)/(Xp_-x2),2)-r3 - 4 * r3 * (Xt_-x2)/hp * r5 * r4/TMath::Power((r4+xi_),2);
//	else {
//		if(TMath::Abs(xi_) > exp(-6.)) {
//			r2_tail = log(1 + 4 * xi_ * r4 * (Xt_-Xp_)/hp)/log(1+2*xi_*(xi_-r4));
//			r2_tail=-r3*(TMath::Power(r2_tail,2));
//		}
//		else{
//			r2_tail	=-4*r3*TMath::Power(((Xt_-Xp_)/hp),2);
//		}
//	}
	//--- Right expo tail
	double norm = 1;
	if(x_ >= Xt_) {
		norm = TMath::Exp(r2_tail);
		fit_result = TMath::Exp( - ( x_ - Xt_) / tau_);
		fit_result *= norm;
	}
//	std::cout<<"WTF: x1 = "<<x1<<" x2 = "<<x2<<" Xt = "<<Xt_<<std::endl;

	return fit_result;
}



