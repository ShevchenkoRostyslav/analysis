#ifndef Analysis_BackgroundModel_RooDoubleBukin_h
#define Analysis_BackgroundModel_RooDoubleBukin_h 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooRealVar;
class RooAbsReal;

namespace analysis {
  namespace backgroundmodel {

    class RooDoubleBukin : public RooAbsPdf {
    public:
      inline RooDoubleBukin() = default;
      RooDoubleBukin(const char *name, const char *title,
                    RooAbsReal& x,
                    RooAbsReal& xg1,
					RooAbsReal& xg2,
					RooAbsReal& ksi,
					RooAbsReal& sigma1,
					RooAbsReal& lambda1,
					RooAbsReal& sigma2,
					RooAbsReal& lambda2
					);
      RooDoubleBukin(const RooDoubleBukin& other, const char* name=0) ;
      virtual TObject* clone(const char* newname) const;
      inline virtual ~RooDoubleBukin() {};

    protected:
      double evaluate() const ;

      RooRealProxy x_;
      RooRealProxy xg1_;
      RooRealProxy xg2_;
      RooRealProxy ksi_;
      RooRealProxy sigma1_;
      RooRealProxy sigma2_;
	  RooRealProxy lambda1_;
	  RooRealProxy lambda2_;

    private:
      ClassDef(RooDoubleBukin,1);
    };

  }
}

#endif  // Analysis_BackgroundModel_RooDoubleBukin_h
