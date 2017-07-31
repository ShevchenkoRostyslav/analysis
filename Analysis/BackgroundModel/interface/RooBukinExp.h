#ifndef Analysis_BackgroundModel_RooBukinExp_h
#define Analysis_BackgroundModel_RooBukinExp_h 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooRealVar;
class RooAbsReal;

namespace analysis {
  namespace backgroundmodel {

    class RooBukinExp : public RooAbsPdf {
    public:
      inline RooBukinExp() = default;
      RooBukinExp(const char *name, const char *title,
                    RooAbsReal& x,
                    RooAbsReal& Xp,
                    RooAbsReal& sigp,
                    RooAbsReal& xi,
                    RooAbsReal& rho1,
					RooAbsReal& rho2,
					RooAbsReal& Xt,
					RooAbsReal& tau);
      RooBukinExp(const RooBukinExp& other, const char* name=0) ;
      virtual TObject* clone(const char* newname) const;
      inline virtual ~RooBukinExp() {};

    protected:
      double evaluate() const ;

      RooRealProxy x_;
      RooRealProxy Xp_;
      RooRealProxy sigp_;
      RooRealProxy xi_;
      RooRealProxy rho1_;
      RooRealProxy rho2_;
      RooRealProxy Xt_;
      RooRealProxy tau_;

    private:
      ClassDef(RooBukinExp,1);
    };

  }
}

#endif  // Analysis_BackgroundModel_RooBukinExp_h
