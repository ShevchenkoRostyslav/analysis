#ifndef Analysis_BackgroundModel_RooExpPolynom_h
#define Analysis_BackgroundModel_RooExpPolynom_h 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooListProxy.h"

class RooRealVar;
class RooAbsReal;

namespace analysis {
  namespace backgroundmodel {

    class RooExpPolynom : public RooAbsPdf {
    public:
      inline RooExpPolynom() = default;
      RooExpPolynom(const char *name, const char *title,
    		  RooAbsReal& x,
			  RooAbsReal& para,
    		  const RooArgList& coefList);
      RooExpPolynom(const RooExpPolynom& other, const char* name=0) ;
      virtual TObject* clone(const char* newname) const;
      inline virtual ~RooExpPolynom() {};

    protected:
      double evaluate() const ;

      RooRealProxy x_;
      RooRealProxy para_;
      RooListProxy coefList_;

    private:
      ClassDef(RooExpPolynom,1);
    };

  }
}

#endif  // Analysis_BackgroundModel_RooExpPolynom_h
