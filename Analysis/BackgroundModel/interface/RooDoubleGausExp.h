#ifndef Analysis_BackgroundModel_RooDoubleGausExp_h
#define Analysis_BackgroundModel_RooDoubleGausExp_h 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"


namespace analysis {
  namespace backgroundmodel {

    class RooDoubleGausExp : public RooAbsPdf {
    public:
      inline RooDoubleGausExp() = default;
      //RooDoubleGausExp() : fValue(-1) {};
      RooDoubleGausExp(const char *name, const char *title,
                  RooAbsReal& x,
                  RooAbsReal& mean,
                  RooAbsReal& sigmaL,
				  RooAbsReal& sigmaR,
				  RooAbsReal& tail_shift,
				  RooAbsReal& tail_sigma);
      RooDoubleGausExp(const RooDoubleGausExp& other, const char* name=0) ;
      virtual TObject* clone(const char* newname) const;
      inline virtual ~RooDoubleGausExp() = default;

    protected:
      double evaluate() const ;

      RooRealProxy x_ ;
      RooRealProxy mean_ ;
      RooRealProxy sigmaL_ ;
      RooRealProxy sigmaR_ ;
      RooRealProxy tail_shift_;
      RooRealProxy tail_sigma_;

    private:
      ClassDef(RooDoubleGausExp,1);
    };

  }
}

#endif  // Analysis_BackgroundModel_RooDoubleGausExp_h
