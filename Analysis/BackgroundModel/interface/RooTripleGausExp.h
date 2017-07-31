#ifndef Analysis_BackgroundModel_RooTripleGausExp_h
#define Analysis_BackgroundModel_RooTripleGausExp_h 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"


namespace analysis {
  namespace backgroundmodel {

    class RooTripleGausExp : public RooAbsPdf {
    public:
      inline RooTripleGausExp() = default;
      //RooTripleGausExp() : fValue(-1) {};
      RooTripleGausExp(const char *name, const char *title,
                  RooAbsReal& x,
                  RooAbsReal& mean,
                  RooAbsReal& sigmaL1,
				  RooAbsReal& sigmaR1,
				  RooAbsReal& sigmaR2,
				  RooAbsReal& tail_shift,
				  RooAbsReal& tail_sigma,
				  RooAbsReal& norm_g2);
      RooTripleGausExp(const RooTripleGausExp& other, const char* name=0) ;
      virtual TObject* clone(const char* newname) const;
      inline virtual ~RooTripleGausExp() = default;

    protected:
      double evaluate() const ;

      RooRealProxy x_ ;
      RooRealProxy mean_ ;
      RooRealProxy sigmaL1_ ;
      RooRealProxy sigmaR1_ ;
      RooRealProxy sigmaR2_ ;
      RooRealProxy tail_shift_;
      RooRealProxy tail_sigma_;
      RooRealProxy norm_g2_;

    private:
      ClassDef(RooTripleGausExp,1)
    };

  }
}

#endif  // Analysis_BackgroundModel_RooTripleGausExp_h
