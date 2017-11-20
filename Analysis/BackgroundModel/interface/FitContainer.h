#ifndef Analysis_BackgroundModel_FitContainer_h
#define Analysis_BackgroundModel_FitContainer_h 1

#include <vector>
#include <string>
#include <memory>
#include <array>
#include <fstream>

#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TString.h"
#include "TIterator.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "RooMinimizer.h"

#include "RooPlot.h"
#include "RooCurve.h"
#include "RooHist.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooList.h"

#include "RooGaussian.h"
#include "RooAddPdf.h"

#include "Analysis/BackgroundModel/interface/HistContainer.h"
#include "Analysis/BackgroundModel/interface/TreeContainer.h"
#include "Analysis/BackgroundModel/interface/ParamModifier.h"
#include "Analysis/BackgroundModel/interface/ProbabilityDensityFunctions.h"
#include "Analysis/BackgroundModel/interface/RooFitQuality.h"
#include "Analysis/Tools/interface/RooFitUtils.h"

namespace analysis {
  namespace backgroundmodel {

    class FitContainer {
    public:
      enum class Type { signal, background };
      inline static std::string toString(const Type& type) {
        switch (type) {
        case Type::signal: return "signal";
        case Type::background: return "background";
        };
        return "";              // to silence compiler
      };

      struct hBands{
      	TH1D *h1sigmaU = nullptr;
      	TH1D *h2sigmaU = nullptr;
      	TH1D *h1sigmaD = nullptr;
      	TH1D *h2sigmaD = nullptr;
      	TH1D *central  = nullptr;
	};

		int getBinsToPlot() const {
			return nBinsToPlot_;
		}

		void setBinsToPlot(int binsToPlot) {
			nBinsToPlot_ = binsToPlot;
		}

      FitContainer(const TH1* data, const TH1* signal, const TH1* background,
		   const std::string& outputDir = defaultOutputDir_);
      FitContainer(const TH1* data, const std::string& outputDir = defaultOutputDir_, const std::string & type = "data");
//      FitContainer(const TH1* data, const std::string& outputDir = defaultOutputDir_);
      FitContainer(TTree& data, const std::string& outputDir = defaultOutputDir_);
      FitContainer(const HistContainer& container,
		   const std::string& outputDir = defaultOutputDir_);
      FitContainer(const TreeContainer& container,
		   const std::string& outputDir = defaultOutputDir_);
      virtual ~FitContainer();
      void initialize();

      FitContainer(const FitContainer&) = default;
      FitContainer& operator=(const FitContainer&) = default;
      FitContainer(FitContainer&&) = default;
      FitContainer& operator=(FitContainer&&) = default;

      FitContainer& verbosity(int level);
      FitContainer& fitRangeMin(float min);
      FitContainer& fitRangeMax(float max);
      FitContainer& setNBins(int nbins);
      RooWorkspace& getWorkspace();

      void setModel(const Type& type, const std::string& model);
      void setModel(const Type& type, const std::string& model,
                    const std::vector<ParamModifier>& modifiers);
      std::unique_ptr<RooFitResult> backgroundOnlyFit(const std::string& model, const bool& plot_params = 0, const bool& control_region = true);
      std::unique_ptr<RooFitResult> FitSignal(const std::string & model, const bool& plot_params = 0);

      void profileModel(const Type& type);
      void showModels() const;
      void Import(const RooAbsArg& inArg);
      void Write();

      void SetChi2CalcLowEdge(const double & val);
      void SetChi2CalcHighEdge(const double & val);

      void SetupTopFrame(RooPlot *frame1, const double& pad_w, const double& pad_h);
      void SetupBottomFrame(RooPlot *frame1, const double& pad_w, const double& pad_h);

    private:

      //Private constructor to avoid code duplication for private members initialisation
      FitContainer(const std::string& outputDir);

      // methods to set the fit model
      double getPeakStart_(const Type& type,const double& max);
      double getPeakStart_(const Type& type);
      double getMaxPosition_(const RooAbsData& data);
      hBands getPullBands_(RooPlot * frame, const std::string& curve_name, TH1D * hData, RooRealVar & x, RooAbsPdf & fit);

      // internal methods
      static void prepareCanvas_(TCanvas& raw);
      static void prepareFrame_(RooPlot& raw);
      std::string getOutputPath_(const std::string& subdirectory = "");
      int getNonZeroBins_(const RooAbsData& data);
      int getBlindedBins_(const RooAbsData& data, double blind_lowEdge, double blind_highEdge);
      bool applyModifiers_(RooAbsPdf& pdf,
                           const std::vector<ParamModifier>& modifiers);
      void makeLog_(const RooFitResult& fitResult);

      // data member
      static const std::string defaultOutputDir_;
      bool initialized_;
      bool written_;
      bool splitrange_;
      std::string outputDir_;
      std::string plotDir_;
      std::string workspaceDir_;
      std::string fullRangeId_;
      std::string fitRangeId_;
      std::string fitRangeLowId_;	// for split range by CA
      //std::string fitRangeMedId_;       // for split range by CA
      std::string fitRangeHighId_;	// for split range by CA
      std::string fitSplRangeId_;	// for split range by CA
      double blind_lowEdge_;
      double blind_highEdge_;
      int verbosity_;
      RooWorkspace workspace_;
      std::string outRootFileName_;
      std::string mbb_;
      std::string weight_;
      std::string data_;
      std::string signal_;
      std::string bkg_;
      float fitRangeMin_;
      float fitRangeMax_;
      TTree bkgOnlyFit_;
      double chi2_lowEdge_;
      double chi2_highEdge_;
      float chi2BkgOnly_;
      float normChi2BkgOnly_;
      int ndfBkgOnly_;
      double covMatrix_[400];
      double eigenVector_[400];
      int nbins_;
      int nBinsToPlot_;
      float lumi_;
      float obs_;
    };

    inline void FitContainer::Import(const RooAbsArg& inArg){ workspace_.import(inArg);}
    inline void FitContainer::Write(){ if(!written_) { workspace_.writeToFile(outRootFileName_.c_str()); written_ = true;}   }
    inline RooWorkspace& FitContainer::getWorkspace() {return workspace_;};
    inline void FitContainer::SetChi2CalcLowEdge(const double& val){chi2_lowEdge_ = val;};
    inline void FitContainer::SetChi2CalcHighEdge(const double& val){chi2_highEdge_ = val;};

  }
}

#endif  // Analysis_BackgroundModel_FitContainer_h
