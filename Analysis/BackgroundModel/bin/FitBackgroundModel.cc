#include <vector>
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/join.hpp>
#include "TSystem.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "RooFitResult.h"
#include "Analysis/BackgroundModel/interface/HistContainer.h"
#include "Analysis/BackgroundModel/interface/TreeContainer.h"
#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/BackgroundModel/interface/ParamModifier.h"
#include "Analysis/BackgroundModel/interface/Tools.h"
#include "Analysis/BackgroundModel/interface/ProbabilityDensityFunctions.h"
#include "Analysis/MssmHbb/interface/HbbStyleClass.h"


namespace po = boost::program_options;
namespace ab = analysis::backgroundmodel;

int backgroundOnlyFit(ab::FitContainer&, po::variables_map&);


int main(int argc, char* argv[]) {
  const auto cmsswBase = static_cast<std::string>(gSystem->Getenv("CMSSW_BASE"));
  HbbStyle style;
  style.setTDRstyle(PRELIMINARY);


  // general command line options
  po::options_description cmdLineOptions("Optional arguments");
  cmdLineOptions.add_options()
    ("help,h", "Produce help message.")
    ("verbose,v", po::value<int>()->default_value(0), "More verbose output.")
    ("profile,p", "Create profile likelihoods of the fit parameters.")
    ("list_parameters,l", "List parameters of the chosen fit models and exit.")
    ("input_file,i", po::value<std::string>()
     ->default_value(cmsswBase+"/src/Analysis/BackgroundModel/"
                     "data/HIG14017_HighMass2012_Packed_M350_inputs.root"),
     "ROOT file from which input histograms are retrieved.")
    ("input_tree_file,t", po::value<std::string>(),
     "ROOT file from which input tree is retrieved. If this parameter is not "
     "given a binned fit is performed using 'input_file'.")
    ("output_directory,o", po::value<std::string>()
     ->default_value(cmsswBase+"/src/Analysis/BackgroundModel/test/"),
     "Directory where the output is stored.")
    ("fit_min", po::value<float>(), "Lower bound of the fit range.")
    ("fit_max", po::value<float>(), "Upper bound of the fit range.")
    ("nbins", po::value<int>(), "Number of bins in the fit range.")
	("control_region", po::value<bool>()->default_value(true), "Number of bins in the fit range.")
    ("modify_param,m", po::value<std::vector<std::string> >()->composing()
     ->default_value(std::vector<std::string>(), ""),
     "Modify parameters as follows: "
     "\"name: [start=<value>,] [min=<value>,] [max=<value>,] "
     "[constant,] [floating]\"")
    ;
  po::variables_map vm;
//  po::store(po::command_line_parser( argc, argv ).options(cmdLineOptions).run(), vm);
  po::store(po::command_line_parser(argc, argv).options(cmdLineOptions)
            .allow_unregistered().run(), vm);

  // now add required options
  po::options_description requiredOptions("Required arguments");
  requiredOptions.add_options()
    ("background_model,b", po::value<std::string>()->required(),
     ("Name of the background model (" +
      boost::algorithm::join(ab::ProbabilityDensityFunctions::getAvailableModels(), ", ") +
      ").").c_str())
    ;
  po::store(po::command_line_parser(argc, argv).options(requiredOptions)
            .allow_unregistered().run(), vm);
//  po::store(po::command_line_parser(argc, argv).options(requiredOptions).run(), vm);


  po::options_description allOptions("Allowed arguments");
  allOptions.add(cmdLineOptions).add(requiredOptions);

  po::notify(vm);

  // check for help flag before checking for required options
  if (vm.count("help")) {
    std::cout << allOptions << std::endl;
    return 0;
  }

  try {
    po::notify(vm);
  } catch (const po::required_option& ex) {
    std::cerr << ex.what() << std::endl;
    return 1;
  }


  // check for required verbosity:
  //
  // ROOT error message levels:
  //    kUnset, kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal
  // RooFit error message levels:
  //    DEBUG=0, INFO=1, PROGRESS=2, WARNING=3, ERROR=4, FATAL=5
  auto verbosity = vm["verbose"].as<int>();
  switch (verbosity) {
  case 0:
    gErrorIgnoreLevel = kError;
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    break;
  case 1:
    gErrorIgnoreLevel = kInfo;
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
    break;
  case 2:
    RooMsgService::instance().setGlobalKillBelow(RooFit::PROGRESS);
    break;
  case 3:
    gErrorIgnoreLevel = kPrint;
    RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
    break;
  default:
      if (verbosity < 0) {
        gErrorIgnoreLevel = kFatal;
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
      } else {  // highest level of verbosity reached
        gErrorIgnoreLevel = kUnset;
        RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);
      }
  }


  std::cout << "Fitting Background Model..." << std::endl;

  if (vm["input_tree_file"].empty()) {
    ab::HistContainer input(vm["input_file"].as<std::string>());
    input.show();
    ab::FitContainer fitter(input, vm["output_directory"].as<std::string>());
    fitter.verbosity(verbosity - 1);
    return backgroundOnlyFit(fitter, vm);
  } else {
    ab::TreeContainer input(vm["input_tree_file"].as<std::string>());
    input.show();
    ab::FitContainer fitter(input, vm["output_directory"].as<std::string>());
    fitter.verbosity(verbosity - 1);
    return backgroundOnlyFit(fitter, vm);
  }
}


int backgroundOnlyFit(ab::FitContainer& fitter, po::variables_map& vm) {
  if (!vm["fit_min"].empty()) fitter.fitRangeMin(vm["fit_min"].as<float>());
  if (!vm["fit_max"].empty()) fitter.fitRangeMax(vm["fit_max"].as<float>());
  if (!vm["nbins"].empty())   fitter.setNBins(vm["nbins"].as<int>());
  std::vector<ab::ParamModifier> bkgModifiers =
    ab::parseModifiers(vm["modify_param"].as<std::vector<std::string> >());
  fitter.setModel(ab::FitContainer::Type::background,
                  vm["background_model"].as<std::string>(),
                  bkgModifiers);

  if (vm.count("list_parameters")) {
    fitter.showModels();
    return 0;
  }

  int returnValue = 0;
  std::unique_ptr<RooFitResult> bkgOnlyFit = fitter.backgroundOnlyFit(vm["background_model"].as<std::string>(),0,vm["control_region"].as<bool>());
  if (bkgOnlyFit) {
    std::cout << "\nCorrelation matrix of background-only fit:" << std::endl;
    if (&(bkgOnlyFit->correlationMatrix()) != nullptr) {
      bkgOnlyFit->correlationMatrix().Print("v");

      TMatrixDSymEigen CM = bkgOnlyFit->correlationMatrix();
      TMatrixD CM_eigen = CM.GetEigenVectors();
      std::cout << "Eigenvectors of Correlation Matrix" << std::endl;
      CM_eigen.Print("v");

    } else {
      std::cout << ">>> no correlation matrix available" << std::endl;
      returnValue = 1;
    }
    if (&(bkgOnlyFit->covarianceMatrix()) != nullptr) {
      std::cout << "Covariance Matrix" << std::endl;
      bkgOnlyFit->covarianceMatrix().Print("v");

      TMatrixDSymEigen VM = bkgOnlyFit->covarianceMatrix();
      TMatrixD VM_eigen = VM.GetEigenVectors();
      std::cout << "Eigenvectors of Covariance Matrix" << std::endl;
      VM_eigen.Print("v");
    } else {
      std::cout << ">>> no covariance matrix available" << std::endl;
      returnValue = 1;
    }
  } else {
    returnValue = 1;
  }

  if (vm.count("profile")) {
    fitter.profileModel(ab::FitContainer::Type::background);
  }

  for (const auto& m : bkgModifiers) m.show();

  return returnValue;
}
