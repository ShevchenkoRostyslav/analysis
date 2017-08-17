/*
 * ProbabilityDensityFunctions.h
 *
 *  Created on: Sep 10, 2016
 *      Author: shevchen
 */

#ifndef BACKGROUNDMODEL_INTERFACE_PROBABILITYDENSITYFUNCTIONS_H_
#define BACKGROUNDMODEL_INTERFACE_PROBABILITYDENSITYFUNCTIONS_H_

#include <sstream>
#include <exception>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

#include <memory>
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooFormulaVar.h"
#include "RooNovosibirsk.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooEffProd.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooBukinPdf.h"
#include "RooBernstein.h"
#include "RooChebychev.h"
#include "RooBreitWigner.h"

#include "Analysis/BackgroundModel/interface/RooDoubleCB.h"
#include "Analysis/BackgroundModel/interface/RooExpGausExp.h"
#include "Analysis/BackgroundModel/interface/RooGausExp.h"
#include "Analysis/BackgroundModel/interface/RooExpBWExp.h"
#include "Analysis/BackgroundModel/interface/RooPhaseSpace.h"
#include "Analysis/BackgroundModel/interface/RooPhaseSpacePol4.h"
#include "Analysis/BackgroundModel/interface/RooRelBreitWigner.h"
#include "Analysis/BackgroundModel/interface/RooQuadGausExp.h"
#include "Analysis/BackgroundModel/interface/RooMyNovosibirsk.h"
#include "Analysis/BackgroundModel/interface/RooSuperNovosibirsk.h"
#include "Analysis/BackgroundModel/interface/RooExtendNovosibirsk.h"
#include "Analysis/BackgroundModel/interface/RooDoubleGausExp.h"
#include "Analysis/BackgroundModel/interface/RooPolyDijet.h"
#include "Analysis/BackgroundModel/interface/RooSuperDiJet.h"
#include "Analysis/BackgroundModel/interface/RooExpPolynom.h"

//#include "Analysis/BackgroundModel/interface/FitContainer.h"
#include "Analysis/BackgroundModel/interface/Tools.h"

namespace analysis {
  namespace backgroundmodel {

  class ProbabilityDensityFunctions {
  public:
  	ProbabilityDensityFunctions(RooWorkspace& workspace, const std::string& var_name, const bool& modify_par_names = false);
  	virtual ~ProbabilityDensityFunctions();

  	void setPdf(const std::string& function, const std::string& name);
  	//PDFs getters

  	void getNovosibirsk(const std::string& name, const std::string& title, RooRealVar& peak, RooRealVar& width, RooRealVar& tail);
  	void getNovosibirsk(const std::string& name);

  	void getNovoPSProd(const std::string& name, const std::string& title, RooRealVar& peak, RooRealVar& width, RooRealVar& tail, RooFormulaVar& PS);
  	void getNovoPSProd(const std::string& name);

  	void getNovoEffProd(const std::string& name);
	void getNovoEfffixProd(const std::string& name);
  	void getNovoPSHighMPol4(const std::string& name);
	void getMyNovosibirsk(const std::string& name);
	void getSuperNovosibirsk(const std::string& name, const int& degree);
	void getMyNovoPSProd(const std::string& name);
	void getExtNovosibirsk(const std::string& name);
	void getExtNovoPSProd(const std::string& name);
	void getExtNovoEffProd(const std::string& name);
	void getSuperNovoEffProd(const std::string& name, const int& degree);
	void getExtNovoEfffixProd(const std::string& name);
	void getExtNovoLogisticProd(const std::string& name);
	void getExtNovoExtLogisticProd(const std::string& name);
	void getExtNovoHypertanProd(const std::string& name);
  	void getCrystalBall(const std::string& name);
  	void getCrystalPSProd(const std::string& name);
  	void getCrystalEffProd(const std::string& name);
  	void getExpEffProd(const std::string& name);
  	void getDoubleCB(const std::string& name);
  	void getDijetv1(const std::string& name);
  	void getDijetv1PSProd(const std::string& name);
  	void getDijetv2(const std::string& name);
  	void getDijetv2PSProd(const std::string& name);
	void getDijetv3(const std::string& name);
	void getDijetv3EffProd(const std::string& name);
	void getDijetv3LogisticProd(const std::string& name);
	void getDijetv4(const std::string& name);
	void getDijetv4LogisticProd(const std::string& name);
  	void getExpGausExp(const std::string& name);
  	void getGausExp(const std::string& name);
  	void getGausExpPSProd(const std::string& name);
	void getGausExpEffProd(const std::string& name);
  	void getExpBWExp(const std::string& name);
  	void getBukin(const std::string& name);
  	void getBukinPSProd(const std::string& name);		
  	void getBernstein(const std::string& name, const int& numCoeffs);
  	void getChebychev(const std::string& name, const int& numCoeffs);
  	void getBernEffProd(const std::string& name, const int& numCoeffs);
	void getBernEfffixProd(const std::string& name, const int& numCoeffs);
  	void getBernPSProd(const std::string& name, const int& numCoeffs);	
  	void getChebEffProd(const std::string& name, const int& numCoeffs);
  	void getDoubleGausExp(const std::string& name);
  	void getTripleGausExp(const std::string& name);
  	void getBreitWigner(const std::string& name);
  	void getRelBreitWigner(const std::string& name);
  	void getRooQuadGausExp(const std::string& name);
  	void getSuperDiJet(const std::string& name, const int& degree);
  	void getExpPolynom(const std::string& name, const int& degree);
  	void getSuperDiJetEffProd(const std::string& name, const int& degree);
  	void getSuperDiJetLinearProd(const std::string& name, const int& degree);

  	//Section of base function that can be used in combination
  	//with more comlicated
  	void getPhaseSpace(const std::string& name);	//PhaseSpace factor
  	void getEfficiency(const std::string& name);	//Efficiency function
	void getEfficiencyFix(const std::string& name);    //Efficiency function fix parameters
	void getLogistic(const std::string& name);
	void getLogistic_ext1(const std::string& name);  
	void getHyperbolictan(const std::string& name);
	void getExp(const std::string& name,const std::string& name1 = "");			//Exponential function
  	void getGaus(const std::string& name,const std::string& mean_name = "", const std::string& sigma_name = "");			//Gaus function


  	void setPeakStart(const double& peak);
  	const double& getPeakStart(const double& max = 250);

  	void showAvailableModels();
  	inline static const std::vector<std::string>& getAvailableModels() {
  	        return availableModels_; };

  private:


  	double peak_;
  	RooWorkspace* workspace_;
  	std::string var_;
  	bool modify_par_names_;
  	static const std::vector<std::string> availableModels_;
  	const int defaultNumberOfCoefficients_ = 9;

  	const char* paramName_(const std::string& name, const std::string& default_name);
  	static std::unique_ptr<RooArgList>
       	getCoefficients_(const int numCoeffs, const std::string& name);

  	void ModifyVarNames_(const std::string& name);


  };

  inline void ProbabilityDensityFunctions::setPeakStart(const double& peak){ peak_ = peak;}
  inline const double& ProbabilityDensityFunctions::getPeakStart(const double& max) {if(peak_ == -100) return max; else return peak_;}
//  inline static const std::vector<std::string>& ProbabilityDensityFunctions::getAvailableModels(){return availableModels_;}
  inline const char* ProbabilityDensityFunctions::paramName_(const std::string& name, const std::string& default_name){
  																if(name != "") return name.c_str(); else return default_name.c_str();}

  }
}
#endif /* BACKGROUNDMODEL_INTERFACE_PROBABILITYDENSITYFUNCTIONS_H_ */
