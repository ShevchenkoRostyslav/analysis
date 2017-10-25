/*
 * ProbabilityDensityFunctions.cpp
 *
 *  Created on: Sep 10, 2016
 *      Author: shevchen
 */

#include "Analysis/BackgroundModel/interface/ProbabilityDensityFunctions.h"

using namespace analysis::backgroundmodel;

ProbabilityDensityFunctions::ProbabilityDensityFunctions(RooWorkspace& workspace, const std::string& var_name, const bool& modify_par_names) : peak_(-100), workspace_(&workspace), var_(var_name), modify_par_names_(modify_par_names) {
	// TODO Auto-generated constructor stub
	//	workspace_.reset(&workspace);// = &workspace;
	//Check whether "X" variable exist in the workspace:
	if(!workspace_->var(var_name.c_str())){
		std::cerr<<"Variable "<<var_name<<" wasn't found in the provided workspace. Interupt"<<std::endl;
		exit(-1);
	}
}

ProbabilityDensityFunctions::~ProbabilityDensityFunctions() {
	// TODO Auto-generated destructor stub
//	delete workspace_;
}

void ProbabilityDensityFunctions::setPdf(const std::string& function, const std::string& name){

	if (workspace_->allPdfs().find(name.c_str())) {
		std::stringstream msg;
	    msg << "Model for " << name << " has already been set!";
	    throw std::runtime_error(msg.str());
	}
	std::string function_name = function;

	int numCoeffs = defaultNumberOfCoefficients_;
	std::vector<std::string> nameSplitted;
	boost::split(nameSplitted, function, boost::is_any_of(","));
	switch (nameSplitted.size()) {
	case 1:
		if(function_name.find("supernovo") !=std::string::npos) numCoeffs = 0;
		break;
	case 2:
		try {
	      numCoeffs = std::stoi(nameSplitted[1]);
		} catch (const std::invalid_argument& ex) {
	      if (std::string(ex.what()) == "stoi") {
	        std::cerr << "Conversion of '" << nameSplitted[1]
	                  << "' to an integral value for the number of coefficients "
	                  << "failed. Using the default value of " << numCoeffs << "."
	                  << std::endl;
	      } else {
	        throw;
	      }
	    }
		function_name = nameSplitted[0];
	    break;
	  default:
	    std::stringstream msg;
	    msg << "Unsupported number of arguments for fit model: "
	        << nameSplitted.size();
	    throw std::runtime_error(msg.str());
	  }

	if (function_name == "novosibirsk") getNovosibirsk(name);			//Novosibirsk (3)
	else if (function_name == "novopsprod") getNovoPSProd(name);		//PS x Novosibirsk (3)
	else if (function_name == "novoeffprod") getNovoEffProd(name);             //Eff x Novosibirsk (5)
	else if (function_name == "novoefffixprod") getNovoEfffixProd(name);
	else if (function_name == "novopshighMpol4") getNovoPSHighMPol4(name);	//PS x Novosibirsk (3) HighM Poly[4]
	else if (function_name == "crystalball") getCrystalBall(name);
	else if (function_name == "crystalpsprod") getCrystalPSProd(name);		//PS x CrystalBall (4)
	else if (function_name == "crystaleffprod") getCrystalEffProd(name);       //Eff x CrystalBall (6)
	else if (function_name == "expeffprod") getExpEffProd(name);
	else if (function_name == "doublecb") getDoubleCB(name);
	else if (function_name == "dijetv1") getDijetv1(name);                     //Dijet (EXO-16-020) (3)
	else if (function_name == "dijetv1psprod") getDijetv1PSProd(name);         //PS x Dijet (3)
	else if (function_name == "dijetv2") getDijetv2(name); 			//Dijet X(750) (3)
	else if (function_name == "dijetv2psprod") getDijetv2PSProd(name);		//PS x Dijet (3)
	else if (function_name == "dijetv3") getDijetv3(name);			//version 2 without ps term
	else if (function_name == "dijetv3effprod") getDijetv3EffProd(name);
	else if (function_name == "dijetv3logprod") getDijetv3LogisticProd(name);
	else if (function_name == "dijetv4") getDijetv4(name);			//version 3 with extended
	else if (function_name == "dijetv4logprod") getDijetv4LogisticProd(name);
 	else if (function_name == "expgausexp") getExpGausExp(name);
	else if (function_name == "gausexp") getGausExp(name);
	else if (function_name == "doublegausexp") getDoubleGausExp(name);
	else if (function_name == "triplegausexp") getTripleGausExp(name);
	else if (function_name == "gausexppsprod") getGausExpPSProd(name);		//PS x GausExp (3)
	else if (function_name == "gausexpeffprod") getGausExpEffProd(name);		//Eff x GausExp (5)
	else if (function_name == "expbwexp") getExpBWExp(name);
	else if (function_name == "bukin") getBukin(name);
	else if (function_name == "bukinpsprod") getBukinPSProd(name);		//PS x Bukin (5)
	else if (function_name == "bernstein") getBernstein(name, numCoeffs);
	else if (function_name == "chebychev") getChebychev(name, numCoeffs);
	else if (function_name == "berneffprod") getBernEffProd(name, numCoeffs);
	else if (function_name == "bernefffixprod") getBernEfffixProd(name, numCoeffs);
	else if (function_name == "bernpsprod") getBernPSProd(name, numCoeffs);
	else if (function_name == "chebeffprod") getChebEffProd(name, numCoeffs);
	else if (function_name == "breitwigner") getBreitWigner(name);
	else if (function_name == "relbreitwigner") getRelBreitWigner(name);
	else if (function_name == "quadgausexp") getRooQuadGausExp(name);
	else if (function_name == "mynovosibirsk") getMyNovosibirsk(name);
	else if (function_name == "mynovopsprod") getMyNovoPSProd(name);
	else if (function_name == "extnovosibirsk") getExtNovosibirsk(name);
	else if (function_name == "extnovopsprod") getExtNovoPSProd(name);
	else if (function_name == "extnovoeffprod") getExtNovoEffProd(name);
	else if (function_name == "extnovoefffixprod") getExtNovoEfffixProd(name);
	else if (function_name == "extnovologprod") getExtNovoLogisticProd(name);
	else if (function_name == "extnovoextlogprod") getExtNovoExtLogisticProd(name);
	else if (function_name == "extnovohypertanprod") getExtNovoHypertanProd(name);
	else if (function_name == "supernovosibirsk") getSuperNovosibirsk(name, numCoeffs);
	else if (function_name == "supernovoeffprod") getSuperNovoEffProd(name, numCoeffs);
	else if (function_name == "superdijet") getSuperDiJet(name, numCoeffs);
	else if (function_name == "superdijeteffprod") getSuperDiJetEffProd(name, numCoeffs);
	else if (function_name == "superdijetlinearprod") getSuperDiJetLinearProd(name, numCoeffs);
	else if (function_name == "exppolynom") getExpPolynom(name, numCoeffs);
	else {
		std::stringstream msg;
		msg << "Model '" << function
				<< "' not implemented! Choose one of the following: "
		<< boost::algorithm::join(availableModels_, ", ");
	    throw std::runtime_error(msg.str());
	}
}

void ProbabilityDensityFunctions::getNovosibirsk(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
//	RooRealVar peak = RooRealVar("peak1", "peak", getPeakStart(), 50.0, 500.0, "GeV");
	//	RooRealVar tail = RooRealVar("tail1", "tail", -0.1, -1.0, 1.0);
	RooRealVar peak = RooRealVar("peak1", "peak", getPeakStart(), 50.0, 1000.0, "GeV");	//only for subrange 3
	RooRealVar width = RooRealVar("width1", "width", 50.0, 5.0, var.getMax()/2.0, "GeV");
	RooRealVar tail("tail", "tail", -0.1, -10.0, 10.0); 					//only for subrange 3

//	RooRealVar peak = RooRealVar("peak1", "peak", 461.819, 50.0, 1000.0, "GeV");	//only for subrange 3
//	RooRealVar width = RooRealVar("width1", "width", 58.3072, 5.0, var.getMax()/2.0, "GeV");
//	RooRealVar tail("tail", "tail", -0.93615, -10.0, 10.0); 					//only for subrange 3
	RooNovosibirsk novo(name.c_str(),(name + "_novosibirsk").c_str(),var,peak,width,tail);

	workspace_->import(novo);
}

void ProbabilityDensityFunctions::getNovosibirsk(const std::string& name, const std::string& title, RooRealVar& peak, RooRealVar& width, RooRealVar& tail){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooNovosibirsk novo2(name.c_str(),title.c_str(),var,peak,width,tail);
	workspace_->import(novo2);
}

void ProbabilityDensityFunctions::getNovoPSProd(const std::string& name,const std::string& title, RooRealVar& peak, RooRealVar& width, RooRealVar& tail, RooFormulaVar& PS){
	getNovosibirsk((name + "_novo").c_str(),(name + "_novo").c_str(),peak,width,tail);
	RooNovosibirsk& novo = (RooNovosibirsk&) *workspace_->pdf((name + "_novo").c_str());
	RooEffProd novopsprod(name.c_str(), title.c_str(), novo, PS);
	workspace_->import(novopsprod);
}

void ProbabilityDensityFunctions::getNovoPSProd(const std::string& name){
	getPhaseSpace(name+"_phasespace");
	getNovosibirsk((name+"_novo"));

	RooNovosibirsk& novo 		= (RooNovosibirsk&) *workspace_->pdf((name + "_novo").c_str());
	RooFormulaVar&  phasespace  = (RooFormulaVar&)  *workspace_->function((name + "_phasespace").c_str());
	RooEffProd novoEffProd(name.c_str(),(name + "_novopsprod").c_str(),novo,phasespace);

	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getNovoEffProd(const std::string& name){
	std::string eff_name = name + "_novosibirskeff";
	getEfficiency(eff_name);
	std::string novo_name = name + "_novosibirsk";
	getNovosibirsk(novo_name);

	RooNovosibirsk& novo 		= (RooNovosibirsk&) *workspace_->pdf((novo_name).c_str());
	RooFormulaVar&  eff  = (RooFormulaVar&)  *workspace_->function((eff_name).c_str());
	RooEffProd novoEffProd(name.c_str(),(name + "_novoeffprod").c_str(),novo,eff);
	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getSuperNovosibirsk(const std::string& name, const int& degree){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooArgList arg_list;
	int npars = 0;
	double arg_val = -0.0005;
	for(int i =0;i < degree; ++i) {
//		arg_val /= (i+1);
		std::string var_name = "par" + std::to_string(i);
		if(modify_par_names_) var_name += "_" + name;
		RooRealVar v(var_name.c_str(),var_name.c_str(),arg_val,-10,10); v.setConstant(false);
		workspace_->import(v);
		arg_list.add(*workspace_->var(var_name.c_str()));
		++npars;
		arg_val /= 1000;
	}
	//Novosibirsk PDF
	std::string peak_name = "peak", width_name = "width", tail_name = "tail";
	if(modify_par_names_){
		peak_name += "_" + name;
		width_name+= "_" + name;
		tail_name += "_" + name;
	}
	RooRealVar peak  = RooRealVar( peak_name.c_str(), "peak", getPeakStart(), 50.0, 1000.0, "GeV");
	RooRealVar width = RooRealVar( width_name.c_str(), "width", 50.0, 5.0, var.getMax()/2.0, "GeV");
	RooRealVar tail  = RooRealVar( tail_name.c_str(), "tail", -0.1, -10.0, 10.0);
	RooSuperNovosibirsk pdf(name.c_str(),name.c_str(),var,peak,width,tail,arg_list);

	workspace_->import(pdf);
}

void ProbabilityDensityFunctions::getSuperDiJet(const std::string& name, const int& degree){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooArgList arg_list;
	int npars = 0;
	double start = 0.1;
	for(int i =0;i < degree; ++i) {
		std::string var_name = "par" + std::to_string(i) + "_dijet";
		if(modify_par_names_) var_name += "_" + name;
		//Hardcoded solution for a single fit: main
		if(degree == 2){
			if(i == 0) start = 0.15;
			else if (i==1) start = 0.028;
		}
		RooRealVar v(var_name.c_str(),var_name.c_str(),start,-5,5); v.setConstant(false);
		workspace_->import(v);
		arg_list.add(*workspace_->var(var_name.c_str()));
		++npars;
	}
	//DiJet PDF
	std::string mean_name = "mean_dijet", par1_name = "para_dijet", par2_name = "parb_dijet";
	if(modify_par_names_){
		mean_name += "_" + name;
		par1_name += "_" + name;
		par2_name += "_" + name;
	}
	RooRealVar mean = RooRealVar( mean_name.c_str(), "mean", 130.0, 0., 199.0, "GeV"); //170
	RooRealVar par1 = RooRealVar( par1_name.c_str(), "para", 2.0, 0.0, 50.0);
	RooRealVar par2 = RooRealVar( par2_name.c_str(), "parb",  0.001, 0.001, 1.0);
	RooSuperDiJet pdf(name.c_str(),name.c_str(),var,mean,par1,par2,arg_list);

	workspace_->import(pdf);
}

void ProbabilityDensityFunctions::getExpPolynom(const std::string& name, const int& degree){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooArgList arg_list;
	int npars = 0;
	for(int i =0;i < degree; ++i) {
		std::string var_name = "par" + std::to_string(i) + "_exppolynom";
		if(modify_par_names_) var_name += "_" + name;
		RooRealVar v(var_name.c_str(),var_name.c_str(),std::pow(0.1,i),-5000,5000); v.setConstant(false);
		workspace_->import(v);
		arg_list.add(*workspace_->var(var_name.c_str()));
		++npars;
	}

	RooRealVar para = RooRealVar( "para", "para", 2.0, 0.0, 50.0);
	RooExpPolynom pdf(name.c_str(),name.c_str(),var,para,arg_list);
//	RooRealVar mean = RooRealVar( "mean", "mean", 130.0, 0., 800.0, "GeV"); //170
//	RooExpPolynom pdf(name.c_str(),name.c_str(),var,mean,arg_list);

	workspace_->import(pdf);
}

void ProbabilityDensityFunctions::getSuperDiJetEffProd(const std::string& name, const int& degree){
	std::string super_name = name + "_superdijet" + std::to_string(degree);
	std::string eff_name = name + "_superdijeteff";
	getEfficiency(eff_name);
	getSuperDiJet(super_name, degree);

//	std::string mean_name = "mean_dijet";
//	std::string para_name = "para_dijet";
//	std::string parb_name = "parb_dijet";
//
//	if(modify_par_names_) {
//		mean_name += "_" + super_name;
//		para_name += "_" + super_name;
//		parb_name += "_" + super_name;
//	}
//	workspace_->var(mean_name.c_str())->setVal(167.0);
//	workspace_->var(para_name.c_str())->setVal(11.0);
//	workspace_->var(parb_name.c_str())->setVal(0.1);
	RooSuperDiJet& super = static_cast<RooSuperDiJet&>(*workspace_->pdf(super_name.c_str()));

	RooFormulaVar&  eff  = static_cast<RooFormulaVar&>(*workspace_->function((eff_name).c_str()));
	RooEffProd dijetEffProd(name.c_str(),(name + "_superdijet" + std::to_string(degree) + "_effprod").c_str(),super,eff);
	workspace_->import(dijetEffProd);
}

void ProbabilityDensityFunctions::getSuperDiJetLinearProd(const std::string& name, const int& degree){
	std::string super_name = name + "_superdijet" + std::to_string(degree);
	std::string linear_name = name + "_superdijetlinear";
	getSuperDiJet(super_name, degree);
	RooSuperDiJet& super = static_cast<RooSuperDiJet&>(*workspace_->pdf(super_name.c_str()));
	//Fix mean value to 0
	std::string mean_name = "mean_dijet";
	std::string para_name = "para_dijet";
	std::string parb_name = "parb_dijet";

	if(modify_par_names_) {
		mean_name += "_" + super_name;
		para_name += "_" + super_name;
		parb_name += "_" + super_name;
	}
	workspace_->var( mean_name.c_str())->setVal(0);
	workspace_->var( mean_name.c_str())->setConstant();

	workspace_->var(para_name.c_str())->setVal(8.0);
	workspace_->var(parb_name.c_str())->setVal(0.06);


	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar par_c("par_c", "par_c", 70., 40., 300.);		//control non-zero values when start mbb > 200 GeV
	RooGenericPdf linear(linear_name.c_str(),(linear_name).c_str(),("(-1+par_c*(" + var_name + "/13000.))").c_str(),
	                     RooArgList(var,par_c));

	RooEffProd dijetLinearProd(name.c_str(),(name + "_superdijet" + std::to_string(degree) + "_lnearprod").c_str(),super,linear); //(-1+cx)e^-alnx(1+blnx)
	workspace_->import(dijetLinearProd);
}

void ProbabilityDensityFunctions::getSuperNovoEffProd(const std::string& name, const int& degree){
	std::string supernovo_name = name + "_supernovo" + std::to_string(degree);
	std::string eff_name = name + "_supernovoeff";
	getEfficiency(eff_name);
	getSuperNovosibirsk(supernovo_name, degree);
	RooSuperNovosibirsk& supernovo = static_cast<RooSuperNovosibirsk&>(*workspace_->pdf(supernovo_name.c_str()));

	RooFormulaVar&  eff  = static_cast<RooFormulaVar&>(*workspace_->function((eff_name).c_str()));
	RooEffProd novoEffProd(name.c_str(),(name + "_supernovo" + std::to_string(degree) + "_effprod").c_str(),supernovo,eff);
	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getNovoEfffixProd(const std::string& name){
        std::string novo_name = name + "_novosibirsk";
        getNovosibirsk(novo_name);
        RooNovosibirsk& novo = (RooNovosibirsk&) *workspace_->pdf(novo_name.c_str());

        std::string ps_name = name + "_ps";
        getEfficiencyFix(ps_name);
        RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

        RooEffProd novoEfffixprod(name.c_str(),(name + "_novoefffixprod").c_str(),novo,ps);
        workspace_->import(novoEfffixprod);
}

void ProbabilityDensityFunctions::getNovoPSHighMPol4(const std::string& name){
	std::string novo_name = name + "_novosibirsk";
	getNovosibirsk(novo_name);
	RooNovosibirsk& novo 		= (RooNovosibirsk&) *workspace_->pdf((novo_name).c_str());

	RooRealVar& var = *workspace_->var(var_.c_str());
	RooPhaseSpacePol4 phasespace((name+"_pspol4").c_str(),(name+"_phasespace").c_str(),var);

	RooEffProd novopsprod(name.c_str(),
	                        (name + "_novopsprod").c_str(), novo, phasespace);
	workspace_->import(novopsprod);
}

void ProbabilityDensityFunctions::getCrystalBall(const std::string& name){
	  RooRealVar& var = *workspace_->var(var_.c_str());

	  RooRealVar m0("m0", "m0",getPeakStart(), 50.0, 500.0, "GeV");
	  RooRealVar sigma("sigma", "sigma", 35.0, 10.0, 100.0, "GeV");
	  //RooRealVar alpha("alpha", "alpha", -1.0, -0.1);
	  RooRealVar alpha("alpha", "alpha", -5.0, -0.001);
	  RooRealVar n("n", "n", 20.0, 3.0, 100.0);
	  RooCBShape cb(name.c_str(), (name + "_crystalball").c_str(),
	                var, m0, sigma, alpha, n);
	  workspace_->import(cb);
}

void ProbabilityDensityFunctions::getCrystalPSProd(const std::string& name){
	std::string cb_name = name + "_cb";
	getCrystalBall(cb_name);
	RooCBShape& cb 		= (RooCBShape&) *workspace_->pdf((cb_name).c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd crystalpsprod(name.c_str(),
	                          (name + "_crystalpsprod").c_str(), cb, ps);
	workspace_->import(crystalpsprod);
}

void ProbabilityDensityFunctions::getCrystalEffProd(const std::string& name){
	std::string cb_name = name + "_cb";
	getCrystalBall(cb_name);
	RooCBShape& cb 		= (RooCBShape&) *workspace_->pdf((cb_name).c_str());

	std::string eff_name = name + "_eff";
	getEfficiency(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd crystalEffprod(name.c_str(),
	                            (name + "_crystaleffprod").c_str(), cb, eff);
	workspace_->import(crystalEffprod);
}

void ProbabilityDensityFunctions::getExpEffProd(const std::string& name){
	std::string exp_name = name + "_exp";
	getExp(exp_name);
	RooExponential& exp = (RooExponential&) *workspace_->pdf(exp_name.c_str());

	std::string eff_name = name + "_eff";
	getEfficiency(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd expEffProd(name.c_str(),(name + "_expeffprod").c_str(), exp, eff);
	workspace_->import(expEffProd);
}

void ProbabilityDensityFunctions::getDoubleCB(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean", "mean", getPeakStart(), 50.0, 500.0, "GeV");
	RooRealVar width("width", "width", 35.0, 5.0, 100.0, "GeV");
	RooRealVar alpha1("alpha1", "alpha1", -1.0, -0.1);
	RooRealVar n1("n1", "n1", 20.0, 0.1, 100.0);
	RooRealVar alpha2("alpha2", "alpha2", 0.1, 1.0);
	RooRealVar n2("n2", "n2", 20.0, 0.1, 100.0);
	RooDoubleCB doubleCB(name.c_str(),
	                       (name + "_doublecb").c_str(),
	                       var, mean, width, alpha1, n1, alpha2, n2);
	workspace_->import(doubleCB);
}

void ProbabilityDensityFunctions::getDijetv1(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar par_b("par_b", "par_b", -13, -1000, 1000.);
	RooRealVar par_c("par_c", "par_c", -1.4, -1000, 1000.);
	std::string formula = "pow(" + std::string(var.GetName()) + "/13000. , par_b + par_c * TMath::Log(" + std::string(var.GetName()) + "/13000.))";
	RooGenericPdf dijet(name.c_str(),
	                     (name + "_dijet").c_str(),
	                     formula.c_str(),	//(x)^(P1+P2log(x)) where x = mbb/sqrt(s)
	                     RooArgList(var,par_b,par_c));
	workspace_->import(dijet);
}

void ProbabilityDensityFunctions::getDijetv1PSProd(const std::string& name){
	std::string di_name = name + "_dijet";
	getDijetv1(di_name);
	RooGenericPdf& dijet = (RooGenericPdf&) *workspace_->pdf(di_name.c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd dijetpsprod(name.c_str(),
	                        (name + "_dijetpsprod").c_str(), dijet, ps);
	workspace_->import(dijetpsprod);
}

void ProbabilityDensityFunctions::getDijetv2(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar par_a("par_a", "par_a", 8.0, 0.0, 50.0);
	RooRealVar par_b("par_b", "par_b", 0.06, 0.001, 0.1);
	RooRealVar par_c("par_c", "par_c", 70., 40., 300.);		//control non-zero values when start mbb > 200 GeV
	RooGenericPdf dijet(name.c_str(),				//(-1+cx)e^-alnx(1+blnx)
	                     (name + "_dijet").c_str(),
	                     ("TMath::Exp(-par_a * TMath::Log(" + var_name + "/13000.) - par_b*par_a*(pow(TMath::Log(" + var_name + "/13000.),2)))*(-1+par_c*(" + var_name + "/13000.))").c_str(),
	                     RooArgList(var,par_a,par_b,par_c));
	workspace_->import(dijet);
}

void ProbabilityDensityFunctions::getDijetv2PSProd(const std::string& name){
	std::string di_name = name + "_dijet";
	getDijetv2(di_name);
	RooGenericPdf& dijet = (RooGenericPdf&) *workspace_->pdf(di_name.c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd dijetpsprod(name.c_str(),
	                        (name + "_dijetpsprod").c_str(), dijet, ps);
	workspace_->import(dijetpsprod);
}

void ProbabilityDensityFunctions::getDijetv3(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar par_a("par_a", "par_a", 10.0, 1.0, 50.0);
	RooRealVar par_b("par_b", "par_b", 0.1, 0.001, 1.0);
	RooRealVar mean("mean", "mean", 150.0, 0., 199.0, "GeV");
	RooGenericPdf dijet(name.c_str(),
			     (name + "_dijet").c_str(),
			     ("TMath::Exp(-par_a * TMath::Log((" + var_name + "- mean)/13000.) - par_b*par_a*(pow(TMath::Log((" + var_name + "- mean)/13000.),2)))").c_str(),
			     RooArgList(var,par_a,par_b,mean));
	workspace_->import(dijet);
}

void ProbabilityDensityFunctions::getDijetv3EffProd(const std::string& name){
	std::string di_name = name + "_dijet";
	getDijetv3(di_name);
	RooGenericPdf& dijet = (RooGenericPdf&) *workspace_->pdf(di_name.c_str());

	std::string eff_name = name + "_eff";
	getEfficiency(eff_name);
        RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd dijetEffProd(name.c_str(),(name + "_dijeteffprod").c_str(), dijet, eff);
	workspace_->import(dijetEffProd);
}	

void ProbabilityDensityFunctions::getDijetv3LogisticProd(const std::string& name){
	std::string di_name = name + "_dijet";
	getDijetv3(di_name);
	RooGenericPdf& dijet = (RooGenericPdf&) *workspace_->pdf(di_name.c_str());

        std::string eff_name = name + "_eff";
        getLogistic(eff_name);
        RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

        RooEffProd dijetEffProd(name.c_str(),(name + "_dijeteffprod").c_str(), dijet, eff);
        workspace_->import(dijetEffProd);
}

void ProbabilityDensityFunctions::getDijetv4(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar par_a("par_a", "par_a", 10.0, 1.0, 50.0);
	RooRealVar par_b("par_b", "par_b", 0.1, 0.001, 1.0);
	RooRealVar par_c("par_c", "par_c", -0.001, -0.1, 0.0);
	RooRealVar mean("mean", "mean", 150.0, 0., 199.0, "GeV");
	//RooGenericPdf dijet(name.c_str(),
	//		     (name + "_dijet").c_str(),
	//("TMath::Exp(-(par_a*TMath::Log(("+var_name+"- mean)/13000.)*(1+par_b*TMath::Log(("+var_name+"- mean)/13000.)+par_c*TMath::Sqrt(TMath::Abs(TMath::Log(("+var_name+"- mean)/13000.))) )))").c_str(),
	//		     RooArgList(var,par_a,par_b,par_c,mean));
	RooPolyDijet dijet(name.c_str(),(name + "_dijet").c_str(),var,mean,par_a,par_b,par_c);

	workspace_->import(dijet);
}

void ProbabilityDensityFunctions::getDijetv4LogisticProd(const std::string& name){
        std::string di_name = name + "_dijet";
        getDijetv4(di_name);
        RooPolyDijet& dijet = (RooPolyDijet&) *workspace_->pdf(di_name.c_str());

        std::string eff_name = name + "_eff";
        getLogistic(eff_name);
        RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

        RooEffProd dijetEffProd(name.c_str(),(name + "_dijeteffprod").c_str(), dijet, eff);
        workspace_->import(dijetEffProd);
}

void ProbabilityDensityFunctions::getExpGausExp(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean", "mean", getPeakStart(), 50.0, 500.0, "GeV");
	RooRealVar sigma("sigma", "sigma", 35.0, 5.0, 150.0, "GeV");
	RooRealVar left("left", "left", 0.1, 15.0);
	RooRealVar right("right", "right", 0.1, 15.0);
	RooExpGausExp expGausExp(name.c_str(),
	                           (name + "_expgausexp").c_str(),
	                           var, mean, sigma, left, right);
	workspace_->import(expGausExp);
}

void ProbabilityDensityFunctions::getGausExp(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean", "mean", getPeakStart(), 50.0, 1500.0, "GeV");
	RooRealVar sigma("sigma", "sigma", 35.0, 5.0, 150.0, "GeV");
	RooRealVar tail_shift("tail_shift", "tail_shift", getPeakStart() * 1.2, 50., 1800);
	RooRealVar tail_sigma("tail_sigma", "tail_sigma", 50., 5., 300);
	RooGausExp gausExp(name.c_str(),
	                           (name + "_gausexp").c_str(),
	                           var, mean, sigma, tail_shift,tail_sigma);
	workspace_->import(gausExp);
}

void ProbabilityDensityFunctions::getDoubleGausExp(const std::string& name){
	double mean_tail_sigma = 50., mean_sigmaL = 35., mean_sigmaR = 35.;
	if(getPeakStart() > 400) {
		mean_sigmaL = 90.;
		mean_tail_sigma = 30.;
		mean_sigmaR = 40.;
	}
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean", "mean", getPeakStart(), 50.0, 1500.0, "GeV");
	RooRealVar sigmaL("sigmaL", "sigmaL", mean_sigmaL, 5.0, 400.0, "GeV");
	RooRealVar sigmaR("sigmaR", "sigmaR", mean_sigmaR, 5.0, 400.0, "GeV");
	RooRealVar tail_shift("tail_shift", "tail_shift", getPeakStart() * 1.2, 50., 1800);
	RooRealVar tail_sigma("tail_sigma", "tail_sigma", mean_tail_sigma, 5., 300);
	RooDoubleGausExp res(name.c_str(),name.c_str(),var,mean,sigmaL,sigmaR,tail_shift,tail_sigma);
	workspace_->import(res);
//	std::string ge_name = name + "_gausexp";
//	getGausExp(ge_name);
//	RooGausExp& gausExp = (RooGausExp&) *workspace_->pdf(ge_name.c_str());
//
//	std::string gaus_name = name + "_gaus2";
//	getGaus(gaus_name,"mean2","sigma2");
//	RooGaussian& gaus2 = (RooGaussian&) *workspace_->pdf(gaus_name.c_str());
//
//	RooRealVar g1frac("g1frac","fraction of gauss1",0.5);
//	RooRealVar g2frac("g2frac","fraction of gauss2",0.5);
//	RooAddPdf sum(name.c_str(),(name + "_doublegausexp").c_str(),RooArgList(gausExp,gaus2),RooArgList(g1frac,g2frac));
//	workspace_->import(sum);
}

void ProbabilityDensityFunctions::getTripleGausExp(const std::string& name){
	std::string ge_name = name + "_gausexp";
	getGausExp(ge_name);
	RooGausExp& gausExp = (RooGausExp&) *workspace_->pdf(ge_name.c_str());

	std::string gaus2_name = name + "_gaus2";
	getGaus(gaus2_name,"mean2","sigma2");
	RooGaussian& gaus2 = (RooGaussian&) *workspace_->pdf(gaus2_name.c_str());

	std::string gaus3_name = name + "_gaus3";
	getGaus(gaus3_name,"mean3","sigma3");
	RooGaussian& gaus3 = (RooGaussian&) *workspace_->pdf(gaus3_name.c_str());

	RooRealVar g1frac("g1frac","fraction of gauss1",0.3);
	RooRealVar g2frac("g2frac","fraction of gauss2",0.3);
	RooRealVar g3frac("g3frac","fraction of gauss3",0.3);
	RooAddPdf sum(name.c_str(),(name + "_triplegausexp").c_str(),RooArgList(gausExp,gaus2,gaus3),RooArgList(g1frac,g2frac,g3frac));
	workspace_->import(sum);

}

void ProbabilityDensityFunctions::getGausExpPSProd(const std::string& name){
	std::string ge_name = name + "_gausexp";
	getGausExp(ge_name);
	RooGausExp& gausExp = (RooGausExp&) *workspace_->pdf(ge_name.c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd gausExppsprod(name.c_str(),
	                          (name + "_gausExppsprod").c_str(), gausExp, ps);
	workspace_->import(gausExppsprod);
}

void ProbabilityDensityFunctions::getGausExpEffProd(const std::string& name){
        std::string ge_name = name + "_gausexp";
        getGausExp(ge_name);
        RooGausExp& gausExp = (RooGausExp&) *workspace_->pdf(ge_name.c_str());

        std::string eff_name = name + "_eff";
        getEfficiency(eff_name);
        RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

        RooEffProd gausExpEffprod(name.c_str(),
                                  (name + "_gausExpeffprod").c_str(), gausExp, eff);
        workspace_->import(gausExpEffprod);
}

void ProbabilityDensityFunctions::getExpBWExp(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar peak("peak", "peak", getPeakStart(), 50.0, 500.0, "GeV");
	RooRealVar width("width", "width", 35.0, 5.0, 150.0, "GeV");
	RooRealVar left("left", "left", 0.1, 15.0);
	RooRealVar right("right", "right", 0.1, 15.0);
	RooExpBWExp expBWExp(name.c_str(),
	                       (name + "_expbwexp").c_str(),
	                       var, peak, width, left, right);
	workspace_->import(expBWExp);
}

void ProbabilityDensityFunctions::getBukin(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	double Xp_ = getPeakStart();
	if(getPeakStart() > 600 && getPeakStart() < 700) Xp_ = getPeakStart()-50;
	else if ( getPeakStart() < 900) Xp_ = getPeakStart() + 50;
	if( getPeakStart() > 400 && getPeakStart() < 500) Xp_ = 450;
	RooRealVar Xp("Xp", "Xp", Xp_, 00.0, 1800.0, "GeV");
	RooRealVar sigp("sigp", "sigp", 100,20.0, 200.0, "GeV");
	RooRealVar xi("xi", "xi", 0.2,-10.0, 10.0);
	RooRealVar rho1("rho1", "rho1", -0.07,-20.,0.5);
	RooRealVar rho2("rho2", "rho2", -0.3,-3.,3.);
	RooBukinPdf bukin(name.c_str(),
	                    (name + "_bukin").c_str(),
	                    var, Xp, sigp, xi, rho1, rho2);
	workspace_->import(bukin);
}

void ProbabilityDensityFunctions::getBukinPSProd(const std::string& name){
	std::string bukin_name = name + "_bukin";
	getBukin(bukin_name);
	RooBukinPdf& bukin = (RooBukinPdf&) *workspace_->pdf(bukin_name.c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd bukinpsprod(name.c_str(),
			(name + "_bukinpsprod").c_str(), bukin, ps);
	workspace_->import(bukinpsprod);
}

void ProbabilityDensityFunctions::getBernstein(const std::string& name,const int& numCoeffs){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::unique_ptr<RooArgList> coeffs(getCoefficients_(numCoeffs, "bernstein"));
	RooBernstein bern(name.c_str(),
			(name+"_bernstein").c_str(), var, *coeffs);
	workspace_->import(bern);
}

void ProbabilityDensityFunctions::getChebychev(const std::string& name, const int& numCoeffs){
	if (numCoeffs > 7) {
		throw std::invalid_argument
		("Chebychev polynomials support not more than 7 coefficients.");
	}
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::unique_ptr<RooArgList> coeffs(getCoefficients_(numCoeffs, "chebychev"));
	RooChebychev cheb(name.c_str(),
			(name + "_chebychev").c_str(), var, *coeffs);
	std::cout<<"WTF"<<std::endl;
	workspace_->import(cheb);
}

void ProbabilityDensityFunctions::getChebEffProd(const std::string& name, const int& numCoeffs){
	std::string cheb_name = name + "_cheb";
	getChebychev(cheb_name,numCoeffs);
	RooChebychev& cheb = (RooChebychev&) *workspace_->pdf(cheb_name.c_str());

	std::string eff_name = name + "_eff";
	getEfficiency(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd chebEffProd(name.c_str(),
			(name+"_chebeffprod").c_str(), cheb, eff);
	workspace_->import(chebEffProd);
}

void ProbabilityDensityFunctions::getBernEffProd(const std::string& name, const int& numCoeffs){
	std::string bern_name = name + "_bern";
	getBernstein(bern_name,numCoeffs);
	RooBernstein& bern = (RooBernstein&) *workspace_->pdf(bern_name.c_str());

	std::string eff_name = name + "_eff";
	getEfficiency(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd bernEffProd(name.c_str(),
			(name + "_berneffprod").c_str(), bern, eff);
	workspace_->import(bernEffProd);
}

void ProbabilityDensityFunctions::getBernEfffixProd(const std::string& name,const int& numCoeffs){
        std::string bern_name = name + "_bern";
        getBernstein(bern_name,numCoeffs);
        RooBernstein& bern = (RooBernstein&) *workspace_->pdf(bern_name.c_str());

        std::string ps_name = name + "_ps";
        getEfficiencyFix(ps_name);
        RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

        RooEffProd bernPSProd(name.c_str(),
                        (name+"_bernPSProd").c_str(), bern, ps);
        workspace_->import(bernPSProd);
}

void ProbabilityDensityFunctions::getBernPSProd(const std::string& name,const int& numCoeffs){
	std::string bern_name = name + "_bern";
	getBernstein(bern_name,numCoeffs);
	RooBernstein& bern = (RooBernstein&) *workspace_->pdf(bern_name.c_str());

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd bernPSProd(name.c_str(),
			(name+"_bernPSProd").c_str(), bern, ps);
	workspace_->import(bernPSProd);
}

void ProbabilityDensityFunctions::getBreitWigner(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean","mean",getPeakStart(),200.,1500.,"GeV");
	RooRealVar width("width","width",35.,5.0,300.,"GeV");
	RooBreitWigner bw(name.c_str(),(name + "_breitwigner").c_str(),var,mean,width);
	workspace_->import(bw);
}

void ProbabilityDensityFunctions::getMyNovosibirsk(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar peak("peak", "peak", getPeakStart(), 50.0, 500.0, "GeV");
	RooRealVar width("width", "width", 50.0, 5.0, var.getMax()/2.0, "GeV");
	RooRealVar tail("tail", "tail", -0.1, -1.0, 1.0);
	RooMyNovosibirsk novo(name.c_str(),(name + "_mynovosibirsk").c_str(),var, peak, width, tail);
	workspace_->import(novo);
}

void ProbabilityDensityFunctions::getMyNovoPSProd(const std::string& name){
	std::string novo_name = name + "_novosibirsk";
	getMyNovosibirsk(novo_name);
	RooMyNovosibirsk& novo = (RooMyNovosibirsk&) *workspace_->pdf(novo_name.c_str());	

	std::string ps_name = name + "_ps";
	getPhaseSpace(ps_name);
	RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

	RooEffProd novopsprod(name.c_str(),(name + "_mynovopsprod").c_str(),novo,ps);
	workspace_->import(novopsprod);
}

void ProbabilityDensityFunctions::getExtNovosibirsk(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar peak("peak", "peak", getPeakStart(), 50.0, 500.0, "GeV");
	RooRealVar width("width", "width", 50.0, 5.0, var.getMax()/2.0, "GeV");
	RooRealVar tail("tail", "tail", -0.1, -1.0, 1.0);
	RooRealVar par4("par4", "par4", -0.0001, -1.0, 1.0);

	//RooRealVar peak("peak", "peak", getPeakStart(), 50.0, 5000.0, "GeV");
	//RooRealVar width("width", "width", 50.0, 0.0, 1000.0, "GeV");
	//RooRealVar tail("tail", "tail", -0.1, -5.0, 5.0);
	//RooRealVar par4("par4", "par4", -0.0001, -100.0, 100.0);

//	RooRealVar peak("peak", "peak", 274.299, 50.0, 500.0, "GeV");
//	RooRealVar width("width", "width", 63.8268, 5.0, var.getMax()/2.0, "GeV");
//	RooRealVar tail("tail", "tail", -0.408301, -1.0, 1.0);
//	RooRealVar par4("par4", "par4", -0.000725092, -1.0, 1.0);

	RooExtendNovosibirsk novo(name.c_str(),(name + "_extnovosibirsk").c_str(),var,peak, width, tail, par4);

	workspace_->import(novo);
}

void ProbabilityDensityFunctions::getExtNovoPSProd(const std::string& name){
        std::string novo_name = name + "_extnovosibirsk";
        getExtNovosibirsk(novo_name);
        RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());

        std::string ps_name = name + "_ps";
        getPhaseSpace(ps_name);
        RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

        RooEffProd novopsprod(name.c_str(),(name + "_extnovopsprod").c_str(),novo,ps);
        workspace_->import(novopsprod);
}

void ProbabilityDensityFunctions::getExtNovoEffProd(const std::string& name){
        std::string novo_name = name + "_extnovosibirsk";
        getExtNovosibirsk(novo_name);
        RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());

        std::string eff_name = name + "_eff";
        getEfficiency(eff_name);
        RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());	

	RooEffProd novoEffProd(name.c_str(),(name + "_extnovoeffprod").c_str(),novo,eff);
	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getExtNovoEfffixProd(const std::string& name){
        std::string novo_name = name + "_extnovosibirsk";
        getExtNovosibirsk(novo_name);
        RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());

        std::string ps_name = name + "_ps";
        getEfficiencyFix(ps_name);
        RooFormulaVar& ps = (RooFormulaVar&) *workspace_->function(ps_name.c_str());

        RooEffProd novoEfffixprod(name.c_str(),(name + "_extnovopsprod").c_str(),novo,ps);
        workspace_->import(novoEfffixprod);
}

void ProbabilityDensityFunctions::getExtNovoLogisticProd(const std::string& name){
	std::string novo_name = name + "_extnovosibirsk";
	getExtNovosibirsk(novo_name);
	RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());

	std::string eff_name = name + "_eff";
	getLogistic(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());
	
	RooEffProd novoEffProd(name.c_str(),(name + "_extnovoeffprod").c_str(),novo,eff);
	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getExtNovoExtLogisticProd(const std::string& name){
	std::string novo_name = name + "_extnovosibirsk";
	getExtNovosibirsk(novo_name);
	RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());
	
	std::string eff_name = name + "_eff";
	getLogistic_ext1(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());

	RooEffProd novoEffProd(name.c_str(),(name + "_extnovoeffprod").c_str(),novo,eff);
	workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getExtNovoHypertanProd(const std::string& name){
	std::string novo_name = name + "_extnovosibirsk";
	getExtNovosibirsk(novo_name);
	RooExtendNovosibirsk& novo = (RooExtendNovosibirsk&) *workspace_->pdf(novo_name.c_str());
	
	std::string eff_name = name + "_eff";
	getHyperbolictan(eff_name);
	RooFormulaVar& eff = (RooFormulaVar&) *workspace_->function(eff_name.c_str());
	
	RooEffProd novoEffProd(name.c_str(),(name + "_extnovoeffprod").c_str(),novo,eff);
        workspace_->import(novoEffProd);
}

void ProbabilityDensityFunctions::getRelBreitWigner(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean("mean","mean",getPeakStart(),200.,1500.,"GeV");
	RooRealVar width("width","width",120,5.0,300.,"GeV");
//	RooRealVar sigma("sigma","sigma",0.1,5.0,300.,"GeV");
	RooRelBreitWigner bw(name.c_str(),(name + "_relbreitwigner").c_str(),var,mean,width);
	workspace_->import(bw);
}

void ProbabilityDensityFunctions::getRooQuadGausExp(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	double peak_pos = getPeakStart() * 0.9;
	RooRealVar mean("mean","mean",peak_pos,200,1500,"GeV");
	RooRealVar sigmaL1("sigmaL1", "sigmaL1", 0.11*getPeakStart(), 20., 800.0, "GeV");
	RooRealVar sigmaL2("sigmaL2", "sigmaL2", 0.21*getPeakStart(), 20., 800.0, "GeV");
	RooRealVar sigmaR1("sigmaR1", "sigmaR1", 0.11*getPeakStart(), 10., 800.0, "GeV");
	RooRealVar sigmaR2("sigmaR2", "sigmaR2", 0.11*getPeakStart(), 20., 800.0, "GeV");
	RooRealVar tail_shift("tail_shift", "tail_shift", 1.21*getPeakStart(), 200.0, 1500.0, "GeV");
	RooRealVar tail_sigma("tail_sigma", "tail_sigma", 0.51*getPeakStart(), 0.5, 1500.0, "GeV");
	RooRealVar norm_g1("norm_g1", "norm_g1", 0.5, 0, 1);
	RooRealVar norm_g2("norm_g2", "norm_g2", 0.5, 0, 1);
	RooQuadGausExp quadgexp(name.c_str(),(name + "_quadgexp").c_str(),var,mean,sigmaL1,sigmaL2,sigmaR1,sigmaR2,tail_shift,tail_sigma,norm_g1,norm_g2);
	workspace_->import(quadgexp);
}

void ProbabilityDensityFunctions::getPhaseSpace(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooFormulaVar phasespace((name).c_str(),
				  ("(1.326 / (1 + (2.287e+03 * TMath::Exp(-3.331e-02 * " + var_name + ")))) - (1.326-1.)").c_str(),RooArgList(var));
	workspace_->import(phasespace);
}

void ProbabilityDensityFunctions::getEfficiency(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	std::string slope_name  = "slope_novoeff";
	std::string turnon_name = "turnon_novoeff";
	if(modify_par_names_) {
		slope_name += "_" + name;
		turnon_name+= "_" + name;
	}
	RooRealVar slope_novoeff ( (slope_name).c_str(), "slope_novoeff", 0.01, 0.0, 0.1);
	RooRealVar turnon_novoeff( (turnon_name).c_str(), "turnon_novoeff",var.getMin()+ 5.0, var.getMin(), var.getMax());
	RooFormulaVar eff((name).c_str(),("0.5*(TMath::Erf(" + slope_name + "*(" + var_name + "-" + turnon_name +")) + 1)").c_str(),
	                    RooArgSet(var, slope_novoeff, turnon_novoeff));
	
	//slope_novoeff.setConstant(kTRUE);
	//turnon_novoeff.setConstant(kTRUE);
	workspace_->import(eff);
}

void ProbabilityDensityFunctions::getEfficiencyFix(const std::string& name){
        RooRealVar& var = *workspace_->var(var_.c_str());
        std::string var_name = var.GetName();
        RooFormulaVar phasespace((name).c_str(),
                                  ("0.5*(TMath::Erf(0.0159765*(" + var_name + "- 223.027)) + 1)").c_str(),RooArgList(var));	// for G4 ExtNovo new WP
				  //("0.5*(TMath::Erf(0.0161701*(" + var_name + "- 218.191)) + 1)").c_str(),RooArgList(var));	// for old WP
				  //("0.5*(TMath::Erf(0.0166748*(" + var_name + "- 245.507)) + 1)").c_str(),RooArgList(var));	// for G7 Novosibirsk
				  //("0.5*(TMath::Erf(0.0127066*(" + var_name + "- 280.87)) + 1)").c_str(),RooArgList(var));	// QCD bEnriched

        workspace_->import(phasespace);
}

void ProbabilityDensityFunctions::getLogistic(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
        std::string var_name = var.GetName();
	RooRealVar slope("slope", "slope", 0.1, 0.0, 0.5);
	RooRealVar turnon("turnon", "turnon", var.getMin()+ 5.0, var.getMin()- 50., var.getMax());
	RooFormulaVar eff((name).c_str(),("1/(1+TMath::Exp(-slope*(" + var_name + "-turnon)))").c_str(),
			    RooArgSet(var, slope, turnon));
	workspace_->import(eff);
}

void ProbabilityDensityFunctions::getLogistic_ext1(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar slope("slope", "slope", 0.1, 0.0, 0.5);
	RooRealVar turnon("turnon", "turnon", var.getMin()+ 5.0, var.getMin()- 50., var.getMax());
	RooRealVar par3("par3", "par3", 0.5, 0.0, 1.0);
	RooFormulaVar eff((name).c_str(),("1/(1 + TMath::Exp(-slope*("+ var_name +"-turnon)) + par3*TMath::Exp(-slope*("+ var_name +"-turnon))*TMath::Exp(-slope*("+ var_name +"-turnon)) )").c_str(),
			    RooArgSet(var, slope, turnon, par3));
	workspace_->import(eff);
}

void ProbabilityDensityFunctions::getHyperbolictan(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
        std::string var_name = var.GetName();
	RooRealVar slope("slope", "slope", 0.1, 0.0, 1.0);		
	RooRealVar turnon("turnon", "turnon", var.getMin()+ 5.0, var.getMin()- 50., var.getMax());
	RooFormulaVar eff((name).c_str(),("0.5*(1 + TMath::TanH(slope*("+ var_name +"-turnon)) )").c_str(),
			    RooArgSet(var, slope, turnon));
	workspace_->import(eff);
}

/*
void ProbabilityDensityFunctions::getHyperbolictan_ext1(const std::string& name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	std::string var_name = var.GetName();
	RooRealVar slope("slope", "slope", 0.1, 0.0, 0.5);
	RooRealVar turnon("turnon", "turnon", var.getMin()+ 5.0, var.getMin()- 50., var.getMax());
	RooRealVar par3("par3", "par3", 0.5, 0.0, 1.0);
	RooFormulaVar eff((name).c_str(),("(1/(2*(1+par3)) * (1 + TMath::TanH(slope*("+ var_name +"-turnon)) + par3*(1+TMath::TanH(slope*("+ var_name +"-turnon)))*(1+TMath::TanH(slope*("+ var_name +"-turnon))) )").c_str(),
			   RooArgSet(var, slope, turnon, par3));
	workspace_->import(eff);
}
*/

void ProbabilityDensityFunctions::getExp(const std::string& name,const std::string& name1){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar tau(paramName_(name1,"tau"), "tau", -0.1);
	tau.setConstant(false);
	RooExponential exp(name.c_str(),
	                     (name+"_exp").c_str(), var, tau);
	workspace_->import(exp);
}

void ProbabilityDensityFunctions::getGaus(const std::string& name,const std::string& mean_name, const std::string& sigma_name){
	RooRealVar& var = *workspace_->var(var_.c_str());
	RooRealVar mean(paramName_(mean_name,"mean"), "mean2", getPeakStart(), 100.0, 1500.0, "GeV");
	RooRealVar sigma(paramName_(sigma_name,"sigma"), "sigma2", 35.0, 5.0, 300.0, "GeV");
	RooGaussian gaus(name.c_str(),"gaussian PDF",var,mean,sigma);
	workspace_->import(gaus);
}

std::unique_ptr<RooArgList> ProbabilityDensityFunctions::getCoefficients_(const int numCoeffs, const std::string& name) {
  std::unique_ptr<RooArgList> coefficients
    (new RooArgList((name+"_coefficients").c_str()));
  //Hardcoded solution for parameters boundaries:
  double max = 15.;	//max value
  if(numCoeffs > 6) max = 25;
  double start = max /5.;	//starting value
  for (int c = 0; c < numCoeffs; ++c) {
	  start = max /5.;
    std::string id(Form((name+"_coefficient_%02d").c_str(), c));
    if(c > 2) start /= 10.;
    std::unique_ptr<RooRealVar> coefficient
      (new RooRealVar(id.c_str(), id.c_str(), start, 0.0, max));
    //hardcoded soluten to set higher limits.
    coefficient->setRange(0, max);
    coefficients->addClone(*coefficient);
  }
  return coefficients;
}

const std::vector<std::string> ProbabilityDensityFunctions::availableModels_ =
  {"novosibirsk",
   "novopsprod",
   "novoeffprod",
   "novoefffixprod",
   "novopshighMpol4",
   "crystalball",
   "crystalpsprod",
   "crystaleffprod",
   "cbeffprod",
   "expeffprod",
   "doublecb",
   "dijetv1",
   "dijetv1psprod",
   "dijetv2",
   "dijetv2psprod",
   "dijetv3",
   "dijetv3effprod",
   "dijetv3logprod",
   "dijetv4",
   "dijetv4logprod",
   "expgausexp",
   "gausexp",
   "doublegausexp",
   "triplegausexp",
   "gausexppsprod",
   "gausexpeffprod",
   "expbwexp",
   "bukin",
   "bukinpsprod",
   "bernstein",
   "chebychev",
   "berneffprod",
   "bernefffixprod",
   "bernpsprod",
   "chebeffprod",
   "breitwigner",
   "mynovosibirsk",
   "mynovopsprod",
   "extnovosibirsk",
   "extnovopsprod",
   "extnovoeffprod",
   "extnovoefffixprod",
   "extnovologprod",
   "extnovoextlogprod",
   "extnovohypertanprod",
   "relbreitwigner",
   "doublegausexp",
   "quadgausexp",
   "supernovosibirsk",
   "supernovoeffprod",
   "superdijet",
   "superdijeteffprod",
   "superdijetlinearprod",
   "exppolynom"};
