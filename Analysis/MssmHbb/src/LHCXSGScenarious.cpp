/*
 * LHCXSGScenarious.cpp
 *
 *  Created on: 18 Aug 2017
 *      Author: shevchen
 */

#include "Analysis/MssmHbb/interface/LHCXSGScenarious.h"

namespace analysis {
namespace mssmhbb {

std::string AvailableScenariosToString(AvailableScenarios scenario){
	/*
	 * Translate scenario code to a string
	 */
	if(scenario == MHMODP_200) return "mhmodp_200";
	else if(scenario == LIGHT_STOP) return "light_stop";
	else if(scenario == LIGHT_STAU) return "light_stau";
	else if(scenario == HMSSM) return "hMSSM";
	else if(scenario == TAU_PHOBIC) return "tau_phobic";
	else if(scenario == TYPE2) return "type2";
	else if(scenario == FLIPPED) return "flipped";
	else throw std::logic_error("ERROR at LHCXSGScenarious::AvailableScenariosToString. Wrong Scenario ");
}

AvailableScenarios AvailableScenariosFromString(const std::string& scenario_string){
	/*
	 * Translate string to scenario
	 */
	 
	if(findStrings(scenario_string,"mhmodp_200")) return MHMODP_200;
	else if(findStrings(scenario_string,"light_stop")) return LIGHT_STOP;
	else if(findStrings(scenario_string,"light_stau")) return LIGHT_STAU;
	else if(findStrings(scenario_string,"hMSSM")) return HMSSM;
	else if(findStrings(scenario_string,"tau_phobic")) return TAU_PHOBIC;
	else if(scenario_string == "type1") return TYPE1;
	else if(findStrings(scenario_string,"type2")) return TYPE2;
	else if(findStrings(scenario_string,"flipped")) return FLIPPED;
	else if(findStrings(scenario_string,"lepton_specific")) return LEPTON_SPECIFIC;
	else throw std::logic_error("ERROR at LHCXSGScenarious::AvailableScenariosFromString. No scenario: " + scenario_string);
}

std::unique_ptr<Scenario> Scenario::Create(AvailableScenarios scenario){
	/*
	* Realisation of a factory pattern
	*/
	if(scenario == MHMODP_200) return std::unique_ptr<mhmodp_200>(new mhmodp_200());
	else if(scenario == LIGHT_STOP) return std::unique_ptr<light_stop>(new light_stop());
	else if(scenario == LIGHT_STAU) return std::unique_ptr<light_stau>(new light_stau());
	else if(scenario == HMSSM) return std::unique_ptr<hMSSM>(new hMSSM());
	else if(scenario == TAU_PHOBIC) return std::unique_ptr<tau_phobic>(new tau_phobic());
	else if(scenario == TYPE2) return std::unique_ptr<type2>(new type2());
	else if(scenario == FLIPPED) return std::unique_ptr<flipped>(new flipped());
	else throw std::logic_error("ERROR at LHCXSGScenarious::Create. Wrong Scenario: " + AvailableScenariosToString(scenario));
}

TGraph mhmodp_200::getPreviousResults() const{
	/*
	* Previous results from 8 TeV according to AN
	*/
	std::array<double,8> x = {{100,140,160,200,300,350,400,500}};
	std::array<double,8> y = {{19.8,19.6,18.0,19.2,24.8,30.0,36.5,53.5}};
	TGraph gr(8,&x[0],&y[0]);
	gr.SetLineStyle(2);
	gr.SetLineColor(kBlack);
	return gr;
}

TText mhmodp_200::getPreviousResultsLabel() const{
	/*
	* Return a TText box with a sign for a previous results
	* when available
	*/
	TText text(400,37.5,"7+8 TeV expected");
	text.SetTextSize(0.02);
	text.SetTextAngle(58);
	return text;
}

TGraph light_stop::getPreviousResults() const{
	/*
	* Previous results from 8 TeV according to AN
	*/
	const int npoints = 7;
	std::array<double,npoints> x = {{100,140,160,200,300,350,400}};
	std::array<double,npoints> y = {{24.7,24.3,22.1,24.2,32.7,41.4,58.6}};
	TGraph gr(npoints,&x[0],&y[0]);
	gr.SetLineStyle(2);
	gr.SetLineColor(kBlack);
	return gr;
}

TText light_stop::getPreviousResultsLabel() const{
	/*
	* Return a TText box with a sign for a previous results
	* when available
	*/
	TText text(350,42.0,"7+8 TeV expected");
	text.SetTextSize(0.02);
	double angle = atan( 10 *(58.6 - 41.4)/(400-350)) * 180. / TMath::Pi();
	text.SetTextAngle(angle);
	return text;
}

TGraph light_stau::getPreviousResults() const{
	/*
	* Previous results from 8 TeV according to AN
	*/
	const int npoints = 7;
	std::array<double,npoints> x = {{100,140,160,200,300,350,400}};
	std::array<double,npoints> y = {{22.2,22.,20.,21.4,27.6,33.1,42.4}};
	TGraph gr(npoints,&x[0],&y[0]);
	gr.SetLineStyle(2);
	gr.SetLineColor(kBlack);
	return gr;
}

TText light_stau::getPreviousResultsLabel() const{
	/*
	* Return a TText box with a sign for a previous results
	* when available
	*/
	TText text(350,33.6,"7+8 TeV expected");
	text.SetTextSize(0.02);
	double angle = atan( 10 *(42.4 - 33.1)/(400-350)) * 180. / TMath::Pi();
	text.SetTextAngle(angle);
	return text;
}

TGraph hMSSM::getPreviousResults() const{
	/*
	* Previous results from 8 TeV according to AN
	*/
	TGraph gr;
	return gr;
}

TText hMSSM::getPreviousResultsLabel() const{
	/*
	* Return a TText box with a sign for a previous results
	* when available
	*/
	TText text;
	return text;
}

TGraph tau_phobic::getPreviousResults() const{
	/*
	* Previous results from 8 TeV according to AN
	*/
	const int npoints = 5;
	std::array<double,npoints> x = {{100,140,160,200,300}};
	std::array<double,npoints> y = {{25.2,33.0,25.3,28.9,43.8}};
	TGraph gr(npoints,&x[0],&y[0]);
	gr.SetLineStyle(2);
	gr.SetLineColor(kBlack);
	return gr;
}

TText tau_phobic::getPreviousResultsLabel() const{
	/*
	* Return a TText box with a sign for a previous results
	* when available
	*/
	TText text(300,44.3,"7+8 TeV expected");
	text.SetTextSize(0.02);
	double angle = atan( 5 *(43.8 - 28.9)/(400-350)) * 180. / TMath::Pi();
	text.SetTextAngle(angle);
	return text;
}

std::vector<TGraph> type1::getAtlasResults(const std::string& var) const{
	/*
	 * ATLAS results according to the scenario.
	 */
	std::vector<TGraph> gr;
	return gr;
}

std::vector<TGraph> type2::getAtlasResults(const std::string& var) const{
	/*
	 * ATLAS results according to the scenario and VAR (dedicated to the specific axis)
	 */
	std::vector<TGraph> grs;
	if(var=="z"||var=="Z"){
		std::vector<std::pair<double,double>> lower_left= {
			std::make_pair(220,0.5),
			std::make_pair(225,1.21),
			std::make_pair(230,1.68),
			std::make_pair(235,1.95),
			std::make_pair(240,2.18),
			std::make_pair(245,2.27),
			std::make_pair(250,2.36),
			std::make_pair(255,2.50),
			std::make_pair(260,2.59),
			std::make_pair(265,2.64),
			std::make_pair(270,2.68),
			std::make_pair(275,2.77),
			std::make_pair(280,2.82),
			std::make_pair(285,2.86),
			std::make_pair(290,2.91),
			std::make_pair(295,2.95),
			std::make_pair(300,3.00),
			std::make_pair(305,3.10),
			std::make_pair(310,3.12),
			std::make_pair(315,3.19),
			std::make_pair(320,3.25),
			std::make_pair(325,3.38),
			std::make_pair(330,3.50),
			std::make_pair(335,3.63),
			std::make_pair(337,3.69),
			std::make_pair(340,3.50),
			std::make_pair(345,3.13),
			std::make_pair(350,0.5)
	};

	std::vector<double> cos_lower_left, tan_lower_left;
	for(const auto& val : lower_left) {
		cos_lower_left.push_back(val.first);
		tan_lower_left.push_back(val.second);
	}

	TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());
	grs.push_back(fl_lower_left);
	}
	else if (var=="x"||var=="X"){
		//LEFT
		std::vector<std::pair<double,double>> curve_left= {
			std::make_pair(-0.9,50), // ADDITIONAL POINT TO DRAW

			std::make_pair(-0.35,50),
			std::make_pair(-0.32,32),
			std::make_pair(-0.35,15),
			std::make_pair(-0.4,10),
			std::make_pair(-0.42,8),
			std::make_pair(-0.4,7),
			std::make_pair(-0.35,6.2),
			std::make_pair(-0.32,5.9),
			std::make_pair(-0.3,5.7),
			std::make_pair(-0.25,5.2),
			std::make_pair(-0.2,4.8),
			std::make_pair(-0.15,4.1),
			std::make_pair(-0.1,3.5),
			std::make_pair(-0.05,2.5),
			std::make_pair(0.0,0.5),

			std::make_pair(-0.9,0.5),	//ADDITIONAL POINT TO DRAW
		};

		std::vector<double> cos_left, tan_left;
		for(const auto& val : curve_left) {
			cos_left.push_back(val.first);
			tan_left.push_back(val.second);
		}

		TGraph fl_left(cos_left.size(),cos_left.data(),tan_left.data());
		grs.push_back(fl_left);
		
		//LOWER RIGHT
			std::vector<std::pair<double,double>> curve_lright = {
			std::make_pair(0.0,0.5),		//ADDITIONAL POINT TO DRAW

			std::make_pair(0.05,2.4),
			std::make_pair(0.1,3.0),
			std::make_pair(0.15,3.1),
			std::make_pair(0.2,2.9),
			std::make_pair(0.25,2.75),
			std::make_pair(0.3,2.4),
			std::make_pair(0.35,2.1),
			std::make_pair(0.4,1.9),
			std::make_pair(0.45,1.7),
			std::make_pair(0.5,1.6),
			std::make_pair(0.55,1.5),
			std::make_pair(0.6,1.32),
			std::make_pair(0.65,1.16),
			std::make_pair(0.7,0.97),
			std::make_pair(0.75,0.85),
			std::make_pair(0.8,0.7),
			std::make_pair(0.85,0.6),
			std::make_pair(0.88,0.5),

		};

		std::vector<double> cos_lright, tan_lright;
		for(const auto& val : curve_lright) {
			cos_lright.push_back(val.first);
			tan_lright.push_back(val.second);
		}

		TGraph fl_lright(cos_lright.size(),cos_lright.data(),tan_lright.data());
		grs.push_back(fl_lright);
		
		//UPPER RIGHT
	std::vector<std::pair<double,double>> curve_hright= {
			std::make_pair(0.9,50),		//ADDITIONAL POINT TO DRAW

			std::make_pair(0.35,50),
			std::make_pair(0.31,32),
			std::make_pair(0.35,15),
			std::make_pair(0.4,12),
			std::make_pair(0.45,7.5),
			std::make_pair(0.4,5.7),
			std::make_pair(0.37,4.5),
			std::make_pair(0.4,3.5),
			std::make_pair(0.45,2.7),
			std::make_pair(0.5,2.3),
			std::make_pair(0.55,1.7),
			std::make_pair(0.6,1.53),
			std::make_pair(0.65,1.26),
			std::make_pair(0.7,1.11),
			std::make_pair(0.75,0.92),
			std::make_pair(0.8,0.8),
			std::make_pair(0.85,0.65),
			std::make_pair(0.9,0.5)


	};

	std::vector<double> cos_hright, tan_hright;
	for(const auto& val : curve_hright) {
		cos_hright.push_back(val.first);
		tan_hright.push_back(val.second);
	}

	TGraph fl_hright(cos_hright.size(),cos_hright.data(),tan_hright.data());
	grs.push_back(fl_hright);
		
	
	}
	return grs;
}

std::vector<TGraph> flipped::getAtlasResults(const std::string& var) const{
	/*
	 * ATLAS results according to the scenario and VAR (dedicated to the specific axis)
	 */
	std::vector<TGraph> grs;
	if(var=="z"||var=="Z"){
	std::vector<std::pair<double,double>> lower_left= {
			std::make_pair(220,0.5),
			std::make_pair(225,1),
			std::make_pair(230,1.8),
			std::make_pair(235,2.),
			std::make_pair(240,2.3),
			std::make_pair(245,2.5),
			std::make_pair(250,2.7),
			std::make_pair(255,2.8),
			std::make_pair(260,2.9),
			std::make_pair(265,2.95),
			std::make_pair(270,3.),
			std::make_pair(275,3.05),
			std::make_pair(280,3.1),
			std::make_pair(285,3.15),
			std::make_pair(290,3.2),
			std::make_pair(295,3.3),
			std::make_pair(300,3.4),
			std::make_pair(305,3.5),
			std::make_pair(310,3.6),
			std::make_pair(315,3.7),
			std::make_pair(320,3.8),
			std::make_pair(325,3.9),
			std::make_pair(330,3.95),
			std::make_pair(335,4),
			std::make_pair(340,4),
			std::make_pair(345,3.7),
			std::make_pair(347.5,3),
			std::make_pair(350,0.5),

				//ADDITIONAL POINT TO DRAW!!!
//				std::make_pair(-0.95,0.5)
		};

		std::vector<double> cos_lower_left, tan_lower_left;
		for(const auto& val : lower_left) {
			cos_lower_left.push_back(val.first);
			tan_lower_left.push_back(val.second);
		}
	TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());

	grs.push_back(fl_lower_left);
	}
	else if(var=="x"||var=="X"){
		//LOWER LEFT
	std::vector<std::pair<double,double>> lower_left= {
			std::make_pair(-0.02,0.5),
			std::make_pair(-0.021,0.6),
			std::make_pair(-0.022,0.7),
			std::make_pair(-0.023,0.8),
				std::make_pair(-0.024,0.9),
				std::make_pair(-0.025,1.),
				std::make_pair(-0.026,1.3),
				std::make_pair(-0.05,2.),
				std::make_pair(-0.054,2.5),
				std::make_pair(-0.1,3.),
				std::make_pair(-0.15,3.5),
				std::make_pair(-0.21,4),
				std::make_pair(-0.3,4.5),
				std::make_pair(-0.45,5),
				std::make_pair(-0.64,5.5),
				std::make_pair(-0.95,6),

				//ADDITIONAL POINT TO DRAW!!!
				std::make_pair(-0.95,0.5)
		};

		std::vector<double> cos_lower_left, tan_lower_left;
		for(const auto& val : lower_left) {
			cos_lower_left.push_back(val.first);
			tan_lower_left.push_back(val.second);
		}

		TGraph fl_lower_left(cos_lower_left.size(),cos_lower_left.data(),tan_lower_left.data());
		grs.push_back(fl_lower_left);
		
		//UPPER LEFT
	std::vector<std::pair<double,double> > upper_left = {
			std::make_pair(-0.95,8.5),
			std::make_pair(-0.75,9.),
			std::make_pair(-0.7,9.5),
			std::make_pair(-0.675,10.),
			std::make_pair(-0.65,11),
			std::make_pair(-0.6,12),
			std::make_pair(-0.55,12.5),
			std::make_pair(-0.5,13.5),
			std::make_pair(-0.45,14),
			std::make_pair(-0.4,18),
			std::make_pair(-0.375,25.),
			std::make_pair(-0.36,40.),
			std::make_pair(-0.35,50.),

			//ADDITIONAL POINT TO DRAW!!!
			std::make_pair(-0.95,50)
	};

	std::vector<double> cos_upper_left, tan_upper_left;
	for(const auto& val : upper_left) {
		cos_upper_left.push_back(val.first);
		tan_upper_left.push_back(val.second);
	}

	TGraph fl_upper_left(cos_upper_left.size(),cos_upper_left.data(),tan_upper_left.data());
	grs.push_back(fl_upper_left);
	
	//LOWER RIGHT
	std::vector<std::pair<double,double> > lower_right= {
			std::make_pair(0.95,6),
			std::make_pair(0.725,5.5),
			std::make_pair(0.575,5),
			std::make_pair(0.5,4.5),
			std::make_pair(0.425,4),
			std::make_pair(0.4,3.7),
			std::make_pair(0.35,3.5),
			std::make_pair(0.3,3.5),
			std::make_pair(0.275,3.5),
			std::make_pair(0.25,3.4),
			std::make_pair(0.15,3.2),
			std::make_pair(0.12,3.),
			std::make_pair(0.1,2.7),
			std::make_pair(0.07,2.5),
			std::make_pair(0.05,2.),
			std::make_pair(0.026,1.3),
			std::make_pair(0.025,1.),
			std::make_pair(0.024,0.9),
			std::make_pair(0.023,0.8),
			std::make_pair(0.022,0.7),
			std::make_pair(0.021,0.6),
			std::make_pair(0.02,0.5),

			//ADDITIONAL POINT TO DRAW!!!
			std::make_pair(0.95,0.5)
	};

	std::vector<double> cos_lower_right, tan_lower_right;
	for(const auto& val : lower_right) {
		cos_lower_right.push_back(val.first);
		tan_lower_right.push_back(val.second);
	}

	TGraph fl_lower_right(cos_lower_right.size(),cos_lower_right.data(),tan_lower_right.data());
	grs.push_back(fl_lower_right);
	
	//UPPER RIGHT
	std::vector<std::pair<double,double> > upper_right = {
			std::make_pair(0.35,50),
			std::make_pair(0.36,40.),
			std::make_pair(0.375,25.),
			std::make_pair(0.4,20),
			std::make_pair(0.45,14),
			std::make_pair(0.5,13.5),
			std::make_pair(0.55,12.5),
			std::make_pair(0.6,12),
			std::make_pair(0.65,11),
			std::make_pair(0.675,10.),
			std::make_pair(0.7,9.7),
			std::make_pair(0.75,9.),
			std::make_pair(0.95,8.5),

			//ADDITIONAL POINT TO DRAW!!!
			std::make_pair(0.95,50)
	};

	std::vector<double> cos_upper_right, tan_upper_right;
	for(const auto& val : upper_right) {
		cos_upper_right.push_back(val.first);
		tan_upper_right.push_back(val.second);
	}

	TGraph fl_upper_right(cos_upper_right.size(),cos_upper_right.data(),tan_upper_right.data());
	grs.push_back(fl_upper_right);
		
	}

	return grs;
}

} /* namespace mssmhbb */
} /* namespace analysis */
