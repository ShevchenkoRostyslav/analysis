/*
 * DrawBiasVsM.cpp
 *
 *  Created on: 29 марта 2017 г.
 *      Author: rostyslav
 */

#include <map>
#include <string>
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"

using namespace std;

class Measurement{
public:
	Measurement(const double& central, const double& err) : central_(central), error_(err) {};
	Measurement(const double& central) : central_(central), error_(0) {};
	Measurement(const Measurement& ) = default;
	double getCentral() const { return central_;}
	void setCentral(double central) {central_ = central;}
	double getError() const {return error_;}
	void setError(double error) {error_ = error;}
private:
	double central_;
	double error_;
};

TGraphErrors MakeTGraphErrors(const map<int,Measurement>& my_map);
void FindMaxBias(const map<int,Measurement>& my_map);

int DrawBiasVsM(){
	/*
	 * Macro to draw Bias values as a function of M
	 *
	 * Bias values are hardcoded.
	 */

	HbbStyle style;
	style.set(PRIVATE);
	gStyle->SetMarkerSize(1.5);

	map<int,Measurement> bias_novo_dijet_sr1_r0 =
	{
			{300, Measurement(0.4)},
			{350, Measurement(1.06)},
			{400, Measurement(-0.86)},
//			{500, Measurement(0.21)}
	};

	map<int,Measurement> bias_novo_dijet_sr1_r1 =
	{
			{300, Measurement(0.4)},
			{350, Measurement(1.08)},
			{400, Measurement(-0.87)},
//			{500, Measurement(0.23)}
	};

	map<int,Measurement> bias_novo_dijet_sr1_r2 =
	{
			{300, Measurement(0.4)},
			{350, Measurement(1.04)},
			{400, Measurement(-0.83)},
//			{500, Measurement(0.23)}
	};

	map<int,Measurement> bias_novo_dijet_sr1_r4 =
	{
			{300, Measurement(0.39)},
			{350, Measurement(1.04)},
			{400, Measurement(-0.85)},
//			{500, Measurement(0.21)}
	};

	map<int,Measurement> bias_novo_bernstein_sr1_r0 =
	{
			{300, Measurement(0.14)},
			{350, Measurement(0.15)},
			{400, Measurement(-0.01)},
//			{500, Measurement(0.68)}
	};

	map<int,Measurement> bias_novo_bernstein_sr1_r1 =
	{
			{300, Measurement(0.11)},
			{350, Measurement(0.14)},
			{400, Measurement(0)},
//			{500, Measurement(0.7)}
	};

	map<int,Measurement> bias_novo_bernstein_sr1_r2 =
	{
			{300, Measurement(0.12)},
			{350, Measurement(0.16)},
			{400, Measurement(0.02)},
//			{500, Measurement(0.64)}
	};

	map<int,Measurement> bias_novo_bernstein_sr1_r4 =
	{
			{300, Measurement(0.1)},
			{350, Measurement(0.13)},
			{400, Measurement(0)},
//			{500, Measurement(0.67)}
	};

	//SR2
	map<int,Measurement> bias_novo_dijet_sr2_r0 =
	{
			{500, Measurement(0.14)},
			{600, Measurement(0.18)},
			{700, Measurement(-0.12)},
			{900, Measurement(-0.17)}
	};

	map<int,Measurement> bias_novo_dijet_sr2_r1 =
	{
			{500, Measurement(0.14)},
			{600, Measurement(0.18)},
			{700, Measurement(-0.15)},
			{900, Measurement(-0.21)}
	};

	map<int,Measurement> bias_novo_dijet_sr2_r2 =
	{
			{500, Measurement(0.15)},
			{600, Measurement(0.16)},
			{700, Measurement(-0.09)},
			{900, Measurement(-0.2)}
	};

	map<int,Measurement> bias_novo_dijet_sr2_r4 =
	{
			{500, Measurement(0.13)},
			{600, Measurement(0.16)},
			{700, Measurement(-0.11)},
			{900, Measurement(-0.18)}
	};

	//sr3
	map<int,Measurement> bias_novo_dijet_sr3_r0 =
	{
			{1100, Measurement(-0.25)},
			{1300, Measurement(-0.17)}
	};

	map<int,Measurement> bias_novo_dijet_sr3_r1 =
	{
			{1100, Measurement(-0.28)},
			{1300, Measurement(-0.195)}
	};

	map<int,Measurement> bias_novo_dijet_sr3_r2 =
	{
			{1100, Measurement(-0.27)},
			{1300, Measurement(-0.194)}
	};

	map<int,Measurement> bias_novo_dijet_sr3_r4 =
	{
			{1100, Measurement(-0.26)},
			{1300, Measurement(-0.175)}
	};

	// Transform maps -> TGraphs
	TGraphErrors gr_bias_novo_dijet_sr1_r0(MakeTGraphErrors(bias_novo_dijet_sr1_r0));
	TGraphErrors gr_bias_novo_dijet_sr1_r1(MakeTGraphErrors(bias_novo_dijet_sr1_r1));
	TGraphErrors gr_bias_novo_dijet_sr1_r2(MakeTGraphErrors(bias_novo_dijet_sr1_r2));
	TGraphErrors gr_bias_novo_dijet_sr1_r4(MakeTGraphErrors(bias_novo_dijet_sr1_r4));

	TGraphErrors gr_bias_novo_bernstein_sr1_r0(MakeTGraphErrors(bias_novo_bernstein_sr1_r0));
	TGraphErrors gr_bias_novo_bernstein_sr1_r1(MakeTGraphErrors(bias_novo_bernstein_sr1_r1));
	TGraphErrors gr_bias_novo_bernstein_sr1_r2(MakeTGraphErrors(bias_novo_bernstein_sr1_r2));
	TGraphErrors gr_bias_novo_bernstein_sr1_r4(MakeTGraphErrors(bias_novo_bernstein_sr1_r4));

	TGraphErrors gr_bias_novo_dijet_sr2_r0(MakeTGraphErrors(bias_novo_dijet_sr2_r0));
	TGraphErrors gr_bias_novo_dijet_sr2_r1(MakeTGraphErrors(bias_novo_dijet_sr2_r1));
	TGraphErrors gr_bias_novo_dijet_sr2_r2(MakeTGraphErrors(bias_novo_dijet_sr2_r2));
	TGraphErrors gr_bias_novo_dijet_sr2_r4(MakeTGraphErrors(bias_novo_dijet_sr2_r4));

	TGraphErrors gr_bias_novo_dijet_sr3_r0(MakeTGraphErrors(bias_novo_dijet_sr3_r0));
	TGraphErrors gr_bias_novo_dijet_sr3_r1(MakeTGraphErrors(bias_novo_dijet_sr3_r1));
	TGraphErrors gr_bias_novo_dijet_sr3_r2(MakeTGraphErrors(bias_novo_dijet_sr3_r2));
	TGraphErrors gr_bias_novo_dijet_sr3_r4(MakeTGraphErrors(bias_novo_dijet_sr3_r4));

	TLegend leg = *style.legend("right,top",10);
	leg.SetX1(0.2);
	leg.SetTextSize(gStyle->GetTextSize()*0.7);
	TH2D frame_sr1("h","h",2,250,450,2,-1.5,5.0);
	frame_sr1.SetTitle("Bias=f(M_{H/A}); M_{H/A} [GeV]; a.u.");
	frame_sr1.SetStats(0);

	TCanvas can_sr1("can_sr1","can_sr1",800,600);
	frame_sr1.Draw();

	gr_bias_novo_dijet_sr1_r0.SetMarkerColor(1);
	gr_bias_novo_dijet_sr1_r0.SetLineColor(1);
	gr_bias_novo_dijet_sr1_r0.SetMarkerStyle(20);
//	leg.AddEntry(&gr_bias_novo_dijet_sr1_r0,"#mu=0, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r0,"#mu=0, Fit: Novo, Gen: DiJet-2","pl");
	gr_bias_novo_dijet_sr1_r0.Draw("PLsame");

	gr_bias_novo_dijet_sr1_r1.SetMarkerColor(2);
	gr_bias_novo_dijet_sr1_r1.SetLineColor(2);
	gr_bias_novo_dijet_sr1_r1.SetMarkerStyle(21);
//	leg.AddEntry(&gr_bias_novo_dijet_sr1_r1,"#mu=1, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r1,"#mu=1, Fit: Novo, Gen: DiJet-2","pl");
	gr_bias_novo_dijet_sr1_r1.Draw("PLsame");

	gr_bias_novo_dijet_sr1_r2.SetMarkerColor(3);
	gr_bias_novo_dijet_sr1_r2.SetLineColor(3);
	gr_bias_novo_dijet_sr1_r2.SetMarkerStyle(22);
//	leg.AddEntry(&gr_bias_novo_dijet_sr1_r2,"#mu=2, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r2,"#mu=2, Fit: Novo, Gen: DiJet-2","pl");
	gr_bias_novo_dijet_sr1_r2.Draw("PLsame");

	gr_bias_novo_dijet_sr1_r4.SetMarkerColor(4);
	gr_bias_novo_dijet_sr1_r4.SetLineColor(4);
	gr_bias_novo_dijet_sr1_r4.SetMarkerStyle(24);
//	leg.AddEntry(&gr_bias_novo_dijet_sr1_r4,"#mu=4, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r4,"#mu=4, Fit: Novo, Gen: DiJet-2","pl");
	gr_bias_novo_dijet_sr1_r4.Draw("PLsame");

	gr_bias_novo_bernstein_sr1_r0.SetMarkerColor(6);
	gr_bias_novo_bernstein_sr1_r0.SetLineColor(6);
	gr_bias_novo_bernstein_sr1_r0.SetMarkerStyle(25);
//	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r0,"#mu=0, Fit: ExtendedNovosibirsk*Turn-on, Gen: Bernstein-8*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r0,"#mu=0, Fit: Novo, Gen: Bernstein-8","pl");
	gr_bias_novo_bernstein_sr1_r0.Draw("PLsame");

	gr_bias_novo_bernstein_sr1_r1.SetMarkerColor(9);
	gr_bias_novo_bernstein_sr1_r1.SetLineColor(9);
	gr_bias_novo_bernstein_sr1_r1.SetMarkerStyle(26);
//	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r1,"#mu=1, Fit: ExtendedNovosibirsk*Turn-on, Gen: Bernstein-8*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r1,"#mu=1, Fit: Novo, Gen: Bernstein-8","pl");
	gr_bias_novo_bernstein_sr1_r1.Draw("PLsame");

	gr_bias_novo_bernstein_sr1_r2.SetMarkerColor(30);
	gr_bias_novo_bernstein_sr1_r2.SetLineColor(30);
	gr_bias_novo_bernstein_sr1_r2.SetMarkerStyle(28);
//	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r2,"#mu=2, Fit: ExtendedNovosibirsk*Turn-on, Gen: Bernstein-8*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r2,"#mu=2, Fit: Novo, Gen: Bernstein-8","pl");
	gr_bias_novo_bernstein_sr1_r2.Draw("PLsame");

	gr_bias_novo_bernstein_sr1_r4.SetMarkerColor(35);
	gr_bias_novo_bernstein_sr1_r4.SetLineColor(35);
	gr_bias_novo_bernstein_sr1_r4.SetMarkerStyle(26);
//	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r4,"#mu=4, Fit: ExtendedNovosibirsk*Turn-on, Gen: Bernstein-8*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_bernstein_sr1_r4,"#mu=4, Fit: Novo, Gen: Bernstein-8","pl");
	gr_bias_novo_bernstein_sr1_r4.Draw("PLsame");

	leg.Draw();
	can_sr1.Print("../pictures/Bias_test/bias_vs_m_sr1.pdf");
	leg.Clear();

	TH2D frame_sr2("h","h",2,400,1000,2,0,1);
	frame_sr2.SetTitle("Bias=f(M_{H/A}); M_{H/A} [GeV]; a.u.");
	frame_sr2.GetYaxis()->SetRangeUser(-1.1,1.1);
	frame_sr2.SetStats(0);

	TCanvas can_sr2("can_sr2","can_sr2",800,600);
	frame_sr2.Draw();

	gr_bias_novo_dijet_sr2_r0.SetMarkerColor(1);
	gr_bias_novo_dijet_sr2_r0.SetLineColor(1);
	gr_bias_novo_dijet_sr2_r0.SetMarkerStyle(20);
//	leg.AddEntry(&gr_bias_novo_dijet_sr2_r0,"#mu=0, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r0,"#mu=0, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr2_r0.Draw("PLsame");

	gr_bias_novo_dijet_sr2_r1.SetMarkerColor(2);
	gr_bias_novo_dijet_sr2_r1.SetLineColor(2);
	gr_bias_novo_dijet_sr2_r1.SetMarkerStyle(21);
//	leg.AddEntry(&gr_bias_novo_dijet_sr2_r1,"#mu=1, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r1,"#mu=1, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr2_r1.Draw("PLsame");

	gr_bias_novo_dijet_sr2_r2.SetMarkerColor(3);
	gr_bias_novo_dijet_sr2_r2.SetLineColor(3);
	gr_bias_novo_dijet_sr2_r2.SetMarkerStyle(22);
//	leg.AddEntry(&gr_bias_novo_dijet_sr2_r2,"#mu=2, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r2,"#mu=2, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr2_r2.Draw("PLsame");

	gr_bias_novo_dijet_sr2_r4.SetMarkerColor(4);
	gr_bias_novo_dijet_sr2_r4.SetLineColor(4);
	gr_bias_novo_dijet_sr2_r4.SetMarkerStyle(24);
//	leg.AddEntry(&gr_bias_novo_dijet_sr2_r4,"#mu=4, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r4,"#mu=4, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr2_r4.Draw("PLsame");

	leg.Draw();
	can_sr2.Print("../pictures/Bias_test/bias_vs_m_sr2.pdf");
	leg.Clear();

	TH2D frame_sr3("h","h",2,1000,1400,2,0,1);
	frame_sr3.SetTitle("Bias=f(M_{H/A}); M_{H/A} [GeV]; a.u.");
	frame_sr3.GetYaxis()->SetRangeUser(-1.1,1.1);
	frame_sr3.SetStats(0);

	TCanvas can_sr3("can_sr3","can_sr3",800,600);
	frame_sr3.Draw();

	gr_bias_novo_dijet_sr3_r0.SetMarkerColor(1);
	gr_bias_novo_dijet_sr3_r0.SetLineColor(1);
	gr_bias_novo_dijet_sr3_r0.SetMarkerStyle(20);
//	leg.AddEntry(&gr_bias_novo_dijet_sr3_r0,"#mu=0, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r0,"#mu=0, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr3_r0.Draw("PLsame");

	gr_bias_novo_dijet_sr3_r1.SetMarkerColor(2);
	gr_bias_novo_dijet_sr3_r1.SetLineColor(2);
	gr_bias_novo_dijet_sr3_r1.SetMarkerStyle(21);
//	leg.AddEntry(&gr_bias_novo_dijet_sr3_r1,"#mu=1, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r1,"#mu=1, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr3_r1.Draw("PLsame");

	gr_bias_novo_dijet_sr3_r2.SetMarkerColor(3);
	gr_bias_novo_dijet_sr3_r2.SetLineColor(3);
	gr_bias_novo_dijet_sr3_r2.SetMarkerStyle(22);
//	leg.AddEntry(&gr_bias_novo_dijet_sr3_r2,"#mu=2, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r2,"#mu=2, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr3_r2.Draw("PLsame");

	gr_bias_novo_dijet_sr3_r4.SetMarkerColor(4);
	gr_bias_novo_dijet_sr3_r4.SetLineColor(4);
	gr_bias_novo_dijet_sr3_r4.SetMarkerStyle(24);
//	leg.AddEntry(&gr_bias_novo_dijet_sr3_r4,"#mu=4, Fit: ExtendedNovosibirsk*Turn-on, Gen: ExtendedDiJet-2*Turn-on","pl");
	leg.AddEntry(&gr_bias_novo_dijet_sr1_r4,"#mu=4, Fit: Novo, Gen: DiJet","pl");
	gr_bias_novo_dijet_sr3_r4.Draw("PLsame");

	leg.Draw();
	can_sr3.Print("../pictures/Bias_test/bias_vs_m_sr3.pdf");
	leg.Clear();




	return 0;
}

TGraphErrors MakeTGraphErrors(const map<int,Measurement>& my_map){
	/*
	 * Fill TGraphError from the map
	 */
	TGraphErrors gr(my_map.size());
	int index =0;
	for(const auto& m : my_map){
		Measurement p = m.second;
		gr.SetPoint(index,m.first,p.getCentral());
		gr.SetPointError(index,0,p.getError());
		++index;
	}

	return gr;
}

void FindMaxBias(const map<int,Measurement>& my_map){
	/*
	 * Function to find a max bias within available
	 */
	vector<int> points = {300,350,400,450,500,600,700,900,1100,1300};
	double max = -9999.;
	cout<<"Maximum Bias: "<<endl;
	for(const auto& p : points){


	}

}
