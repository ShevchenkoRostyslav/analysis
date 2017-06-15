#include "Analysis/MssmHbb/macros/Drawer/HbbStyle.cc"
#include "Analysis/MssmHbb/macros/Drawer/RatioPlots.cpp"
#include "TROOT.h"

using namespace std;

void renormalise(TH1*);
void cosmetics(TH1* h);
void blind(TH1* h);

const bool findStrings(const std::string & input, const std::string & needful)
{
	std::string input1 = input;
	std::string input2 = needful;
	std::transform(input1.begin(),input1.end(),input1.begin(),::tolower);
	std::transform(input2.begin(),input2.end(),input2.begin(),::tolower);
	if(input1.find(input2) != std::string::npos) return true;
	else return false;
}

int data_vs_mc(){
	gROOT->Reset();

	std::string selection = "DataMC_ReReco";
	std::string lowM = "lowM";
	std::string MCsample = "bbbvsbbx";//"bEnriched_Tot";//"bbbvsbbx";//"QCD_Pt";//
	std::string lumi = "35673fb";
//	TFile *data_low_m = new TFile( ("../output/" + selection + "_" + lumi + "_" + lowM + "_Run2016.root").c_str() );
//	TFile *mc_low_m 	 = new TFile( ("../output/" + selection + "_" + lumi + "_" + lowM + "_" + MCsample + ".root").c_str() );

	//Input files
	TFile *data_low_m = new TFile( "../output/DataMC_3b_ReReco_35673fb_lowM_QCD.root","READ" );
	TFile *mc_low_m 	 = new TFile( "../output/bbx_ReReco_35673fb_lowM_QCD.root","READ" );
	
	TFile *fCorrection;

	//Ratio plots
	RatioPlots style;
	style.set(PRELIMINARY);

	//Declare list of Legends
	map<string,TLegend*> leg;

	//variables vector:
	vector<string> Var={//"diJet_b_m","jet_b_Ht","NumberOfJets_b","jet_b_pt1","jet_b_pt2","jet_b_eta1","jet_b_eta2","jet_b_phi1","jet_b_phi2","jet12_b_assym","jet_b_btag_csv1",
			//			"jet_b_btag_csv2","jet_b_btag_cmva1","jet_b_btag_cmva2","jet_b_deta12","jet_b_dR12", "jet_b_dphi12",
			//			"diJet_b_pt","diJet_b_eta","diJet_b_phi"};
						"diJet_b_m"};

	//histo list:
	map<string,TH1D*> h_data_low_m;
	map<string,TH1D*> h_mc_low_m;
	map<string,TH1D*> h_mc_stat;
	map<string,TH1D*> h_ratio_low_m;

	//TCanvas list
	map<string,TCanvas*> can_low_m;

	//Pt1
	for(const string &v : Var){
		string name = (v+"/"+v);//v;//
		cout<<"Current histo: "<<name<<endl;

		//Create TCavnas
		can_low_m[v] 	= new TCanvas((v+"_lowM").c_str(),(v+"_lowM").c_str(),1000,800);

		//Get histograms
		h_data_low_m[v]	= (TH1D*) data_low_m->Get(v.c_str());
		h_mc_low_m[v]	= (TH1D*) mc_low_m->Get(name.c_str());

		//Normalise MC to Data
		std::cout<<"Data/MC, Scale = "<<h_data_low_m[v]->Integral()/h_mc_low_m[v]->Integral()<<std::endl;
		h_mc_low_m[v] ->Scale(h_data_low_m[v]->Integral()/h_mc_low_m[v]->Integral());

		//Adjust style
		h_mc_low_m[v]->SetFillColor(kRed-10);
		h_data_low_m[v] -> SetMinimum(0.);
		h_data_low_m[v] ->SetTitle("");
		h_data_low_m[v]->GetYaxis()->SetTitle("a.u.");

		if(findStrings(v,"pt")) h_data_low_m[v]->GetXaxis()->SetRangeUser(0.,h_data_low_m[v]->GetXaxis()->GetXmax());
//		if( findStrings(v,"diJet") && findStrings(v,"m")) blind(h_data_low_m[v]);

		//Renormalise
		renormalise(h_data_low_m[v]);
		renormalise(h_mc_low_m[v]);

		//Stat histo for MC statistics
		h_mc_stat[v]	= (TH1D*) h_mc_low_m[v]->Clone();


		//TLegend position
		string leg_position = "top,right";
		if(findStrings(v,"NumberOfJets_b")) leg_position = "top,left";
		if(findStrings(v,"eta2") || findStrings(v,"deta") ||findStrings(v,"dR") || findStrings(v,"dphi")) leg_position = "top,left";
		if(findStrings(v,"btag")) leg_position = "top,left";

		leg[v] = (TLegend*) style.legend(leg_position,3,0.3);
		leg[v]->AddEntry(h_data_low_m[v],"Data","p");
		leg[v]->AddEntry(h_mc_low_m[v],"QCD","f");
		leg[v]->AddEntry(h_mc_stat[v],"MC Stat.","lf");

		h_ratio_low_m[v] = (TH1D*) style.DrawRatio(h_data_low_m[v],h_mc_low_m[v],can_low_m[v],leg[v],nullptr,h_mc_stat[v]);
		h_ratio_low_m[v]->GetYaxis()->SetTitle("Data/MC");
		h_ratio_low_m[v] -> SetAxisRange(0.5,1.5,"Y");
		gPad->RedrawAxis();
		style.GetTopPad();		style.drawStandardTitle();

		//Draw specific
		if(findStrings(v,"phi")) h_data_low_m[v]->SetMaximum(h_data_low_m[v]->GetMaximum()*1.2);
		if(findStrings(v,"assym")) h_ratio_low_m[v]->GetXaxis()->SetTitle("(p^{1}_{T} - p^{2}_{T})/(p^{1}_{T} + p^{2}_{T})");
		if(findStrings(v,"pt") || findStrings(v,"dphi") || findStrings(v,"Ht")) {
			h_data_low_m[v]->SetMinimum(0.01);
			gPad->SetLogy();
			h_data_low_m[v]->SetMaximum(h_data_low_m[v]->GetMaximum()*2.);
			gPad->RedrawAxis();
		}

		cout<<"Marker: "<<h_ratio_low_m[v]->GetMarkerSize()<<endl;
		can_low_m[v]->Print(("pictures/DataVsMC/" + v + "_lumi_" + lowM + "_" + MCsample + ".pdf").c_str());
		/**/
//		if(findStrings(v,"diJet") && findStrings(v,"m")){
////		if(findStrings(v,"Ht")){
//			/*
//			 * Write correction function
//			 */
//			fCorrection = new TFile("bbbvsbbnb_M12_weight.root","RECREATE");
//			h_ratio_low_m[v]->Write();
//			fCorrection->Close();
//		}
	}
//	Ht_correction->Close();
	return 0;
}

void renormalise(TH1 * h){
	for(int i=1; i<= h->GetNbinsX();++i){
		h->SetBinContent(i,h->GetBinContent(i)/h->GetBinWidth(i));
		h->SetBinError(i,h->GetBinError(i)/h->GetBinWidth(i));
	}
//	h->Scale(1./h->Integral());
}

void cosmetics(TH1* h){
	if( string(h->GetName()).find("pt") != std::string::npos){
		h->SetMinimum(0.01);
		gPad->SetLogy();
	}
}

void blind(TH1* h){
	for(int i=1;i<h->GetNbinsX();++i){
		if(h->GetXaxis()->GetBinUpEdge(i) >= 650 &&
		   h->GetXaxis()->GetBinLowEdge(i) <= 850){
			h->SetBinContent(i,0);
			h->SetBinError(i,0);
		}
	}
}
