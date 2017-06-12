/*
 * Histograms.cpp
 *
 *  Created on: 16 апр. 2016 г.
 *      Author: rostyslav
 */

#include "Analysis/MssmHbb/interface/Histograms.h"

Histograms::Histograms() : size_(0) {
	// TODO Auto-generated constructor stub

}

Histograms::~Histograms() {
	std::cout<<"I'm at ~Histograms"<<std::endl;
	// TODO Auto-generated destructor stub
	//Clean memory:
	histo_.clear();
	histo2D_.clear();
}

void Histograms::Make(const int &size, const bool & lowM) {
	size_ = size;
	TH1::SetDefaultSumw2();

	//Declare All Basic Histograms
	histo_["TotalNumberOfGenEvents"] 				= new TH1D("TotalNumberOfGenEvents","Total number of generated events",1,0,5.e+08);
	histo_["NumberOfGenEvents_afterMHat"] 			= new TH1D("NumberOfGenEvents_afterMHat","Total number of generated events after mHat cut",1,0,5.e+08);
	histo_["NumberOfGenEvents_afterMHat_rewPU"] 	= new TH1D("NumberOfGenEvents_afterMHat_rewPU","Total number of generated events after mHat cut and PU reweighting",1,0,5.e+08);

	histo_["lumi_weight"]							= new TH1D("lumi_weight","Lumi weight = (Data Luminosity / MC luminositry) ",1,0,5.e+08);
	histo_["xsection"]								= new TH1D("xsection","Cross section in pb",1,0,5.e+08);
	histo_["nTruePileup"]							= new TH1D("nTruePileup","True Pileup",75,0,75);

	histo_["NumberOfJets"]					= new TH1D("NumberOfJets","Total Number Of jets; Jet multiplicity",14,0,14);

	histo_["jet_pt1"]			=  new TH1D("jet_pt1","p_{T} of the first Leading Jet; Leading Jet p_{T}, [GeV]",size,0.,1500.);
	histo_["jet_pt2"]			=  new TH1D("jet_pt2","p_{T} of the second Leading Jet; sub-Leading Jet p_{T}, [GeV]",size,0.,1500.);
	histo_["jet_pt3"]			=  new TH1D("jet_pt3","p_{T} of the third Leading Jet; third-Leading Jet p_{T}, [GeV]",size,0.,1500.);

	histo_["jet_Ht"]			= new TH1D("jet_Ht","H_{T} of the jets; H_{T}, [GeV]",size,0.,2500.);

	histo_["jerResolution1"]	=	new TH1D("jerResolution1","JER for first Jet, [GeV]",size,0,1);
	histo_["jerResolution2"]	=	new TH1D("jerResolution2","JER for second Jet, [GeV]",size,0,1);

	histo_["jet_turnOn_pt1"]			=  new TH1D("jet_turnOn_pt1","p_{T} of the first Leading Jet; Leading Jet p_{T}, [GeV]",(int)size/2,0.,500.);
	histo_["jet_turnOn_pt2"]			=  new TH1D("jet_turnOn_pt2","p_{T} of the second Leading Jet; sub-Leading Jet p_{T}, [GeV]",(int)size/2,0.,500.);
	histo_["jet_turnOn_pt3"]			=  new TH1D("jet_turnOn_pt3","p_{T} of the third Leading Jet; third-Leading Jet p_{T}, [GeV]",(int)size/2,0.,500.);

	histo_["jet12_assym"]			= new TH1D("jet12_assym","p_{T} Leading - sub-Leading jet assymetry",size,0.,1.);

	histo_["jet_eta1"]			=  new TH1D("jet_eta1","#eta of the first Leading Jet; Leading Jet #eta ",(int)size/2,-2.5,2.5);
	histo_["jet_eta2"]			=  new TH1D("jet_eta2","#eta of the second Leading Jet; sub-Leading Jet #eta ",(int)size/2,-2.5,2.5);
	histo_["jet_eta3"]			=  new TH1D("jet_eta3","#eta of the third Leading Jet; third-Leading Jet #eta ",(int)size/2,-2.5,2.5);

	histo_["jet_phi1"]			=  new TH1D("jet_phi1","#phi of the first Leading Jet; Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_phi2"]			=  new TH1D("jet_phi2","#phi of the second Leading Jet; sub-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_phi3"]			=  new TH1D("jet_phi3","#phi of the third Leading Jet; third-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());

    histo_["jet_flavor1"]        	 =  new TH1D("jet_flavor1","True Flavor, first Leading Jet; Leading Jet btag_{csv} discr.",5,0.,6.);
    histo_["jet_flavor2"]         =  new TH1D("jet_flavor2","True Flavor, second Leading Jet; Leading Jet btag_{csv} discr.",5,0.,6.);
    histo_["jet_flavor3"]         =  new TH1D("jet_flavor3","True Flavor, third Leading Jet; Leading Jet btag_{csv} discr.",5,0.,6.);

	histo_["jet_btag_csv1"]		=  new TH1D("jet_btag_csv1","btag_{csv} discr. first Leading Jet; Leading Jet btag_{csv} discr.",size,0.5,1.);
	histo_["jet_btag_csv2"]		=  new TH1D("jet_btag_csv2","btag_{csv} discr. second Leading Jet; Leading Jet btag_{csv} discr.",size,0.5,1.);
	histo_["jet_btag_csv3"]		=  new TH1D("jet_btag_csv3","btag_{csv} discr. third Leading Jet; Leading Jet btag_{csv} discr.",size,0.,1.);

	histo_["jet_btag_cmva1"]	=  new TH1D("jet_btag_cmva1","btag_{cmva} discr. first Leading Jet; Leading Jet btag_{cmva} discr.",size,0.3,1.);
	histo_["jet_btag_cmva2"]	=  new TH1D("jet_btag_cmva2","btag_{cmva} discr. second Leading Jet; Leading Jet btag_{cmva} discr.",size,0.3,1.);
	histo_["jet_btag_cmva3"]	=  new TH1D("jet_btag_cmva3","btag_{cmva} discr. third Leading Jet; Leading Jet btag_{cmva} discr.",size,0.3,1.);

	histo_["jet_deta12"]		=  new TH1D("jet_deta12","#Delta #eta between Leading and sub-Leading jets;#Delta #eta_{12}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_dphi12"]		=  new TH1D("jet_dphi12","#Delta #phi between Leading and sub-Leading jets;#Delta #phi_{12}",(int)size/2,0,TMath::Pi());
	histo_["jet_deta13"]		=  new TH1D("jet_deta13","#Delta #eta between Leading and third-Leading jets;#Delta #eta_{13}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_dphi13"]		=  new TH1D("jet_dphi13","#Delta #phi between Leading and third-Leading jets;#Delta #phi_{13}",(int)size/2,0,TMath::Pi());
	histo_["jet_deta23"]		=  new TH1D("jet_deta23","#Delta #eta between sub-Leading and third-Leading jets;#Delta #eta_{23}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_dphi23"]		=  new TH1D("jet_dphi23","#Delta #phi between sub-Leading and third-Leading jets;#Delta #phi_{23}",(int)size/2,0,TMath::Pi());

	histo_["jet_dR12"]			=  new TH1D("jet_dR12","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);
	histo_["jet_dR13"]			=  new TH1D("jet_dR13","#Delta R between Leading and third-Leading jets;#Delta R_{13}",size,0,5);
	histo_["jet_dR23"]			=  new TH1D("jet_dR23","#Delta R between sub-Leading and third-Leading jets;#Delta R_{23}",size,0,5);

	histo_["jet_dR14"]			=  new TH1D("jet_dR14","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);
	histo_["jet_dR15"]			=  new TH1D("jet_dR15","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);
	histo_["jet_dR16"]			=  new TH1D("jet_dR16","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);
	histo_["jet_dR17"]			=  new TH1D("jet_dR17","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);
	histo_["jet_dR18"]			=  new TH1D("jet_dR18","#Delta R between Leading and sub-Leading jets;#Delta R_{12}",size,0,5);

	histo_["diJet_pt"]			=  new TH1D("diJet_pt","p_{T} of the di-Jet object;di-Jet p_{T}, [GeV]",size,0.,1500.);
	histo_["diJet_eta"]			=  new TH1D("diJet_eta","#eta of the di-Jet object; di-Jet #eta",(int)size/2,-2.5,2.5);
	histo_["diJet_phi"]			=  new TH1D("diJet_phi","#phi of the di-Jet object; di-Jet #phi",(int)size/2,0,2.*TMath::Pi());
	histo_["diJet_m"]			=  new TH1D("diJet_m","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_5GeV"]			=  new TH1D("diJet_m_5GeV","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",340,0.,1700);
	histo_["diJet_m_20GeV"]			=  new TH1D("diJet_m_20GeV","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",85,0.,1700);
	histo_["diJet_m_bb"]		=  new TH1D("diJet_m_bb","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_cc"]		=  new TH1D("diJet_m_cc","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_ll"]		=  new TH1D("diJet_m_ll","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_bc"]		=  new TH1D("diJet_m_bc","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_cl"]		=  new TH1D("diJet_m_cl","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);
	histo_["diJet_m_bl"]		=  new TH1D("diJet_m_bl","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",170,0.,1700);

	/***********************Trigger Object Parameters*********************/
	histo_["onlobject_L1_pt1"]	= new TH1D("onlobject_L1_pt1","p_{T} of the L1 Obj. matched to the first Leading Jet; L1-Obj. p^{(1)}_{T}, [GeV]",size,0.,1500.);
	histo_["onlobject_L1_pt2"]	= new TH1D("onlobject_L1_pt2","p_{T} of the L1 Obj. matched to the second Leading Jet; L1-Obj. p^{(2)}_{T}, [GeV]",size,0.,1500.);

	histo_["onlobject_L2_pt1"]	= new TH1D("onlobject_L2_pt1","p_{T} of the L2 Obj. matched to the first Leading Jet; L2-Obj. p^{(1)}_{T}, [GeV]",size,0.,1500.);
	histo_["onlobject_L2_pt2"]	= new TH1D("onlobject_L2_pt2","p_{T} of the L2 Obj. matched to the second Leading Jet; L2-Obj. p^{(2)}_{T}, [GeV]",size,0.,1500.);

	histo_["onlobject_L3_pt1"]	= new TH1D("onlobject_L3_pt1","p_{T} of the L3 Obj. matched to the first Leading Jet; L3-Obj. p^{(1)}_{T}, [GeV]",size,0.,1500.);
	histo_["onlobject_L3_pt2"]	= new TH1D("onlobject_L3_pt2","p_{T} of the L3 Obj. matched to the second Leading Jet; L3-Obj. p^{(2)}_{T}, [GeV]",size,0.,1500.);

	histo_["onlobject_turnOn_L1_pt1"]	= new TH1D("onlobject_turnOn_L1_pt1","p_{T} of the L1 Obj. matched to the first Leading Jet; L1-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,0.,500.);
	histo_["onlobject_turnOn_L1_pt2"]	= new TH1D("onlobject_turnOn_L1_pt2","p_{T} of the L1 Obj. matched to the second Leading Jet; L1-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,0.,500.);

	histo_["onlobject_turnOn_L2_pt1"]	= new TH1D("onlobject_turnOn_L2_pt1","p_{T} of the L2 Obj. matched to the first Leading Jet; L2-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,0.,500.);
	histo_["onlobject_turnOn_L2_pt2"]	= new TH1D("onlobject_turnOn_L2_pt2","p_{T} of the L2 Obj. matched to the second Leading Jet; L2-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,0.,500.);

	histo_["onlobject_turnOn_L3_pt1"]	= new TH1D("onlobject_turnOn_L3_pt1","p_{T} of the L3 Obj. matched to the first Leading Jet; L3-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,0.,500.);
	histo_["onlobject_turnOn_L3_pt2"]	= new TH1D("onlobject_turnOn_L3_pt2","p_{T} of the L3 Obj. matched to the second Leading Jet; L3-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,0.,500.);

	histo_["onlobject_turnOn_L1ovOff_pt1"] = new TH1D("onlobject_turnOn_L1ovOff_pt1","p_{T} of the L1 Obj. matched to the first Leading Jet Over Offline p^{(1)}_{T}; L1-Obj. p^{(1)}_{T} / (Offline p^{(1)}_{T}), [GeV]",(int)size/2,0.,3.);
	histo_["onlobject_turnOn_L2ovOff_pt1"] = new TH1D("onlobject_turnOn_L2ovOff_pt1","p_{T} of the L2 Obj. matched to the first Leading Jet Over Offline p^{(1)}_{T}; L1-Obj. p^{(1)}_{T} / (Offline p^{(1)}_{T}), [GeV]",(int)size/2,0.,3.);
	histo_["onlobject_turnOn_L3ovOff_pt1"] = new TH1D("onlobject_turnOn_L3ovOff_pt1","p_{T} of the L3 Obj. matched to the first Leading Jet Over Offline p^{(1)}_{T}; L1-Obj. p^{(1)}_{T} / (Offline p^{(1)}_{T}), [GeV]",(int)size/2,0.,3.);

	histo_["onlobject_L1_eta1"]	= new TH1D("onlobject_L1_eta1","#eta of the L1 Obj. matched to the first Leading Jet; L1-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,-2.5,2.5);
	histo_["onlobject_L1_eta2"]	= new TH1D("onlobject_L1_eta2","#eta of the L1 Obj. matched to the second Leading Jet; L1-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,-2.5,2.5);

	histo_["onlobject_L2_eta1"]	= new TH1D("onlobject_L2_eta1","#eta of the L2 Obj. matched to the first Leading Jet; L2-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,-2.5,2.5);
	histo_["onlobject_L2_eta2"]	= new TH1D("onlobject_L2_eta2","#eta of the L2 Obj. matched to the second Leading Jet; L2-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,-2.5,2.5);

	histo_["onlobject_L3_eta1"]	= new TH1D("onlobject_L3_eta1","#eta of the L3 Obj. matched to the first Leading Jet; L3-Obj. p^{(1)}_{T}, [GeV]",(int)size/2,-2.5,2.5);
	histo_["onlobject_L3_eta2"]	= new TH1D("onlobject_L3_eta2","#eta of the L3 Obj. matched to the second Leading Jet; L3-Obj. p^{(2)}_{T}, [GeV]",(int)size/2,-2.5,2.5);

	histo2D_["onlobject_L1_pt1_vs_off"] = new TH2D("onlobject_L1_pt1_vs_off","onlobject_L1_pt1_vs_off",(int)size/2,0.,1500,(int)size/2,0.,1500);
	histo2D_["onlobject_L2_pt1_vs_off"] = new TH2D("onlobject_L2_pt1_vs_off","onlobject_L2_pt1_vs_off",(int)size/2,0.,1500,(int)size/2,0.,1500);
	histo2D_["onlobject_L3_pt1_vs_off"] = new TH2D("onlobject_L3_pt1_vs_off","onlobject_L3_pt1_vs_off",(int)size/2,0.,1500,(int)size/2,0.,1500);

	histo2D_["onlobject_turnOn_L1_pt1_vs_off"] = new TH2D("onlobject_turnOn_L1_pt1_vs_off","onlobject_L1_pt1_vs_off",(int)size/4,0.,500,(int)size/4,0.,500);
	histo2D_["onlobject_turnOn_L2_pt1_vs_off"] = new TH2D("onlobject_turnOn_L2_pt1_vs_off","onlobject_L2_pt1_vs_off",(int)size/4,0.,500,(int)size/4,0.,500);
	histo2D_["onlobject_turnOn_L3_pt1_vs_off"] = new TH2D("onlobject_turnOn_L3_pt1_vs_off","onlobject_L3_pt1_vs_off",(int)size/4,0.,500,(int)size/4,0.,500);
	/*********************************************************************/
	
	/***********************Generated Jets Parameters*********************/
	histo_["Genjet_pt1"]		= new TH1D("Genjet_pt1","p_{T} of the first Gen Leading Jet; Gen Leading Jet p_{T}, [GeV]",size,0.,1500.);
	histo_["Genjet_pt2"]            = new TH1D("Genjet_pt2","p_{T} of the second Gen Leading Jet; Gen Leading Jet p_{T}, [GeV]",size,0.,1500.);
	histo_["Genjet_pt3"]            = new TH1D("Genjet_pt3","p_{T} of the third Gen Leading Jet; Gen Leading Jet p_{T}, [GeV]",size,0.,1500.);

	histo_["Genjet_turnOn_pt1"]		= new TH1D("Genjet_turnOn_pt1","p_{T} of the first Gen Leading Jet; Gen Leading Jet p_{T}, [GeV]",(int)size/2,0.,500.);

    histo_["Genjet12_assym"]        = new TH1D("Genjet12_assym","Gen p_{T} Leading - sub-Leading jet assymetry",size,0.,1.);

    histo_["Genjet_eta1"]           =  new TH1D("Genjet_eta1","#eta of the first Leading Jet; Leading Jet #eta ",(int)size/2,-2.5,2.5);
    histo_["Genjet_eta2"]           =  new TH1D("Genjet_eta2","#eta of the second Leading Jet; sub-Leading Jet #eta ",(int)size/2,-2.5,2.5);
    histo_["Genjet_eta3"]           =  new TH1D("Genjet_eta3","#eta of the third Leading Jet; third-Leading Jet #eta ",(int)size/2,-2.5,2.5);

    histo_["Genjet_phi1"]           =  new TH1D("jet_phi1","#phi of the first Leading Jet; Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
    histo_["Genjet_phi2"]           =  new TH1D("jet_phi2","#phi of the second Leading Jet; sub-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
    histo_["Genjet_phi3"]           =  new TH1D("jet_phi3","#phi of the third Leading Jet; third-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());

    histo_["Genjet_flavor1"]         =  new TH1D("Genjet_flavor1","True Flavor, first Leading Gen Jet; Leading Jet btag_{csv} discr.",5,0.,6.);
    histo_["Genjet_flavor2"]         =  new TH1D("Genjet_flavor2","True Flavor, second Leading Gen Jet; Leading Jet btag_{csv} discr.",5,0.,6.);
    histo_["Genjet_flavor3"]         =  new TH1D("Genjet_flavor3","True Flavor, third Leading Gen Jet; Leading Jet btag_{csv} discr.",5,0.,6.);

    histo2D_["Genjet_vs_off_pt1"] = new TH2D("Genjet_vs_off_pt1","Genjet_vs_off_pt1",(int)size/2,0.,1500,(int)size/2,0.,1500);
    histo2D_["Genjet_turnOn_vs_off_pt1"] = new TH2D("Genjet_vs_off_pt1","Genjet_vs_off_pt1",(int)size/4,0.,500,(int)size/4,0.,500);
	
	/*********************************************************************/

	//Signal templates
	//SFb
    /*
	int temp_bins;
	double temp_bin_lowM[] = {240.,260.,280.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,780.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.};
	double temp_bin_highM[] = {350.,375.,400.,425.,450.,475.,500.,525.,550.,575.,600.,625.,650.,675.,700.,725.,750.,775.,800.,850.,900.,950.,1000.,1050.,1100.,1150.,1200.,1250.,1300.,1350.,1400.,1450.,1500.,1550.,1600.,1650.,1700.,1750.,1800.,1850.,1900.,1950.,2000.,2050.,2100.,2150.,2200.,2250.,2300.,2350.,2400.,2450.,2500.,2550.,2600.,2650.,2700.,2750.,2800.,2850.,2900.,2950.,3000.};
	
	if (lowM) temp_bins = sizeof(temp_bin_lowM)/sizeof(double) - 1;
	else temp_bins = sizeof(temp_bin_highM)/sizeof(double) - 1;	
	double temp_bin[100];
	if(lowM) {
		std::copy(std::begin(temp_bin_lowM),std::end(temp_bin_lowM),std::begin(temp_bin));
	}
	else{
		 std::copy(std::begin(temp_bin_highM),std::end(temp_bin_highM),std::begin(temp_bin));
	}
	histo_["template_Mbb"]			=  new TH1D("template_Mbb","M_{12} of the di_jet object; di-Jet M_{12}, [GeV]",temp_bins,temp_bin);
	histo_["template_SFb_down"]			= new TH1D("template_SFb_down","",temp_bins,temp_bin);
	histo_["template_SFb_up"]			= new TH1D("template_SFb_up","",temp_bins,temp_bin);

	//SFl
	histo_["template_SFl_down"]			= new TH1D("template_SFl_down","",temp_bins,temp_bin);
	histo_["template_SFl_up"]			= new TH1D("template_SFl_up","",temp_bins,temp_bin);

	//JES
	histo_["template_JES_down"]			= new TH1D("template_JES_down","",temp_bins,temp_bin);
	histo_["template_JES_up"]			= new TH1D("template_JES_up","",temp_bins,temp_bin);

	//JER
	histo_["template_JER_down"]			= new TH1D("template_JER_down","",temp_bins,temp_bin);
	histo_["template_JER_up"]			= new TH1D("template_JER_up","",temp_bins,temp_bin);

	//Pileup
	histo_["template_PU_down"]			= new TH1D("template_PU_down","",temp_bins,temp_bin);
	histo_["template_PU_up"]				= new TH1D("template_PU_up","",temp_bins,temp_bin);

	//Kinematic trigger efificiency
	histo_["template_PtEff_down"]		= new TH1D("template_PtEff_down","",temp_bins,temp_bin);
	histo_["template_PtEff_up"]			= new TH1D("template_PtEff_up","",temp_bins,temp_bin);
	*/
	int vis_bins = 80;
	double xMax = 1700;
//	vis_bins = temp_bins*(1 + xMin/(xMax-xMin));

	histo_["template_Mbb_VIS"]				= new TH1D("template_Mbb_VIS","M_{12} of the di_jet object; di-Jet M_{12}, [GeV]",vis_bins,0.,xMax);
	histo_["template_SFb_VIS_down"]			= new TH1D("template_SFb_VIS_down","",vis_bins,0,xMax);
	histo_["template_SFb_VIS_up"]			= new TH1D("template_SFb_VIS_up","",vis_bins,0,xMax);
	histo_["template_SFl_VIS_down"]			= new TH1D("template_SFl_VIS_down","",vis_bins,0,xMax);
	histo_["template_SFl_VIS_up"]			= new TH1D("template_SFl_VIS_up","",vis_bins,0,xMax);
	histo_["template_JES_VIS_down"]			= new TH1D("template_JES_VIS_down","",vis_bins,0,xMax);
	histo_["template_JES_VIS_up"]			= new TH1D("template_JES_VIS_up","",vis_bins,0,xMax);
	histo_["template_JER_VIS_down"]			= new TH1D("template_JER_VIS_down","",vis_bins,0,xMax);
	histo_["template_JER_VIS_up"]			= new TH1D("template_JER_VIS_up","",vis_bins,0,xMax);
	histo_["template_PU_VIS_down"]			= new TH1D("template_PU_VIS_down","",vis_bins,0,xMax);
	histo_["template_PU_VIS_up"]			= new TH1D("template_PU_VIS_up","",vis_bins,0,xMax);
	histo_["template_PtEff_VIS_down"]		= new TH1D("template_PtEff_VIS_down","",vis_bins,0,xMax);
	histo_["template_PtEff_VIS_up"]			= new TH1D("template_PtEff_VIS_up","",vis_bins,0,xMax);
	histo_["template_BTagEff_VIS_down"]		= new TH1D("template_BTagEff_VIS_down","",vis_bins,0,xMax);
	histo_["template_BTagEff_VIS_up"]		= new TH1D("template_BTagEff_VIS_up","",vis_bins,0,xMax);
	histo_["template_PDF_VIS_down"]		        = new TH1D("template_PDF_VIS_down","",vis_bins,0,xMax);
	histo_["template_PDF_VIS_up"]		        = new TH1D("template_PDF_VIS_up","",vis_bins,0,xMax);

//	syst_["Pileup"]

	// Trigger Efficiency Histograms:
	histo_["KinTrigEff_pt1_distr_matched_PFJet80"]	= new TH1D("KinTrigEff_pt1_distr_matched_PFJet80","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_pt1_distr_not_matched_PFJet80"]	= new TH1D("KinTrigEff_pt1_distr_not_matched_PFJet80","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);

	histo_["KinTrigEff_eta1_distr_matched_PFJet80"]	= new TH1D("KinTrigEff_eta1_distr_matched_PFJet80","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_eta1_distr_not_matched_PFJet80"]	= new TH1D("KinTrigEff_eta1_distr_not_matched_PFJet80","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Num_1D_PF80_PF100_pt1_central_eta"] = new TH1D("KinTrigEff_Num_1D_PF80_PF100_pt1_central_eta","1D Kinematic Trigger Efficiency with PF80 f(pt1), central Eta",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF80_PF100_pt1_central_eta"] = new TH1D("KinTrigEff_Den_1D_PF80_PF100_pt1_central_eta","1D Kinematic Trigger Efficiency with PF80 f(pt1), central Eta",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF80_PF100_pt1_central_eta"] = new TH1D("KinTrigEff_1D_PF80_PF100_pt1_central_eta","1D Kinematic Trigger Efficiency with PF80 f(pt1), central Eta",(int)size/2,0.,500.);

	//Latest version - f(pt,eta) 2D histograms;
	double Bins2D_Pt_PF100[] = {50.,100.,105.,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,210,220,230,240,250,260,270,280,290,300,330,360,390,420,450,500};
	double Bins2D_Eta_PF100[]= {0.,0.6,1.7,2.2};
	//PF80 for PF100
	histo2D_["KinTrigEff_2D_PF80_PF100_pt1eta1"]		= new TH2D("KinTrigEff_2D_PF80_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_pt1eta1"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_pt1eta1"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_2D_PF80_PF100_pt2eta2"]		= new TH2D("KinTrigEff_2D_PF80_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_pt2eta2"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_pt2eta2"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);

	//PF60 for PF100
	histo2D_["KinTrigEff_2D_PF60_PF100_pt1eta1"]		= new TH2D("KinTrigEff_2D_PF60_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Num_2D_PF60_PF100_pt1eta1"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Den_2D_PF60_PF100_pt1eta1"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF100_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_2D_PF60_PF100_pt2eta2"]		= new TH2D("KinTrigEff_2D_PF60_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Num_2D_PF60_PF100_pt2eta2"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);
	histo2D_["KinTrigEff_Den_2D_PF60_PF100_pt2eta2"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF100_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF100)/sizeof(double)-1,Bins2D_Pt_PF100,sizeof(Bins2D_Eta_PF100)/sizeof(double)-1,Bins2D_Eta_PF100);

	double Bins2D_Pt_PF160[] = {160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,260,270,280,290,300,310,320,330,340,350,360,390,420,450,500};
	double Bins2D_Eta_PF160[]= {0.,0.6,1.7,2.2};
	//PF80 for PF160
	histo2D_["KinTrigEff_2D_PF80_PF160_pt1eta1"]		= new TH2D("KinTrigEff_2D_PF80_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Num_2D_PF80_PF160_pt1eta1"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Den_2D_PF80_PF160_pt1eta1"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_2D_PF80_PF160_pt2eta2"]		= new TH2D("KinTrigEff_2D_PF80_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Num_2D_PF80_PF160_pt2eta2"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Den_2D_PF80_PF160_pt2eta2"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);

	//PF60 for PF160
	histo2D_["KinTrigEff_2D_PF60_PF160_pt1eta1"]		= new TH2D("KinTrigEff_2D_PF60_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Num_2D_PF60_PF160_pt1eta1"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Den_2D_PF60_PF160_pt1eta1"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF160_pt1eta1","2D Kinematic Trigger Efficiency with PF80 f(pt1,eta1)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_2D_PF60_PF160_pt2eta2"]		= new TH2D("KinTrigEff_2D_PF60_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Num_2D_PF60_PF160_pt2eta2"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);
	histo2D_["KinTrigEff_Den_2D_PF60_PF160_pt2eta2"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF160_pt2eta2","2D Kinematic Trigger Efficiency with PF80 f(pt2,eta2)",sizeof(Bins2D_Pt_PF160)/sizeof(double)-1,Bins2D_Pt_PF160,sizeof(Bins2D_Eta_PF160)/sizeof(double)-1,Bins2D_Eta_PF160);

	//..........................
	//PF60 for PF100
	histo_["KinTrigEff_1D_PF60_PF100_pt1"]		= new TH1D("KinTrigEff_1D_PF60_PF100_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF60_PF100_pt2"]		= new TH1D("KinTrigEff_1D_PF60_PF100_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF60_PF100_eta1"]		= new TH1D("KinTrigEff_1D_PF60_PF100_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_1D_PF60_PF100_eta2"]		= new TH1D("KinTrigEff_1D_PF60_PF100_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Num_1D_PF60_PF100_pt1"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF100_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF60_PF100_pt2"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF100_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF60_PF100_eta1"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF100_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Num_1D_PF60_PF100_eta2"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF100_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Den_1D_PF60_PF100_pt1"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF100_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF60_PF100_pt2"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF100_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF60_PF100_eta1"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF100_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Den_1D_PF60_PF100_eta2"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF100_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	//dEta PF60 for PF100
	histo_["KinTrigEff_1D_PF60_PF100_dEta"]		= new TH1D("KinTrigEff_1D_PF60_PF100_dEta","1D Kinematic Trigger Efficiency with PF60 f(dEta)",(int)size/2,0.,2.5);
	histo_["KinTrigEff_Num_1D_PF60_PF100_dEta"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF100_dEta","1D Kinematic Trigger Efficiency with PF60 f(dEta)",(int)size/2,0.,2.5);
	histo_["KinTrigEff_Den_1D_PF60_PF100_dEta"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF100_dEta","1D Kinematic Trigger Efficiency with PF60 f(dEta)",(int)size/2,0.,2.5);

	//2D PF60 for PF100
	double Bins2DPF100[] = {50.,100.,112.,118.,124,130.,135.,140.,145,150,170,200,250,360,500};
	histo2D_["KinTrigEff_2D_PF60_PF100_pt1vspt2"]		= new TH2D("KinTrigEff_2D_PF60_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF60_PF100_pt1vspt2"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF60_PF100_pt1vspt2"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	histo2D_["KinTrigEff_2D_PF60_PF100_pt2vspt1"]		= new TH2D("KinTrigEff_2D_PF60_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF60_PF100_pt2vspt1"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF60_PF100_pt2vspt1"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	//PF80 for PF100
	histo_["KinTrigEff_1D_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_1D_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF80_PF100_pt2"]		= new TH1D("KinTrigEff_1D_PF80_PF100_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF80_PF100_eta1"]		= new TH1D("KinTrigEff_1D_PF80_PF100_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_1D_PF80_PF100_eta2"]		= new TH1D("KinTrigEff_1D_PF80_PF100_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Num_1D_PF80_PF100_pt1"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF80_PF100_pt2"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF80_PF100_eta1"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Num_1D_PF80_PF100_eta2"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Den_1D_PF80_PF100_pt1"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF80_PF100_pt2"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF80_PF100_eta1"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Den_1D_PF80_PF100_eta2"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	//dEta for PF100
	histo_["KinTrigEff_1D_PF80_PF100_dEta"]		= new TH1D("KinTrigEff_1D_PF80_PF100_dEta","1D Kinematic Trigger Efficiency with PF80 f(dEta)",(int)size/2,0.,2.5);
	histo_["KinTrigEff_Num_1D_PF80_PF100_dEta"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_dEta","1D Kinematic Trigger Efficiency with PF80 f(dEta)",(int)size/2,0.,2.5);
	histo_["KinTrigEff_Den_1D_PF80_PF100_dEta"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_dEta","1D Kinematic Trigger Efficiency with PF80 f(dEta)",(int)size/2,0.,2.5);

	//2D PF80 for PF100
	histo2D_["KinTrigEff_2D_PF80_PF100_pt1vspt2"]		= new TH2D("KinTrigEff_2D_PF80_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_pt1vspt2"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_pt1vspt2"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	histo2D_["KinTrigEff_2D_PF80_PF100_pt2vspt1"]		= new TH2D("KinTrigEff_2D_PF80_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_pt2vspt1"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_pt2vspt1"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	//the same binning as for 2D
	histo_["KinTrigEff_1D_PF80_PF100_finebins_pt1"]		= new TH1D("KinTrigEff_1D_PF80_PF100_finebins_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo_["KinTrigEff_Num_1D_PF80_PF100_finebins_pt1"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF100_finebins_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo_["KinTrigEff_Den_1D_PF80_PF100_finebins_pt1"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF100_finebins_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	//Eta binned
	histo2D_["KinTrigEff_2D_PF80_PF100_eta_l1_pt1vspt2"]		= new TH2D("KinTrigEff_2D_PF80_PF100_eta_l1_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_eta_l1_pt1vspt2"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_eta_l1_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_eta_l1_pt1vspt2"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_eta_l1_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	histo2D_["KinTrigEff_2D_PF80_PF100_eta_l1_pt2vspt1"]		= new TH2D("KinTrigEff_2D_PF80_PF100_eta_l1_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Num_2D_PF80_PF100_eta_l1_pt2vspt1"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF100_eta_l1_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);
	histo2D_["KinTrigEff_Den_2D_PF80_PF100_eta_l1_pt2vspt1"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF100_eta_l1_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100,sizeof(Bins2DPF100)/sizeof(double)-1,Bins2DPF100);

	//PF60 for PF160
	histo_["KinTrigEff_1D_PF60_PF160_pt1"]		= new TH1D("KinTrigEff_1D_PF60_PF160_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF60_PF160_pt2"]		= new TH1D("KinTrigEff_1D_PF60_PF160_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF60_PF160_eta1"]		= new TH1D("KinTrigEff_1D_PF60_PF160_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_1D_PF60_PF160_eta2"]		= new TH1D("KinTrigEff_1D_PF60_PF160_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Num_1D_PF60_PF160_pt1"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF160_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF60_PF160_pt2"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF160_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF60_PF160_eta1"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF160_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Num_1D_PF60_PF160_eta2"]	= new TH1D("KinTrigEff_Num_1D_PF60_PF160_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Den_1D_PF60_PF160_pt1"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF160_pt1","1D Kinematic Trigger Efficiency with PF60 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF60_PF160_pt2"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF160_pt2","1D Kinematic Trigger Efficiency with PF60 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF60_PF160_eta1"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF160_eta1","1D Kinematic Trigger Efficiency with PF60 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Den_1D_PF60_PF160_eta2"]	= new TH1D("KinTrigEff_Den_1D_PF60_PF160_eta2","1D Kinematic Trigger Efficiency with PF60 f(eta2)",(int)size/2,-2.5,2.5);

	//2D PF60 for PF160
	double Bins2DPF160[] = {110,160,165,170,180,190,200,210,220,240,260,300,350,400,500};
	histo2D_["KinTrigEff_2D_PF60_PF160_pt1vspt2"]		= new TH2D("KinTrigEff_2D_PF60_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Num_2D_PF60_PF160_pt1vspt2"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Den_2D_PF60_PF160_pt1vspt2"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF60 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);

	histo2D_["KinTrigEff_2D_PF60_PF160_pt2vspt1"]		= new TH2D("KinTrigEff_2D_PF60_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Num_2D_PF60_PF160_pt2vspt1"]	= new TH2D("KinTrigEff_Num_2D_PF60_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Den_2D_PF60_PF160_pt2vspt1"]	= new TH2D("KinTrigEff_Den_2D_PF60_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF60 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	//PF80 for PF160
	histo_["KinTrigEff_1D_PF80_PF160_pt1"]		= new TH1D("KinTrigEff_1D_PF80_PF160_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF80_PF160_pt2"]		= new TH1D("KinTrigEff_1D_PF80_PF160_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_1D_PF80_PF160_eta1"]		= new TH1D("KinTrigEff_1D_PF80_PF160_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_1D_PF80_PF160_eta2"]		= new TH1D("KinTrigEff_1D_PF80_PF160_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Num_1D_PF80_PF160_pt1"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF160_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF80_PF160_pt2"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF160_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_PF80_PF160_eta1"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF160_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Num_1D_PF80_PF160_eta2"]	= new TH1D("KinTrigEff_Num_1D_PF80_PF160_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	histo_["KinTrigEff_Den_1D_PF80_PF160_pt1"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF160_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF80_PF160_pt2"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF160_pt2","1D Kinematic Trigger Efficiency with PF80 f(pt2)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_PF80_PF160_eta1"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF160_eta1","1D Kinematic Trigger Efficiency with PF80 f(eta1)",(int)size/2,-2.5,2.5);
	histo_["KinTrigEff_Den_1D_PF80_PF160_eta2"]	= new TH1D("KinTrigEff_Den_1D_PF80_PF160_eta2","1D Kinematic Trigger Efficiency with PF80 f(eta2)",(int)size/2,-2.5,2.5);

	//2D PF80 for PF160
	histo2D_["KinTrigEff_2D_PF80_PF160_pt1vspt2"]		= new TH2D("KinTrigEff_2D_PF80_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Num_2D_PF80_PF160_pt1vspt2"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Den_2D_PF80_PF160_pt1vspt2"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF160_pt1vspt2","2D Kinematic Trigger Efficiency with PF80 f(pt1,pt2)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);

	histo2D_["KinTrigEff_2D_PF80_PF160_pt2vspt1"]		= new TH2D("KinTrigEff_2D_PF80_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Num_2D_PF80_PF160_pt2vspt1"]	= new TH2D("KinTrigEff_Num_2D_PF80_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);
	histo2D_["KinTrigEff_Den_2D_PF80_PF160_pt2vspt1"]	= new TH2D("KinTrigEff_Den_2D_PF80_PF160_pt2vspt1","2D Kinematic Trigger Efficiency with PF80 f(pt2,pt1)",sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160,sizeof(Bins2DPF160)/sizeof(double)-1,Bins2DPF160);

	// WITH TRUE jets
	histo_["KinTrigEff_1D_True_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_1D_True_PF80_PF100_pt1"," True 1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Num_1D_True_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Num_1D_True_PF80_PF100_pt1"," True 1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);
	histo_["KinTrigEff_Den_1D_True_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Den_1D_True_PF80_PF100_pt1"," True 1D Kinematic Trigger Efficiency with PF80 f(pt1)",(int)size/2,0.,500.);

	histo_["KinTrigEff_Num_1D_L1_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Num_1D_L1_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L1 matching",(int)size/2,0,500);
	histo_["KinTrigEff_Den_1D_L1_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Den_1D_L1_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L1 matching",(int)size/2,0,500);
	histo_["KinTrigEff_Num_1D_L2_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Num_1D_L2_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L2 matching",(int)size/2,0,500);
	histo_["KinTrigEff_Den_1D_L2_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Den_1D_L2_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L2 matching",(int)size/2,0,500);
	histo_["KinTrigEff_Num_1D_L3_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Num_1D_L3_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L3 matching",(int)size/2,0,500);
	histo_["KinTrigEff_Den_1D_L3_PF80_PF100_pt1"]			= new TH1D("KinTrigEff_Den_1D_L3_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with only L3 matching",(int)size/2,0,500);

	histo_["KinTrigEff_Num_1D_noMuons_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Num_1D_noMuons_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) without MUONS",(int)size/2,0,500);
	histo_["KinTrigEff_Den_1D_noMuons_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Den_1D_noMuons_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) without MUONS",(int)size/2,0,500);

	histo_["KinTrigEff_Num_1D_csvv2T_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Num_1D_csvv2T_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with csvv2T jets",(int)size/2,0,500);
	histo_["KinTrigEff_Den_1D_csvv2T_PF80_PF100_pt1"]		= new TH1D("KinTrigEff_Den_1D_csvv2T_PF80_PF100_pt1","1D Kinematic Trigger Efficiency with PF80 f(pt1) with csvv2T jets",(int)size/2,0,500);

	//For Chayanit
	histo_["data_Mbb_230"]		= new TH1D("data_Mbb_230","Mbb ",76,230,1750);
	histo_["data_Mbb_220"]		= new TH1D("data_Mbb_220","Mbb ",88,0,1760);
	histo_["data_Mbb_20"]		= new TH1D("data_Mbb_20","Mbb ",100,0,2000);
	histo_["data_Mbb_30"]		= new TH1D("data_Mbb_30","Mbb ",70,0,2100);
	histo_["data_Mbb_50"]		= new TH1D("data_Mbb_50","Mbb ",40,0,2000);

	//TRUE Bs
	for(int i = 0 ; i < 10; ++i){
		std::string name = "bJet" + std::to_string(i) + "_pt";
		histo_[name] = new TH1D(name.c_str(),name.c_str(),100,0,2000);
		name = "bJet" + std::to_string(i) + "_eta";
		histo_[name] = new TH1D(name.c_str(),name.c_str(),30,-5,5);
		name = "bJet" + std::to_string(i) + "_btag";
		histo_[name] = new TH1D(name.c_str(),name.c_str(),100,0,1);
		name = "bJet" + std::to_string(i) + "_order";
		histo_[name] = new TH1D(name.c_str(),name.c_str(),15,0,15);
	}
	histo_["nTrue_b"]			= new TH1D("nTrue_b","nTrue_b",15,0,15);

	//Cut optimisation
	std::vector<std::string> pt_cuts = {"40","50","60","70","80","90","100","120","140","160","180","200"};
	std::vector<std::string> d_csv	 = {"0.46","0.8","0.935"};
	std::vector<std::string> d_cmva	 = {"-0.715","0.185","0.875"};
	for(const auto & s : pt_cuts) {
		histo_[("jet1_m12_pt>"+s)] = new TH1D( ("jet1_m12_pt>"+s).c_str(),("jet1_pt>"+s).c_str(),100,0,2000 );
		histo_[("jet2_m12_pt>"+s)] = new TH1D( ("jet2_m12_pt>"+s).c_str(),("jet2_pt>"+s).c_str(),100,0,2000 );
		histo_[("jet3_m12_pt>"+s)] = new TH1D( ("jet3_m12_pt>"+s).c_str(),("jet3_pt>"+s).c_str(),100,0,2000 );
	}
	for(const auto & s1 : d_csv) {
		for(const auto & s2 : d_csv) {
			for(const auto & s3 : d_csv) {
		histo_[("jet1_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet1_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet1_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
		histo_[("jet2_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet2_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet2_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
		histo_[("jet3_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet3_m12_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet3_d_csv>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
			}
		}
	}
	for(const auto & s1 : d_cmva) {
		for(const auto & s2 : d_cmva) {
			for(const auto & s3 : d_cmva) {
		histo_[("jet1_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet1_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet1_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
		histo_[("jet2_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet2_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet2_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
		histo_[("jet3_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3)] = new TH1D( ("jet3_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),("jet3_m12_d_cmva>"+s1+"_j2>"+s2+"_j3>"+s3).c_str(),100,0,2000 );
			}
		}
	}

	//Data vs MC comparison
	DeclareDataMCHistograms(size);
}

void Histograms::DeclareDataMCHistograms(const int &size){
	TH1::SetDefaultSumw2();

	double PtBins[]={100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,375,400,450,500,575,650,750,1000};
	double Pt3Bins[]={40,60,80,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,375,400,450,500,575,650,750,1000};
	double PtAssymBins[] = {0.,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.24,0.28,0.32,0.36,0.4,0.45,0.5,0.6};
	double M12Bins[] = {200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,700,750,800,850,940,1100};
	double NJetsBins[] = {2,3,4,5,6,7,8,9,14};
	double HtBins[]={200,205,210,215,220,240,260,280,300,330,360,390,420,450,480,510,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2200,2600};

	//Central
	histo_["jet_b_pt1"]			=  new TH1D( ("jet_b_pt1"),"p_{T} of the first Leading Jet; Leading Jet p_{T}, [GeV]",sizeof(PtBins)/sizeof(double)-1,PtBins);
	histo_["jet_b_pt2"]			=  new TH1D( ("jet_b_pt2"),"p_{T} of the second Leading Jet; sub-Leading Jet p_{T}, [GeV]",sizeof(PtBins)/sizeof(double)-1,PtBins);
	histo_["jet_b_pt3"]			=  new TH1D( ("jet_b_pt3"),"p_{T} of the third Leading Jet; third-Leading Jet p_{T}, [GeV]",sizeof(Pt3Bins)/sizeof(double)-1,Pt3Bins);

	histo_["jet_b_Ht"]			= new TH1D("jet_b_Ht","H_{T} of the jets; H_{T}, [GeV]",sizeof(HtBins)/sizeof(double)-1,HtBins);

	histo_["jet12_b_assym"]			= new TH1D( ("jet12_b_assym"),"p_{T} Leading - sub-Leading jet assymetry; #frac{p^{(1)}_{T} - p^{(2)}_{T}}{p^{(1)}_{T} + p^{(2)}_{T}}",sizeof(PtAssymBins)/sizeof(double)-1,PtAssymBins);

	histo_["jet_b_eta1"]			=  new TH1D( ("jet_b_eta1"),"#eta of the first Leading Jet; Leading Jet #eta ",30,-2.2,2.2);
	histo_["jet_b_eta2"]			=  new TH1D( ("jet_b_eta2"),"#eta of the second Leading Jet; sub-Leading Jet #eta ",30,-2.2,2.2);
	histo_["jet_b_eta3"]			=  new TH1D( ("jet_b_eta3"),"#eta of the third Leading Jet; third-Leading Jet #eta ",30,-2.2,2.2);

	histo_["jet_b_phi1"]			=  new TH1D( ("jet_b_phi1"),"#phi of the first Leading Jet; Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_b_phi2"]			=  new TH1D( ("jet_b_phi2"),"#phi of the second Leading Jet; sub-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_b_phi3"]			=  new TH1D( ("jet_b_phi3"),"#phi of the third Leading Jet; third-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());

	histo_["jet_b_btag_csv1"]		=  new TH1D( ("jet_b_btag_csv1"),"btag_{csv} discr. first Leading Jet; Leading Jet btag_{csv} discr.",30,0.935,1.);
	histo_["jet_b_btag_csv2"]		=  new TH1D( ("jet_b_btag_csv2"),"btag_{csv} discr. second Leading Jet; sub-Leading Jet btag_{csv} discr.",30,0.935,1.);
	histo_["jet_b_btag_csv3"]		=  new TH1D( ("jet_b_btag_csv3"),"btag_{csv} discr. third Leading Jet; third-Leading Jet btag_{csv} discr.",60,0.,1.);

	histo_["jet_b_btag_cmva1"]	=  new TH1D( ("jet_b_btag_cmva1"),"btag_{cmva} discr. first Leading Jet; Leading Jet btag_{cmva} discr.",30,0.5,1.);
	histo_["jet_b_btag_cmva2"]	=  new TH1D( ("jet_b_btag_cmva2"),"btag_{cmva} discr. second Leading Jet; sub-Leading Jet btag_{cmva} discr.",30,0.5,1.);
	histo_["jet_b_btag_cmva3"]	=  new TH1D( ("jet_b_btag_cmva3"),"btag_{cmva} discr. third Leading Jet; third-Leading Jet btag_{cmva} discr.",60,0.,1.);

	histo_["jet_b_deta12"]		=  new TH1D( ("jet_b_deta12"),"#Delta #eta between Leading and sub-Leading jets;#Delta #eta_{12}",20,-1.6,1.6);
	histo_["jet_b_dphi12"]		=  new TH1D("jet_b_dphi12","#Delta #phi between Leading and sub-Leading jets;#Delta #phi_{12}",20,0,TMath::Pi());
	histo_["jet_b_deta13"]		=  new TH1D("jet_b_deta13","#Delta #eta between Leading and third-Leading jets;#Delta #eta_{13}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_b_dphi13"]		=  new TH1D("jet_b_dphi13","#Delta #phi between Leading and third-Leading jets;#Delta #phi_{13}",(int)size/2,0,TMath::Pi());
	histo_["jet_b_deta23"]		=  new TH1D("jet_b_deta23","#Delta #eta between sub-Leading and third-Leading jets;#Delta #eta_{23}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
	histo_["jet_b_dphi23"]		=  new TH1D("jet_b_dphi23","#Delta #phi between sub-Leading and third-Leading jets;#Delta #phi_{23}",(int)size/2,0,TMath::Pi());

	histo_["jet_b_dR12"]			=  new TH1D( ("jet_b_dR12"),"#Delta R between Leading and sub-Leading jets;#Delta R_{12}",30,1,3.6);
	histo_["jet_b_dR13"]			=  new TH1D( ("jet_b_dR13"),"#Delta R between Leading and third-Leading jets;#Delta R_{13}",30,1,3.6);
	histo_["jet_b_dR23"]			=  new TH1D( ("jet_b_dR23"),"#Delta R between sub-Leading and third-Leading jets;#Delta R_{23}",30,1,3.6);

	histo_["diJet_b_pt"]			=  new TH1D( ("diJet_b_pt"),"p_{T} of the di-Jet object;di-Jet p_{T}, [GeV]",30,0.,600.);
	histo_["diJet_b_eta"]			=  new TH1D( ("diJet_b_eta"),"#eta of the di-Jet object; di-Jet #eta",(int)size/4,-2.5,2.5);
	histo_["diJet_b_phi"]			=  new TH1D( ("diJet_b_phi"),"#phi of the di-Jet object; di-Jet #phi",(int)size/2,-2.*TMath::Pi(),2.*TMath::Pi());
	histo_["diJet_b_m"]			=  new TH1D( ("diJet_b_m"),"M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_bb"]		=  new TH1D("diJet_b_m_bb","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_cc"]		=  new TH1D("diJet_b_m_cc","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_ll"]		=  new TH1D("diJet_b_m_ll","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_bc"]		=  new TH1D("diJet_b_m_bc","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_cl"]		=  new TH1D("diJet_b_m_cl","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);
	histo_["diJet_b_m_bl"]		=  new TH1D("diJet_b_m_bl","M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);

	histo_["NumberOfJets_b"]		= new TH1D( "NumberOfJets_b","Total Number Of jets; Jet multiplicity",sizeof(NJetsBins)/sizeof(double)-1,NJetsBins);


	//Systematic list:
	std::vector<std::string> Syst = {"_PtEff_","_PU_","_SFb_","_SFl_","_JES_","_JER_","_OnlBTag_", "_PDF_"};
	//Variation list:
	std::vector<std::string> Variation = {"up","down"};

	for(const std::string v : Variation){
		for(const std::string s : Syst){
			histo_["jet_b_pt1" + s + v]			=  new TH1D( ("jet_b_pt1" + s + v).c_str(),"p_{T} of the first Leading Jet; Leading Jet p_{T}, [GeV]",sizeof(PtBins)/sizeof(double)-1,PtBins);
			histo_["jet_b_pt2" + s + v]			=  new TH1D( ("jet_b_pt2" + s + v).c_str(),"p_{T} of the second Leading Jet; sub-Leading Jet p_{T}, [GeV]",sizeof(PtBins)/sizeof(double)-1,PtBins);
			histo_["jet_b_pt3" + s + v]			=  new TH1D( ("jet_b_pt3" + s + v).c_str(),"p_{T} of the third Leading Jet; third-Leading Jet p_{T}, [GeV]",sizeof(Pt3Bins)/sizeof(double)-1,Pt3Bins);

			histo_["jet_b_Ht" + s + v]			= new TH1D( ("jet_b_Ht" + s + v).c_str(),"H_{T} of the jets; H_{T}, [GeV]",sizeof(HtBins)/sizeof(double)-1,HtBins);

			histo_["jet12_b_assym" + s + v]			= new TH1D( ("jet12_b_assym" + s + v).c_str(),"p_{T} Leading - sub-Leading jet assymetry; #frac{p^{(1)}_{T} - p^{(2)}_{T}}{p^{(1)}_{T} + p^{(2)}_{T}}",sizeof(PtAssymBins)/sizeof(double)-1,PtAssymBins);

			histo_["jet_b_eta1" + s + v]			=  new TH1D( ("jet_b_eta1" + s + v).c_str(),"#eta of the first Leading Jet; Leading Jet #eta ",30,-2.2,2.2);
			histo_["jet_b_eta2" + s + v]			=  new TH1D( ("jet_b_eta2" + s + v).c_str(),"#eta of the second Leading Jet; sub-Leading Jet #eta ",30,-2.2,2.2);
			histo_["jet_b_eta3" + s + v]			=  new TH1D( ("jet_b_eta3" + s + v).c_str(),"#eta of the third Leading Jet; third-Leading Jet #eta ",30,-2.2,2.2);

			histo_["jet_b_phi1" + s + v]			=  new TH1D( ("jet_b_phi1" + s + v).c_str(),"#phi of the first Leading Jet; Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
			histo_["jet_b_phi2" + s + v]			=  new TH1D( ("jet_b_phi2" + s + v).c_str(),"#phi of the second Leading Jet; sub-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
			histo_["jet_b_phi3" + s + v]			=  new TH1D( ("jet_b_phi3" + s + v).c_str(),"#phi of the third Leading Jet; third-Leading Jet #phi ",(int)size/2,-1.* TMath::Pi(),TMath::Pi());

			histo_["jet_b_btag_csv1" + s + v]		=  new TH1D( ("jet_b_btag_csv1" + s + v).c_str(),"btag_{csv} discr. first Leading Jet; Leading Jet btag_{csv} discr.",30,0.935,1.);
			histo_["jet_b_btag_csv2" + s + v]		=  new TH1D( ("jet_b_btag_csv2" + s + v).c_str(),"btag_{csv} discr. second Leading Jet; Leading Jet btag_{csv} discr.",30,0.935,1.);
			histo_["jet_b_btag_csv3" + s + v]		=  new TH1D( ("jet_b_btag_csv3" + s + v).c_str(),"btag_{csv} discr. third Leading Jet; Leading Jet btag_{csv} discr.",60,0.,1.);

			histo_["jet_b_btag_cmva1" + s + v]	=  new TH1D( ("jet_b_btag_cmva1" + s + v).c_str(),"btag_{cmva} discr. first Leading Jet; Leading Jet btag_{cmva} discr.",30,0.5,1.);
			histo_["jet_b_btag_cmva2" + s + v]	=  new TH1D( ("jet_b_btag_cmva2" + s + v).c_str(),"btag_{cmva} discr. second Leading Jet; Leading Jet btag_{cmva} discr.",30,0.5,1.);
			histo_["jet_b_btag_cmva3" + s + v]	=  new TH1D( ("jet_b_btag_cmva3" + s + v).c_str(),"btag_{cmva} discr. third Leading Jet; Leading Jet btag_{cmva} discr.",60,0.,1.);

			histo_["jet_b_deta12" + s + v]		=  new TH1D( ("jet_b_deta12" + s + v).c_str(),"#Delta #eta between Leading and sub-Leading jets;#Delta #eta_{12}",20,-1.6,1.6);
			histo_["jet_b_dphi12" + s + v]		=  new TH1D( ("jet_b_dphi12" + s + v).c_str(),"#Delta #phi between Leading and sub-Leading jets;#Delta #phi_{12}",20,0,TMath::Pi());
			histo_["jet_b_deta13" + s + v]		=  new TH1D( ("jet_b_deta13" + s + v).c_str(),"#Delta #eta between Leading and third-Leading jets;#Delta #eta_{13}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
			histo_["jet_b_dphi13" + s + v]		=  new TH1D( ("jet_b_dphi13" + s + v).c_str(),"#Delta #phi between Leading and third-Leading jets;#Delta #phi_{13}",(int)size/2,0,TMath::Pi());
			histo_["jet_b_deta23" + s + v]		=  new TH1D( ("jet_b_deta23" + s + v).c_str(),"#Delta #eta between sub-Leading and third-Leading jets;#Delta #eta_{23}",(int)size/2,-1.* TMath::Pi(),TMath::Pi());
			histo_["jet_b_dphi23" + s + v]		=  new TH1D( ("jet_b_dphi23" + s + v).c_str(),"#Delta #phi between sub-Leading and third-Leading jets;#Delta #phi_{23}",(int)size/2,0,TMath::Pi());

			histo_["jet_b_dR12" + s + v]		=  new TH1D( ("jet_b_dR12" + s + v).c_str(),"#Delta R between Leading and sub-Leading jets;#Delta R_{12}",30,1,3.6);
			histo_["jet_b_dR13" + s + v]		=  new TH1D( ("jet_b_dR13" + s + v).c_str(),"#Delta R between Leading and sub-sub-Leading jets;#Delta R_{13}",30,1,3.6);
			histo_["jet_b_dR23" + s + v]		=  new TH1D( ("jet_b_dR23" + s + v).c_str(),"#Delta R between second-Leading and third-Leading jets;#Delta R_{23}",30,1,3.6);

			histo_["diJet_b_pt" + s + v]		=  new TH1D( ("diJet_b_pt" + s + v).c_str(),"p_{T} of the di-Jet object;di-Jet p_{T}, [GeV]",30,0.,600.);
			histo_["diJet_b_eta" + s + v]		=  new TH1D( ("diJet_b_eta" + s + v).c_str(),"#eta of the di-Jet object; di-Jet #eta",(int)size/4,-2.5,2.5);
			histo_["diJet_b_phi" + s + v]		=  new TH1D( ("diJet_b_phi" + s + v).c_str(),"#phi of the di-Jet object; di-Jet #phi",(int)size/2,-2.*TMath::Pi(),2.*TMath::Pi());
			histo_["diJet_b_m" + s + v]			=  new TH1D( ("diJet_b_m" + s + v).c_str(),"M_{12} of the di-Jet object; di-Jet M_{12}, [GeV]",sizeof(M12Bins)/sizeof(double)-1,M12Bins);

			histo_["NumberOfJets_b" + s + v]		= new TH1D( ("NumberOfJets_b" + s + v).c_str(),"Total Number Of jets; Jet multiplicity",sizeof(NJetsBins)/sizeof(double)-1,NJetsBins);
		}
	}
}

std::map<std::string, TH1* >& Histograms::getHisto() {
	return histo_;
}

std::map<std::string, TH2*>& Histograms::getHisto2D() {
	return histo2D_;
}

void Histograms::MakeM12Templates(const std::size_t& nbins, const double& xmin, const double& xmax){

	std::array<std::string,17> syst = {{"Mbb","SFb_down","SFb_up","SFl_down","SFl_up","JES_down","JES_up","JER_down","JER_up","PU_down","PU_up","PtEff_down","PtEff_up","BTagEff_up","BTagEff_down", "PDF_up", "PDF_down"}};
	for(const auto& h: syst){
		std::string name = "template_" + h;
		if(histo_[name] != nullptr) throw std::invalid_argument("Histogram: " + name + " already exists. Check spelling");
		histo_[name] = new TH1D( name.c_str(),"",nbins,xmin,xmax);
	}
}
