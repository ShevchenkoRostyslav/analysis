

int compareJEC(){

	TFile *miniAODjecTree = new TFile("output_miniaod.root");
	TFile *newJecTree     = new TFile("output_new_jec.root");

	TTree *miniAOD, *reJec;

	miniAODjecTree -> GetObject("MssmHbb/Events/slimmedJetsPuppi",miniAOD);
	newJecTree     -> GetObject("MssmHbb/Events/patJetsReapplyJEC",reJec);

	miniAOD->Draw("e[0]>>oldJecTH");
	reJec->Draw("e[0]>>newJecTH","","E same");

	return 0;
}
