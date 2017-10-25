using namespace std;

const string cmsswBase = getenv("CMSSW_BASE");
const pair<string,string> GetInputFileNames(const int& mass);
void createAsimov(const int& mass, const string& out, const int& n_mu = 1);

int SpBg_asimov(){
	/*
	 * Macro to create S+Bg Asimov data manually
	 */
	const int mass = 500;
	const string output = cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_bg_fit/sr1/background_SpBg_Asimov.root";
	createAsimov(500,output,2);


	return 0;
}

void createAsimov(const int& mass, const string& out, const int& n_mu){
	/*
	 * Method to create workspace with Asimov S+Bg data
	 */
	pair<string,string> path = GetInputFileNames(mass); //first - signal, second - BG

	//Read current data_obs
	TFile *fBg = new TFile(path.second.c_str(),"READ");
	RooWorkspace *wBg = (RooWorkspace*) fBg->Get("workspace");
	TH1 *h_data_old = wBg->pdf("background")->createHistogram("data_obs",*wBg->var("mbb"),RooFit::Binning(5000,wBg->var("mbb")->getMin(),wBg->var("mbb")->getMax()));
	h_data_old->Scale(wBg->var("background_norm")->getValV()/h_data_old->Integral());

	//Read signal pdf and generate histogram
	TFile *fSg = new TFile(path.first.c_str(),"READ");
	RooWorkspace *wSg = (RooWorkspace*) fSg->Get("workspace");
	TH1 *hSg = wSg->pdf("signal")->createHistogram("signal_obs",*wSg->var("mbb"),RooFit::Binning(5000,wSg->var("mbb")->getMin(),wSg->var("mbb")->getMax()));
	hSg->Scale(n_mu *( wSg->function("signal_norm")->getValV()) / hSg->Integral());
	cout<<"norm3"<<endl;

	//Sum two this histos
	h_data_old->Add(hSg);
	cout<<"norm4"<<endl;

	//New data_obs
	RooDataHist *data_obs = new RooDataHist("data_obs","data_obs",(RooRealVar&) *wBg->var("mbb"),RooFit::Import(*h_data_old));
	RooRealVar bg_norm("background_norm","background_norm",h_data_old->Integral());

	//Write to the file
	RooWorkspace wOut("workspace");
	wOut.import(*(RooAbsPdf*)wBg->pdf("background"));
	wOut.import(*data_obs);
	wOut.import(bg_norm);
	wOut.Print("v");
	wOut.writeToFile(out.c_str());

	//Checker:
	cout<<"Sanity check:\n";
	cout<<"Asimov data: old = "<<wBg->var("background_norm")->getValV()<<" signal = "<<n_mu * wSg->function("signal_norm")->getValV()<<" new = "<<h_data_old->Integral()<<endl;


}

const pair<string,string> GetInputFileNames(const int& mass){
	/*
	 * Method to get input Bg root file name and signal root file name
	 * according to provided mass
	 */
	string signal_temp	= cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_signal_M-";
	string bg_temp		= cmsswBase + "/src/Analysis/MssmHbb/output/ReReco_bg_fit/";

	//Construct signal smaple:
	signal_temp += to_string(mass) + "/workspace/signal_workspace.root";
	if(mass < 600) 							bg_temp += "/sr1/background_workspace_TurnOnFix_5000bins.root";
	else if (mass >= 600 && mass <= 900)	bg_temp += "/sr2/background_workspace_5000bins.root";
	else 									bg_temp += "/sr2/background_workspace_5000bins.root";

	if(gSystem->AccessPathName(signal_temp.c_str())) throw invalid_argument("ERROR: at SpBg_asimov::GetInputFileNames - no signal file at: " + signal_temp);
	if(gSystem->AccessPathName(bg_temp.c_str())) throw invalid_argument("ERROR: at SpBg_asimov::GetInputFileNames - no background file at: " + bg_temp);

	return make_pair(signal_temp,bg_temp);
}
