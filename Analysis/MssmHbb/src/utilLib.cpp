#include "Analysis/MssmHbb/interface/utilLib.h"

bool is_empty(std::ifstream& pFile)
{
	/*
	 * Check whether file is empty or not
	 */
    return pFile.peek() == std::ifstream::traits_type::eof();
}

bool findStrings(const std::string & input, const std::string & needful)
{
	std::string input1 = input;
	std::string input2 = needful;
	std::transform(input1.begin(),input1.end(),input1.begin(),::tolower);
	std::transform(input2.begin(),input2.end(),input2.begin(),::tolower);
	if(input1.find(input2) != std::string::npos) return true;
	else return false;
}

int returnMassPoint(const std::string& name){
	int Mpoint = 0;
	std::string MassPos = "_M-";
	auto p1 = name.std::string::find(MassPos) + 3;
	if(p1 == 3) {
		std::cerr<<"FileNames were changed!!!!"<<std::endl;
		exit(1);
	}
	auto p2 = name.std::string::find("_",p1);
	std::string MpointString = name.std::string::substr(p1,size_t(p2-p1));
	Mpoint = std::stoi(MpointString);
	return Mpoint;
}

void CheckZombie(const TFile& file){
	if(file.IsZombie()){
		throw std::invalid_argument("ERROR: TFile " + std::string(file.GetName()) + " is Zombie");
	}
}

void CheckZombieObjectInTFile(const TFile& file, const std::string& name){
	if(! file.GetListOfKeys()->Contains(name.c_str()))
		throw std::invalid_argument("ERROR: Object: " + name + " is not in TFile " + std::string(file.GetName()));
}

RooWorkspace* GetRooWorkspace(const std::string& path_to_file, const std::string& workspace_name){
	/*
	 * Function to get RooWorkspace from TFile
	 */
	TFile f(path_to_file.c_str(),"READ");
	CheckZombie(f);
	CheckZombieObjectInTFile(f,workspace_name);
	RooWorkspace *w = static_cast<RooWorkspace*>( f.Get(workspace_name.c_str()));
	return w;
}

void CheckOutputDir(const std::string& oDir){
	/*
	 * Function to check an output path
	 * and create a dir if doesn't exist
	 */
	try {
		boost::filesystem::path opath(oDir);
		if(!boost::filesystem::exists(opath)) {
			boost::filesystem::create_directory(opath);
			std::cout << "Creating output directory : " << oDir << std::endl;
		}
	}
	catch (const boost::filesystem::filesystem_error& e){
		//If folder doesn't exist - go one higher and try it
		std::string parent_folder = oDir.std::string::substr(0,oDir.find_last_of('/'));
		if(oDir[oDir.length()] == '/')  parent_folder = oDir.std::string::substr(0,oDir.find_last_of("/"));
		CheckOutputDir(parent_folder);
		CheckOutputDir(oDir);
	}
}
