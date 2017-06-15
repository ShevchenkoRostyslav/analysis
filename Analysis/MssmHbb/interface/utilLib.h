
#ifndef MSSMHBB_INTERFACE_UTILLIB_H_
#define MSSMHBB_INTERFACE_UTILLIB_H_

#include <sys/stat.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cctype>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"

#include "RooFit.h"
#include "RooWorkspace.h"


//include boost
#include <boost/filesystem.hpp>

template <typename T1, typename T2>
struct less_second {
    typedef std::pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second < b.second;
    }
};

template <typename T1, typename T2>
struct great_second {
    typedef std::pair<T1, T2> type;
    bool operator ()(type const& a, type const& b) const {
        return a.second > b.second;
    }
};

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

// Function to check whether file exists or not
inline bool file_exists(const std::string&name);

bool findStrings(const std::string & input, const std::string & needful);
//Function to check whether root file isSombie or not
void CheckZombie(const TFile& name);
void CheckZombieObjectInTFile(const TFile& file, const std::string& name);

//Function to get a Pointer to the workspace in the TFile
RooWorkspace* GetRooWorkspace(const std::string& path_to_file, const std::string& workspace_name = "workspace");

int returnMassPoint(const std::string& name);

//Function to check an output path and create a dir if doesn't exist
void CheckOutputDir(const std::string& oDir);

inline bool file_exists(const std::string&name){
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

//Function to Get something from the TFile. Exception safety:
template <typename T> T* GetFromTFile(TFile& file, const std::string& obj_name){
	auto *obj = static_cast<T*>(file.Get(obj_name.c_str()));
	if(!obj) throw std::invalid_argument("Invalid TObject name: " + obj_name);
	return obj;
}
template <typename T> T* GetFromTFile(const std::string& file_name, const std::string& obj_name){
	TFile f(file_name.c_str(),"READ");
	CheckZombie(f);
	auto *obj =  static_cast<T*>(GetFromTFile<T>(f,obj_name));
	return obj;
}

#endif
