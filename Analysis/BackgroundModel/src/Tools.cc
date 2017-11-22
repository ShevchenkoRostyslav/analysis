#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.
#include "Analysis/BackgroundModel/interface/Tools.h"


namespace analysis {
  namespace backgroundmodel {

    std::string uniquify(const std::string& original) {
      boost::uuids::uuid uuid = boost::uuids::random_generator()();
      std::stringstream uniqueName;
      uniqueName << original << "_" << uuid;
      return uniqueName.str();
    }


    TObject* uniqueClone(const TObject* original) {
      if (original == nullptr) return nullptr;
      return original->Clone(uniquify(original->GetName()).c_str());
    }


    std::vector<ParamModifier> parseModifiers(const std::vector<std::string>& input) {
      std::vector<ParamModifier> output;
      for (const auto& i: input) {
        std::vector<std::string> nameAndModifiers;
        boost::split(nameAndModifiers, i, boost::is_any_of(":"));
        const std::string name(nameAndModifiers.at(0));
        ParamModifier paramModifier(name);
        std::vector<std::string> modifiers;
        boost::split(modifiers, nameAndModifiers.at(1), boost::is_any_of(","));
        for (const auto& m : modifiers) {
          if (m.find("start") != std::string::npos)
            paramModifier.start(std::stof(m.substr(m.find("=") + 1)));
          if (m.find("min") != std::string::npos)
            paramModifier.min(std::stof(m.substr(m.find("=") + 1)));
          if (m.find("max") != std::string::npos)
            paramModifier.max(std::stof(m.substr(m.find("=") + 1)));
          if (m.find("constant") != std::string::npos) paramModifier.constant();
          if (m.find("floating") != std::string::npos) paramModifier.floating();
        }
        paramModifier.show();
        output.push_back(paramModifier);
      }
      return output;
    }

  }
}
