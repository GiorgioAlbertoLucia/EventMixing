#include <string>
#include <TString.h>
#include <TSystem.h>

void build(TString myopt = "fast") {
    
    gSystem->AddIncludePath((std::string("-I ")+"build").c_str());
    
    //gSystem->AddIncludePath((std::string("-I ")+"/home/galucia/local/include").c_str());
    //gSystem->Load("/home/galucia/local/lib/libyaml-cpp.so");
    
    //gSystem->AddIncludePath((std::string("-I ")+"/opt/homebrew/opt/yaml-cpp/include").c_str());
    //gSystem->Load("/opt/homebrew/opt/yaml-cpp/lib/libyaml-cpp.0.8.dylib");
    
    TString opt;
    if(myopt.Contains("force"))   opt = "kfg";
    else                          opt = "kg";
  
    gSystem->CompileMacro("src/eventMixingLi4.cxx", opt.Data(), "", "build");
}