#include <iostream>
#include <vector>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>

#include <TRandom3.h>

#include "../include/core/treeUtils.hh"
#include "../include/li4/li4candidates.hh"
#include "../include/li4/mixing.hh"

#include <yaml-cpp/yaml.h>

void mergeTrees(const char * inputFileName,
                const char * candidatesFileName = "inputCands.root", 
                const char * collisionsFileName = "inputColls.root",
                const char * treeNameCands = "O2he3hadtable", const char * treeNameColls = "O2he3hadmult")
{

    TFile * inputCandsFile = TFile::Open(candidatesFileName, "RECREATE");
    TFile * inputCollsFile = TFile::Open(collisionsFileName, "RECREATE");

    treeUtils::treeMerging(inputFileName, treeNameCands, inputCandsFile);
    treeUtils::treeMerging(inputFileName, treeNameColls, inputCollsFile);
    
    inputCandsFile->Close();
    inputCollsFile->Close();
}

void mixingLi4(const char * configFileName = "config/configMixingLi4.yml")
{   
    TStopwatch timer;
    HistogramsQA histQA;

    const char * candidatesFileName = "/home/galucia/EventMixing/output/inputCands.root";
    const char * collisionsFileName = "/home/galucia/EventMixing/output/inputColls.root";
    const char * candidatesTreeName = "O2he3hadtable";
    const char * collisionsTreeName = "O2he3hadmult";

    YAML::Node config = YAML::LoadFile(configFileName);
    const bool doMerge = config["doMerge"].as<bool>();
    const int mixingStrategy = config["mixingStrategy"].as<int>();
    const int mixingDepth = config["mixingDepth"].as<int>();
    const bool is23 = config["is23"].as<bool>();
    const bool applyCuts = config["applyCuts"].as<bool>();
    const int randomSeed = config["randomSeed"].as<int>();
    gRandom->SetSeed(randomSeed);

    if (doMerge) {
        std::string inputFileName = config["inputFileName"].as<std::string>();
        mergeTrees(inputFileName.c_str(), candidatesFileName, collisionsFileName, candidatesTreeName, collisionsTreeName);
    }

    TFile *inputCandsFile = TFile::Open(candidatesFileName);
    TTree *inputCandidateTree = (TTree *)inputCandsFile->Get(candidatesTreeName);
    TFile *inputCollsFile = TFile::Open(collisionsFileName);
    TTree *inputCollisionTree = (TTree *)inputCollsFile->Get(collisionsTreeName);

    std::vector<He3Candidate> he3Candidates;
    std::vector<HadCandidate> hadCandidates;
    std::vector<CollisionCandidate> collisionCandidates;
    auto collisionBrackets = mixing::fillParticlesFromTree(inputCollisionTree, inputCandidateTree, hadCandidates,
                                                           he3Candidates, collisionCandidates, histQA, applyCuts);

    inputCandsFile->Close();
    inputCollsFile->Close();

    std::string outputFileName = config["outputFileName"].as<std::string>();
    auto outputFile = TFile::Open(outputFileName.c_str(), "RECREATE"); 
    auto outputTree = new TTree("MixedTree", "MixedTree");

    timer.Start();
    Mixer mixer(hadCandidates, he3Candidates, collisionCandidates, collisionBrackets, mixingDepth, is23);
    if (mixingStrategy == mixing::MixingStrategy::kEvent) {
        mixer.performEventMixing(outputTree, histQA);
    } else if (mixingStrategy == mixing::MixingStrategy::kRotation) {
        mixer.performAngleMixing(outputTree, histQA);
    } else {
        std::cout << "Unknown mixing strategy." << std::endl;
        return;
    }
    timer.Stop();
    std::cout << "Event mixing completed in " << timer.RealTime() << " seconds." << std::endl;

    outputFile->cd();
    outputTree->Write();

    auto qaDirectory = outputFile->mkdir("HistogramsQA");
    histQA.saveHistograms(qaDirectory);
    outputFile->Close();
    
}
