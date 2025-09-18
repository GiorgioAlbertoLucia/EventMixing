#include <iostream>
#include <vector>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TStopwatch.h>

#include <TRandom3.h>

#include "../include/treeUtils.hh"
#include "../include/li4/li4candidates.hh"
#include "../include/li4/mixing.hh"



void mergeTrees(const char * candidatesFileName = "inputCands.root", 
                const char * collisionsFileName = "inputColls.root",
                const char * treeNameCands = "O2he3hadtable", const char * treeNameColls = "O2he3hadmult")
{
    const char * inputFileName = "/data/galucia/lithium_local/same/LHC23_PbPb_pass4_long_same_lsus.root";
    //const char * inputFileName = "/data/galucia/lithium_local/same/LHC24as_pass1_same.root";
    //const char * inputFileName = "/Users/glucia/Projects/ALICE/data/lithium/same/LHC24as_pass1_same.root";
    //const char * inputFileName = "/Users/glucia/Projects/ALICE/data/lithium/same/LHC24ag_pass1_skimmed_same.root";

    TFile * inputCandsFile = TFile::Open(candidatesFileName, "RECREATE");
    TFile * inputCollsFile = TFile::Open(collisionsFileName, "RECREATE");

    treeUtils::treeMerging(inputFileName, treeNameCands, inputCandsFile);
    treeUtils::treeMerging(inputFileName, treeNameColls, inputCollsFile);
    
    inputCandsFile->Close();
    inputCollsFile->Close();
}

void eventMixingLi4(const bool doMerge = false, const int mixingStrategy = mixing::MixingStrategy::kEvent,
                    const int mixingDepth = 2)
{   
    TStopwatch timer;
    HistogramsQA histQA;
    gRandom->SetSeed(1995);

    const char * candidatesFileName = "/home/galucia/EventMixing/output/inputCands.root";
    const char * collisionsFileName = "/home/galucia/EventMixing/output/inputColls.root";
    const char * candidatesTreeName = "O2he3hadtable";
    const char * collisionsTreeName = "O2he3hadmult";

    if (doMerge) {
        mergeTrees(candidatesFileName, collisionsFileName, candidatesTreeName, collisionsTreeName);
    }

    TFile *inputCandsFile = TFile::Open(candidatesFileName);
    TTree *inputCandidateTree = (TTree *)inputCandsFile->Get(candidatesTreeName);
    TFile *inputCollsFile = TFile::Open(collisionsFileName);
    TTree *inputCollisionTree = (TTree *)inputCollsFile->Get(collisionsTreeName);

    std::vector<He3Candidate> he3Candidates;
    std::vector<HadCandidate> hadCandidates;
    std::vector<CollisionCandidate> collisionCandidates;
    auto collisionBrackets = mixing::fillParticlesFromTree(inputCollisionTree, inputCandidateTree, hadCandidates,
                                                           he3Candidates, collisionCandidates, histQA);

    std::string strategyName = (mixingStrategy == mixing::MixingStrategy::kEvent) ? "event" : "rotation";
    auto outputFile = TFile::Open(Form("/data/galucia/lithium_local/mixing/LHC23_PbPb_pass4_%s_mixing_lsus_small.root", strategyName.c_str()), "RECREATE"); 
    //auto outputFile = TFile::Open("/data/galucia/lithium_local/mixing/LHC24as_pass1_mixing_lsus_new.root", "RECREATE");
    //auto outputFile = TFile::Open("/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24ar_pass1_mixing.root", "RECREATE");
    //auto outputFile = TFile::Open("/Users/glucia/Projects/ALICE/data/lithium/mixing/LHC24as_pass1_mixing_small.root", "RECREATE");
    //auto outputFile = TFile::Open("/Users/glucia/Projects/ALICE/data/lithium/same/LHC24ag_pass1_skimmed_mixing.root", "RECREATE");
    auto outputTree = new TTree("MixedTree", "MixedTree");

    timer.Start();
    Mixer mixer(hadCandidates, he3Candidates, collisionCandidates, collisionBrackets, mixingDepth);
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
