#pragma once

#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TKey.h>
#include <TDirectory.h>

namespace treeUtils {

    void treeMerging(const char *inputFileName, const char *treeName, TFile *outputFile)
    {

        TFile *inputFile = TFile::Open(inputFileName, "READ");
        TList *treeList = new TList();
        TIter nextDir(inputFile->GetListOfKeys());
        TKey *key;
        while ((key = (TKey *)nextDir()))
        {
            std::cout << "Reading directory: " << key->GetName() << std::endl;
            TObject *obj = key->ReadObj();

            if (obj->InheritsFrom(TDirectory::Class()))
            {
                TDirectory *dir = (TDirectory *)obj;
                TTree *tmpTree = (TTree *)dir->Get(treeName);
                treeList->Add(tmpTree);
            }
            else
            {
                std::cerr << "Missing trees in directory: " << key->GetName() << std::endl;
            }
        }

        outputFile->cd();
        TTree *tree = TTree::MergeTrees(treeList);
        tree->Write();
        inputFile->Close();

        std::cout << "Merged tree written to file: " << outputFile->GetName() << std::endl;
    }

}   // namespace treeHandling
