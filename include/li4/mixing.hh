#pragma once

#include <vector>
#include <Riostream.h>
#include <TTree.h>

#include "histograms.hh"
#include "../candidates.hh"
#include "../indexTableUtils.hh"
#include "li4candidates.hh"

namespace mixing 
{
    enum MixingStrategy {
        kEvent = 0,
        kRotation = 1
    };

    std::vector<std::vector<CollHadBracket>> fillParticlesFromTree(TTree* inputCollisionTree, TTree* inputCandidateTree, 
                                                                   std::vector<HadCandidate>& hadrons, std::vector<He3Candidate>& he3s,
                                                                   std::vector<CollisionCandidate>& collisions, HistogramsQA& histQA) {

        HistVertexMultiplicity hVertexMultiplicity;
        std::vector<std::vector<CollHadBracket>> collisionBracket;
        collisionBracket.resize(hVertexMultiplicity.mZetaBins * hVertexMultiplicity.mMultBins + 1);

        CollisionCandidate collCand;
        CollisionCandidate collCandPrev;
        CollHadBracket collBracket;

        He3Candidate he3Cand;
        HadCandidate hadCand;

        collCand.setBranchAddress(inputCollisionTree);
        he3Cand.setBranchAddress(inputCandidateTree);
        hadCand.setBranchAddress(inputCandidateTree);


        for (int iEntry = 0; iEntry < inputCollisionTree->GetEntries(); iEntry++)
        {
            inputCollisionTree->GetEntry(iEntry);
            inputCandidateTree->GetEntry(iEntry);

            hadCand.fZHad = collCand.fZVertex;
            hadCand.fCentralityFT0C = collCand.fCentralityFT0C;
            hadrons.emplace_back(hadCand);

            if (he3Cand.fPtHe3 < 0.) {
                if (hadCand.fPtHad < 0.) {
                    histQA.hInvMassBeforeEMLikeSign->Fill(Li4Candidate::li4InvMass(he3Cand, hadCand));
                } else {
                    histQA.hInvMassBeforeEMUnlikeSign->Fill(Li4Candidate::li4InvMass(he3Cand, hadCand));
                }
            }
            histQA.hHe3BeforeEMAll->Fill(he3Cand.fPtHe3);

            if (abs(collCandPrev.fZVertex - collCand.fZVertex) < 1e-5)
            {
                collBracket.SetMax(hadrons.size() - 1);
                continue;
            }

            if (collBracket.GetMin() != -1)
            {
                int iBin = hVertexMultiplicity.getBinIndex(collCandPrev.fZVertex, collCandPrev.fCentralityFT0C);
                collBracket.CollID = iEntry - 1;
                collisionBracket[iBin].push_back(collBracket);
            }

            // a new collision has been found, dumping collision and he3 candidates

            he3s.emplace_back(he3Cand); 
            collisions.emplace_back(collCand);
            histQA.hHe3BeforeEM->Fill(he3Cand.fPtHe3);

            collBracket.SetMin(hadrons.size() - 1);
            collBracket.SetMax(hadrons.size() - 1);
            collCandPrev = collCand;
        }

        std::cout << "--------------------------------" << std::endl;
        std::cout << "Filled candidates!" << std::endl;
        std::cout << "Size of Hadron Candidates to be mixed: " << hadrons.size() << std::endl;
        std::cout << "Size of He3 Candidates to be mixed: " << he3s.size() << std::endl;
        std::cout << "Size of Coll Candidates: " << collisions.size() << std::endl;
        std::cout << "--------------------------------" << std::endl;

        return collisionBracket;
    }

}   // namespace mixing

class Mixer
{
    public:
        Mixer() = default;
        Mixer(const std::vector<HadCandidate>& hadrons, const std::vector<He3Candidate>& he3s, 
              const std::vector<CollisionCandidate>& collisions, 
              const std::vector<std::vector<CollHadBracket>>& collisionBrackets,
              const int mixingDepth = 5)
            : fHadrons(hadrons), fHe3s(he3s), fCollisions(collisions), fCollisionBrackets(collisionBrackets),
             fMixingDepth(mixingDepth) {}
        ~Mixer() = default;

        void performEventMixing(TTree* outputTree, HistogramsQA& histQA);
        void performAngleMixing(TTree* outputTree, HistogramsQA& histQA);

    private:
        std::vector<HadCandidate> fHadrons;
        std::vector<He3Candidate> fHe3s;
        std::vector<CollisionCandidate> fCollisions;
        std::vector<std::vector<CollHadBracket>> fCollisionBrackets;
        int fMixingDepth = 5;
};

void Mixer::performEventMixing(TTree* outputTree, HistogramsQA& histQA)
{
    Li4Candidate li4Candidate;

    std::vector<int> hadronProcessTimes(fHadrons.size(), 0);
    const int maxProcessTimes = 10;

    HistVertexMultiplicity hVertexMultiplicity;
    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Starting event mixing with " << fHadrons.size() << " hadrons and " 
              << fHe3s.size() << " He3 candidates." << std::endl;

    for (size_t iHe3 = 0; iHe3 < fHe3s.size(); iHe3++)
    {
        if (iHe3 % (fHe3s.size()/100) == 0) {
            std::cout << "Processing He3 candidate " << iHe3 << " / " << fHe3s.size() << " (" 
                      << static_cast<float>(iHe3)/fHe3s.size()*100 << "%)" << std::endl;
        }

        const He3Candidate& he3Cand = fHe3s[iHe3];
        const CollisionCandidate& collCand = fCollisions[iHe3];
        li4Candidate.setHe3(he3Cand);
        int iBin = hVertexMultiplicity.getBinIndex(collCand.fZVertex, collCand.fCentralityFT0C);
        histQA.hHe3Unique->Fill(he3Cand.fPtHe3);

        for (size_t iDepth = 0; static_cast<int>(iDepth) < fMixingDepth; iDepth++)
        {

            if (fCollisionBrackets[iBin].size() == 0 || iDepth >= fCollisionBrackets[iBin].size())
            {
                break;
            }

            int iCollEM;
            iCollEM = gRandom->Integer(fCollisionBrackets[iBin].size());
            int collIDHad = fCollisionBrackets[iBin][iCollEM].CollID;
            if (collIDHad == static_cast<int>(iHe3))
            {
                continue;
            }
            for (int iHad = fCollisionBrackets[iBin][iCollEM].GetMin(); iHad <= fCollisionBrackets[iBin][iCollEM].GetMax(); iHad++)
            {
                hadronProcessTimes[iHad]++;
                if (hadronProcessTimes[iHad] > maxProcessTimes)
                {
                    continue;
                }

                const auto& hadCand = fHadrons[iHad];
                li4Candidate.setHad(hadCand);
                li4Candidate.setZVertex(collCand.fZVertex);
                li4Candidate.setCentralityFT0C(collCand.fCentralityFT0C);
                if (!true){ // conditions on the li4 pair
                    continue;
                }

                if (li4Candidate.getPtHe3() < 0) {
                    if (li4Candidate.getPtHad() < 0.) {
                        histQA.hInvMassAfterEMLikeSign->Fill(li4Candidate.calcInvMass());
                    } else {
                        histQA.hInvMassAfterEMUnlikeSign->Fill(li4Candidate.calcInvMass());
                    }
                }
                histQA.hHe3AfterEM->Fill(he3Cand.fPtHe3);
                outputTree->Fill();
            }
        }
    }
}

void Mixer::performAngleMixing(TTree* outputTree, HistogramsQA& histQA)
{
    Li4Candidate li4Candidate;

    HistVertexMultiplicity hVertexMultiplicity;
    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Starting angle mixing with " << fHadrons.size() << " hadrons and " 
              << fHe3s.size() << " He3 candidates." << std::endl;

    for (size_t iHe3 = 0; iHe3 < fHe3s.size(); iHe3++)
    {
        if (iHe3 % (fHe3s.size()/100) == 0) {
            std::cout << "Processing He3 candidate " << iHe3 << " / " << fHe3s.size() << " (" 
                      << iHe3/fHe3s.size()*100 << "%)" << std::endl;
        }

        const He3Candidate& he3Cand = fHe3s[iHe3];
        const CollisionCandidate& collCand = fCollisions[iHe3];
        li4Candidate.setHe3(he3Cand);
        int iBin = hVertexMultiplicity.getBinIndex(collCand.fZVertex, collCand.fCentralityFT0C);
        histQA.hHe3Unique->Fill(he3Cand.fPtHe3);

        int iBracketIdx = 0;
        for (size_t iBracketIdx = 0; iBracketIdx < fCollisionBrackets[iBin].size(); ++iBracketIdx) {
            if (fCollisionBrackets[iBin][iBracketIdx].CollID == he3Cand.CollID) {
                break;
            }
        }

        for (int iHad = fCollisionBrackets[iBin][iBracketIdx].GetMin(); iHad <= fCollisionBrackets[iBin][iBracketIdx].GetMax(); iHad++) {
            
            li4Candidate.setHad(fHadrons[iHad]);
            if (!true){ // conditions on the li4 pair
                continue;
            }
            if (li4Candidate.getPtHe3() < 0) {
                if (li4Candidate.getPtHad() < 0.) {
                    histQA.hInvMassAfterEMLikeSign->Fill(li4Candidate.calcInvMass());
                } else {
                    histQA.hInvMassAfterEMUnlikeSign->Fill(li4Candidate.calcInvMass());
                }
            }
            histQA.hHe3AfterEM->Fill(he3Cand.fPtHe3);
            outputTree->Fill();
        }
    }
}
