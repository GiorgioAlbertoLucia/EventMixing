#pragma once

#include <vector>
#include <unordered_set>

#include <Riostream.h>
#include <TTree.h>

#include "histograms.hh"
#include "../core/candidates.hh"
#include "../core/indexTableUtils.hh"
#include "li4candidates.hh"
#include "selections.h"

namespace mixing 
{
    enum MixingStrategy {
        kEvent = 0,
        kRotation = 1
    };

    bool preliminaryCuts(const He3Candidate& he3, const HadCandidate& had, const CollisionCandidate& collision, const bool is23) {
        
        const double pthe3 = (he3.fPIDtrkHe3 == 7) || (he3.fPIDtrkHe3 == 8) || (std::abs(he3.fPtHe3) > 2.5) ? std::abs(he3.fPtHe3) : CorrectPidTrkHe(std::abs(he3.fPtHe3));
        const double fNSigmaDCAxyHe3 = ComputeNsigmaDCAxyHe(std::abs(pthe3), he3.fDCAxyHe3);
        const double fNSigmaDCAzHe3 = ComputeNsigmaDCAzHe(std::abs(pthe3), he3.fDCAzHe3);
        const double fNSigmaDCAxyHad = ComputeNsigmaDCAxyPr(std::abs(had.fPtHad), had.fDCAxyHad);
        const double fNSigmaDCAzHad = ComputeNsigmaDCAzPr(std::abs(had.fPtHad), had.fDCAzHad);

        if ((std::abs(pthe3) < 2.5 || he3.fPIDtrkHe3 == 7) &&
            std::abs(had.fNSigmaTPCHad) < 2.5 &&
            std::abs(fNSigmaDCAxyHe3) < 4 &&
            std::abs(fNSigmaDCAzHe3) < 4 &&
            std::abs(fNSigmaDCAxyHad) < 4 &&
            std::abs(fNSigmaDCAzHad) < 4 &&
            std::abs(he3.fEtaHe3) < 0.9 &&
            std::abs(had.fEtaHad) < 0.9 &&
            ((he3.fChi2TPCHe3 > 0.5) || (is23 == false)) && 
            (he3.fChi2TPCHe3 < 4) &&
            (had.fChi2TPCHad < 4))

            return true;
        
        return false;
    }

    std::vector<std::vector<CollHadBracket>> fillParticlesFromTree(TTree* inputCollisionTree, TTree* inputCandidateTree, 
                                                                   std::vector<HadCandidate>& hadrons, std::vector<He3Candidate>& he3s,
                                                                   std::vector<CollisionCandidate>& collisions, HistogramsQA& histQA,
                                                                   const bool applyCuts = false, const bool is23 = false) {

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

            if (applyCuts)
                if (!preliminaryCuts(he3Cand, hadCand, collCand, is23))
                    continue;

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
              const int mixingDepth = 5, const bool  is23 = false)
            : fHadrons(hadrons), fHe3s(he3s), fCollisions(collisions), fCollisionBrackets(collisionBrackets),
             fMixingDepth(mixingDepth), fIs23(is23) {}
        ~Mixer() = default;

        void performEventMixing(TTree* outputTree, HistogramsQA& histQA);
        void performAngleMixing(TTree* outputTree, HistogramsQA& histQA);

    private:
        std::vector<HadCandidate> fHadrons;
        std::vector<He3Candidate> fHe3s;
        std::vector<CollisionCandidate> fCollisions;
        std::vector<std::vector<CollHadBracket>> fCollisionBrackets;
        int fMixingDepth = 5;
        bool fIs23 = false;
};

void Mixer::performEventMixing(TTree* outputTree, HistogramsQA& histQA)
{
    Li4Candidate li4Candidate;
    li4Candidate.setBranch(outputTree);

    std::vector<int> hadronProcessTimes(fHadrons.size(), 0);
    const int maxProcessTimes = 10;

    HistVertexMultiplicity hVertexMultiplicity;
    const float tolerance = 1e-7;
    
    //std::vector<std::tuple<float, float, float, float, float, float>> processedCandidates; // pt, eta, phi for both He3 and hadron
    
    auto floatBitsHash = [](float v) -> size_t {
        uint32_t bits;
        std::memcpy(&bits, &v, sizeof(bits));
        return std::hash<uint32_t>{}(bits);
    };
    auto roundFloat = [](float v, float tol) -> float {
        return std::round(v / tol) * tol;
    };
    auto pairKey = [&](float pt3, float eta3, float phi3, float ptH, float etaH, float phiH) -> size_t {
        size_t seed = 0;
        for (float v : {roundFloat(pt3, tolerance), roundFloat(eta3, tolerance), roundFloat(phi3, tolerance),
                        roundFloat(ptH, tolerance), roundFloat(etaH, tolerance), roundFloat(phiH, tolerance)}) {
            seed ^= floatBitsHash(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    };
    std::unordered_set<size_t> processedCandidates;

    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Starting event mixing with " << fHadrons.size() << " hadrons and " 
              << fHe3s.size() << " He3 candidates." << std::endl;

    for (size_t iHe3 = 0; iHe3 < fHe3s.size(); iHe3++)
    {
        if (iHe3 % (fHe3s.size()/100) == 0) {
            std::cout << "Processing He3 candidate " << iHe3 << " / " << fHe3s.size() << " (" 
                      << static_cast<int>(static_cast<float>(iHe3)/fHe3s.size()*100) << "%)" << std::endl;
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
                li4Candidate.setIs23(fIs23);
                if (!true){ // conditions on the li4 pair
                    continue;
                }

                // check if the pair was already filled
                //bool alreadyProcessed = false;
                //for (const auto& procCand : processedCandidates) {
                //    if (std::abs(std::get<0>(procCand) - he3Cand.fPtHe3) < tolerance &&
                //        std::abs(std::get<1>(procCand) - he3Cand.fEtaHe3) < tolerance &&
                //        std::abs(std::get<2>(procCand) - he3Cand.fPhiHe3) < tolerance &&
                //        std::abs(std::get<3>(procCand) - hadCand.fPtHad) < tolerance &&
                //        std::abs(std::get<4>(procCand) - hadCand.fEtaHad) < tolerance &&
                //        std::abs(std::get<5>(procCand) - hadCand.fPhiHad) < tolerance) {
                //        alreadyProcessed = true;
                //        break;
                //    }
                //} if (alreadyProcessed) {
                //    continue;
                //}
                //processedCandidates.emplace_back(he3Cand.fPtHe3, he3Cand.fEtaHe3, he3Cand.fPhiHe3, hadCand.fPtHad, hadCand.fEtaHad, hadCand.fPhiHad);
                
                // check if the pair was already filled
                size_t key = pairKey(he3Cand.fPtHe3, he3Cand.fEtaHe3, he3Cand.fPhiHe3, hadCand.fPtHad, hadCand.fEtaHad, hadCand.fPhiHad);
                if (!processedCandidates.insert(key).second) {
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
    li4Candidate.setBranch(outputTree);

    HistVertexMultiplicity hVertexMultiplicity;
    
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Starting angle mixing with " << fHadrons.size() << " hadrons and " 
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

        int iBracketIdx = 0;
        for (size_t iBracketIdx = 0; iBracketIdx < fCollisionBrackets[iBin].size(); ++iBracketIdx) {
            if (fCollisionBrackets[iBin][iBracketIdx].CollID == he3Cand.CollID) {
                break;
            }
        }

        for (int iHad = fCollisionBrackets[iBin][iBracketIdx].GetMin(); iHad <= fCollisionBrackets[iBin][iBracketIdx].GetMax(); iHad++) {
            
            li4Candidate.setHad(fHadrons[iHad]);
            li4Candidate.setIs23(fIs23);
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
