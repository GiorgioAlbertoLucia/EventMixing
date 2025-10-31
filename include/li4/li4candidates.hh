#pragma once

#include <TTree.h>
#include "../core/candidates.hh"
#include "../core/physics.hh"
#include "histograms.hh"

class HadCandidate: public Candidate
{
    public:

        HadCandidate() = default;
        ~HadCandidate() = default;
        HadCandidate(const HadCandidate& other) = default;
        HadCandidate& operator= (const HadCandidate& other) = default;

        float fPtHad, fEtaHad, fPhiHad, fDCAxyHad, fDCAzHad, fSignalTPCHad, fInnerParamTPCHad, fMassTOFHad;
        unsigned int fItsClusterSizeHad, fPIDtrkHad;
        unsigned char fSharedClustersHad;
        float fNSigmaTPCHad, fNSigmaTOFHad, fChi2TPCHad;
        float fZHad, fCentralityFT0C;
        int CollID = -1;

        void setBranchAddress(TTree * tree) override 
        {
            tree->SetBranchAddress("fPtHad", &fPtHad);
            tree->SetBranchAddress("fEtaHad", &fEtaHad);
            tree->SetBranchAddress("fPhiHad", &fPhiHad);
            tree->SetBranchAddress("fDCAxyHad", &fDCAxyHad);
            tree->SetBranchAddress("fDCAzHad", &fDCAzHad);
            tree->SetBranchAddress("fSignalTPCHad", &fSignalTPCHad);
            tree->SetBranchAddress("fInnerParamTPCHad", &fInnerParamTPCHad);
            tree->SetBranchAddress("fMassTOFHad", &fMassTOFHad);
            tree->SetBranchAddress("fItsClusterSizeHad", &fItsClusterSizeHad);
            tree->SetBranchAddress("fPIDtrkHad", &fPIDtrkHad);
            tree->SetBranchAddress("fSharedClustersHad", &fSharedClustersHad);
            
            //tree->SetBranchAddress("fNSigmaTPCHad", &fNSigmaTPCHad);
            tree->SetBranchAddress("fNSigmaTPCHadPr", &fNSigmaTPCHad);
            tree->SetBranchAddress("fNSigmaTOFHadPr", &fNSigmaTOFHad);
            
            tree->SetBranchAddress("fChi2TPCHad", &fChi2TPCHad);
        }
};

class He3Candidate: public Candidate
{
    public:

        He3Candidate() = default;
        ~He3Candidate() = default;
        He3Candidate(const He3Candidate& other) = default;
        He3Candidate& operator= (const He3Candidate& other) = default;

        float fPtHe3, fEtaHe3, fPhiHe3, fDCAxyHe3, fDCAzHe3, fSignalTPCHe3, fInnerParamTPCHe3, fMassTOFHe3;
        unsigned int fItsClusterSizeHe3, fPIDtrkHe3;
        unsigned char fNClsTPCHe3, fSharedClustersHe3;
        float fNSigmaTPCHe3, fChi2TPCHe3;
        float fZHe3, fCentralityFT0C;
        int CollID = -1;
        
        void setBranchAddress(TTree * tree) override 
        {
            tree->SetBranchAddress("fPtHe3", &fPtHe3);
            tree->SetBranchAddress("fEtaHe3", &fEtaHe3);
            tree->SetBranchAddress("fPhiHe3", &fPhiHe3);
            tree->SetBranchAddress("fDCAxyHe3", &fDCAxyHe3);
            tree->SetBranchAddress("fDCAzHe3", &fDCAzHe3);
            tree->SetBranchAddress("fSignalTPCHe3", &fSignalTPCHe3);
            tree->SetBranchAddress("fInnerParamTPCHe3", &fInnerParamTPCHe3);
            tree->SetBranchAddress("fMassTOFHe3", &fMassTOFHe3);
            tree->SetBranchAddress("fNClsTPCHe3", &fNClsTPCHe3);
            tree->SetBranchAddress("fItsClusterSizeHe3", &fItsClusterSizeHe3);
            tree->SetBranchAddress("fPIDtrkHe3", &fPIDtrkHe3);
            tree->SetBranchAddress("fSharedClustersHe3", &fSharedClustersHe3);
            tree->SetBranchAddress("fNSigmaTPCHe3", &fNSigmaTPCHe3);
            tree->SetBranchAddress("fChi2TPCHe3", &fChi2TPCHe3);
        }
};

class CollisionCandidate: public Candidate
{

    public:
        CollisionCandidate() = default;
        ~CollisionCandidate() = default;
        CollisionCandidate(const CollisionCandidate& other) = default;
        CollisionCandidate& operator= (const CollisionCandidate& other) = default;

        float fZVertex = -99., fCentralityFT0C;
        int CollID = -1;
        bool fIs23 = false;

        void setBranchAddress(TTree * tree) override {
            tree->SetBranchAddress("fZVertex", &fZVertex);
            tree->SetBranchAddress("fCentralityFT0C", &fCentralityFT0C);
        }

};

class Li4Candidate 
{
    public:
        
        Li4Candidate() = default;
        ~Li4Candidate() = default;
        Li4Candidate(const Li4Candidate& other) = default;
        Li4Candidate& operator= (const Li4Candidate& other) = default;

        inline void setHe3(const He3Candidate& he3);
        inline void setHad(const HadCandidate& had);
        inline void setColl(const CollisionCandidate& coll) { fColl = coll; }
        inline void setZVertex(const float z) { fColl.fZVertex = z; }
        inline void setCentralityFT0C(const float cent) { fColl.fCentralityFT0C = cent; }
        inline void setIs23(const bool is23) {  fColl.fIs23 = is23; }
        
        void setBranch(TTree* tree);
        float calcInvMass() const;
        static float li4InvMass(const He3Candidate& he3, const HadCandidate& had);
        float calcPt() const;
        
        inline float getPtHe3() const { return fHe3.fPtHe3; }
        inline float getPtHad() const { return fHad.fPtHad; }
        inline float getZVertex() const { return fColl.fZVertex; }
        inline float getCentralityFT0C() const { return fColl.fCentralityFT0C; }

    private:

        He3Candidate fHe3;
        HadCandidate fHad;
        CollisionCandidate fColl;
        bool fIsUnlikeSign = false;
        bool fIsHadSet = false;
        bool fIsHe3Set = false;
};

void Li4Candidate::setHe3(const He3Candidate& he3) 
{
    fHe3 = he3;

    if (fIsHadSet) {
        if ( (fHe3.fPtHe3 > 0 && fHad.fPtHad < 0) || (fHe3.fPtHe3 < 0 && fHad.fPtHad > 0) ) {
            fIsUnlikeSign = true;
        } else {
            fIsUnlikeSign = false;
        }
    }

    fIsHe3Set = true;
}

void Li4Candidate::setHad(const HadCandidate& had) 
{
    fHad = had;

    if (fIsHe3Set) {
        if ( (fHe3.fPtHe3 > 0 && fHad.fPtHad < 0) || (fHe3.fPtHe3 < 0 && fHad.fPtHad > 0) ) {
            fIsUnlikeSign = true;
        } else {
            fIsUnlikeSign = false;
        }
    }

    fIsHadSet = true;
}

void Li4Candidate::setBranch(TTree* tree)
{
    tree->Branch("fPtHe3", &fHe3.fPtHe3);
    tree->Branch("fEtaHe3", &fHe3.fEtaHe3);
    tree->Branch("fPhiHe3", &fHe3.fPhiHe3);
    tree->Branch("fPtHad", &fHad.fPtHad);
    tree->Branch("fEtaHad", &fHad.fEtaHad);
    tree->Branch("fPhiHad", &fHad.fPhiHad);
    tree->Branch("fDCAxyHe3", &fHe3.fDCAxyHe3);
    tree->Branch("fDCAzHe3", &fHe3.fDCAzHe3);
    tree->Branch("fDCAxyHad", &fHad.fDCAxyHad);
    tree->Branch("fDCAzHad", &fHad.fDCAzHad);
    tree->Branch("fSignalTPCHe3", &fHe3.fSignalTPCHe3);
    tree->Branch("fInnerParamTPCHe3", &fHe3.fInnerParamTPCHe3);
    tree->Branch("fSignalTPCHad", &fHad.fSignalTPCHad);
    tree->Branch("fInnerParamTPCHad", &fHad.fInnerParamTPCHad);
    tree->Branch("fMassTOFHe3", &fHe3.fMassTOFHe3);
    tree->Branch("fMassTOFHad", &fHad.fMassTOFHad);
    tree->Branch("fItsClusterSizeHe3", &fHe3.fItsClusterSizeHe3);
    tree->Branch("fItsClusterSizeHad", &fHad.fItsClusterSizeHad);
    tree->Branch("fPIDtrkHe3", &fHe3.fPIDtrkHe3);
    tree->Branch("fPIDtrkHad", &fHad.fPIDtrkHad);
    tree->Branch("fSharedClustersHe3", &fHe3.fSharedClustersHe3);
    tree->Branch("fSharedClustersHad", &fHad.fSharedClustersHad);
    tree->Branch("fNSigmaTPCHe3", &fHe3.fNSigmaTPCHe3);
    
    //tree->Branch("fNSigmaTPCHad", &fHad.fNSigmaTPCHad);
    tree->Branch("fNSigmaTPCHadPr", &fHad.fNSigmaTPCHad);
    tree->Branch("fNSigmaTOFHadPr", &fHad.fNSigmaTOFHad);

    tree->Branch("fChi2TPCHe3", &fHe3.fChi2TPCHe3);
    tree->Branch("fChi2TPCHad", &fHad.fChi2TPCHad);
    tree->Branch("fZVertex", &fColl.fZVertex);
    tree->Branch("fCentralityFT0C", &fColl.fCentralityFT0C);
    tree->Branch("fIs23", &fColl.fIs23);
}

float Li4Candidate::li4InvMass(const He3Candidate& he3, const HadCandidate& had)
{
    std::array<float, 3> pHe3 = {std::abs(he3.fPtHe3) * std::cos(he3.fPhiHe3), 
                                std::abs(he3.fPtHe3) * std::sin(he3.fPhiHe3),
                                std::abs(he3.fPtHe3) * std::sinh(he3.fEtaHe3)};
    std::array<float, 3> pHad = {std::abs(had.fPtHad) * std::cos(had.fPhiHad), 
                                std::abs(had.fPtHad) * std::sin(had.fPhiHad),
                                std::abs(had.fPtHad) * std::sinh(had.fEtaHad)};
    return physics::invariantMass(pHe3, pHad, physics::mass::kHelium3, physics::mass::kProton);
}

float Li4Candidate::calcInvMass() const
{
    return Li4Candidate::li4InvMass(fHe3, fHad);
}

float Li4Candidate::calcPt() const
{
    std::array<float, 3> pHe3 = {std::abs(fHe3.fPtHe3) * std::cos(fHe3.fPhiHe3), 
                                std::abs(fHe3.fPtHe3) * std::sin(fHe3.fPhiHe3),
                                std::abs(fHe3.fPtHe3) * std::sinh(fHe3.fEtaHe3)};
    std::array<float, 3> pHad = {std::abs(fHad.fPtHad) * std::cos(fHad.fPhiHad), 
                                std::abs(fHad.fPtHad) * std::sin(fHad.fPhiHad),
                                std::abs(fHad.fPtHad) * std::sinh(fHad.fEtaHad)};
    return physics::momentumMother(pHe3, pHad);
}
