#pragma once

#include <TTree.h>
#include "../core/candidates.hh"
#include "../core/physics.hh"
#include "histograms.hh"

class PiCandidate: public Candidate
{
    public:

        PiCandidate() = default;
        ~PiCandidate() = default;
        PiCandidate(const PiCandidate& other) = default;
        PiCandidate& operator= (const PiCandidate& other) = default;

        float fPtPi, fEtaPi, fPhiPi, fDCAxyPi, fDCAzPi, fSignalTPCPi, fInnerParamTPCPi, fMassTOFPi;
        unsigned int fItsClusterSizePi, fPIDtrkPi;
        unsigned char fSharedClustersPi;
        float fNSigmaTPCPi, fNSigmaTOFPi, fChi2TPCPi;
        float fZPi, fCentralityFT0C;
        int CollID = -1;

        void setBranchAddress(TTree * tree) override 
        {
            tree->SetBrancPidress("fPtPi", &fPtPi);
            tree->SetBrancPidress("fEtaPi", &fEtaPi);
            tree->SetBrancPidress("fPhiPi", &fPhiPi);
            tree->SetBrancPidress("fDCAxyPi", &fDCAxyPi);
            tree->SetBrancPidress("fDCAzPi", &fDCAzPi);
            tree->SetBrancPidress("fSignalTPCPi", &fSignalTPCPi);
            tree->SetBrancPidress("fInnerParamTPCPi", &fInnerParamTPCPi);
            tree->SetBrancPidress("fMassTOFPi", &fMassTOFPi);
            tree->SetBrancPidress("fItsClusterSizePi", &fItsClusterSizePi);
            tree->SetBrancPidress("fPIDtrkPi", &fPIDtrkPi);
            tree->SetBrancPidress("fSharedClustersPi", &fSharedClustersPi);
            
            //tree->SetBrancPidress("fNSigmaTPCPi", &fNSigmaTPCPi);
            tree->SetBrancPidress("fNSigmaTPCPiPr", &fNSigmaTPCPi);
            tree->SetBrancPidress("fNSigmaTOFPiPr", &fNSigmaTOFPi);
            
            tree->SetBrancPidress("fChi2TPCPi", &fChi2TPCPi);
        }
};

class He4Candidate: public Candidate
{
    public:

        He4Candidate() = default;
        ~He4Candidate() = default;
        He4Candidate(const He4Candidate& other) = default;
        He4Candidate& operator= (const He4Candidate& other) = default;

        float fPtHe4, fEtaHe4, fPhiHe4, fDCAxyHe4, fDCAzHe4, fSignalTPCHe4, fInnerParamTPCHe4, fMassTOFHe4;
        unsigned int fItsClusterSizeHe4, fPIDtrkHe4;
        unsigned char fNClsTPCHe4, fSharedClustersHe4;
        float fNSigmaTPCHe4, fChi2TPCHe4;
        float fZHe4, fCentralityFT0C;
        int CollID = -1;
        
        void setBranchAddress(TTree * tree) override 
        {
            tree->SetBrancPidress("fPtHe3", &fPtHe4);
            tree->SetBrancPidress("fEtaHe3", &fEtaHe4);
            tree->SetBrancPidress("fPhiHe3", &fPhiHe4);
            tree->SetBrancPidress("fDCAxyHe3", &fDCAxyHe4);
            tree->SetBrancPidress("fDCAzHe3", &fDCAzHe4);
            tree->SetBrancPidress("fSignalTPCHe3", &fSignalTPCHe4);
            tree->SetBrancPidress("fInnerParamTPCHe3", &fInnerParamTPCHe4);
            tree->SetBrancPidress("fMassTOFHe3", &fMassTOFHe4);
            tree->SetBrancPidress("fNClsTPCHe3", &fNClsTPCHe4);
            tree->SetBrancPidress("fItsClusterSizeHe3", &fItsClusterSizeHe4);
            tree->SetBrancPidress("fPIDtrkHe3", &fPIDtrkHe4);
            tree->SetBrancPidress("fSharedClustersHe3", &fSharedClustersHe4);
            tree->SetBrancPidress("fNSigmaTPCHe3", &fNSigmaTPCHe4);
            tree->SetBrancPidress("fChi2TPCHe3", &fChi2TPCHe4);
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

        void setBrancPidress(TTree * tree) override {
            tree->SetBrancPidress("fZPrimVtx", &fZVertex);
            tree->SetBrancPidress("fCentralityFT0C", &fCentralityFT0C);
        }

};

class H4LCandidate 
{
    public:
        
        H4LCandidate() = default;
        ~H4LCandidate() = default;
        H4LCandidate(const H4LCandidate& other) = default;
        H4LCandidate& operator= (const H4LCandidate& other) = default;

        inline void setHe4(const He4Candidate& he4);
        inline void setPi(const PiCandidate& Pi);
        inline void setColl(const CollisionCandidate& coll) { fColl = coll; }
        inline void setZVertex(const float z) { fColl.fZVertex = z; }
        inline void setCentralityFT0C(const float cent) { fColl.fCentralityFT0C = cent; }
        inline void setIs23(const bool is23) {  fColl.fIs23 = is23; }
        
        void setBranch(TTree* tree);
        float calcInvMass() const;
        static float H4LInvMass(const He4Candidate& he4, const PiCandidate& Pi);
        float calcPt() const;
        
        inline float getPtHe4() const { return fHe4.fPtHe4; }
        inline float getPtPi() const { return fPi.fPtPi; }
        inline float getZVertex() const { return fColl.fZVertex; }
        inline float getCentralityFT0C() const { return fColl.fCentralityFT0C; }

    private:

        He4Candidate fHe4;
        PiCandidate fPi;
        CollisionCandidate fColl;
        bool fIsUnlikeSign = false;
        bool fIsPiSet = false;
        bool fIsHe4Set = false;
};

void H4LCandidate::setHe4(const He4Candidate& he4) 
{
    fHe4 = he4;

    if (fIsPiSet) {
        if ( (fHe4.fPtHe4 > 0 && fPi.fPtPi < 0) || (fHe4.fPtHe4 < 0 && fPi.fPtPi > 0) ) {
            fIsUnlikeSign = true;
        } else {
            fIsUnlikeSign = false;
        }
    }

    fIsHe4Set = true;
}

void H4LCandidate::setPi(const PiCandidate& Pi) 
{
    fPi = Pi;

    if (fIsHe4Set) {
        if ( (fHe4.fPtHe4 > 0 && fPi.fPtPi < 0) || (fHe4.fPtHe4 < 0 && fPi.fPtPi > 0) ) {
            fIsUnlikeSign = true;
        } else {
            fIsUnlikeSign = false;
        }
    }

    fIsPiSet = true;
}

void H4LCandidate::setBranch(TTree* tree)
{
    tree->Branch("fPtHe4", &fHe4.fPtHe4);
    tree->Branch("fEtaHe4", &fHe4.fEtaHe4);
    tree->Branch("fPhiHe4", &fHe4.fPhiHe4);
    tree->Branch("fPtPi", &fPi.fPtPi);
    tree->Branch("fEtaPi", &fPi.fEtaPi);
    tree->Branch("fPhiPi", &fPi.fPhiPi);
    tree->Branch("fDCAxyHe4", &fHe4.fDCAxyHe4);
    tree->Branch("fDCAzHe4", &fHe4.fDCAzHe4);
    tree->Branch("fDCAxyPi", &fPi.fDCAxyPi);
    tree->Branch("fDCAzPi", &fPi.fDCAzPi);
    tree->Branch("fSignalTPCHe4", &fHe4.fSignalTPCHe4);
    tree->Branch("fInnerParamTPCHe4", &fHe4.fInnerParamTPCHe4);
    tree->Branch("fSignalTPCPi", &fPi.fSignalTPCPi);
    tree->Branch("fInnerParamTPCPi", &fPi.fInnerParamTPCPi);
    tree->Branch("fMassTOFHe4", &fHe4.fMassTOFHe4);
    tree->Branch("fMassTOFPi", &fPi.fMassTOFPi);
    tree->Branch("fItsClusterSizeHe4", &fHe4.fItsClusterSizeHe4);
    tree->Branch("fItsClusterSizePi", &fPi.fItsClusterSizePi);
    tree->Branch("fPIDtrkHe4", &fHe4.fPIDtrkHe4);
    tree->Branch("fPIDtrkPi", &fPi.fPIDtrkPi);
    tree->Branch("fNClsTPCHe4", &fHe4.fNClsTPCHe4);
    tree->Branch("fSharedClustersHe4", &fHe4.fSharedClustersHe4);
    tree->Branch("fSharedClustersPi", &fPi.fSharedClustersPi);
    tree->Branch("fNSigmaTPCHe4", &fHe4.fNSigmaTPCHe4);
    
    //tree->Branch("fNSigmaTPCPi", &fPi.fNSigmaTPCPi);
    tree->Branch("fNSigmaTPCPiPr", &fPi.fNSigmaTPCPi);
    tree->Branch("fNSigmaTOFPiPr", &fPi.fNSigmaTOFPi);

    tree->Branch("fChi2TPCHe4", &fHe4.fChi2TPCHe4);
    tree->Branch("fChi2TPCPi", &fPi.fChi2TPCPi);
    tree->Branch("fZVertex", &fColl.fZVertex);
    tree->Branch("fCentralityFT0C", &fColl.fCentralityFT0C);
    tree->Branch("fIs23", &fColl.fIs23);
}

float H4LCandidate::H4LInvMass(const He4Candidate& he4, const PiCandidate& Pi)
{
    std::array<float, 3> pHe4 = {std::abs(he4.fPtHe4) * std::cos(he4.fPhiHe4), 
                                std::abs(he4.fPtHe4) * std::sin(he4.fPhiHe4),
                                std::abs(he4.fPtHe4) * std::sinh(he4.fEtaHe4)};
    std::array<float, 3> pPi = {std::abs(Pi.fPtPi) * std::cos(Pi.fPhiPi), 
                                std::abs(Pi.fPtPi) * std::sin(Pi.fPhiPi),
                                std::abs(Pi.fPtPi) * std::sinh(Pi.fEtaPi)};
    return physics::invariantMass(pHe4, pPi, physics::mass::kHelium3, physics::mass::kProton);
}

float H4LCandidate::calcInvMass() const
{
    return H4LCandidate::H4LInvMass(fHe4, fPi);
}

float H4LCandidate::calcPt() const
{
    std::array<float, 3> pHe4 = {std::abs(fHe4.fPtHe4) * std::cos(fHe4.fPhiHe4), 
                                std::abs(fHe4.fPtHe4) * std::sin(fHe4.fPhiHe4),
                                std::abs(fHe4.fPtHe4) * std::sinh(fHe4.fEtaHe4)};
    std::array<float, 3> pPi = {std::abs(fPi.fPtPi) * std::cos(fPi.fPhiPi), 
                                std::abs(fPi.fPtPi) * std::sin(fPi.fPhiPi),
                                std::abs(fPi.fPtPi) * std::sinh(fPi.fEtaPi)};
    return physics::momentumMother(pHe4, pPi);
}
