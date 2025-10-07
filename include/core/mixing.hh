#include 

void fillParticlesFromTree(TTree* tree, std::vector<Particle>& pions, std::vector<Particle>& protons) {
    
    std::cout << "Filling particles from tree with " << tree->GetEntries() << " entries." << std::endl; 

    Particle particle;
    particle.setBranchAddress(tree);

    const int nEntries = tree->GetEntries();
    for (int i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (particle.partId == 2) { // pion
            pions.push_back(particle);
        } else if (particle.partId == 4) { // proton
            protons.push_back(particle);
        }
    }
}

void perforMixingRotation(const std::vector<Particle>& protons, const std::vector<Particle>& pions, 
                    TH2F* h2PInvariantMass, const int mixingDepth = 2) {

    std::cout << "Starting mixing with " << protons.size() << " protons and " << pions.size() << " pions." << std::endl;

    const int nPions = pions.size();
    for (int iproton = 0; iproton < protons.size(); ++iproton) {

        if (iproton % (protons.size()/100) == 0) {
            std::cout << "Processing proton " << iproton << " / " << protons.size() << " (" 
                      << iproton/protons.size()*100 << "%)" << std::endl;
        }
        
        const Particle& proton = protons[iproton];
        std::array<float, 3> p1 = {std::abs(proton.p), proton.eta, proton.phi};
        for (int i = 0; i < mixingDepth; ++i) {
            Particle pion = pions[gRandom->Integer(nPions)];
            while (pion.p * proton.p > 0) { // ensure unlike-sign
                pion = pions[gRandom->Integer(nPions)];
            }
            const float rotatedPhiPion = physics::randomAngleRotation(pion.phi);
            std::array<float, 3> p2 = {std::abs(pion.p), pion.eta, rotatedPhiPion};
            const float invMass = physics::invariantMass(p1, p2, physics::kMassProton, physics::kMassPion);
            const float p = physics::momentumMother(p1, p2);
            const float charge = proton.p > 0 ? 1 : -1;
            h2PInvariantMass->Fill(charge * p, invMass);
        }
    }
}

void perforMixingLikeSign(const std::vector<Particle>& protons, const std::vector<Particle>& pions, 
                    TH2F* h2PInvariantMass, const int mixingDepth = 2) {

    std::cout << "Starting mixing with " << protons.size() << " protons and " << pions.size() << " pions." << std::endl;

    const int nPions = pions.size();
    for (int iproton = 0; iproton < protons.size(); ++iproton) {

        if (iproton % (protons.size()/100) == 0) {
            std::cout << "Processing proton " << iproton << " / " << protons.size() << " (" 
                      << iproton/protons.size()*100 << "%)" << std::endl;
        }
        
        const Particle& proton = protons[iproton];
        std::array<float, 3> p1 = {std::abs(proton.p), proton.eta, proton.phi};
        for (int i = 0; i < mixingDepth; ++i) {
            Particle pion = pions[gRandom->Integer(nPions)];
            while (pion.p * proton.p < 0) { // ensure like-sign
                pion = pions[gRandom->Integer(nPions)];
            }
            std::array<float, 3> p2 = {std::abs(pion.p), pion.eta, pion.phi};
            const float invMass = physics::invariantMass(p1, p2, physics::kMassProton, physics::kMassPion);
            const float p = physics::momentumMother(p1, p2);
            const float charge = proton.p > 0 ? 1 : -1;
            h2PInvariantMass->Fill(charge * p, invMass);
        }
    }
}

void run_lambda_mixing(const bool mergeTrees = false) {

    const char * mergedInputFileName = "merged_trees/data_04_08_2025_merged.root";
    const char * treeName = "O2clsttable";
    
    if (mergeTrees) {
        const char * unmergedInputFileName = "/data/galucia/its_pid/LHC24_pass1_skimmed/data_04_08_2025.root";
        TFile * mergedTreeFile = TFile::Open(mergedInputFileName, "RECREATE");
        treeHandling::treeMerging(unmergedInputFileName, treeName, mergedTreeFile);
        mergedTreeFile->Close();
    }

    TFile * inputFile = TFile::Open(mergedInputFileName, "READ");
    TTree * tree = (TTree *)inputFile->Get(treeName);
    std::vector<Particle> pions, protons;

    fillParticlesFromTree(tree, pions, protons);

    auto h2PInvariantMassRotation = new TH2F("h2PInvariantMassRotation", "Invariant Mass Distribution - Rotation strategy; #it{p}_{p#pi} (GeV/#it{c}); M_{p#pi} (GeV/c^{2});", 100, -5, 5, 50, 1.08, 1.18);
    perforMixingRotation(protons, pions, h2PInvariantMassRotation, 2);

    auto h2PInvariantMassLikeSign = new TH2F("h2PInvariantMassLikeSign", "Invariant Mass Distribution - Like-sign strategy; #it{p}_{p#pi} (GeV/#it{c}); M_{p#pi} (GeV/c^{2});", 100, -5, 5, 50, 1.08, 1.18);
    perforMixingLikeSign(protons, pions, h2PInvariantMassLikeSign, 2);

    const char * outputFileName = "output/v0_cascade_mixing.root";
    TFile * outputFile = TFile::Open(outputFileName, "RECREATE");
    outputFile->cd();
    h2PInvariantMassRotation->Write();
    h2PInvariantMassLikeSign->Write();
    outputFile->Close();

}