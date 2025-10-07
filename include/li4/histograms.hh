#pragma once

#include <TH1F.h>
#include <TDirectory.h>
#include <TCanvas.h>

struct HistogramsQA
{
    TH1F* hHe3BeforeEMAll = new TH1F("hHe3BeforeEMAll", "; #it{p}_{T} (GeV/#it{c}); Entries", 200, -10, 0);
    TH1F* hHe3BeforeEM = new TH1F("hHe3BeforeEM", "; #it{p}_{T} (GeV/#it{c}); Entries", 200, -10, 0);
    TH1F* hHe3Unique = new TH1F("hHe3Unique", "; #it{p}_{T} (GeV/#it{c}); Entries", 200, -10, 0);
    TH1F* hHe3AfterEM = new TH1F("hHe3AfterEM", "; #it{p}_{T} (GeV/#it{c}); Entries", 200, -10, 0);
    TH1F* hInvMassBeforeEMUnlikeSign = new TH1F("hInvMassBeforeEMUnlikeSign", "; m (p+^{3}He) (GeV/#it{c}^{2}); Entries", 600, 3.743, 4.343);
    TH1F* hInvMassAfterEMUnlikeSign = new TH1F("hInvMassAfterEMUnlikeSign", "; m (p+^{3}He) (GeV/#it{c}^{2}); Entries", 600, 3.743, 4.343);
    TH1F* hInvMassBeforeEMLikeSign = new TH1F("hInvMassBeforeEMLikeSign", "; m (p+^{3}He) (GeV/#it{c}^{2}); Entries", 600, 3.743, 4.343);
    TH1F* hInvMassAfterEMLikeSign = new TH1F("hInvMassAfterEMLikeSign", "; m (p+^{3}He) (GeV/#it{c}^{2}); Entries", 600, 3.743, 4.343);

    ~HistogramsQA() {
        delete hHe3BeforeEMAll;
        delete hHe3BeforeEM;
        delete hHe3Unique;
        delete hHe3AfterEM;

        delete hInvMassBeforeEMUnlikeSign;
        delete hInvMassAfterEMUnlikeSign;
        delete hInvMassBeforeEMLikeSign;
        delete hInvMassAfterEMLikeSign;
    }

    void saveHistograms(TDirectory* output)
    {
        output->cd();

        hHe3BeforeEMAll->Write();
        hHe3BeforeEM->Write();
        hHe3Unique->Write();
        hHe3AfterEM->Write();

        hInvMassBeforeEMUnlikeSign->Write();
        hInvMassAfterEMUnlikeSign->Write();
        hInvMassBeforeEMLikeSign->Write();
        hInvMassAfterEMLikeSign->Write();

        // plot in the same canvas
        hHe3BeforeEMAll->SetLineColor(kBlack);
        hHe3BeforeEM->SetLineColor(kRed);
        hHe3AfterEM->SetLineColor(kBlue);
        TCanvas* cHe3PtComparison = new TCanvas("c_he3PtComparison", "", 800, 600);
        hHe3BeforeEMAll->DrawNormalized();
        hHe3BeforeEM->DrawNormalized("SAME");
        hHe3AfterEM->DrawNormalized("SAME");
        cHe3PtComparison->BuildLegend();
        cHe3PtComparison->Write();

        hInvMassBeforeEMLikeSign->SetLineColor(kRed);
        hInvMassAfterEMLikeSign->SetLineColor(kBlue);
        TCanvas* cInvariantMassComparisonLikeSign = new TCanvas("c_invariantMassComparisonLikeSign", "", 800, 600);
        hInvMassBeforeEMLikeSign->DrawNormalized();
        hInvMassAfterEMLikeSign->DrawNormalized("SAME");
        cInvariantMassComparisonLikeSign->BuildLegend();
        cInvariantMassComparisonLikeSign->Write();

        hInvMassBeforeEMUnlikeSign->SetLineColor(kRed);
        hInvMassAfterEMUnlikeSign->SetLineColor(kBlue);
        TCanvas* cInvariantMassComparisonUnlikeSign = new TCanvas("c_invariantMassComparisonUnlikeSign", "", 800, 600);
        hInvMassBeforeEMUnlikeSign->DrawNormalized();
        hInvMassAfterEMUnlikeSign->DrawNormalized("SAME");
        cInvariantMassComparisonUnlikeSign->BuildLegend();
        cInvariantMassComparisonUnlikeSign->Write();

        delete cHe3PtComparison;
        delete cInvariantMassComparisonLikeSign;
        delete cInvariantMassComparisonUnlikeSign;
    }
};
