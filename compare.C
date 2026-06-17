#include <algorithm>
#include <iostream>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"

namespace {

TH1* GetHist(TFile* file, const char* name)
{
    if(file == nullptr) return nullptr;

    TH1* hist = dynamic_cast<TH1*>(file->Get(name));
    if(hist == nullptr) {
        std::cerr << "Cannot find histogram: " << name << std::endl;
    }
    return hist;
}

void DrawOverlay(TH1* hSim, TH1* hData,
                 const char* title,
                 const char* simLabel,
                 const char* dataLabel)
{
    if(hSim == nullptr || hData == nullptr) return;

    hSim = dynamic_cast<TH1*>(hSim->Clone(Form("%s_norm", hSim->GetName())));
    hData = dynamic_cast<TH1*>(hData->Clone(Form("%s_norm", hData->GetName())));
    if(hSim == nullptr || hData == nullptr) return;

    hSim->SetDirectory(nullptr);
    hData->SetDirectory(nullptr);

    if(hSim->Integral() > 0.) hSim->Scale(1. / hSim->Integral());
    if(hData->Integral() > 0.) hData->Scale(1. / hData->Integral());

    hSim->SetLineColor(kBlue + 1);
    hSim->SetLineWidth(2);
    hSim->SetStats(0);
    hSim->SetTitle(title);

    hData->SetLineColor(kRed + 1);
    hData->SetLineWidth(2);
    hData->SetStats(0);

    const Double_t maxY = std::max(hSim->GetMaximum(), hData->GetMaximum());
    hSim->SetMaximum(maxY * 1.15);

    hSim->Draw("hist");
    hData->Draw("hist same");

    TLegend* legend = new TLegend(0.56, 0.72, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(hSim, simLabel, "l");
    legend->AddEntry(hData, dataLabel, "l");
    legend->Draw();
}

} // namespace

void compare()
{
    TFile* file_data = TFile::Open("output/out_24017_30_60.root");
    if(file_data == nullptr || file_data->IsZombie()) {
        std::cerr << "Cannot open data file." << std::endl;
        return;
    }

    TFile* file_sim = TFile::Open("output/GEMdig_output_test.root");
    if(file_sim == nullptr || file_sim->IsZombie()) {
        std::cerr << "Cannot open sim file." << std::endl;
        return;
    }

    TH1* hDataChargeY = GetHist(file_data, "h_matched_gem_y_charge_hycal_det3");
    TH1* hDataPeakY = GetHist(file_data, "h_matched_gem_y_peak_hycal_det3");
    TH1* hDataSizeY = GetHist(file_data, "h_matched_gem_y_size_hycal_det3");

    TH1* hSimChargeY = GetHist(file_sim, "chargeY_signal");
    TH1* hSimPeakY = GetHist(file_sim, "peakY_signal");
    TH1* hSimSizeY = GetHist(file_sim, "clusterSizeY_signal");

    if(hDataChargeY == nullptr || hDataPeakY == nullptr ||
       hDataSizeY == nullptr || hSimChargeY == nullptr ||
       hSimPeakY == nullptr || hSimSizeY == nullptr) {
        return;
    }

    TCanvas* cYCompare = new TCanvas("cYCompare",
                                     "Data/Sim Y distributions", 1200, 450);
    cYCompare->Divide(3, 1);

    cYCompare->cd(1);
    DrawOverlay(hSimChargeY, hDataChargeY,
                "Y total charge;total charge;Entries",
                "sim", "data");

    cYCompare->cd(2);
    DrawOverlay(hSimPeakY, hDataPeakY,
                "Y peak;peak ADC;Entries",
                "sim", "data");

    cYCompare->cd(3);
    DrawOverlay(hSimSizeY, hDataSizeY,
                "Y cluster size;cluster size;Entries",
                "sim", "data");

    cYCompare->Update();
}
