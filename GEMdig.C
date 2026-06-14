#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"

#include "GEMdigVariables.h"
#include "GEMIonModel.h"
#include "GEMStripMap.h"

using std::cout;
using std::endl;
using std::string;

std::vector<string> FindRootFiles(const string& folder,
                                  const string& filenameKeyword)
{
    std::vector<string> files;
    std::error_code error;
    const std::filesystem::path folderPath(folder);

    if(!std::filesystem::is_directory(folderPath, error)) {
        std::cerr << "ROOT file folder not found: " << folder << endl;
        return files;
    }

    std::filesystem::recursive_directory_iterator iterator(folderPath, error);
    const std::filesystem::recursive_directory_iterator end;
    while(iterator != end) {
        if(error) {
            error.clear();
            iterator.increment(error);
            continue;
        }

        const std::filesystem::directory_entry& entry = *iterator;
        const string filename = entry.path().filename().string();
        if(entry.is_regular_file(error) &&
           entry.path().extension() == ".root" &&
           filename.find(filenameKeyword) != string::npos) {
            files.push_back(entry.path().string());
        }
        iterator.increment(error);
    }

    std::sort(files.begin(), files.end());
    return files;
}

Double_t PulseShape(Double_t t,
                    Double_t C,  // normalization factor, charge number on one strip
				    Double_t Tp) // shaping time 56ns
{   
    if(t < 0.) return 0.;
    else{
        Double_t v;
        Double_t x;
        x = t/Tp;
        v = C/Tp * x * TMath::Exp(-x);
  
        return ( v>0. ) ? v : 0.;
    }
}

//add pedestal noise
void GetPedestal(std::array<Double_t, samplingNum>& pulseXorY){
    for(int n = 0; n < samplingNum; n++){
        pulseXorY[n] += gRandom->Gaus(0., 15.); // pedestal noise, 15 ADC
    }
}

void ADCConvert(Double_t off, Double_t gain, Int_t bits,
                std::array<Double_t, samplingNum>& pulseXorY)
{
  // Convert analog value 'val' to integer ADC reading with 'bits' resolution
    for(int n=0; n<samplingNum; n++){
        Double_t val = pulseXorY[n];
        if( val < 0.) val = 0.;
        Double_t vvv = (val - off)/gain;
        Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
        //saturation = 1400.;
        if( vvv > saturation ) vvv = saturation;
        Double_t dval =
            static_cast<Double_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));
        if( dval < 0 ) dval = 0;
        pulseXorY[n] = dval;
    }
}

void SetTH2D(TH2D* hist) {
    hist->GetXaxis()->SetTitle("sample (25ns)");
    hist->GetYaxis()->SetTitle("strip number");
    hist->GetZaxis()->SetTitle("ADC");
    //hist->GetZaxis()->SetRangeUser(0., 2000.);
    hist->GetXaxis()->SetTitleOffset(1.7);
    hist->GetYaxis()->SetTitleOffset(1.7);
    hist->GetZaxis()->SetTitleOffset(1.5);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
}

void loadFile(TChain& file, string filename) {
    file.Add(filename.c_str());
    file.SetBranchAddress("EventID", &EventID);
    file.SetBranchAddress("GEM.N", &GasN);
    file.SetBranchAddress("GEM.X", &GasX[0]);
    file.SetBranchAddress("GEM.Y", &GasY[0]);
    file.SetBranchAddress("GEM.Z", &GasZ[0]);
    file.SetBranchAddress("GEM.Time", &GasTime[0]);
    file.SetBranchAddress("GEM.Xout", &GasXout[0]);
    file.SetBranchAddress("GEM.Yout", &GasYout[0]);
    file.SetBranchAddress("GEM.Zout", &GasZout[0]);
    file.SetBranchAddress("GEM.PID", &GasPID[0]);
    file.SetBranchAddress("GEM.PTID", &GasPTID[0]);
    file.SetBranchAddress("GEM.TID", &GasTID[0]);
    file.SetBranchAddress("GEM.DID", &GasDID[0]);
    file.SetBranchAddress("GEM.Theta", &GasTheta[0]);
    file.SetBranchAddress("GEM.P", &GasP[0]);
    file.SetBranchAddress("GEM.Edep", &GasEdep[0]);
}

void GEMdig() {
    auto start = std::chrono::high_resolution_clock::now();
    gRandom->SetSeed(0);//1232 // 0 表示用系统时间自动生成种子


    std::vector<string> signalFiles =
        FindRootFiles("data", "signal");
    std::vector<string> backgroundFiles =
        FindRootFiles("2105", "simrun");
    TChain signal("T");
    TChain background("T");

    for(const string& filename : signalFiles) {
        loadFile(signal, filename);
    }
    for(const string& filename : backgroundFiles) {
        loadFile(background, filename);
    }
    
    const int Ndata = signal.GetEntries();
    const int Nbg = background.GetEntries();
    cout << "Number of entries in all signal data: " << Ndata << endl;
    cout << "Number of entries in all background data: " << Nbg << endl;
    background.GetEntry(Nbg-1);
    double fluxN_bg = EventID * backgroundFiles.size() * 10.;
    cout << "background flux (Hz): " << fluxN_bg << endl;
    int bg_entry = int(double(Nbg) / fluxN_bg * 312.5/4. * 350.);
    cout << "background entry numbers for each event: " << bg_entry << endl;
    bg_entry = 2.;
    TH1D *peak = new TH1D("peak", "peak of the pulse", 100, 0, 5000);
    TH1D *peakY = new TH1D("peakY", "peakY of the pulse", 100, 0, 5000);
    TH2D *adc_Num = new TH2D("adc_Num", "Pulse Amplitude vs Fired Strip Number", 40, 0, 4000, 11, 0-0.5, 11-0.5);

    TH1D *Occupancy_X[4];
    TH1D *Occupancy_Y[4];
    for(int did = 0; did < 4; did++){
        Occupancy_X[did] = new TH1D(Form("Occupancy of chamber%d(X strips)", did), Form("Occupancy of chamber%d (X strips)", did), map_x_bins, 0, map_x_bins);
        Occupancy_Y[did] = new TH1D(Form("Occupancy of chamber%d(Y strips)", did), Form("Occupancy of chamber%d (Y strips)", did), map_y_bins, 0, map_y_bins);
    }

    TFile *FiredStripsF = new TFile("output/FiredStrips_test.root", "RECREATE");
    FiredStripsF->cd();
    TTree *FiredStripsT = new TTree("T", "T");

    Int_t outputEventNum = -1;
    Int_t fireNumX[4] = {};
    Int_t fireNumY[4] = {};

    Int_t stripNumberX_0[map_x_bins];
    Double_t pulseSamplesX_0[map_x_bins][samplingNum];
    Double_t peakX_0[map_x_bins];
    Double_t totalChargeX_0[map_x_bins];
    Int_t peakSampleX_0[map_x_bins];

    Int_t stripNumberY_0[map_y_bins];
    Double_t pulseSamplesY_0[map_y_bins][samplingNum];
    Double_t peakY_0[map_y_bins];
    Double_t totalChargeY_0[map_y_bins];
    Int_t peakSampleY_0[map_y_bins];

    Int_t stripNumberX_1[map_x_bins];
    Double_t pulseSamplesX_1[map_x_bins][samplingNum];
    Double_t peakX_1[map_x_bins];
    Double_t totalChargeX_1[map_x_bins];
    Int_t peakSampleX_1[map_x_bins];

    Int_t stripNumberY_1[map_y_bins];
    Double_t pulseSamplesY_1[map_y_bins][samplingNum];
    Double_t peakY_1[map_y_bins];
    Double_t totalChargeY_1[map_y_bins];
    Int_t peakSampleY_1[map_y_bins];

    Int_t stripNumberX_2[map_x_bins];
    Double_t pulseSamplesX_2[map_x_bins][samplingNum];
    Double_t peakX_2[map_x_bins];
    Double_t totalChargeX_2[map_x_bins];
    Int_t peakSampleX_2[map_x_bins];

    Int_t stripNumberY_2[map_y_bins];
    Double_t pulseSamplesY_2[map_y_bins][samplingNum];
    Double_t peakY_2[map_y_bins];
    Double_t totalChargeY_2[map_y_bins];
    Int_t peakSampleY_2[map_y_bins];

    Int_t stripNumberX_3[map_x_bins];
    Double_t pulseSamplesX_3[map_x_bins][samplingNum];
    Double_t peakX_3[map_x_bins];
    Double_t totalChargeX_3[map_x_bins];
    Int_t peakSampleX_3[map_x_bins];

    Int_t stripNumberY_3[map_y_bins];
    Double_t pulseSamplesY_3[map_y_bins][samplingNum];
    Double_t peakY_3[map_y_bins];
    Double_t totalChargeY_3[map_y_bins];
    Int_t peakSampleY_3[map_y_bins];

    FiredStripsT->Branch("eventNum", &outputEventNum, "eventNum/I");

    FiredStripsT->Branch("fireNumX_0", &fireNumX[0], "fireNumX_0/I");
    FiredStripsT->Branch("fireNumY_0", &fireNumY[0], "fireNumY_0/I");
    FiredStripsT->Branch("stripNumberX_0", stripNumberX_0, "stripNumberX_0[fireNumX_0]/I");
    FiredStripsT->Branch("pulseSamplesX_0", pulseSamplesX_0, "pulseSamplesX_0[fireNumX_0][6]/D");
    FiredStripsT->Branch("peakX_0", peakX_0, "peakX_0[fireNumX_0]/D");
    FiredStripsT->Branch("totalChargeX_0", totalChargeX_0, "totalChargeX_0[fireNumX_0]/D");
    FiredStripsT->Branch("peakSampleX_0", peakSampleX_0, "peakSampleX_0[fireNumX_0]/I");
    FiredStripsT->Branch("stripNumberY_0", stripNumberY_0, "stripNumberY_0[fireNumY_0]/I");
    FiredStripsT->Branch("pulseSamplesY_0", pulseSamplesY_0, "pulseSamplesY_0[fireNumY_0][6]/D");
    FiredStripsT->Branch("peakY_0", peakY_0, "peakY_0[fireNumY_0]/D");
    FiredStripsT->Branch("totalChargeY_0", totalChargeY_0, "totalChargeY_0[fireNumY_0]/D");
    FiredStripsT->Branch("peakSampleY_0", peakSampleY_0, "peakSampleY_0[fireNumY_0]/I");

    FiredStripsT->Branch("fireNumX_1", &fireNumX[1], "fireNumX_1/I");
    FiredStripsT->Branch("fireNumY_1", &fireNumY[1], "fireNumY_1/I");
    FiredStripsT->Branch("stripNumberX_1", stripNumberX_1, "stripNumberX_1[fireNumX_1]/I");
    FiredStripsT->Branch("pulseSamplesX_1", pulseSamplesX_1, "pulseSamplesX_1[fireNumX_1][6]/D");
    FiredStripsT->Branch("peakX_1", peakX_1, "peakX_1[fireNumX_1]/D");
    FiredStripsT->Branch("totalChargeX_1", totalChargeX_1, "totalChargeX_1[fireNumX_1]/D");
    FiredStripsT->Branch("peakSampleX_1", peakSampleX_1, "peakSampleX_1[fireNumX_1]/I");
    FiredStripsT->Branch("stripNumberY_1", stripNumberY_1, "stripNumberY_1[fireNumY_1]/I");
    FiredStripsT->Branch("pulseSamplesY_1", pulseSamplesY_1, "pulseSamplesY_1[fireNumY_1][6]/D");
    FiredStripsT->Branch("peakY_1", peakY_1, "peakY_1[fireNumY_1]/D");
    FiredStripsT->Branch("totalChargeY_1", totalChargeY_1, "totalChargeY_1[fireNumY_1]/D");
    FiredStripsT->Branch("peakSampleY_1", peakSampleY_1, "peakSampleY_1[fireNumY_1]/I");

    FiredStripsT->Branch("fireNumX_2", &fireNumX[2], "fireNumX_2/I");
    FiredStripsT->Branch("fireNumY_2", &fireNumY[2], "fireNumY_2/I");
    FiredStripsT->Branch("stripNumberX_2", stripNumberX_2, "stripNumberX_2[fireNumX_2]/I");
    FiredStripsT->Branch("pulseSamplesX_2", pulseSamplesX_2, "pulseSamplesX_2[fireNumX_2][6]/D");
    FiredStripsT->Branch("peakX_2", peakX_2, "peakX_2[fireNumX_2]/D");
    FiredStripsT->Branch("totalChargeX_2", totalChargeX_2, "totalChargeX_2[fireNumX_2]/D");
    FiredStripsT->Branch("peakSampleX_2", peakSampleX_2, "peakSampleX_2[fireNumX_2]/I");
    FiredStripsT->Branch("stripNumberY_2", stripNumberY_2, "stripNumberY_2[fireNumY_2]/I");
    FiredStripsT->Branch("pulseSamplesY_2", pulseSamplesY_2, "pulseSamplesY_2[fireNumY_2][6]/D");
    FiredStripsT->Branch("peakY_2", peakY_2, "peakY_2[fireNumY_2]/D");
    FiredStripsT->Branch("totalChargeY_2", totalChargeY_2, "totalChargeY_2[fireNumY_2]/D");
    FiredStripsT->Branch("peakSampleY_2", peakSampleY_2, "peakSampleY_2[fireNumY_2]/I");

    FiredStripsT->Branch("fireNumX_3", &fireNumX[3], "fireNumX_3/I");
    FiredStripsT->Branch("fireNumY_3", &fireNumY[3], "fireNumY_3/I");
    FiredStripsT->Branch("stripNumberX_3", stripNumberX_3, "stripNumberX_3[fireNumX_3]/I");
    FiredStripsT->Branch("pulseSamplesX_3", pulseSamplesX_3, "pulseSamplesX_3[fireNumX_3][6]/D");
    FiredStripsT->Branch("peakX_3", peakX_3, "peakX_3[fireNumX_3]/D");
    FiredStripsT->Branch("totalChargeX_3", totalChargeX_3, "totalChargeX_3[fireNumX_3]/D");
    FiredStripsT->Branch("peakSampleX_3", peakSampleX_3, "peakSampleX_3[fireNumX_3]/I");
    FiredStripsT->Branch("stripNumberY_3", stripNumberY_3, "stripNumberY_3[fireNumY_3]/I");
    FiredStripsT->Branch("pulseSamplesY_3", pulseSamplesY_3, "pulseSamplesY_3[fireNumY_3][6]/D");
    FiredStripsT->Branch("peakY_3", peakY_3, "peakY_3[fireNumY_3]/D");
    FiredStripsT->Branch("totalChargeY_3", totalChargeY_3, "totalChargeY_3[fireNumY_3]/D");
    FiredStripsT->Branch("peakSampleY_3", peakSampleY_3, "peakSampleY_3[fireNumY_3]/I");

    GEMStripMap gemMap[4]; // strip map for the four chambers, indexed by detector ID
    for(int did = 0; did < 4; did++) gemMap[did] = GEMStripMap(did);

    double validEntries = 0;
    for(int i=0;i<1000;i++) { //max 9000 for 4um
        if (i%100 == 0) cout << "start " << i << endl;        

        //clear strip ADC information for each event
        for(int did = 0; did < 4; did++){
            gemMap[did].Clear();
        }

        /*
        signal.GetEntry(i);
        
        int sig = 0;
        for(int n = 0; n < GasN; n++){
            //for the one beam particle
            if( GasPTID[n] != 0 || GasZout[n] - GasZ[n] < 2.999) continue;
            int did = GasDID[n];
            double t_hit = GasTime[n];
            double hitX[2], hitY[2], hitZ[2];
            if(did == 0 || did == 2){ 
                hitX[0] = GasX[n]; hitY[0] = GasY[n]; hitZ[0] = GasZ[n];
                hitX[1] = GasXout[n]; hitY[1] = GasYout[n]; hitZ[1] = GasZout[n];
            }
            if(did == 1 || did == 3){ 
                hitX[0] = -GasX[n]; hitY[0] = GasY[n]; hitZ[0] = GasZ[n];
                hitX[1] = -GasXout[n]; hitY[1] = GasYout[n]; hitZ[1] = GasZout[n];
            }

            TVector3 in(hitX[0], hitY[0], hitZ[0]);
            TVector3 out(hitX[1], hitY[1], hitZ[1]);
            IonModel(in, out, GasEdep[n]);
            if(fRIon.empty()) continue;

            //Simulate the DAQ and trigger jitter
            Double_t trigger_jitter = gRandom->Gaus(0, fTriggerJitter);
            Double_t apvJitter = (gRandom->Uniform(fAPVTimeJitter) - fAPVTimeJitter/2.);
            trigger_jitter += apvJitter;
            trigger_jitter=0.;

            // to get the region of interest for the ion, only sum the charge on these strips to save time
            //get the mean position of the ion distribution on the readout plane as the hit position
            double hitX_mean = 0, hitY_mean = 0;
            for(auto& ion : fRIon) {
                hitX_mean += ion.Xp;
                hitY_mean += ion.Yp;
            }
            hitX_mean /= fRIon.size();
            hitY_mean /= fRIon.size();
            GEMStripRegion region = findRegion(gemMap[did], hitX_mean, hitY_mean, activeStripNhalf);
            if(!region.IsValid()) continue;

            //Simulate the avalanche and get the total charge on each strip
            for(auto& ion : fRIon) {
                //calculate the charge distribution on the strip plane
                int stripStep = 6;
                int xStripN = region.upX_strip + 1 - region.lowX_strip;
                int yStripN = region.upY_strip + 1 - region.lowY_strip;
                int xBinN = xStripN * stripStep;
                int yBinN = yStripN * stripStep;
                Double_t area = 0.4*0.4/stripStep/stripStep;
                TH2D *integrateMap = new TH2D("inegrate map", "integrate map",
                    xBinN, region.xLow, region.xHigh,
                    yBinN, region.yLow, region.yHigh);

                Double_t eff_sigma = ion.eff_sigma;
                for(int ix = 1; ix <= xBinN; ix++){
                    Double_t dx = ion.Xp - integrateMap->GetXaxis()->GetBinCenter(ix);
                    Double_t dx2 = dx * dx;
                    for(int iy = 1; iy <= yBinN; iy++){
                        Double_t dy = ion.Yp - integrateMap->GetYaxis()->GetBinCenter(iy);
                        Double_t dy2 = dy * dy;
                        Double_t result = fAvaGain * ion.ggnorm * (1./(TMath::Pi()*eff_sigma)) * (eff_sigma*eff_sigma) / ((dx2+dy2) + eff_sigma*eff_sigma);
                        integrateMap->SetBinContent(ix, iy, result); //the charge deposit in each small area, charge / mm^2
                    }
                }

                //Get integrated charge for each blocks then rebin to the strip width
                integrateMap->RebinX(stripStep);
                integrateMap->RebinY(stripStep);

                for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                    GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                    Double_t Charge = 0.;
                    for(int yy = 1; yy <= yBinN/stripStep; yy++){
                        Charge += integrateMap->GetBinContent(kx - region.lowX_strip + 1, yy);
                    }
                    Charge *= area;
                    xStrip->charge += Charge;
                }
                for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                    GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                    Double_t Charge = 0.;
                    for(int xx = 1; xx <= xBinN/stripStep; xx++){
                        Charge += integrateMap->GetBinContent(xx, ky - region.lowY_strip + 1);
                    }
                    Charge *= area;
                    yStrip->charge += Charge;
                }
                delete integrateMap;
            }
            //Simulate the pulse shape and get the ADC value for each strip this particle hits
            double startTime = t_hit + minTravelT + trigger_jitter - 5.; // -5ns (trigger lattency) to move the pulse to the left 
            for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                for(int t = 0; t < samplingNum; t++){
                    double time = double(t+0.5) * samplingTime - startTime;
                    xStrip->pulse[t] += PulseShape(time, xStrip->charge, pulseTimeLength);
                }
                xStrip->charge = 0; // clear charge after getting the pulse
            }
            for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                for(int t = 0; t < samplingNum; t++){
                    double time = double(t+0.5) * samplingTime - startTime;
                    yStrip->pulse[t] += PulseShape(time, yStrip->charge, pulseTimeLength);
                }
                yStrip->charge = 0; // clear charge after getting the pulse
            }
            sig++;
        }
        if(sig < 1) continue;
        */
        validEntries++;
        //for the background particles
        int bg = 0;
        for(int ii=i*bg_entry; ii<(i+1)*bg_entry; ii++) {
            background.GetEntry(ii);
            Double_t t0 = gRandom->Uniform(-200., 6.*25.);

            for(int n = 0; n < GasN; n++){
                int did = GasDID[n];
                double t_hit = GasTime[n] + t0;
                double hitX[2], hitY[2], hitZ[2];
                if(did == 0 || did == 2){ 
                    hitX[0] = GasX[n]; hitY[0] = GasY[n]; hitZ[0] = GasZ[n];
                    hitX[1] = GasXout[n]; hitY[1] = GasYout[n]; hitZ[1] = GasZout[n];
                }
                if(did == 1 || did == 3){ 
                    hitX[0] = -GasX[n]; hitY[0] = GasY[n]; hitZ[0] = GasZ[n];
                    hitX[1] = -GasXout[n]; hitY[1] = GasYout[n]; hitZ[1] = GasZout[n];
                }

                TVector3 in(hitX[0], hitY[0], hitZ[0]);
                TVector3 out(hitX[1], hitY[1], hitZ[1]);
                IonModel(in, out, GasEdep[n]);
                if(fRIon.empty()) continue;

                //Simulate the DAQ and trigger jitter
                Double_t trigger_jitter = gRandom->Gaus(0, fTriggerJitter);
                Double_t apvJitter = (gRandom->Uniform(fAPVTimeJitter) - fAPVTimeJitter/2.);
                trigger_jitter += apvJitter;
                trigger_jitter=0.;

                // to get the region of interest for the ion, only sum the charge on these strips to save time
                GEMStripRegion region = findRegion(gemMap[did], fRIon[0].Xp, fRIon[0].Yp, activeStripNhalf);
                if(!region.IsValid()) continue;

                //Simulate the avalanche and get the total charge on each strip
                for(auto& ion : fRIon) {
                    //calculate the charge distribution on the strip plane
                    int stripStep = 6;
                    int xStripN = region.upX_strip + 1 - region.lowX_strip;
                    int yStripN = region.upY_strip + 1 - region.lowY_strip;
                    int xBinN = xStripN * stripStep;
                    int yBinN = yStripN * stripStep;
                    Double_t area = 0.4*0.4/stripStep/stripStep;
                    TH2D *integrateMap = new TH2D("inegrate map", "integrate map",
                        xBinN, region.xLow, region.xHigh,
                        yBinN, region.yLow, region.yHigh);
                    
                    Double_t eff_sigma = ion.eff_sigma;
                    for(int ix = 1; ix <= xBinN; ix++){
                        Double_t dx = ion.Xp - integrateMap->GetXaxis()->GetBinCenter(ix);
                        Double_t dx2 = dx * dx;
                        for(int iy = 1; iy <= yBinN; iy++){
                            Double_t dy = ion.Yp - integrateMap->GetYaxis()->GetBinCenter(iy);
                            Double_t dy2 = dy * dy;
                            Double_t result = fAvaGain * ion.ggnorm * (1./(TMath::Pi()*eff_sigma)) * (eff_sigma*eff_sigma) / ((dx2+dy2) + eff_sigma*eff_sigma);
                            integrateMap->SetBinContent(ix, iy, result); //the charge deposit in each small area, charge / mm^2
                        }
                    }
                    //Get integrated charge for each blocks then rebin to the strip width
                    integrateMap->RebinX(stripStep);
                    integrateMap->RebinY(stripStep);

                    for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                        GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                        Double_t Charge = 0.;
                        for(int yy = 1; yy <= yBinN/stripStep; yy++){
                            Charge += integrateMap->GetBinContent(kx - region.lowX_strip + 1, yy);
                        }
                        Charge *= area;
                        xStrip->charge += Charge;
                    }
                    for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                        GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                        Double_t Charge = 0.;
                        for(int xx = 1; xx <= xBinN/stripStep; xx++){
                            Charge += integrateMap->GetBinContent(xx, ky - region.lowY_strip + 1);
                        }
                        Charge *= area;
                        yStrip->charge += Charge;
                    }
                    delete integrateMap;
                }
                //Simulate the pulse shape and get the ADC value for each strip this particle hits
                double startTime = t_hit + minTravelT + trigger_jitter - 5.; // -5ns (trigger lattency) to move the pulse to the left
                for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                    GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                    for(int t = 0; t < samplingNum; t++){
                        double time = double(t+0.5) * samplingTime - startTime;
                        xStrip->pulse[t] += PulseShape(time, xStrip->charge, pulseTimeLength);
                    }
                    xStrip->charge = 0; // clear charge after getting the pulse
                }
                for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                    GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                    for(int t = 0; t < samplingNum; t++){
                        double time = double(t+0.5) * samplingTime - startTime;
                        yStrip->pulse[t] += PulseShape(time, yStrip->charge, pulseTimeLength);
                    }
                    yStrip->charge = 0; // clear charge after getting the pulse
                }
                bg++;
            }
        }
        // Already have the pulse shape for each strip, now convert to ADC value, add pedestal noise, threshold cut, and save the results

        Int_t* stripNumberX[4] = {
            stripNumberX_0, stripNumberX_1, stripNumberX_2, stripNumberX_3};
        Double_t (*pulseSamplesX[4])[samplingNum] = {
            pulseSamplesX_0, pulseSamplesX_1, pulseSamplesX_2, pulseSamplesX_3};
        Double_t* peakXOut[4] = {peakX_0, peakX_1, peakX_2, peakX_3};
        Double_t* totalChargeX[4] = {
            totalChargeX_0, totalChargeX_1, totalChargeX_2, totalChargeX_3};
        Int_t* peakSampleX[4] = {
            peakSampleX_0, peakSampleX_1, peakSampleX_2, peakSampleX_3};

        Int_t* stripNumberY[4] = {
            stripNumberY_0, stripNumberY_1, stripNumberY_2, stripNumberY_3};
        Double_t (*pulseSamplesY[4])[samplingNum] = {
            pulseSamplesY_0, pulseSamplesY_1, pulseSamplesY_2, pulseSamplesY_3};
        Double_t* peakYOut[4] = {peakY_0, peakY_1, peakY_2, peakY_3};
        Double_t* totalChargeY[4] = {
            totalChargeY_0, totalChargeY_1, totalChargeY_2, totalChargeY_3};
        Int_t* peakSampleY[4] = {
            peakSampleY_0, peakSampleY_1, peakSampleY_2, peakSampleY_3};

        for(int did = 0; did < 4; did++){
            fireNumX[did] = 0;
            fireNumY[did] = 0;

            for(int kx = 0; kx < map_x_bins; kx++) {
                GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                GetPedestal(xStrip->pulse);
                ADCConvert(0., 1., 12, xStrip->pulse);

                Double_t peakVal = 0.;
                Double_t totalCharge = 0.;
                Int_t peakSample = -1;
                for(int t = 0; t < samplingNum; t++){
                    totalCharge += xStrip->pulse[t];
                    if(xStrip->pulse[t] > peakVal){
                        peakVal = xStrip->pulse[t];
                        peakSample = t;
                    }
                }

                if(peakVal > fADCThreshold &&
                   totalCharge > fChargeThreshold) {
                    const Int_t firedIndex = fireNumX[did];
                    stripNumberX[did][firedIndex] = kx;
                    for(int t = 0; t < samplingNum; t++) {
                        pulseSamplesX[did][firedIndex][t] = xStrip->pulse[t];
                    }
                    peakXOut[did][firedIndex] = peakVal;
                    totalChargeX[did][firedIndex] = totalCharge;
                    peakSampleX[did][firedIndex] = peakSample;
                    fireNumX[did]++;
                }
            }

            for(int ky = 0; ky < map_y_bins; ky++) {
                GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                GetPedestal(yStrip->pulse);
                ADCConvert(0., 1., 12, yStrip->pulse);

                Double_t peakVal = 0.;
                Double_t totalCharge = 0.;
                Int_t peakSample = -1;
                for(int t = 0; t < samplingNum; t++){
                    totalCharge += yStrip->pulse[t];
                    if(yStrip->pulse[t] > peakVal){
                        peakVal = yStrip->pulse[t];
                        peakSample = t;
                    }
                }

                if(peakVal > fADCThreshold &&
                   totalCharge > fChargeThreshold) {
                    const Int_t firedIndex = fireNumY[did];
                    stripNumberY[did][firedIndex] = ky;
                    for(int t = 0; t < samplingNum; t++) {
                        pulseSamplesY[did][firedIndex][t] = yStrip->pulse[t];
                    }
                    peakYOut[did][firedIndex] = peakVal;
                    totalChargeY[did][firedIndex] = totalCharge;
                    peakSampleY[did][firedIndex] = peakSample;
                    fireNumY[did]++;
                }
            }
        }

        outputEventNum = i;
        FiredStripsT->Fill();
    }

    FiredStripsT->Write();
    FiredStripsF->Close();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "运行时间: " << diff.count() << " 秒" << std::endl;
}
