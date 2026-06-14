#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

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
#include "GEMAvalancheModel.h"
#include "GEMIonModel.h"
#include "GEMStripMap.h"

using std::cout;
using std::endl;
using std::string;

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
void GetPedestal(Double_t pulseXorY[samplingNum]){
    for(int n = 0; n < samplingNum; n++){
        pulseXorY[n] += gRandom->Gaus(0., 20.); // pedestal noise, 20 ADC
    }
}

void ADCConvert(Double_t off, Double_t gain, Int_t bits,
                Double_t pulseXorY[samplingNum])
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
    string folder = "data/";
    filename = folder + filename + ".root";
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
    int Num = 10; //number of files
    TChain signal("T");
    TChain background("T");
    loadFile(signal, "Uniform_signal_660-1540MeV");
    loadFile(background, "x17_He_bot_2.2GeV_4um_1");
    
    const int Ndata = signal.GetEntries();
    const int Nbg = background.GetEntries();
    cout << "Number of entries in the signal data: " << Ndata << endl;
    cout << "Number of entries in one background file: " << Nbg << endl;
    background.GetEntry(Nbg-1);
    int bg_entry = int(double(Nbg) / double(EventID) * 312.5/4. * 350.);
    cout << "background entry numbers for each event: " << bg_entry << endl;
    //bg_entry = 0; //550
    bg_entry *= 2.;

    loadFile(background, "x17_He_bot_2.2GeV_4um_2");
    loadFile(background, "x17_He_bot_2.2GeV_4um_3");
    loadFile(background, "x17_He_bot_2.2GeV_4um_4");
    loadFile(background, "x17_He_bot_2.2GeV_4um_5");
    loadFile(background, "x17_He_bot_2.2GeV_4um_6");
    loadFile(background, "x17_He_bot_2.2GeV_4um_7");
    loadFile(background, "x17_He_bot_2.2GeV_4um_8");
    loadFile(background, "x17_He_bot_2.2GeV_4um_9");
    loadFile(background, "x17_He_bot_2.2GeV_4um_10");

    TH1D *peak = new TH1D("peak", "peak of the pulse", 100, 0, 5000);
    TH1D *peakY = new TH1D("peakY", "peakY of the pulse", 100, 0, 5000);
    TH2D *adc_Num = new TH2D("adc_Num", "Pulse Amplitude vs Fired Strip Number", 40, 0, 4000, 11, 0-0.5, 11-0.5);

    TH1D *Occupancy_X[4];
    TH1D *Occupancy_Y[4];
    for(int did = 0; did < 4; did++){
        Occupancy_X[did] = new TH1D(Form("Occupancy of chamber%d(X strips)", did), Form("Occupancy of chamber%d (X strips)", did), map_x_bins, 0, map_x_bins);
        Occupancy_Y[did] = new TH1D(Form("Occupancy of chamber%d(Y strips)", did), Form("Occupancy of chamber%d (Y strips)", did), map_y_bins, 0, map_y_bins);
    }

    TFile *outputFile = new TFile("output_4mm_200percent.root", "RECREATE");

    TFile *FiredStripsF = new TFile("FiredStrips_4mm_200percent.root", "RECREATE");
    TTree *FiredStripsT = new TTree("T", "T");
    int NumX_f[4], NumY_f[4], Ientry;
    float StripX_f[4][800][6+3], StripY_f[4][1200][6+3]; // assuming a maximum of 4000 strips can fire per chamber
    FiredStripsT->Branch("Ientry", &Ientry, "Ientry/I");
    FiredStripsT->Branch("NumX", &NumX_f[0], "NumX[4]/I");
    FiredStripsT->Branch("NumY", &NumY_f[0], "NumY[4]/I");
    FiredStripsT->Branch("StripX", &StripX_f[0], "StripX[4][800][9]/F");
    FiredStripsT->Branch("StripY", &StripY_f[0], "StripY[4][1200][9]/F");

    GEMStripMap gemMap[4]; // strip map for the four chambers, indexed by detector ID
    for(int did = 0; did < 4; did++) gemMap[did] = GEMStripMap(did);

    double validEntries = 0;
    for(int i=0;i<9000;i++) { //max 9000 for 4um
        if (i%100 == 0) cout << "start " << i << endl;        
        signal.GetEntry(i);

        //clear strip ADC information for each event
        for(int did = 0; did < 4; did++){
            gemMap[did].Clear();
        }
        
        int sig = 0;
        for(int n = 0; n < GasN; n++){
            //for the one beam particle
            if( GasPTID[n] != 0 || GasZout[n] - GasZ[n] < 2.999) continue;
            int did = GasDID[n];
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
            }
            //Simulate the pulse shape and get the ADC value for each strip this particle hits
            double startTime = minTravelT + trigger_jitter - 5.; // -5ns (trigger lattency) to move the pulse to the left 
            for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                for(int t = 0; t < samplingNum; t++){
                    double time = double(t+0.5) * samplingTime + startTime;
                    xStrip->pulse[t] += PulseShape(time, xStrip->charge, pulseTimeLength);
                }
                xStrip->charge = 0; // clear charge after getting the pulse
            }
            for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                for(int t = 0; t < samplingNum; t++){
                    double time = double(t+0.5) * samplingTime + startTime;
                    yStrip->pulse[t] += PulseShape(time, yStrip->charge, pulseTimeLength);
                }
                yStrip->charge = 0; // clear charge after getting the pulse
            }
            sig++;
        }
        if(sig < 1) continue;
        validEntries++;
        //for the background particles
        int bg = 0;
        for(int ii=i*bg_entry; ii<(i+1)*bg_entry; ii++) {
            background.GetEntry(ii);
            Double_t t0 = gRandom->Uniform(-200., 6.*25.);

            for(int n = 0; n < GasN; n++){
                int did = GasDID[n];
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
                }
                //Simulate the pulse shape and get the ADC value for each strip this particle hits
                double startTime = minTravelT + trigger_jitter - 5. + t0; // -5ns (trigger lattency) to move the pulse to the left
                for(int kx = region.lowX_strip; kx <= region.upX_strip; kx++) {
                    GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                    for(int t = 0; t < samplingNum; t++){
                        double time = double(t+0.5) * samplingTime + startTime;
                        xStrip->pulse[t] += PulseShape(time, xStrip->charge, pulseTimeLength);
                    }
                    xStrip->charge = 0; // clear charge after getting the pulse
                }
                for(int ky = region.lowY_strip; ky <= region.upY_strip; ky++) {
                    GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                    for(int t = 0; t < samplingNum; t++){
                        double time = double(t+0.5) * samplingTime + startTime;
                        yStrip->pulse[t] += PulseShape(time, yStrip->charge, pulseTimeLength);
                    }
                    yStrip->charge = 0; // clear charge after getting the pulse
                }
                bg++;
            }
        }
        // Already have the pulse shape for each strip, now convert to ADC value, add pedestal noise, threshold cut, and save the results
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "运行时间: " << diff.count() << " 秒" << std::endl;
}
