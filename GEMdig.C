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

struct EntryCounts {
    Long64_t ndata = 0;
    Long64_t nbg = 0;
    Int_t bgEntry = 0;
};

EntryCounts Rsolve_bgEntry(TChain* signal, TChain* bg, int bgFiles){
    // 8.0175e-04 inverse picobarn for 1e6 events
    double luminosity = 8.0175e-06; // pb^-1, 1e4 events
    double C_density = 2.267; // g/cm3, density of carbon
    double C_targetThickness = 1.e-4; // cm, 1 um carbon target
    double C_molarMass = 12.011; // g/mol
    double AvogadroConstant = 6.02214076e23; // mol^-1
    double C_arealDensity =
        C_density * C_targetThickness / C_molarMass * AvogadroConstant; // nuclei/cm2
    double beamN_signal = luminosity * 1.e36 / C_arealDensity;

    Long64_t Ndata = signal->GetEntries();
    const Long64_t availableNbg = bg->GetEntries();
    cout << "Number of entries in all signal data: " << Ndata << endl;
    cout << "Number of entries in all background data: "
         << availableNbg << endl;

    if(Ndata <= 0 || availableNbg <= 0 || bgFiles <= 0) return {};

    bg->GetEntry(availableNbg - 1);
    const double beamN_bg =
        static_cast<double>(EventID) * bgFiles;
    cout << "signal beamN: " << beamN_signal << endl;
    cout << "background beamN: " << beamN_bg << endl;
    if(beamN_bg <= 0.) return {};

    const double beam_current = 20.; // nA
    const double beam_rate =
        beam_current * 1.e-9 / 1.602176634e-19; // Hz
    const double beamN_window =
        400. * 1.e-9 * beam_rate; // beam particles in a 400 ns window
    const Int_t bgEntry = static_cast<Int_t>(
        static_cast<double>(availableNbg) / beamN_bg * beamN_window / 10.)+1;

    Long64_t Nbg = Ndata * static_cast<Long64_t>(bgEntry);
    if(Nbg > availableNbg) {
        cout << "Not enough background entries." << endl;
        Ndata = availableNbg / bgEntry;
        Nbg = Ndata * static_cast<Long64_t>(bgEntry);
        cout << "Use " << Ndata << " signal entries and " << Nbg << " background entries." << endl;
    }

    return {Ndata, Nbg, bgEntry};
}

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
        Double_t vvv = (val + off)/gain;
        Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
        if( vvv > saturation ) vvv = saturation;
        vvv -= off;
        Double_t dval =
            static_cast<Double_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));
        if( dval < 0 ) dval = 0;
        pulseXorY[n] = dval;
    }
}

void AccumulateIonCharge(GEMStripMap& gemStripMap,
                         const GEMStripRegion& region,
                         const IonPar_t& ion)
{
    constexpr Int_t stripStep = 6;
    const Int_t xStripN = region.upX_strip + 1 - region.lowX_strip;
    const Int_t yStripN = region.upY_strip + 1 - region.lowY_strip;
    const Int_t xBinN = xStripN * stripStep;
    const Int_t yBinN = yStripN * stripStep;

    constexpr Int_t maxRegionStripN = 2 * activeStripNhalf + 1;
    if(xStripN > maxRegionStripN || yStripN > maxRegionStripN) return;

    std::array<Double_t, maxRegionStripN> stripChargeX = {};
    std::array<Double_t, maxRegionStripN> stripChargeY = {};

    const Double_t xBinWidth = (region.xHigh - region.xLow) / xBinN;
    const Double_t yBinWidth = (region.yHigh - region.yLow) / yBinN;
    const Double_t area = xBinWidth * yBinWidth;
    constexpr Int_t edgeDeadSteps = 0; //stripStep / 4;
    const Double_t effSigma2 = ion.eff_sigma * ion.eff_sigma;
    const Double_t normalization =
        fAvaGain * ion.ggnorm / (TMath::Pi() * ion.eff_sigma);

    for(Int_t ix = 0; ix < xBinN; ++ix) {
        const Double_t xCenter =
            region.xLow + (static_cast<Double_t>(ix) + 0.5) * xBinWidth;
        const Double_t dx = ion.Xp - xCenter;
        const Double_t dx2 = dx * dx;
        const Int_t xStripIndex = ix / stripStep;
        const Int_t xStepInStrip = ix % stripStep;
        const bool xStepActive =
            xStepInStrip >= edgeDeadSteps &&
            xStepInStrip < stripStep - edgeDeadSteps;

        for(Int_t iy = 0; iy < yBinN; ++iy) {
            const Double_t yCenter =
                region.yLow + (static_cast<Double_t>(iy) + 0.5) * yBinWidth;
            const Double_t dy = ion.Yp - yCenter;
            if(dx2 + dy * dy > 1.7 * ion.R2) continue; // only consider bins within 3 sigma
            const Int_t yStepInStrip = iy % stripStep;
            const bool yStepActive =
                yStepInStrip >= edgeDeadSteps &&
                yStepInStrip < stripStep - edgeDeadSteps;
            const Double_t density =
                normalization * effSigma2 / (dx2 + dy * dy + effSigma2);

            if(xStepActive) stripChargeX[xStripIndex] += density;
            if(yStepActive) stripChargeY[iy / stripStep] += density;
        }
    }

    for(Int_t localX = 0; localX < xStripN; ++localX) {
        GEMStrip* strip =
            gemStripMap.GetXStrip(region.lowX_strip + localX);
        strip->charge += stripChargeX[localX] * area;
    }
    for(Int_t localY = 0; localY < yStripN; ++localY) {
        GEMStrip* strip =
            gemStripMap.GetYStrip(region.lowY_strip + localY);
        strip->charge += stripChargeY[localY] * area;
    }
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

void GEMdig(Long64_t maxEvent = -1) {

    auto start = std::chrono::high_resolution_clock::now();
    gRandom->SetSeed(0);//1232 // 0 表示用系统时间自动生成种子


    std::vector<string> signalFiles =
        FindRootFiles("2105", "signal");
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
    
    const EntryCounts entryCounts =
        Rsolve_bgEntry(
            &signal, &background, static_cast<int>(backgroundFiles.size()));
    Long64_t Ndata = entryCounts.ndata;
    Long64_t Nbg = entryCounts.nbg;
    Int_t bg_entry = entryCounts.bgEntry;

    bg_entry = 0;

    if(maxEvent > 0 && maxEvent < Ndata) {
        Ndata = maxEvent;
        Nbg = Ndata * static_cast<Long64_t>(bg_entry);
    }
    
    cout << "Background entries per signal event: " << bg_entry << endl;
    cout << "Use " << Ndata << " signal entries and "
         << Nbg << " background entries." << endl;
    TH1D *peakX_signal = new TH1D("peakX_signal", "peakX of the pulse", 250/2, 0, 2500);
    TH1D *peakY_signal = new TH1D("peakY_signal", "peakY of the pulse", 250/2, 0, 2500);
    TH1D *chargeX_signal = new TH1D("chargeX_signal", "total charge on X strips", 500/2, 0, 5000);
    TH1D *chargeY_signal = new TH1D("chargeY_signal", "total charge on Y strips", 500/2, 0, 5000);
    TH1D *maxBinX_signal = new TH1D("maxBinX_signal", "max sample bin of the pulse X", samplingNum, 0, samplingNum);
    TH1D *maxBinY_signal = new TH1D("maxBinY_signal", "max sample bin of the pulse Y", samplingNum, 0, samplingNum);
    TH1D *clusterSizeX_signal = new TH1D("clusterSizeX_signal", "cluster size of X strips", 15, 0.5, 15.5);
    TH1D *clusterSizeY_signal = new TH1D("clusterSizeY_signal", "cluster size of Y strips", 15, 0.5, 15.5);

    TH1D *Occupancy_X[4];
    TH1D *Occupancy_Y[4];
    for(int did = 0; did < 4; did++){
        Occupancy_X[did] = new TH1D(Form("Occupancy of chamber%d(X strips)", did), Form("Occupancy of chamber%d (X strips)", did), map_x_bins, 0-0.5, map_x_bins-0.5);
        Occupancy_Y[did] = new TH1D(Form("Occupancy of chamber%d(Y strips)", did), Form("Occupancy of chamber%d (Y strips)", did), map_y_bins, 0-0.5, map_y_bins-0.5);
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
    int totalEvent[4] = {0, 0, 0, 0}, validEvent[4] = {0, 0, 0, 0}; 

    double validEntries = 0;
    double totalOccupancy_X[4] = {0, 0, 0, 0};
    double totalOccupancy_Y[4] = {0, 0, 0, 0};

    for(int i=0;i<Ndata;i++) { //max 9000 for 4um
        if (i%5 == 0) cout << "start " << i << "\r" << flush;        

        //clear strip ADC information for each event
        for(int did = 0; did < 4; did++){
            gemMap[did].Clear();
        }

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
                AccumulateIonCharge(gemMap[did], region, ion);
            }
            //Simulate the pulse shape and get the ADC value for each strip this particle hits
            double startTime = t_hit + minTravelT + trigger_jitter - 120.; // -5ns (trigger lattency) to move the pulse to the left
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
            totalEvent[did]++;
        }
        if(sig < 1) continue;

        //fill the histograms for signal particles
        for(int did = 0; did < 4; did++){
            vector<double> StripPeaksX;
            vector<double> StripChargesX;
            int fireX = 0;
            for(int kx = 0; kx < map_x_bins; kx++){
                GEMStrip* xStrip = gemMap[did].GetXStrip(kx);
                if(xStrip->pulse[2] < 5.) continue;
                GetPedestal(xStrip->pulse);
                ADCConvert(xStrip->pedestal_line, 1., 11, xStrip->pulse);
                double peak = *std::max_element(xStrip->pulse.begin(), xStrip->pulse.end());
                double totalCharge = std::accumulate(xStrip->pulse.begin(), xStrip->pulse.end(), 0.);
                int peakSample = std::distance(xStrip->pulse.begin(), std::max_element(xStrip->pulse.begin(), xStrip->pulse.end()));
                if(peak > fADCThreshold && totalCharge > fChargeThreshold) {
                    StripPeaksX.push_back(peak);
                    StripChargesX.push_back(totalCharge);
                    maxBinX_signal->Fill(peakSample);
                    fireX++;
                }
            }
            double eventPeakX = StripPeaksX.empty() ? 0. : *std::max_element(StripPeaksX.begin(), StripPeaksX.end());
            double eventChargeX = StripChargesX.empty() ? 0. : *std::max_element(StripChargesX.begin(), StripChargesX.end());
            if(eventPeakX > fADCThreshold && eventChargeX > fChargeThreshold) {
                peakX_signal->Fill(eventPeakX);
                chargeX_signal->Fill(eventChargeX);
            }
            clusterSizeX_signal->Fill(fireX);


            vector<double> StripPeaksY;
            vector<double> StripChargesY;
            int fireY = 0;
            for(int ky = 0; ky < map_y_bins; ky++){
                GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                if(yStrip->pulse[2] < 5.) continue;
                GetPedestal(yStrip->pulse);
                ADCConvert(yStrip->pedestal_line, 1., 11, yStrip->pulse);
                double peak = *std::max_element(yStrip->pulse.begin(), yStrip->pulse.end());
                double totalCharge = std::accumulate(yStrip->pulse.begin(), yStrip->pulse.end(), 0.);
                int peakSample = std::distance(yStrip->pulse.begin(), std::max_element(yStrip->pulse.begin(), yStrip->pulse.end()));
                if(peak > fADCThreshold && totalCharge > fChargeThreshold) {
                    StripPeaksY.push_back(peak);
                    StripChargesY.push_back(totalCharge);
                    maxBinY_signal->Fill(peakSample);
                    fireY++;
                }
            }
            double eventPeakY = StripPeaksY.empty() ? 0. : *std::max_element(StripPeaksY.begin(), StripPeaksY.end());
            double eventChargeY = StripChargesY.empty() ? 0. : *std::max_element(StripChargesY.begin(), StripChargesY.end());
            if(eventPeakY > fADCThreshold && eventChargeY > fChargeThreshold) {
                peakY_signal->Fill(eventPeakY);
                chargeY_signal->Fill(eventChargeY);
            }
            clusterSizeY_signal->Fill(fireY);

            if(fireX > 0 && fireY > 0)
                validEvent[did]++;
        }

        validEntries++;
        //for the background particles
        int bg = 0;
        for(int ii=i*bg_entry; ii<(i+1)*bg_entry; ii++) {
            background.GetEntry(ii);
            Double_t t0 = gRandom->Uniform(-25.*8., 8.*25);

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
                    AccumulateIonCharge(gemMap[did], region, ion);
                }
                //Simulate the pulse shape and get the ADC value for each strip this particle hits
                double startTime = t_hit + minTravelT + trigger_jitter - 50.; // -5ns (trigger lattency) to move the pulse to the left
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
                if(xStrip->pulse[2] < 5.) continue; // skip strips with very low signal to save time
                GetPedestal(xStrip->pulse);
                ADCConvert(xStrip->pedestal_line, 1., 11, xStrip->pulse);

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
                    Occupancy_X[did]->Fill(kx);
                    totalOccupancy_X[did] += 1;
                }
            }

            for(int ky = 0; ky < map_y_bins; ky++) {
                GEMStrip* yStrip = gemMap[did].GetYStrip(ky);
                if(yStrip->pulse[2] < 5.) continue; // skip strips with very low signal to save time
                GetPedestal(yStrip->pulse);
                ADCConvert(yStrip->pedestal_line, 1., 11, yStrip->pulse);

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
                    Occupancy_Y[did]->Fill(ky);
                    totalOccupancy_Y[did] += 1;
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

    for(int did = 0; did < 4; did++) {
        std::cout << "Detector " << did << ": Total Events = " << totalEvent[did] << ", Valid Events = " << validEvent[did] << ", Efficiency: " << (totalEvent[did] > 0 ? static_cast<double>(validEvent[did]) / totalEvent[did] : 0.) << std::endl;
    }

    for(int did = 0; did < 4; did++) {
        totalOccupancy_X[did] = totalOccupancy_X[did] / map_x_bins / validEntries;
        totalOccupancy_Y[did] = totalOccupancy_Y[did] / map_y_bins / validEntries;
        std::cout << "Detector " << did << ": Total Occupancy X = " << totalOccupancy_X[did]*100. << "%, Total Occupancy Y = " << totalOccupancy_Y[did]*100. << "%" << std::endl;
    }

    TCanvas *c1 = new TCanvas("c1", "Peak Distributions", 800, 600);
    c1->Divide(2, 2);
    c1->cd(1);
    peakX_signal->Draw();
    c1->cd(2);
    peakY_signal->Draw();
    c1->cd(3);
    chargeX_signal->Draw();
    c1->cd(4);
    chargeY_signal->Draw();

    TCanvas *c2 = new TCanvas("c2", "Max Sample Bin Distributions", 800, 600);
    c2->Divide(2);
    c2->cd(1);
    maxBinX_signal->Draw();
    c2->cd(2);
    maxBinY_signal->Draw();

    TCanvas *c3 = new TCanvas("c3", "Cluster Size Distributions", 800, 600);
    c3->Divide(2, 1);
    c3->cd(1);
    clusterSizeX_signal->Draw();
    c3->cd(2);
    clusterSizeY_signal->Draw();

    for(int did = 0; did < 4; did++) {
        Occupancy_X[did]->Scale(100.0 / validEntries);
        Occupancy_Y[did]->Scale(100.0 / validEntries);
    }
    TCanvas *c4 = new TCanvas("c4", "Occupancy Distributions X", 800, 600);
    c4->Divide(2, 2);
    for(int did = 0; did < 4; did++) {
        c4->cd(did+1);
        Occupancy_X[did]->Draw();
    }
    TCanvas *c5 = new TCanvas("c5", "Occupancy Distributions Y", 800, 600);
    c5->Divide(2, 2);
    for(int did = 0; did < 4; did++) {
        c5->cd(did+1);
        Occupancy_Y[did]->Draw();
    }
    gSystem->ProcessEvents(); // Update the canvases to display the plots

    TFile *outputF = new TFile("output/GEMdig_output_test.root", "RECREATE");
    outputF->cd();
    peakX_signal->Write();
    peakY_signal->Write();
    chargeX_signal->Write();
    chargeY_signal->Write();
    maxBinX_signal->Write();
    maxBinY_signal->Write();
    clusterSizeX_signal->Write();
    clusterSizeY_signal->Write();
    outputF->Close();
}
