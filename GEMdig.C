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

void findRegion(Double_t x, Double_t y, Int_t Num_half) {
    int centerX = int( (x - map_x_low) / stripWidth );
    int centerY = int( (y - map_y_low) / stripWidth );
    lowX = centerX - Num_half;
    lowY = centerY - Num_half;
    upX = centerX + Num_half;
    upY = centerY + Num_half;
    if(lowX < 0) lowX = 0; 
    if(lowY < 0) lowY = 0;
    if(upX > map_x_bins - 1) upX = map_x_bins - 1;
    if(upY > map_y_bins - 1) upY = map_y_bins - 1;
}

//convert charge number to the pulse on each strip
void GetPulse(Double_t* stripX, Double_t* stripY,
              Double_t Xpos, Double_t Ypos, Double_t startTime,
              Double_t pulseX[map_x_bins][samplingNum], Double_t pulseY[map_y_bins][samplingNum]){
    
    findRegion(Xpos, Ypos, activeStripNhalf);
    for(int i = lowX; i < upX+1; i++) {
        for(int n = 0; n < samplingNum; n++){
            Double_t t = n * samplingTime - startTime;
            pulseX[i][n] += PulseShape(t, stripX[i], pulseTimeLength);
        }
    }
    for(int i = lowY; i < upY+1; i++) {
        for(int n = 0; n < samplingNum; n++){
            Double_t t = n * samplingTime - startTime;
            pulseY[i][n] += PulseShape(t, stripY[i], pulseTimeLength);
        }
    }
}

void IonModel(TVector3& in, TVector3& out, const Double_t Edep) {
    //all the unit for position is mm
    fRIon.clear();
    if(out.Z() < in.Z()) {
        TVector3 buffer(in);
        in = out;
        out = buffer;
    }
    TVector3 dir = out - in;
    TVector3 track(dir);
    dir.SetMag(1.); // Normalize the direction vector
    //beam particle's position on the readout plane
    TVector3 readoutPos(out);
    readoutPos = readoutPos + dir * Zreadout;

    //number of primary ions
    if(Edep*1.e+6 > GasIonE){
        Num_ions = Edep*1.e+6 / GasIonE;
        Num_ions = gRandom->Poisson(Num_ions);
    }
    else Num_ions = 0;
    if(Num_ions > 10000) Num_ions = 10000;
    //cout << "Number of ions: " << Num_ions << endl;

    minTravelT = 1000000.;
    
    for(int j = 0; j < Num_ions; j++) {
        // Generate a random start position along the direction vector
        Double_t fraction = gRandom->Uniform(0., 1.0); // Random fraction along the direction vector
        TVector3 ionPos = in + fraction * track;

        //Get the traveling time of the ion
        Double_t Length = Zreadout + out.Z() - ionPos.Z();
        Double_t Time = Length / GasDriftVelocity; // in us

        // spatial sigma on readout [mm] after avalanche
        Double_t Spatial_sigma = sqrt(2. * GasDiffusion * Time);

        //Add the lateral uncertainty(get to the readout plane)
        ionPos.SetX(ionPos.X() + gRandom->Gaus(0., GasLateralUncertainty));
        ionPos.SetY(ionPos.Y() + gRandom->Gaus(0., GasLateralUncertainty));

        //Get total charge of one primary ion pair after avalanche
        //Double_t Charge = gRandom->Exp(gain_mean); // exponential distribution of the charge
        Double_t Gain0 = 20.;
        Double_t Charge = gRandom->Gaus(gain_mean, gain_mean/sqrt(Gain0)); // Gaussian distribution of the charge
        if(Charge < 0.) Charge = 0.;
        IonPar_t ip;
        ip.Xp = ionPos.X();
        ip.Yp = ionPos.Y();
        ip.Charge = Charge;
        ip.S_sigma = Spatial_sigma;
        ip.eff_sigma = Spatial_sigma * Spatial_sigma;
        ip.Ttime = Time * 1.e+3; // convert to ns
        minTravelT = TMath::Min(minTravelT, ip.Ttime);
        ip.ggnorm = ip.Charge * TMath::InvPi() / pow(Spatial_sigma*fSNormNsigma, 2.);

        fRIon.push_back(ip);
        //cout << "X: " << ionPos.X() << " Y: " << ionPos.Y() << " Charge: " << Charge << " time: " << ip.Ttime << " Spatial sigma: " << ip.S_sigma << endl;
    }
}

void avalanche_model(Double_t pulseX[map_x_bins][samplingNum], Double_t pulseY[map_y_bins][samplingNum], Double_t latencyT){

    double total_charge = 0;
    for(int i=0; i<Num_ions;i++){
        //find the active strip, only sum these strips to save time
        findRegion(fRIon[i].Xp, fRIon[i].Yp, activeStripNhalf);

        Int_t xStripN = upX + 1 - lowX;
        Int_t yStripN = upY + 1 - lowY;
        int stripStep = 6;
        int xBinN = xStripN * stripStep;
        int yBinN = yStripN * stripStep;
        Double_t area = 0.4*0.4/stripStep/stripStep;
        TH2D *integrateMap = new TH2D("inegrate map", "integrate map",
            xBinN, stripMap->GetXaxis()->GetBinLowEdge(lowX + 1), stripMap->GetXaxis()->GetBinUpEdge(upX + 1),
            yBinN, stripMap->GetYaxis()->GetBinLowEdge(lowY + 1), stripMap->GetYaxis()->GetBinUpEdge(upY + 1));
        
        Double_t eff_sigma = fRIon[i].eff_sigma;
        for(int ix = 0; ix < xBinN; ix++){
            Double_t dx = fRIon[i].Xp - integrateMap->GetXaxis()->GetBinCenter(ix+1);
            Double_t dx2 = dx * dx;
            for(int iy = 0; iy < yBinN; iy++){
                if(ix % stripStep == 0 || ix % stripStep == 1 || ix % stripStep == 4 || ix % stripStep == 5){
                    if(iy % stripStep == 0 || iy % stripStep == 5){
                        if( int(floor(ix/(stripStep/2.))) % 2 > 0.5 )
                            continue;
                    }
                }
                Double_t dy = fRIon[i].Yp - integrateMap->GetYaxis()->GetBinCenter(iy+1);
                Double_t dy2 = dy * dy;
                Double_t result = fAvaGain * fRIon[i].ggnorm * (1./(TMath::Pi()*eff_sigma)) * (eff_sigma*eff_sigma) / ((dx2+dy2) + eff_sigma*eff_sigma);
                integrateMap->SetBinContent(ix+1, iy+1, result); //the charge deposit in each small area
            }
        }

        //Get integrated charge on each strip
        Double_t stripX[map_x_bins] = {};
        Double_t stripY[map_y_bins] = {};
        int nxx = -1;
        for(int kx = lowX; kx < upX+1; kx++) {
            for(int xx = 0; xx < stripStep; xx++){
                nxx++;
                if(xx == 0 || xx == 1 || xx == 4 || xx == 5) continue;
                for(int yy = 0; yy < yBinN; yy++){
                    stripX[kx] += integrateMap->GetBinContent(nxx+1, yy+1);
                }
            }
            stripX[kx] *= area;
        }
        int nyy = -1;
        for(int ky = lowY; ky < upY+1; ky++) {
            for(int yy = 0; yy < stripStep; yy++){
                nyy++;
                if(yy == 0 || yy == 5) continue;
                for(int xx = 0; xx < xBinN; xx++){
                    if(int(floor(xx/(stripStep/2.))) % 2 > 0.5) continue;
                    stripY[ky] += integrateMap->GetBinContent(xx+1, nyy+1);
                }
            }
            stripY[ky] *= area;
        }
        
        //convert charge number to pulse on each strip
        GetPulse(stripX, stripY, fRIon[i].Xp, fRIon[i].Yp, /*fRIon[i].Ttime*/minTravelT + latencyT, pulseX, pulseY);
        delete integrateMap;
    }
    //cout << total_charge << " total charge in the avalanche model" << endl;
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
    file.SetBranchAddress("GEM_Gas.N", &GasN);
    file.SetBranchAddress("GEM_Gas.X", &GasX[0]);
    file.SetBranchAddress("GEM_Gas.Y", &GasY[0]);
    file.SetBranchAddress("GEM_Gas.Z", &GasZ[0]);
    file.SetBranchAddress("GEM_Gas.Time", &GasTime[0]);
    file.SetBranchAddress("GEM_Gas.Xout", &GasXout[0]);
    file.SetBranchAddress("GEM_Gas.Yout", &GasYout[0]);
    file.SetBranchAddress("GEM_Gas.Zout", &GasZout[0]);
    file.SetBranchAddress("GEM_Gas.PID", &GasPID[0]);
    file.SetBranchAddress("GEM_Gas.PTID", &GasPTID[0]);
    file.SetBranchAddress("GEM_Gas.TID", &GasTID[0]);
    file.SetBranchAddress("GEM_Gas.DID", &GasDID[0]);
    file.SetBranchAddress("GEM_Gas.Theta", &GasTheta[0]);
    file.SetBranchAddress("GEM_Gas.P", &GasP[0]);
    file.SetBranchAddress("GEM_Gas.Edep", &GasEdep[0]);
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

    double validEntries = 0;
    for(int i=0;i<9000;i++) { //max 9000 for 4um
        if (i%100 == 0) cout << "start " << i << endl;        
        signal.GetEntry(i);
        Double_t pulseX[4][map_x_bins][samplingNum] = {}; // DID = 0, 1, 2, 3
        Double_t pulseY[4][map_y_bins][samplingNum] = {};
        Double_t hitX, hitY;
        std::vector<double> bg_hitX;
        std::vector<double> bg_hitY;
        std::vector<double> bg_hitDID;
        int sig = 0;
        for(int n = 0; n < GasN; n++){
            //for the one beam particle
            if( GasPTID[n] == 0 && GasZout[n] - GasZ[n] > 2.999) {
                if(GasDID[n] == 0 || GasDID[n] == 1) {
                    hitX = GasX[n];
                    hitY = GasY[n];
                }
                if(GasDID[n] == 1 || GasDID[n] == 3) {
                    GasX[n] = -GasX[n];
                    GasXout[n] = -GasXout[n];
                }
                TVector3 in(GasX[n], GasY[n], GasZ[n]);
                TVector3 out(GasXout[n], GasYout[n], GasZout[n]);
                IonModel(in, out, GasEdep[n]);
                if(Num_ions < 1) continue;

                Double_t trigger_jitter = gRandom->Gaus(0, fTriggerJitter);
                Double_t apvJitter = (gRandom->Uniform(fAPVTimeJitter) - fAPVTimeJitter/2.);
                trigger_jitter += apvJitter;
                trigger_jitter=0.;
                avalanche_model(pulseX[GasDID[n]], pulseY[GasDID[n]], -minTravelT+trigger_jitter - 5.); //-time to left move
                double max = 0.;
                double maxY = 0.;
                int fireNum = 0;
                for(int k=0; k<map_x_bins; k++){
                    double mean = TMath::Mean(samplingNum, &pulseX[GasDID[n]][k][0]);
                    if(mean > 20.*5. / sqrt(6.)) fireNum++;
                    for(int s=0; s<samplingNum; s++){
                        if(pulseX[GasDID[n]][k][s] > max) max = pulseX[GasDID[n]][k][s];
                    }
                }
                for(int k=0; k<map_y_bins; k++){
                    double mean = TMath::Mean(samplingNum, &pulseY[GasDID[n]][k][0]);
                    for(int s=0; s<samplingNum; s++){
                        if(pulseY[GasDID[n]][k][s] > maxY) maxY = pulseY[GasDID[n]][k][s];
                    }
                }
                peak->Fill(max + gRandom->Gaus(0., 20.));
                peakY->Fill(maxY + gRandom->Gaus(0., 20.));
                adc_Num->Fill(max, fireNum);
                sig++;
            }
        }
        if(sig < 1) continue;
        validEntries++;
        //for the background particles
        int bg = 0;
        for(int ii=i*bg_entry; ii<(i+1)*bg_entry; ii++) {
            background.GetEntry(ii);
            Double_t t0 = gRandom->Uniform(-200., 6.*25.);
            for(int n = 0; n < GasN; n++){
                if(GasDID[n] == 1 || GasDID[n] == 3) {
                    GasX[n] = -GasX[n];
                    GasXout[n] = -GasXout[n];
                }
                TVector3 in(GasX[n], GasY[n], GasZ[n]);
                TVector3 out(GasXout[n], GasYout[n], GasZout[n]);
                IonModel(in, out, GasEdep[n]);
                if(Num_ions < 1) continue;
                if(GasDID[n] == 1 || GasDID[n] == 3) {
                    GasX[n] = -GasX[n];
                    GasXout[n] = -GasXout[n];
                }
                bg_hitX.push_back(GasX[n]);
                bg_hitY.push_back(GasY[n]);
                bg_hitDID.push_back(GasDID[n]);
                Double_t trigger_jitter = gRandom->Gaus(0, fTriggerJitter);
                Double_t apvJitter = (gRandom->Uniform(fAPVTimeJitter) - fAPVTimeJitter/2.);
                trigger_jitter += apvJitter;
                trigger_jitter=0.;
                Double_t bgT = t0 + GasTime[n]; // random time for the background
                avalanche_model(pulseX[GasDID[n]], pulseY[GasDID[n]], -minTravelT+trigger_jitter - 5. + bgT); //-time to left move
                bg++;
            }
        }

        TH2D* shape_pulse_x[4];
        TH2D* shape_pulse_y[4];
        for(int did=0; did<4; did++){
            shape_pulse_x[did] = new TH2D(Form("event%d_chamber%d_x_strips_signal", i, did), Form("event%d_chamber%d_x_strips_signal", i, did),
                                          samplingNum, 0.+0.5, samplingNum+0.5, map_x_bins, 0, map_x_bins);
            shape_pulse_y[did] = new TH2D(Form("event%d_chamber%d_y_strips_signal", i, did), Form("event%d_chamber%d_y_strips_signal", i, did),
                                          samplingNum, 0.+0.5, samplingNum+0.5, map_y_bins, 0, map_y_bins);
        }
        for(int did = 0; did < 4; did++){
            NumX_f[did] = 0;
            for(int k=0; k < map_x_bins; k++){
                double mean;
                int maxIndex = TMath::LocMax(samplingNum, &pulseX[did][k][0]);
                if(pulseX[did][k][maxIndex] > 1.) {
                    GetPedestal(pulseX[did][k]);
                    ADCConvert(0., 1., 12, pulseX[did][k]);
                    mean = TMath::Mean(samplingNum, &pulseX[did][k][0]);
                    maxIndex = TMath::LocMax(samplingNum, &pulseX[did][k][0]);
                }
                if( mean > 20. * 5. / sqrt(6.) && (maxIndex == 1 || maxIndex == 2 || maxIndex == 3) ){
                    Occupancy_X[did]->Fill(k);
                    for(int s=0; s<samplingNum; s++){
                        shape_pulse_x[did]->SetBinContent(s+1, k+1, pulseX[did][k][s]);
                    }
                    StripX_f[did][NumX_f[did]][0] = float(k); //strip index
                    StripX_f[did][NumX_f[did]][1] = float(pulseX[did][k][maxIndex]); //peak
                    StripX_f[did][NumX_f[did]][2] = float(maxIndex); //peak time sample
                    for(int n=3;n<samplingNum+3;n++){
                        StripX_f[did][NumX_f[did]][n] = float(pulseX[did][k][n-3]);
                    }
                    NumX_f[did]++;
                }
                else{
                    for(int s=0; s<samplingNum; s++){
                        pulseX[did][k][s] = 0.;
                    }
                }
            }
        }
        for(int did = 0; did < 4; did++){
            NumY_f[did] = 0;
            for(int k=0; k < map_y_bins; k++){
                double mean;
                int maxIndex = TMath::LocMax(samplingNum, &pulseY[did][k][0]);
                if(pulseY[did][k][maxIndex] > 1.) {
                    GetPedestal(pulseY[did][k]);
                    ADCConvert(0., 1., 12, pulseY[did][k]);
                    mean = TMath::Mean(samplingNum, &pulseY[did][k][0]);
                    maxIndex = TMath::LocMax(samplingNum, &pulseY[did][k][0]);
                }
                if( mean > 20. * 5. / sqrt(6.) && (maxIndex == 1 || maxIndex == 2 || maxIndex == 3) ){
                    Occupancy_Y[did]->Fill(k);
                    for(int s=0; s<samplingNum; s++){
                        shape_pulse_y[did]->SetBinContent(s+1, k+1, pulseY[did][k][s]);
                    }
                    StripY_f[did][NumY_f[did]][0] = float(k); //strip index
                    StripY_f[did][NumY_f[did]][1] = float(pulseY[did][k][maxIndex]); //peak
                    StripY_f[did][NumY_f[did]][2] = float(maxIndex); //peak time sample
                    for(int n=3;n<samplingNum+3;n++){
                        StripY_f[did][NumY_f[did]][n] = float(pulseY[did][k][n-3]);
                    }
                    NumY_f[did]++;
                }
                else{
                    for(int s=0; s<samplingNum; s++){
                        pulseY[did][k][s] = 0;
                    }  
                }
            }
        }
        bool over = false;
        for(int did = 0; did < 4; did++){
            if(NumX_f[did] > 800 || NumY_f[did] > 1200){
                cout << "event " << i << " chamber " << did << " X fired strips: " << NumX_f[did] << " Y fired strips: " << NumY_f[did] << endl;
                cout << "Fired strips number over maximum!!!" << endl;
                over = true;
            }
        }
        if(over) continue;
        Ientry = i;
        FiredStripsF->cd();
        FiredStripsT->Fill();

        outputFile->cd();
        for(int did =0; did < 4; did++){
            SetTH2D(shape_pulse_x[did]);
            SetTH2D(shape_pulse_y[did]);
            shape_pulse_x[did]->Write();
            shape_pulse_y[did]->Write();
            delete shape_pulse_x[did];
            delete shape_pulse_y[did];
        }

    }
    TCanvas *c3 = new TCanvas();
    peak->GetXaxis()->SetTitle("peak ADC count");
    peak->GetYaxis()->SetTitle("counts");
    peak->Draw();
    peakY->SetLineColor(kRed);
    peakY->Draw("SAME");

    TCanvas *c4 = new TCanvas();
    adc_Num->GetXaxis()->SetTitle("max ADC value");
    adc_Num->GetYaxis()->SetTitle("fired strip number");
    adc_Num->GetZaxis()->SetTitle("counts");
    adc_Num->GetXaxis()->SetTitleOffset(1.3);
    adc_Num->GetYaxis()->SetTitleOffset(1.3);
    adc_Num->GetZaxis()->SetTitleOffset(1.3);
    adc_Num->GetXaxis()->CenterTitle();
    adc_Num->GetYaxis()->CenterTitle();
    adc_Num->Draw("colz");

    outputFile->Close();
    FiredStripsF->Write();
    FiredStripsF->Close();

    /*TFile *OccupancyFile = new TFile("Occupancy.root", "RECREATE");
    OccupancyFile->cd();
    for(int did = 0; did < 4; did++){
        Occupancy_X[did]->GetXaxis()->SetTitle("X strip channel");
        Occupancy_X[did]->GetYaxis()->SetTitle("Occupancy[%]");
        Occupancy_X[did]->GetXaxis()->SetTitleOffset(1.3);
        Occupancy_X[did]->GetYaxis()->SetTitleOffset(1.3);
        Occupancy_X[did]->GetXaxis()->CenterTitle();
        Occupancy_X[did]->GetYaxis()->CenterTitle();
        Occupancy_X[did]->Scale(100./validEntries);
        Occupancy_X[did]->Write();
        Occupancy_Y[did]->GetXaxis()->SetTitle("Y strip channel");
        Occupancy_Y[did]->GetYaxis()->SetTitle("Occupancy[%]");
        Occupancy_Y[did]->GetXaxis()->SetTitleOffset(1.3);
        Occupancy_Y[did]->GetYaxis()->SetTitleOffset(1.3);
        Occupancy_Y[did]->GetXaxis()->CenterTitle();
        Occupancy_Y[did]->GetYaxis()->CenterTitle();
        Occupancy_Y[did]->Scale(100./validEntries);
        Occupancy_Y[did]->Write();
    }
    OccupancyFile->Close();*/

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "运行时间: " << diff.count() << " 秒" << std::endl;
}
