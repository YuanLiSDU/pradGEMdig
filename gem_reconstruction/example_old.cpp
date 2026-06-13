#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TStyle.h"

#include "GEMDetectorLayer.h"
#include "GEMDetector.h"
#include "GEMPlane.h"
#include "GEMCluster.h"

int main(int argc, char* argv[])
{
    ////////////////////////////////////////////////////////////////////////////
    // GEM Detector setup

    // declare first GEM layer 
    GEMDetectorLayer *layer0 = new GEMDetectorLayer(0); // layer id: 0
    // declare second GEM layer
    GEMDetectorLayer *layer1 = new GEMDetectorLayer(1); // layer id: 1

    // declare first GEM detector in first layer 
    GEMDetector *gem00 = new GEMDetector(layer0, "GEM00", 0, 0, 0); // layer, detector_name, detector_unique_id, layer_id, position_index
    gem00 -> AddPlane(0, "GEM00X", 550, 1); // axis(x_axis = 0, y_axis=1), plane_unique_name, plane_size, plane_direction
    gem00 -> AddPlane(1, "GEM00Y", 1228.8, 1);
    gem00 -> SetOffset(556.2-20.+0.2, 0); // set detector offset after adding planes
    layer0 -> AddGEMDetector(gem00);
    // declear second GEM detector in first layer
    GEMDetector *gem01 = new GEMDetector(layer0, "GEM01", 1, 0, 1);
    gem01 -> AddPlane(0, "GEM01X", 550, -1); // gem1 detector plane will be flipped
    gem01 -> AddPlane(1, "GEM01Y", 1228.8, 1); //
    gem01 -> SetOffset(-556.2+20.-0.2, 0);
    layer0 -> AddGEMDetector(gem01);

    // declare first GEM detector in second layer
    GEMDetector *gem10 = new GEMDetector(layer1, "GEM10", 2, 1, 0);
    gem10 -> AddPlane(0, "GEM10X", 550, 1);
    gem10 -> AddPlane(1, "GEM10Y", 1228.8, 1);
    gem10 -> SetOffset(556.2-20.+0.2, 0);
    layer1 -> AddGEMDetector(gem10);
    // declare second GEM detector in second layer
    GEMDetector *gem11 = new GEMDetector(layer1, "GEM11", 3, 1, 1);
    gem11 -> AddPlane(0, "GEM11X", 550, -1);
    gem11 -> AddPlane(1, "GEM11Y", 1228.8, 1);
    gem11 -> SetOffset(-556.2+20.-0.2, 0);
    layer1 -> AddGEMDetector(gem11);

    layer0 -> PrintStatus();
    layer1 -> PrintStatus();

    ////////////////////////////////////////////////////////////////////////////
    // initialize GEM clustering methods and cuts for filtering clusters
    GEMCluster *clustering_method = new GEMCluster("config/gem.conf");


    ////////////////////////////////////////////////////////////////////////////
    // usage: example

    // a pointer array for looping through PRad-II GEM detector setup
    GEMDetector *det[2][2]; // det[layer_index][detector_index]
    det[0][0] = gem00;
    det[0][1] = gem01;
    det[1][0] = gem10;
    det[1][1] = gem11;

    ////////////////////////////////////////////////////////////////////////////

    TString RootFile = "/home/liyuan/PRad2sim/gem_digitization/FiredStrips_onlySignal.root";//The pwd of root file
    auto File = new TFile(RootFile,"READ");//只读, RECREATE to write
    TTree* tree = (TTree*)File->Get("T");//Tree name
    int NumX_f[4], NumY_f[4], Ientry;
    float StripX_f[4][800][6+3], StripY_f[4][1200][6+3];
    tree->SetBranchAddress("Ientry", &Ientry);
    tree->SetBranchAddress("NumX", &NumX_f);
    tree->SetBranchAddress("NumY", &NumY_f);
    tree->SetBranchAddress("StripX", &StripX_f);
    tree->SetBranchAddress("StripY", &StripY_f);

    RootFile = "/home/liyuan/PRad2sim/gem_digitization/data/Uniform_signal_660-1540MeV2.root";
    auto File2 = new TFile(RootFile,"READ");//只读, RECREATE to write
    TTree* tree2 = (TTree*)File2->Get("T");//Tree
    Int_t GasN, EventID;
    Double_t GasX[2000], GasY[2000], GasZ[2000];
    Double_t GasXout[2000], GasYout[2000], GasZout[2000];
    Int_t GasPID[2000], GasPTID[2000], GasTID[2000], GasDID[2000];
    Int_t VDN, VDPTID[2000], VDPID[2000];
    Double_t VDX[2000], VDY[2000], VDZ[2000], VDP[2000];
    tree2->SetBranchAddress("GEM_Gas.N", &GasN);
    tree2->SetBranchAddress("GEM_Gas.X", &GasX[0]);
    tree2->SetBranchAddress("GEM_Gas.Y", &GasY[0]);
    tree2->SetBranchAddress("GEM_Gas.Z", &GasZ[0]);
    tree2->SetBranchAddress("GEM_Gas.Xout", &GasXout[0]);
    tree2->SetBranchAddress("GEM_Gas.Yout", &GasYout[0]);
    tree2->SetBranchAddress("GEM_Gas.Zout", &GasZout[0]);
    tree2->SetBranchAddress("GEM_Gas.PID", &GasPID[0]);
    tree2->SetBranchAddress("GEM_Gas.PTID", &GasPTID[0]);
    tree2->SetBranchAddress("GEM_Gas.TID", &GasTID[0]);
    tree2->SetBranchAddress("GEM_Gas.DID", &GasDID[0]);
    tree2->SetBranchAddress("VD.N", &VDN);
    tree2->SetBranchAddress("VD.X", &VDX[0]);
    tree2->SetBranchAddress("VD.Y", &VDY[0]);
    tree2->SetBranchAddress("VD.Z", &VDZ[0]);
    tree2->SetBranchAddress("VD.P", &VDP[0]);
    tree2->SetBranchAddress("VD.PTID", &VDPTID[0]);
    tree2->SetBranchAddress("VD.PID", &VDPID[0]);

    TH1D *deltaX = new TH1D("deltaX", "deltaX", 200, -0.3, 0.3);
    TH1D *deltaY = new TH1D("deltaY", "deltaY", 200, -0.3, 0.3);
    TH1D *deltaR = new TH1D("deltaR", "deltaR", 200, 0, 0.3);
    TH2D *theta_deltaY = new TH2D("theta_deltaY", "theta_deltaY", 200, -0.3, 0.3, 50, 0.3, 5.);


    // usage: fill simulated strip signals

    int count_r = 0, count_f = 0;
    
    int layer_index = -1, detector_index = -1, plane_index = -1, strip_index = -1;
    int Ndata = tree->GetEntries();
    std::cout << "Number of entries in the data: " << Ndata << std::endl;
    for(int event_id = 0; event_id<Ndata; ++event_id) 
    {
        std::cout<<"event id: "<<event_id<<std::endl;
        // for each event, reset each layer to clear contents from previous event
        layer0 -> Reset();
        layer1 -> Reset();

        // fill all fired strips in simulation
        tree->GetEntry(event_id);
        std::cout<<"entry_event_id: "<<Ientry<<std::endl;
        for(int layer_index = 0; layer_index < 2; layer_index++){
            for(int detector_index = 0; detector_index < 2; detector_index++){
                GEMPlane *plnX = det[layer_index][detector_index] -> GetPlane(0); // 0 for X strips, 1 for Y strips
                GEMPlane *plnY = det[layer_index][detector_index] -> GetPlane(1);
                int did = layer_index * 2 + detector_index;
                for(int i=0; i<NumX_f[did]; i++){
                    plnX -> AddStripHit(int(StripX_f[did][i][0]), StripX_f[did][i][1], int(StripX_f[did][i][2]),
                        {StripX_f[did][i][3], StripX_f[did][i][4], StripX_f[did][i][5], StripX_f[did][i][6], StripX_f[did][i][7], StripX_f[did][i][8]});
                }
                for(int i=0; i<NumY_f[did]; i++){
                    plnY -> AddStripHit(int(StripY_f[did][i][0]), StripY_f[did][i][1], int(StripY_f[did][i][2]),
                        {StripY_f[did][i][3], StripY_f[did][i][4], StripY_f[did][i][5], StripY_f[did][i][6], StripY_f[did][i][7], StripY_f[did][i][8]});
                }
            }
        }
        // suppose you have strip# 543 and strip#544 fired in simulation with the following 6 time sample ADCs on this plane
            //std::vector<float> time_sample_adc{50.5, 150.4, 375.0, 302.4, 201.0, 104.0};
            //pln -> AddStripHit(543, 375.0, 3, {50.5, 150.4, 375.0, 302.4, 201.0, 104.0}); // strip_index, strip_charge, max_time_bin, time_sample_adc


        // after done with fired strips for this event, reconstruct cluster positions
        tree2->GetEntry(Ientry);

        for(int ii_layer=0; ii_layer<2; ii_layer++)
        {
            for(int ii_detector=0; ii_detector<2; ii_detector++)
            {
                tree2->GetEntry(Ientry);
                bool haveHit = false;
                double E = 0.;
                for(int i = 0; i < VDN; i++){
                    if(VDPTID[i] == 0 && VDPID[i] == 11){
                        haveHit = true;
                        E = VDP[i];
                        break;
                    }
                }
                if(!haveHit) {
                    continue;
                }

                // do cluster reconstruction for each detector
                det[ii_layer][ii_detector] -> Reconstruct(clustering_method);

                // after reconstruction, you should be able to get reconstructed clusters
                const std::vector<StripCluster> &clusters_x = det[ii_layer][ii_detector] -> GetPlane(0) -> GetStripClusters();
                const std::vector<StripCluster> &clusters_y = det[ii_layer][ii_detector] -> GetPlane(1) -> GetStripClusters();
                float dx = 0., dy = 0.;
                double hitX = 0., hitY = 0., theta = 0.;
                bool hit_found = false;
                for(int n = 0; n < GasN; n++){
                    if(GasPID[n] == 11 && GasPTID[n] == 0 && GasDID[n] == ii_layer*2+ii_detector){
                        hitX = (GasX[n] + GasXout[n]) / 2.;
                        hitY = (GasY[n] + GasYout[n]) / 2.;
                        double hitZ = (GasZ[n] + GasZout[n]) / 2.;
                        theta = atan(sqrt(hitX*hitX + hitY*hitY) / (hitZ + 3000. - 89.)) * 180. / TMath::Pi();
                        hit_found = true;
                        break;
                    }
                }
                if(hit_found){
                    bool recon_hit_x = false, recon_hit_y = false;
                    for(auto &c: clusters_x) {
                        std::cout<<"layer: "<<ii_layer<<", detector: "<<ii_detector<<", planeX"<<", position: "<<c.position<<", charge: "<<c.total_charge<<std::endl;
                        if(fabs(c.position - hitX) < 0.4*2.) {
                            recon_hit_x = true;
                            dx = c.position - hitX;
                            break;
                        }
                    }

                    for(auto &c: clusters_y) {
                        std::cout<<"layer: "<<ii_layer<<", detector: "<<ii_detector<<", planeY"<<", position: "<<c.position<<", charge: "<<c.total_charge<<std::endl;
                        if(fabs(c.position - hitY) < 0.4*2.) {
                            recon_hit_y = true;
                            dy = c.position - hitY;
                            break;
                        }
                    }
                    if(recon_hit_x && recon_hit_y) {
                        //std::cout<<"Reconstructed hit found in layer "<<ii_layer<<", detector "<<ii_detector<<std::endl;
                        deltaX->Fill(dx);
                        deltaY->Fill(dy);
                        deltaR->Fill(sqrt(dx*dx + dy*dy));
                        theta_deltaY->Fill(dy, theta);
                        //count_r++;
                    }
                    else continue;
                }
            }
        }
    }
    gStyle->SetOptFit(1111);
    TCanvas *c1 = new TCanvas();
    //c1->SetLogy();
    deltaX->Draw();
    deltaX->SetTitle("#deltaX = Reconstructed position - MC position");
    deltaX->GetXaxis()->SetTitle("#deltax [mm]");
    deltaX->GetYaxis()->SetTitle("Counts");
    deltaX->GetXaxis()->CenterTitle();
    deltaX->Fit("gaus","Q","",-0.15,0.15);
    c1->SaveAs("deltaX.png");
    c1->SetLogy();
    c1->Update();
    c1->SaveAs("deltaX_log.png");
    TCanvas *c2 = new TCanvas();
    //c2->SetLogy();
    deltaY->Draw();
    deltaY->SetTitle("#deltaY = Reconstructed position - MC position");
    deltaY->GetXaxis()->SetTitle("#deltay [mm]");
    deltaY->GetYaxis()->SetTitle("Counts");
    deltaY->GetXaxis()->CenterTitle();
    deltaY->Fit("gaus","Q","",-0.15,0.15);
    c2->SaveAs("deltaY.png");
    c2->SetLogy();
    c2->Update();
    c2->SaveAs("deltaY_log.png");
    TCanvas *c3 = new TCanvas();
    c3->SetLogy();
    deltaR->Draw();
    deltaR->SetTitle("#deltaR = sqrt(#deltaX^{2} + #deltaY^{2})");
    deltaR->GetXaxis()->SetTitle("#deltar [mm]");
    deltaR->GetYaxis()->SetTitle("Counts");
    deltaR->GetXaxis()->CenterTitle();
    c3->SaveAs("deltaR.png");

    TCanvas *c4 = new TCanvas();
    theta_deltaY->Draw("COLZ");
    theta_deltaY->SetTitle("#deltaY vs #theta");
    theta_deltaY->GetXaxis()->SetTitle("#deltay [mm]");
    theta_deltaY->GetYaxis()->SetTitle("#theta [deg]");
    theta_deltaY->GetXaxis()->CenterTitle();
    theta_deltaY->GetYaxis()->CenterTitle();
    c4->SetLogz();
    c4->SaveAs("theta_deltaY.png");


    std::cout << "Total number of reconstructed hits: " << count_r + count_f << std::endl;
    std::cout << "Probability of right: " << double(count_r) / double(count_f + count_r) << std::endl;
    std::cout << "Probability of false: " << double(count_f) / double(count_f + count_r) << std::endl;

    return 0;
}
