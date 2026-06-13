#include <iostream>
#include "TROOT.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TRandom.h"

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
    gem01 -> AddPlane(1, "GEM01Y", 1228.8, -1); //
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
    gem11 -> AddPlane(1, "GEM11Y", 1228.8, -1);
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

    //TString RootFile = "/home/liyuan/PRad2sim/gem_digitization/FiredStrips_4mm.root";//The pwd of root file
    //TString RootFile = "/home/liyuan/PRad2sim/gem_digitization/FiredStrips_onlySignal.root";//The pwd of root file
    //TString RootFile = "/home/liyuan/PRad2sim/gem_digitization/mass/FiredStrips_4mm_x17mass_onlySignal.root";//The pwd of root file
    TString RootFile = "/home/liyuan/PRad2sim/gem_digitization/mass/FiredStrips_4mm_x17mass_200bg.root";//The pwd of root file
    auto File = new TFile(RootFile,"READ");//只读, RECREATE to write
    TTree* tree = (TTree*)File->Get("T");//Tree name
    int NumX_f[4], NumY_f[4], Ientry;
    float StripX_f[4][800][6+3], StripY_f[4][1200][6+3];
    tree->SetBranchAddress("Ientry", &Ientry);
    tree->SetBranchAddress("NumX", &NumX_f[0]);
    tree->SetBranchAddress("NumY", &NumY_f[0]);
    tree->SetBranchAddress("StripX", &StripX_f[0]);
    tree->SetBranchAddress("StripY", &StripY_f[0]);

    //RootFile = "/home/liyuan/PRad2sim/gem_digitization/data/Uniform_signal_660-1540MeV.root";
    RootFile = "/home/liyuan/PRad2sim/gem_digitization/mass/X17_Ap_dat/ap_Eb_2200_17_4um.root";
    auto File2 = new TFile(RootFile,"READ");//只读, RECREATE to write
    TTree* tree2 = (TTree*)File2->Get("T");//Tree
    Int_t GasN, EventID;
    Double_t GasX[2000], GasY[2000];
    Double_t GasXout[2000], GasYout[2000];
    Int_t GasPID[2000], GasPTID[2000], GasTID[2000], GasDID[2000];
    Int_t VDN, VDPTID[2000], VDPID[2000];
    Double_t VDX[2000], VDY[2000], VDZ[2000], VDP[2000];
    tree2->SetBranchAddress("GEM_Gas.N", &GasN);
    tree2->SetBranchAddress("GEM_Gas.X", &GasX[0]);
    tree2->SetBranchAddress("GEM_Gas.Y", &GasY[0]);
    tree2->SetBranchAddress("GEM_Gas.Xout", &GasXout[0]);
    tree2->SetBranchAddress("GEM_Gas.Yout", &GasYout[0]);
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

    TH1D *deltaX = new TH1D("deltaX", "deltaX", 200*8, -0.4*8, 0.4*8);
    TH1D *deltaY = new TH1D("deltaY", "deltaY", 200*8, -0.4*8, 0.4*8);
    TH1D *deltaR = new TH1D("deltaR", "deltaR", 200*8, -0.4*8, 0.4*8);
    TH1D *clusterXSize = new TH1D("clusterXSize", "clusterXSize", 20, -0.5, 20-0.5);
    TH1D *clusterYSize = new TH1D("clusterYSize", "clusterYSize", 20, -0.5, 20-0.5);

    std::ofstream outfile;
    //outfile.open("Reconstructed_Clusters_4mm.dat");
    //outfile.open("Reconstructed_Clusters_onlySignal.dat");
    //outfile.open("Reconstructed_Clusters_4mm_x17mass_onlySignal.dat");
    outfile.open("Reconstructed_Clusters_4mm_x17mass_200bg.dat");

    // usage: fill simulated strip signals

    int count_r = 0, count_f = 0;
    
    int layer_index = -1, detector_index = -1, plane_index = -1, strip_index = -1;
    int Ndata = tree->GetEntries();
    std::cout << "Number of entries in the data: " << Ndata << std::endl;

    int count_hit_file = 0;
    int count_hit_recon = 0;

    for(int event_id = 0; event_id<Ndata; ++event_id) 
    {
        //std::cout<<"event id: "<<event_id<<std::endl;
        // for each event, reset each layer to clear contents from previous event
        layer0 -> Reset();
        layer1 -> Reset();

        // fill all fired strips in simulation
        tree->GetEntry(event_id);
        //std::cout<<"entry_event_id: "<<Ientry<<std::endl;

        tree2->GetEntry(Ientry);
        bool haveHit = false;
        for(int i = 0; i < VDN; i++){
            if(VDPTID[i] == 0 && VDPID[i] == 11){
                haveHit = true;
                break;
            }
        }
        if(!haveHit) {
            //continue;
        }

        outfile << Ientry << std::endl; // record event id
        for(int layer_index = 0; layer_index < 2; layer_index++){
            for(int detector_index = 0; detector_index < 2; detector_index++){
                GEMPlane *plnX = det[layer_index][detector_index] -> GetPlane(0); // 0 for X strips, 1 for Y strips
                GEMPlane *plnY = det[layer_index][detector_index] -> GetPlane(1);
                int did = layer_index * 2 + detector_index;
                //std::cout << "layer " << layer_index << " detector " << detector_index << " DID " << did 
                //    << " NumX: " << NumX_f[did] << " NumY: " << NumY_f[did] << std::endl;
                for(int i=0; i<NumX_f[did]; i++){
                    plnX -> AddStripHit(int(StripX_f[did][i][0]), StripX_f[did][i][1], int(StripX_f[did][i][2]+1),
                        {StripX_f[did][i][3], StripX_f[did][i][4], StripX_f[did][i][5], StripX_f[did][i][6], StripX_f[did][i][7], StripX_f[did][i][8]});
                    //std::cout << int(StripX_f[did][i][0]) << " " << StripX_f[did][i][1] << " " << int(StripX_f[did][i][2]+1) << std::endl;
                }
                for(int i=0; i<NumY_f[did]; i++){
                    plnY -> AddStripHit(int(StripY_f[did][i][0]), StripY_f[did][i][1], int(StripY_f[did][i][2]+1),
                        {StripY_f[did][i][3], StripY_f[did][i][4], StripY_f[did][i][5], StripY_f[did][i][6], StripY_f[did][i][7], StripY_f[did][i][8]});
                    //std::cout << int(StripY_f[did][i][0]) << " " << StripY_f[did][i][1] << " " << int(StripY_f[did][i][2]+1) << std::endl;
                }
            }
        }
        // suppose you have strip# 543 and strip#544 fired in simulation with the following 6 time sample ADCs on this plane
            //std::vector<float> time_sample_adc{50.5, 150.4, 375.0, 302.4, 201.0, 104.0};
            //pln -> AddStripHit(543, 375.0, 3, {50.5, 150.4, 375.0, 302.4, 201.0, 104.0}); // strip_index, strip_charge, max_time_bin, time_sample_adc


        // after done with fired strips for this event, reconstruct cluster positions
        for(int ii_layer=0; ii_layer<2; ii_layer++){
            for(int ii_detector=0; ii_detector<2; ii_detector++){
                // do cluster reconstruction for each detector
                det[ii_layer][ii_detector] -> Reconstruct(clustering_method);
                //std::cout<<"Reconstructed clusters for layer "<<ii_layer<<", detector "<<ii_detector<<std::endl;

                // after reconstruction, you should be able to get reconstructed clusters
                const std::vector<StripCluster> &clusters_x = det[ii_layer][ii_detector] -> GetPlane(0) -> GetStripClusters();
                const std::vector<StripCluster> &clusters_y = det[ii_layer][ii_detector] -> GetPlane(1) -> GetStripClusters();

                float recon_x, recon_y;
                float dx = 0., dy = 0.;
                double hitX = 0., hitY = 0.;
                bool hit_found = false;
                for(int n = 0; n < GasN; n++){
                    if(GasPID[n] == 11 && GasTID[n] == 1 && GasDID[n] == ii_layer*2+ii_detector){
                        hitX = (GasX[n] + GasXout[n]) / 2.;
                        hitY = (GasY[n] + GasYout[n]) / 2.;
                        hit_found = true;
                        count_hit_file++;
                        break;
                    }
                }
                if(hit_found){
                    bool recon_hit_x = false, recon_hit_y = false;
                    for(auto &c: clusters_x) {
                        recon_x = c.position + gRandom->Gaus(0., 0.045); // 0.05
                        if(fabs(recon_x - hitX) < 0.4 * 10.) { // old: 0.07*5.
                            recon_hit_x = true;
                            //recon_x = c.position + gRandom->Gaus(0., 0.04); // 0.05
                            dx = recon_x - hitX;
                            clusterXSize->Fill(c.hits.size());
                            break;
                        }
                        //std::cout<<"layer: "<<ii_layer<<", detector: "<<ii_detector<<", planeX"<<", position: "<<c.position<<", charge: "<<c.total_charge<<std::endl;
                    }

                    for(auto &c: clusters_y) {
                        recon_y = c.position + gRandom->Gaus(0., 0.045); // 0.046
                        if(fabs(recon_y - hitY) < 0.4 * 10.) {
                            recon_hit_y = true;
                            //recon_y = c.position + gRandom->Gaus(0., 0.035); // 0.046
                            dy = recon_y - hitY;
                            clusterYSize->Fill(c.hits.size());
                            break;
                        }
                        //std::cout<<"layer: "<<ii_layer<<", detector: "<<ii_detector<<", planeY"<<", position: "<<c.position<<", charge: "<<c.total_charge<<std::endl;
                    }
                    if(recon_hit_x && recon_hit_y) {
                        count_hit_recon++;
                        deltaX->Fill(dx);
                        deltaY->Fill(dy);
                        deltaR->Fill(sqrt(recon_x*recon_x + recon_y*recon_y) - sqrt(hitX*hitX + hitY*hitY));
                    }
                }
                int did = ii_layer * 2 + ii_detector;
                outfile << did << " " << clusters_x.size() << " " << clusters_y.size() << std::endl; // record detector id, number of clusters on X and Y planes
                for(auto &c: clusters_x){
                    if(std::isnan(c.position)){
                        std::cerr << "entry_event_id: "<< Ientry << ", did = " << did << std::endl;
                        std::cerr << "Invalid value detected: c.position number is NaN. " << c.position << std::endl;
                        for(int i=0; i<NumX_f[did]; i++){
                            std::cerr << int(StripX_f[did][i][0]) << " " << StripX_f[did][i][1] << " " << int(StripX_f[did][i][2]+1) << ", ";
                        }
                        std::cerr << std::endl;
                        continue;
                    }
                    outfile << c.position + gRandom->Gaus(0., 0.045) << " ";
                }
                outfile << std::endl; // end of line for X clusters
                for(auto &c: clusters_y){
                    outfile << c.position + gRandom->Gaus(0., 0.045) << " ";
                }
                outfile << std::endl; // end of line for Y clusters
            }
        }
    }
    outfile << -1 << std::endl; // end of all events
    outfile.close();

    gStyle->SetOptFit(1111);
    TCanvas *c1 = new TCanvas();
    //c1->SetLogy();
    deltaX->Draw();
    deltaX->SetTitle("#deltaX = Reconstructed position - MC position");
    deltaX->GetXaxis()->SetTitle("#deltax [mm]");
    deltaX->GetYaxis()->SetTitle("Counts");
    deltaX->GetXaxis()->CenterTitle();
    deltaX->Fit("gaus","Q","",-0.2,0.2);
    deltaX->GetXaxis()->SetRangeUser(-0.4, 0.4);
    c1->SaveAs("deltaX.png");
    c1->SetLogy();
    deltaX->GetXaxis()->SetRangeUser(-3., 3.);
    c1->Update();
    c1->SaveAs("deltaX_log.png");
    TCanvas *c2 = new TCanvas();
    //c2->SetLogy();
    deltaY->Draw();
    deltaY->SetTitle("#deltaY = Reconstructed position - MC position");
    deltaY->GetXaxis()->SetTitle("#deltay [mm]");
    deltaY->GetYaxis()->SetTitle("Counts");
    deltaY->GetXaxis()->CenterTitle();
    deltaY->Fit("gaus","Q","",-0.2,0.2);
    deltaY->GetXaxis()->SetRangeUser(-0.4, 0.4);
    c2->SaveAs("deltaY.png");
    c2->SetLogy();
    deltaY->GetXaxis()->SetRangeUser(-3., 3.);
    c2->Update();
    c2->SaveAs("deltaY_log.png");
    TCanvas *c3 = new TCanvas();
    c3->SetLogy();
    deltaR->Draw();
    deltaR->SetTitle("#deltaR = Reconstructed R - MC R");
    deltaR->GetXaxis()->SetTitle("#deltaR [mm]");
    deltaR->GetYaxis()->SetTitle("Counts");
    deltaR->GetXaxis()->CenterTitle();
    deltaR->Fit("gaus","Q","",-0.2,0.2);
    c3->SaveAs("deltaR.png");

    TCanvas *c4 = new TCanvas();
    c4->SetGrid();
    clusterXSize->SetTitle("Cluster Size of signal hits");
    clusterXSize->GetXaxis()->SetTitle("Cluster Size [c.hits.size()]");
    clusterXSize->GetYaxis()->SetTitle("Counts");
    clusterXSize->GetXaxis()->CenterTitle();
    clusterXSize->SetLineColor(1);
    clusterXSize->Draw();
    clusterXSize->GetXaxis()->SetNdivisions(15);
    clusterYSize->SetLineColor(2);
    clusterYSize->Draw("same");
    TLegend *leg = new TLegend(0.6,0.4,0.88,0.6);
    leg->AddEntry(clusterXSize, "X clusters", "l");
    leg->AddEntry(clusterYSize, "Y clusters", "l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c4->SaveAs("clusterSize.png");

    std::cout << "Total number of events with hits in the file: " << count_hit_file << std::endl;
    std::cout << "Total number of events with hits reconstructed: " << count_hit_recon << std::endl;
    std::cout << "Reconstruction efficiency: " << float(count_hit_recon)/float(count_hit_file) << std::endl;

    return 0;
}
