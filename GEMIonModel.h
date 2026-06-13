
#ifndef GEM_ION_MODEL_H
#define GEM_ION_MODEL_H

#include <cmath>

#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"

#include "GEMdigVariables.h"

inline void IonModel(TVector3& in, TVector3& out, const Double_t Edep) {
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

#endif
