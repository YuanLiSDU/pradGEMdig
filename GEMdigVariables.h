
#include <vector>

#include "Rtypes.h"
#include "TH2D.h"

// Database: gas
Double_t GasIonE = 26.0; // eV, ionization energy of gas
Double_t GasDriftVelocity = 55.; // from db: mm/us, drift velocity of gas
Double_t GasLateralUncertainty = 0.; // mm, lateral uncertainty of avalanche center
Double_t gain_mean = 41000.; // from db: mean gain of the GEM, 4500
Double_t GasDiffusion = 0.1; // mm2/us, gas diffusion coefficient after avalanche
Double_t fSNormNsigma = 18.0;
Double_t fAvaGain = 54.;

// GEM
const Int_t activeStripNhalf = 6;
Double_t Zreadout = 2.0 + 2.0 + 2.0 + 180. * 1.e-3; // mm, distance from the first GEM foil Cu
const Double_t stripWidth = 0.4; // mm
Double_t pulseTimeLength = 56.; // ns, 50, 56
Double_t map_x_low = -20.; // mm
Double_t map_x_high = 528.4; // mm
// const Int_t map_x_bins = int((map_x_high - map_x_low) / stripWidth); // 1371
const Int_t map_x_bins = 1371;
Double_t map_y_low = -614.4; // mm
Double_t map_y_high = 614.4; // mm
// const Int_t map_y_bins = int((map_y_high - map_y_low) / stripWidth); // 3072
const Int_t map_y_bins = 3072;
Double_t hole_length = 74.; // mm, length of the hole
TH2D* stripMap = new TH2D("strip map DID=0", "strip map DID=0",
                          map_x_bins, map_x_low, map_x_high,
                          map_y_bins, map_y_low, map_y_high); // mm, strip width 0.4 mm

// Signal
const Int_t samplingNum = 6;
const Double_t samplingTime = 25.; // ns
const Double_t fTriggerJitter = 10.; // ns
const Double_t fAPVTimeJitter = 25.; // ns

// Variables
struct IonPar_t {
    Double_t Xp; // position of the point on the projection
    Double_t Yp;
    Double_t Charge; // charge deposited by this ion
    Double_t S_sigma; // sigma of spatial radius of ion diffusion area
    Double_t Ttime; // traveling time of the ion
    Double_t ggnorm;
    Double_t eff_sigma;
};

std::vector<IonPar_t> fRIon;
Double_t Num_ions; // number of ions
Double_t minTravelT = 1000000.;
Int_t lowX, lowY, upX, upY;

// Parameters in ROOT file
Int_t GasN, EventID;
Double_t GasX[2000], GasY[2000], GasZ[2000], GasTime[2000];
Double_t GasXout[2000], GasYout[2000], GasZout[2000];
Double_t GasTheta[2000], GasP[2000], GasEdep[2000];
Int_t GasPID[2000], GasPTID[2000], GasTID[2000], GasDID[2000];
