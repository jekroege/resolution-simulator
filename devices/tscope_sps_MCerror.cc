// Simon Spannagel (DESY) January 2016
// Modified by Jens Kroeger (CERN) March 2021

#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TRandom.h"
#include "TRandom3.h"

#include "assembly.h"
#include "propagate.h"
#include "materials.h"
#include "constants.h"
#include "log.h"

#include "/home/jens/CLICdpStyle.C"

using namespace std;
using namespace gblsim;
using namespace unilog;

int main(int argc, char* argv[]) {

    CLICdpStyle();
    gStyle->SetOptFit(1111);

    /*
    * Telescope resolution simulation for the CLICdp Timepix3 telescope at the SPS H6B beam line
    * Seven Timepix3 planes with different spacing, intrinsic sensor resolution 4.0um
    * DUT with variable thickness (scan)
    */


    Log::ReportingLevel() = Log::FromString("INFO");

    int mode;
    if(argc == 1) {
        std::cout << "Please choose a mode!" << std::endl;
        return 0;
    }

    for (int i = 1; i < argc; i++) {
        // Setting verbosity:
        if (std::string(argv[i]) == "-v") {
            Log::ReportingLevel() = Log::FromString(std::string(argv[++i]));
            continue;
        } else {
            mode = atoi(argv[i]);
            if(mode == 0) {
                std::cout << "Please choose your mode:" << std::endl;
                std::cout << "\t1: 7 Timepix3 planes, DUT = APX" << std::endl;
                std::cout << "\t2: 6 Timepix3 planes, DUT = APX" << std::endl;
                std::cout << "\t3: 7 Timepix3 planes, DUT = CP2" << std::endl;
                std::cout << "\t4: 6 Timepix3 planes, DUT = CP2" << std::endl;
                return 0;
            }
            std::cout << "You chose mode = " << mode << std::endl;
        }
    }

    TFile * out;
    if(mode == 1){
        out = TFile::Open("output/sps-resolution-nov2018_apx_7planes.root","RECREATE");
    } else if(mode == 2){
        out = TFile::Open("output/sps-resolution-nov2018_apx_6planes.root","RECREATE");
    } else if(mode == 3){
        out = TFile::Open("output/sps-resolution-nov2018_cp2_7planes.root","RECREATE");
    } else if(mode == 4){
        out = TFile::Open("output/sps-resolution-nov2018_cp2_6planes.root","RECREATE");
    } else {
        std::cout << "Invalid mode...try again..." << std::endl;
    }
    gDirectory->pwd();

    TCanvas *c1 = new TCanvas();
    TH1F* hResolution = new TH1F("hResolution","hResolution",1000,0,10);
    hResolution->GetXaxis()->SetTitle("resolution at DUT [#mum]");
    hResolution->GetYaxis()->SetTitle("# entries");

    //----------------------------------------------------------------------------
    // Preparation of the telescope and beam properties:

    // Timepix3 telescope planes have X/X_0
    // (see PhD thesis Niloufar Tehrani)
    double X_TPX3 = 4.0e-2; // X/X_0
    double ERR_X_TPX3 = 0.5e-2; // same as used by Morag

    // The intrinsic resolution has been measured to be around 4.0um
    // (also see PhD thesis Niloufar Tehrani):
    double RES = 4e-3; // in mm
    double ERR_RES = 0.2e-3; // in mm

    // ATLASpix_Simple radiation length (see PhD thesis Jens Kroeger)
    double X_DUT;
    double ERR_X_DUT;
    if(mode == 1 || mode == 2) {
        // X_DUT = 0.985e-2; // thickness = 62um (my calculation)
        X_DUT = 1.025e-2; // thickness = 100um (my calculation)
        ERR_X_DUT = 0.5e-2;
    }
    if(mode == 3 || mode == 4 ) {
        X_DUT = 2.4e-2; // CLICpix2 -> See PhD thesis Morag Williams
        ERR_X_DUT = 0.5e-2;
    }

    //      D04  E04  G02      DUT      G03  J05  L09  F09
    // beam  |    |    |        |        |    |    |    |
    // --->  |    |    |        |        |    |    |    |
    //       |    |    |        |        |    |    |    |
    //      0.0                                        336.5 mm

    // Positions of telescope planes in mm:
    std::vector<double> Z_TEL;
    Z_TEL.emplace_back(0.0);
    Z_TEL.emplace_back(21.5);
    Z_TEL.emplace_back(43.5);
    Z_TEL.emplace_back(186.5);
    Z_TEL.emplace_back(208.0);
    Z_TEL.emplace_back(231.5);
    if(mode == 1 || mode == 3) {
        Z_TEL.emplace_back(336.5);
    }
    // Position of the DUT in mm:
    double Z_DUT = 105.0;
    // double ERR_Z = 0.5;
    double ERR_Z = 1.0;

    // Beam energy 120 GeV pions at SPS:
    double EBEAM = 120.0;


    //----------------------------------------------------------------------------
    // Build the trajectory through the telescope device:

    TRandom3* rand = new TRandom3();
    // rand->SetSeed(33333);
    for(int it=0; it<1e4; it++) {

        // Prepare the DUT (no measurement, just scatterer
        plane dut(rand->Gaus(Z_DUT,ERR_Z),
                  rand->Gaus(X_DUT,ERR_X_DUT),
                  false);

        // Telescope setup:
        std::vector<plane> tpx3_tel;
        // Build a vector of all telescope planes:
        for(int i = 0; i < Z_TEL.size(); i++) {
            tpx3_tel.emplace_back(plane(rand->Gaus(Z_TEL.at(i), ERR_Z),
                                        rand->Gaus(X_TPX3,ERR_X_TPX3),
                                        true,
                                        rand->Gaus(RES,ERR_RES)));
        }

        // Duplicate the planes vector and add the current DUT:
        std::vector<plane> planes = tpx3_tel;
        planes.emplace_back(dut);

        // Build the telescope:
        telescope mytel(planes, EBEAM);

        // Get the resolution at plane-vector position (x):
        LOG(logRESULT) << "Track resolution at DUT in iteration it = " << it << ": " << mytel.getResolution(3) << "% X0";
        hResolution->Fill(mytel.getResolution(3));
    }

    LOG(logRESULT) << "Histgram has " << hResolution->GetEntries() << " entries.";
    c1->cd();
    hResolution->Draw();
    double mean = hResolution->GetMean();
    double sigma = hResolution->GetRMS();
    std::cout << "Set range to " << mean-5*sigma << ", " << mean+5*sigma << std::endl;
    hResolution->GetXaxis()->SetRangeUser(mean-5*sigma,mean+5*sigma);

    hResolution->Fit("gaus");
    TF1* func = hResolution->GetFunction("gaus");
    double fMean = func->GetParameter(1);
    double fSigma = func->GetParameter(2);

    std::cout << std::setprecision(4) << "Track poiinting resolution (mean+/-sigma): " << fMean << "+/-" << fSigma << "um" << std::endl;


    if(mode == 1) {
        hResolution->SetTitle("APX: Res. at DUT (7 planes)");
    }
    if(mode == 2) {
        hResolution->SetTitle("APX: Res. at DUT (6 planes)");
    }
    if(mode == 3) {
        hResolution->SetTitle("CP2: Res. at DUT (7 planes)");
    }
    if(mode == 4) {
        hResolution->SetTitle("CP2: Res. at DUT (6 planes)");
    }

    c1->Write();

    // Write result to file
    out->Write();
    out->Close();
    return 0;
}
