// Simon Spannagel (DESY) January 2016
// Modified by Jens Kroeger (CERN) March 2021

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TString.h"
#include "TFile.h"

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

    /*
    * Telescope resolution simulation for the CLICdp Timepix3 telescope at the SPS H6B beam line
    * Seven Timepix3 planes with different spacing, intrinsic sensor resolution 4.0um
    * DUT with variable thickness (scan)
    */


    Log::ReportingLevel() = Log::FromString("INFO");

    int mode;
    for (int i = 1; i < argc; i++) {
        // Setting verbosity:
        if (std::string(argv[i]) == "-v") {
            Log::ReportingLevel() = Log::FromString(std::string(argv[++i]));
            continue;
        } else {
            mode = atoi(argv[i]);
            if(mode == 0) {
                std::cout << "Please choose your mode:" << std::endl;
                std::cout << "\t0: 7 Timepix3 planes, DUT = APX" << std::endl;
                std::cout << "\t1: 6 Timepix3 planes, DUT = APX" << std::endl;
                std::cout << "\t0: 7 Timepix3 planes, DUT = CP2" << std::endl;
                std::cout << "\t1: 6 Timepix3 planes, DUT = CP2" << std::endl;
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
    TGraphErrors *resolution = new TGraphErrors("resolution","resolution");

    //----------------------------------------------------------------------------
    // Preparation of the telescope and beam properties:

    // Timepix3 telescope planes have X/X_0
    // (see PhD thesis Niloufar Tehrani)
    double X_TPX3 = 4.0e-2; // X/X_0
    // The intrinsic resolution has been measured to be around 4.0um
    // (also see PhD thesis Niloufar Tehrani):
    double RES = 4e-3; // in mm

    // ATLASpix_Simple radiation length (see PhD thesis Jens Kroeger)
    double X_DUT;
    if(mode == 1 || mode == 2) {
        X_DUT = 0.985e-2; // thickness = 62um
    }
    if(mode == 3 || mode == 4 ) {
        X_DUT = 2.4e-2; // CLICpix2 -> See PhD thesis Morag Williams
    }
    // double X_DUT = 1.025e-2; // thickness = 100um
    double vary_x = 0.1;

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

    // Beam energy 120 GeV pions at SPS:
    double EBEAM = 120.0;


    //----------------------------------------------------------------------------
    // Build the trajectory through the telescope device:

    // Build a vector of all telescope planes:
    std::vector<plane> tpx3_tel;

    // Telescope setup:
    for(int i = 0; i < Z_TEL.size(); i++) {
        tpx3_tel.emplace_back(plane(Z_TEL.at(i),X_TPX3,true,RES));
    }

    // Calculate for 6 telescope planes:
    int j=0;
    for(double dut_x0 = X_DUT*(1-vary_x); dut_x0 < X_DUT*(1+vary_x); dut_x0 += 0.0001) {

        // Prepare the DUT (no measurement, just scatterer
        plane dut(Z_DUT, dut_x0, false);

        // Duplicate the planes vector and add the current DUT:
        std::vector<plane> planes = tpx3_tel;
        planes.emplace_back(dut);

        // Build the telescope:
        telescope mytel(planes, EBEAM);

        // Get the resolution at plane-vector position (x):
        LOG(logRESULT) << "Track resolution at DUT with " << dut_x0 << "% X0: " << mytel.getResolution(3);
            resolution->SetPoint(j,dut_x0*100,mytel.getResolution(3));
        j++;
    }

    c1->cd();
    if(mode == 1) {
        resolution->SetTitle("APX: Res. at DUT (7 planes);DUT material budget X/X_{0} [%];resolution at DUT [#mum]");
    }
    if(mode == 2) {
        resolution->SetTitle("APX: Res. at DUT (6 planes);DUT material budget X/X_{0} [%];resolution at DUT [#mum]");
    }
    if(mode == 3) {
        resolution->SetTitle("CP2: Res. at DUT (7 planes);DUT material budget X/X_{0} [%];resolution at DUT [#mum]");
    }
    if(mode == 4) {
        resolution->SetTitle("CP2: Res. at DUT (6 planes);DUT material budget X/X_{0} [%];resolution at DUT [#mum]");
    }
    // resolution->GetYaxis()->SetRangeUser(1.65,1.85);
    // resolution->GetXaxis()->SetRangeUser(X_DUT*(1-vary_x)*100,X_DUT*(1+vary_x)*100);

    resolution->Draw();
    resolution->Write();
    c1->Write();

    // Write result to file
    out->Write();
    return 0;
}
