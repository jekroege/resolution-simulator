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
    * Telescope resolution simulation for the Mimosa26 telescopes at the DESY-II testbeam facility
    * Six Mimosa26 planes with different spacing, intrinsic sensor resolution 3.2um,
    * ATLASpix as DUT
    * Timepix3 as additional timing plane (downstream)
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
                std::cout << "\t1: (June 2019) 6 Mimosa26, DUT = APX" << std::endl;
                std::cout << "\t2: (June 2019) 6 Mimosa26 + Timepix3, DUT = APX" << std::endl;
                std::cout << "\t3: (June 2019) 6 Mimosa26, DUT = CP2" << std::endl;
                std::cout << "\t4: (June 2019) 6 Mimosa26 + Timepix3, DUT = CP2" << std::endl;
                std::cout << "\t5: (July 2019) 6 Mimosa26, DUT = APX" << std::endl;
                std::cout << "\t6: (July 2019) 6 Mimosa26 + Timepix3, DUT = APX" << std::endl;
                return 0;
            }
            std::cout << "You chose mode = " << mode << std::endl;
        }
    }

    TFile * out;
    if(mode == 1){
        out = TFile::Open("output/desy-resolution-june2019_apx_M26.root","RECREATE");
    } else if(mode == 2){
        out = TFile::Open("output/desy-resolution-june2019_apx_M26+TPX3.root","RECREATE");
    } else if(mode == 3){
        out = TFile::Open("output/desy-resolution-june2019_cp2_M26.root","RECREATE");
    } else if(mode == 4){
        out = TFile::Open("output/desy-resolution-june2019_cp2_M26+TPX3.root","RECREATE");
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

    // Mimosa26 telescope planes have X/X_0
    // (see doi:10.1140/epjti/s40485-016-0033-2)
    double X_M26 = 0.075e-2;
    double ERR_X_M26 = 0.01e-2;

    // Timepix3 telescope planes have X/X_0
    // (see PhD thesis Niloufar Tehrani)
    double X_TPX3 = 3.8e-2; // X/X_0
    double ERR_X_TPX3 = 0.5e-2; // same as used by Morag

    // The intrinsic resolution has been measured to be around 3.2um
    // (see doi:10.1140/epjti/s40485-016-0033-2)
    double RES_M26 = 3.2e-3; // in mm
    double ERR_RES_M26 = 0.1e-3; // in mm

    // The intrinsic resolution has been measured to be around 12.75um
    // (also see PhD thesis Niloufar Tehrani):
    double RES_TPX3 = 12.75e-3; // in mm
    double ERR_RES_TPX3 = 0.01e-3; // in mm

    // ATLASpix_Simple radiation length (see PhD thesis Jens Kroeger)
    double X_DUT;
    double ERR_X_DUT;
    if(mode == 1 || mode == 2) {
        // X_DUT = 0.985e-2; // thickness = 62um (my calculation)
        X_DUT = 1.025e-2; // thickness = 100um (my calculation)
        ERR_X_DUT = 0.5e-2;
        // ERR_X_DUT = 0.1e-3;
    }
    if(mode == 3 || mode == 4 ) {
        X_DUT = 2.4e-2; // CLICpix2 -> See PhD thesis Morag Williams
        ERR_X_DUT = 0.5e-2;
    }
    // double X_DUT = 1.025e-2; // thickness = 100um

    //      M26_0  M26_1  M26_2  DUT  M26_3  M26_4  M26_5  TPX3_0
    // beam  |      |       |     |     |     |      |       |
    // --->  |      |       |     |     |     |      |       |
    //       |      |       |     |     |     |      |       |
    //      0.0                                             666 mm

    // Positions of telescope planes in mm:
    std::vector<double> Z_TEL;
    double Z_TPX;

    // Position of the DUT in mm:
    double Z_DUT;
    if(mode == 1 || mode == 2 || mode == 3 || mode == 4) {
        // June 2019
        Z_TEL.emplace_back(0);
        Z_TEL.emplace_back(153);
        Z_TEL.emplace_back(305);
        Z_DUT = 333;
        Z_TEL.emplace_back(344);
        Z_TEL.emplace_back(456);
        Z_TEL.emplace_back(576);
        Z_TPX3 = 666;
    } else if (mode == 5 || mode == 6) {
        // July 2019
        Z_TEL.emplace_back(0);
        Z_TEL.emplace_back(153);
        Z_TEL.emplace_back(305);
        Z_DUT = 331;
        Z_TEL.emplace_back(345);
        Z_TEL.emplace_back(455);
        Z_TEL.emplace_back(565);
        Z_TPX3 = 629;
    }
    double ERR_Z = 1.0;


    // Beam energy 5.4 GeV pions at DESY-II:
    double EBEAM = 5.42;
    double ERR_EBEAM = EBEAM*0.02; // 2% as specified in DESY-II paper


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
        std::vector<plane> m26_tel;
        // Build a vector of all telescope planes:
        for(int i = 0; i < Z_TEL.size(); i++) {
            m26_tel.emplace_back(plane(rand->Gaus(Z_TEL.at(i), ERR_Z),
                                        rand->Gaus(X_M26,ERR_X_M26),
                                        true,
                                        rand->Gaus(RES_M26,ERR_RES_M26)));
        }
        if(mode == 1 || mode == 3) {
            m26_tel.emplace_back(plane(rand->Gaus(Z_TPX3, ERR_Z),
                                        rand->Gaus(X_TPX3,ERR_X_TPX3),
                                        true,
                                        rand->Gaus(RES_TPX3,ERR_RES_TPX3)));
        }


        // Duplicate the planes vector and add the current DUT:
        std::vector<plane> planes = m26_tel;
        planes.emplace_back(dut);

        // Build the telescope:
        telescope mytel(planes, rand->Gaus(EBEAM,ERR_EBEAM));

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

    std::cout << std::setprecision(4) << "Track pointing resolution (mean+/-sigma): " << fMean << "+/-" << fSigma << "um" << std::endl;


    if(mode == 1) {
        hResolution->SetTitle("APX: Res. at DUT (M26)");
    }
    if(mode == 2) {
        hResolution->SetTitle("APX: Res. at DUT (M26+TPX3)");
    }
    if(mode == 3) {
        hResolution->SetTitle("CP2: Res. at DUT (M26)");
    }
    if(mode == 4) {
        hResolution->SetTitle("CP2: Res. at DUT (M26+TPX3)");
    }

    c1->Write();

    // Write result to file
    out->Write();
    out->Close();
    return 0;
}
