#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include <fstream>
#include <iomanip>

#include <TMatrixDSym.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>

// Levy reader
#include "Levy_proj_reader.h"

Levy_reader* myLevy_reader;

using namespace RooFit;

using namespace std;


// ===============================
// Particle container
// ===============================
struct Particle {
    Int_t event_id;
    Int_t pid;
    Float_t px, py, pz, E;
    Float_t x, y, z, t;
    
    Float_t pT() const { return TMath::Sqrt(px*px + py*py); }
    Float_t p() const { return TMath::Sqrt(px*px + py*py + pz*pz); }
    
    Float_t eta() const {
        Float_t p_abs = p();
        if (p_abs <= fabs(pz)) return (pz > 0) ? 999.0 : -999.0;
        return 0.5 * TMath::Log((p_abs + pz) / (p_abs - pz));
    }
    
};

// =========================================
// Levy fit function
// =========================================
double LevySource3D(double *x, double *par)
{
  double alpha = par[0];
  double R     = par[1];
  double N     = par[2]; // This will now stay near 1.0
  double Rcc   = R * pow(2.0, 1.0/alpha);
   
    return  (N / (Rcc*Rcc*Rcc)) * myLevy_reader->getValue_3d(alpha, x[0]/Rcc);
}

// ===============================
// Load OSCAR file
// ===============================
vector<Particle> LoadOSCAR(const char* fname) {
    vector<Particle> particles;
    
    ifstream fin(fname);
    if (!fin.is_open()) {
        cout << "ERROR: Cannot open OSCAR file: " << fname << endl;
        return particles;
    }

    string line;
    int event_id = -1;

    while (getline(fin, line)) {
        if (line.empty()) continue;

        if (line[0] == '#') {
            if (line.find("# event") != string::npos) {
                string tmp;
                stringstream ss(line);
                ss >> tmp >> tmp >> event_id;
            }
            continue;
        }

        Particle p;
        float mass;
        int charge, dummy;
        
       

        stringstream ss(line);
        ss >> p.t >> p.x >> p.y >> p.z >> mass
           >> p.E >> p.px >> p.py >> p.pz
           >> p.pid >> dummy >> charge;

        p.event_id = event_id;
        particles.push_back(p);
    }

    fin.close();
    cout << "Loaded " << particles.size() << " particles" << endl;
    return particles;
}

// =========================================
    // Configuration
    // =========================================
    const char* collision_system = "Au Au";
    const char* centrality = "0-5%";
    double sqrt_sNN = 200;
    
    // Kinematic cuts applied to all particles
    double pT_min = 0.15;     // GeV/c
    double pT_max = 1.0;      // GeV/c
    double eta_cut = 1.0;     


// ===============================
// Analyze pairs function
// ===============================
TH1F* AnalyzePairs(const vector<Particle>& particles,
                   double kT_min,
                   double kT_max,
                   bool pions_only,
                   const char* hist_name,
                   const char* hist_title,
                   int &pairs_count_out)
{
    
    // Group particles by event
    map<int, vector<Particle>> events;
    for (auto& p : particles) {
        if (pions_only && abs(p.pid) != 211) continue;
        events[p.event_id].push_back(p);
    }
    
// =======================
// Histogram binning
// =======================

const int n_bins = 150;
double rho_min = 0.7;
double rho_max = 40.0;

 double bins[n_bins + 1];
    
    // constant rho space ratio
    double r = pow(rho_max / rho_min, 1.0 / n_bins);

    // build bin sequence
     bins[0] = rho_min;
    for (int i = 1; i <= n_bins; i++) {
    bins[i] = bins[i-1] * r;
     }
     
TH1F* hRho = new TH1F(hist_name, hist_title, n_bins, bins);
hRho->Sumw2();
    
    // Pair analysis
    int pairs_count = 0;
    
    for (auto& ev : events) {
        const auto& v = ev.second;

        for (size_t i = 0; i < v.size(); ++i) {
            const Particle& p1 = v[i];
            
            // Single particle cuts
            if (p1.pT() < pT_min || p1.pT() > pT_max) continue;
            if (fabs(p1.eta()) > eta_cut) continue;

            
            for (size_t j = i + 1; j < v.size(); ++j) {
                const Particle& p2 = v[j];
                
                // Single particle cuts for second particle
                if (p2.pT() < pT_min || p2.pT() > pT_max) continue;
                if (fabs(p2.eta()) > eta_cut) continue;

                
                // For pions only
                if (pions_only && p1.pid != p2.pid) continue;
                
                // Calculate pair kinematics
                double kx = 0.5 * (p1.px + p2.px);
                double ky = 0.5 * (p1.py + p2.py);
                double kz = 0.5 * (p1.pz + p2.pz);
                double kT = TMath::Sqrt(kx * kx + ky * ky); 
                
                // pair transverse mass
                double mpi = 0.13957;
                double mT = sqrt(kT*kT + mpi*mpi);

                 // mT dependent cut
                 double c = 0.15 ; 
                 double qmCut = sqrt(c * mT);
           
                  // relative momentum  in the LCMS
                 double qx = p1.px - p2.px;
                 double qy = p1.py - p2.py;
                 double qzL = (4*pow((p1.pz * p2.E - p2.pz * p1.E),2))/(pow((p1.E+p2.E),2)-pow((p1.pz+p2.pz),2));

                 double qLCMS = TMath::Sqrt(qx*qx + qy*qy + qzL);

               // Pair cuts due to kT
                if (kT < kT_min || kT > kT_max) continue;
                 
                //pair cuts due to QLCMS frame
                if (qLCMS > qmCut) continue;

                
                pairs_count++;
                
                // Spatial separation
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                double dt = p1.t - p2.t;
                
                // LCMS transformation
                double K0 = 0.5 * (p1.E + p2.E); //pair energy
                double kp = K0*K0 - kz*kz;
                 double phi = atan2(ky, kx); // azimuthal angle
                
                // Bertsch-Pratt coordinates
                double r_out  = TMath::Cos(phi) * dx + TMath::Sin(phi) * dy - (kT/kp) * (K0*dt - kz*dz);
                double r_side = - TMath::Sin(phi) * dx + TMath::Cos(phi) * dy;
                double r_long = (K0*dz - kz*dt) / sqrt(kp);
                
                double rho = sqrt(r_out * r_out + r_side * r_side + r_long * r_long);
                
                hRho->Fill(rho);
  

            }
        }
    }
// =========================================
// Convert Counts to 3D Density D(rho)
// =========================================
double total_pairs = hRho->Integral(1, hRho->GetNbinsX()+1);

for (int i = 1; i <= hRho->GetNbinsX(); i++)
{
    double rho = hRho->GetBinCenter(i);
    double bw  = hRho->GetBinWidth(i);

    // Spherical phase-space element
    double volume_element = 4.0 * TMath::Pi() * rho * rho * bw;

    if (volume_element > 0 && total_pairs > 0)
    {
        double content = hRho->GetBinContent(i);
        double error   = hRho->GetBinError(i);

        // Normalize to probability density
        hRho->SetBinContent(i, content / (volume_element * total_pairs));
        hRho->SetBinError(i,   error   / (volume_element * total_pairs));
    }
}

     // Safety for log plots
     //hRho->SetMinimum(1e-12);
     
     pairs_count_out = pairs_count;
    cout << (pions_only ? "Pion pairs: " : "All pairs: ") << pairs_count << " pairs analyzed" << endl;

    return hRho;
}

// ===============================
// Main function
// ===============================
void auauok() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1::AddDirectory(kFALSE);
    
    // =========================================
    // Configuration
    // =========================================
    const char* oscar_file = "/home/zeinab/Documents/vhlle-smash/hybrid/sampler.out/test_run/particle_lists.oscar";
    
    // =========================================
    // Load particles
    // =========================================
    vector<Particle> particles = LoadOSCAR(oscar_file);
    if (particles.empty()) {
        cout << "ERROR: No particles loaded!" << endl;
        return;
    }
    
    // =========================================
    // Create output files once
    // =========================================
    TFile* outFile = new TFile("hbt_comparison_kt3237results.root","RECREATE");
    ofstream fout("levy_fit_resultsWref30fmkt3237.txt");
     fout << "\n==================================================================================================\n";
fout << "                 Levy Fit Results ,rho range( rhomin(0.1,0.7,1.0)-rhomax (due B param) fm)  \n" ;
fout << "=========================================================================================================\n";

   fout << setw(10) << "kTmin"
     << setw(10) << "kTmax"
     << setw(10) << "rho_min"
     << setw(10) << "B"
     << setw(12) << "rho_max"
     << setw(18) << "alpha"
     << setw(16) << "R [fm]"
     << setw(16) << "lambda"
     << setw(15) << "chi2"
     << setw(15)  << "NDF"
     << setw(15) << "CL"
     << "\n";

    // =========================================
    // Analyze pion pairs
    // =========================================
    vector<double> kT_bins;

double kT_min_val = 0.32;
double kT_max_val = 0.37;
double step = 0.05;

for(double k = kT_min_val; k <= kT_max_val + 1e-6; k += step)
{
    kT_bins.push_back(k);
}

int nBins = kT_bins.size()-1;
myLevy_reader = new Levy_reader("levy_proj3D_values.dat");

//storing map of the mt vs fitting paramters 

map<double, map<double, vector<double>>> rho_map;
map<double, map<double, vector<double>>> R_map, R_map_err;
map<double, map<double, vector<double>>> lambda_map, lambda_map_err;
map<double, map<double, vector<double>>> R_vs_mT, Rerr_vs_mT;
map<double, map<double, vector<double>>> mT_map;

// ===============================
// LOOP OVER kT BINS
// ===============================
for(int ibin = 0; ibin < nBins; ibin++)
{
    double kT_min = kT_bins[ibin];
    double kT_max = kT_bins[ibin+1];

    cout << "Processing kT: " << kT_min << " - " << kT_max << endl;

    int pion_pairs_count = 0;

    TH1F* hRho_pions = AnalyzePairs(
        particles,
        kT_min,
        kT_max,
        true,
        Form("hRho_pions_%d", ibin),
        "Pion pairs",
        pion_pairs_count
    );

    // ===============================
    // Canvas
    // ===============================
    TCanvas* c1 = new TCanvas(Form("c1_%d", ibin), "HBT Comparison", 1400, 900);
    c1->SetLogx();
    c1->SetLogy();
    c1->SetGrid();

    hRho_pions->SetMarkerStyle(20);
    hRho_pions->SetMarkerSize(0.7);
    hRho_pions->Draw("P");
    
    // ===============================
    // SCAN fit PARAMETERS
    // ===============================
vector<double> rho_min_list = {0.5, 0.7, 1.0};
vector<double> B_values     = {1600, 2500, 3600};

TCanvas* c9 = new TCanvas("c9", "Levy Scan 3x3 pads", 2000, 1800);

c9->Divide(3,3); 
int pad_idx = 1;
int color = 1;
double kT_mean = 0.5 * (kT_min + kT_max);
double mpi = 0.13957;
double mT = sqrt(kT_mean*kT_mean + mpi*mpi);

for(double rho_min_scan : rho_min_list)
{
    for(double B : B_values)
    {
        double rho_max_val = sqrt(B / mT);

        c9->cd(pad_idx);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetGridx();
        gPad->SetGridy();
        
       

// Clone histogram
TH1F* h_clone = (TH1F*)hRho_pions->Clone(Form("h_clone_%d", pad_idx));
h_clone->SetDirectory(0);


h_clone->SetMarkerStyle(20);
h_clone->GetXaxis()->SetTitle("#rho [fm]");
h_clone->GetYaxis()->SetTitle("D(#rho)");
h_clone->Draw("P");


// Define fit
TF1* levy_scan = new TF1(
    Form("levy_r%.1f_B%.0f_pad%d", rho_min_scan, B, pad_idx),
    LevySource3D,
    rho_min_scan,
    rho_max_val,
    3
);

levy_scan->SetParameters(1.7, 4.5, 1.0);

  levy_scan->SetParLimits(0, 0.5, 2.0);
    levy_scan->SetParLimits(1, 2.0, 15.0);
    levy_scan->SetParLimits(2, 0.5, 2.0);


// Fit
h_clone->Fit(levy_scan, "QRN", "", rho_min_scan, rho_max_val);

// Draw fit 
levy_scan->SetLineColor(kRed);
levy_scan->SetLineWidth(2);
levy_scan->Draw("SAME");


        // Info box with fit parameters
        double alpha = levy_scan->GetParameter(0);
        double R_val = levy_scan->GetParameter(1);
        double N_val = levy_scan->GetParameter(2);
        double alpha_err = levy_scan->GetParError(0);
        double R_err_val = levy_scan->GetParError(1);
        double N_err_val = levy_scan->GetParError(2);
        double chi2 = levy_scan->GetChisquare();
        double ndf  = levy_scan->GetNDF();
        double CL   = TMath::Prob(chi2, ndf);
        
        // ----- Fit results box -----

        TPaveText* info = new TPaveText(0.55,0.55,0.95,0.9,"NDC");
        info->SetFillColor(0);
        info->SetBorderSize(1);
        info->AddText(Form("rho_min = %.2f", rho_min_scan));
        info->AddText(Form("B       = %.0f", B));
        info->AddText(Form("#alpha = %.3f #pm %.3f", levy_scan->GetParameter(0), levy_scan->GetParError(0)));
        info->AddText(Form("R = %.3f #pm %.3f fm", levy_scan->GetParameter(1), levy_scan->GetParError(1)));
        info->AddText(Form("#lambda = %.3f #pm %.3f", levy_scan->GetParameter(2), levy_scan->GetParError(2)));
        info->AddText(Form("chi2/ndf= %.2f/%d", chi2, int(ndf)));
        info->AddText(Form("C.L. = %.5f%%", CL * 100));
        info->Draw();

        pad_idx++;
        
        // ===== Save to file ===== 
        
        fout << setw(10) << Form("%.2f", kT_min) << setw(10) << Form("%.2f", kT_max) << setw(10) << rho_min_scan << setw(10) << B << setw(12) << Form("%.3f", rho_max_val) << setw(19) << Form("%.4f ± %.4f", alpha, alpha_err) << setw(19) << Form("%.4f ± %.4f", R_val, R_err_val) << setw(19) << Form("%.4f ± %.4f", N_val, N_err_val) << setw(15) << Form("%.2f", chi2) << setw(12) << ndf << setw(12) << Form("%.6f", CL) << "\n";
    }
}

// Save canvas
c9->SaveAs(Form("Levy_scan_3x3_bin_%d.png", ibin));
c9->Write();
}
fout.close();
    outFile->Close();
}
