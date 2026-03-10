#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

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
  double R = par[1];
  double N = par[2];
  double Rcc = (R*pow(2.,1./alpha));
  return (2.*N/Rcc)*(myLevy_reader->getValue_1d(alpha, x[0]/Rcc));
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
    double eta_cut = 1.0;     // |η| < 1.0
    double kT_min = 0.325;    // GeV/c
    double kT_max = 0.375;    // GeV/c
    
// ===============================
// Analyze pairs function
// ===============================
TH1F* AnalyzePairs(const vector<Particle>& particles, 
                   bool pions_only = false,
                   const char* hist_name = "hRho",
                   const char* hist_title = "Distance distribution") {
    
    // Group particles by event
    map<int, vector<Particle>> events;
    for (auto& p : particles) {
        if (pions_only && abs(p.pid) != 211) continue;
        events[p.event_id].push_back(p);
    }
    
// =======================
// Histogram binning
// =======================

const int n_bins = 200;
double rho_min = 0.2;
double rho_max = 100.0;

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
    
     // optimize counts to bin width and 4*pi space
    for (int i = 1; i <= hRho->GetNbinsX(); i++)
    {
    double rho = hRho->GetBinCenter(i);
    double bw  = hRho->GetBinWidth(i);
    
    double source_volume = 4.0 * TMath::Pi() * rho * rho * bw;

    if (source_volume > 0){
        hRho->SetBinContent(i, hRho->GetBinContent(i) / source_volume);
        hRho->SetBinError(i,  hRho->GetBinError(i) / source_volume);
       }
    }
    
    // Normalizing histogram
    double total_integral = hRho->Integral("width");
    if (total_integral > 0) {
        hRho->Scale(1.0 / total_integral);
        
    }
    
    cout << (pions_only ? "Pion pairs: " : "All pairs: ") << pairs_count << " pairs analyzed" << endl;
    
    return hRho;
}



// ===============================
// Main function
// ===============================
void auauspital() {
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
    // Analyze pion pairs
    // =========================================
    cout << "\n=== Analyzing Pion Pairs ===" << endl;
    TH1F* hRho_pions = AnalyzePairs(particles, true, 
                                   "hRho_pions", 
                                   "Pion pairs distance distribution; #rho [fm]; D(#rho)");
    // =========================================
    // Analyze all particle pairs
    // =========================================
    cout << "\n=== Analyzing All Particle Pairs ===" << endl;
    TH1F* hRho_all = AnalyzePairs(particles, false, 
                                 "hRho_all", 
                                 "All particle pairs distance distribution; #rho [fm]; D(#rho)");          
            
     TFile* outFile = new TFile("hbt_comparison_results.root","RECREATE");
   // =========================================
// Create canvas
// =========================================
TCanvas* c1 = new TCanvas("c1","HBT Comparison",1000,800);
c1->SetLogx();
c1->SetLogy();
c1->SetGridx();
c1->SetGridy();

// Histogram styling
hRho_pions->SetLineColor(kRed);
hRho_pions->SetMarkerColor(kRed);
hRho_pions->SetMarkerStyle(20);
hRho_pions->SetMarkerSize(0.8);

hRho_all->SetLineColor(kBlue);
hRho_all->SetMarkerColor(kBlue);
hRho_all->SetMarkerStyle(21);
hRho_all->SetMarkerSize(0.8);

// Axis labels
hRho_pions->GetXaxis()->SetTitle("#rho [fm]");
hRho_pions->GetYaxis()->SetTitle("D(#rho)");

// Find range
double max_val = std::max(hRho_pions->GetMaximum(), hRho_all->GetMaximum());

// Draw histogram first
hRho_pions->SetMinimum(1e-10);
hRho_pions->SetMaximum(max_val*2);
hRho_pions->Draw("P");

// Overlay second histogram
hRho_all->Draw("P SAME");

// =========================================
// Perform levy fit
// =========================================

// Initialize Levy reader object
myLevy_reader = new Levy_reader("levy_proj3D_values.dat");

TF1 *levy = new TF1("levy", LevySource3D, 5, 50, 3);
levy->SetParNames("alpha","R","N");
levy->SetParameters(2.0, 4.38, 0.85); 
levy->FixParameter(0,2.0);
levy->SetParLimits(1, 0.1, 20.0); 
levy->SetParLimits(2, 0.0, 2.0); 
levy->SetLineColor(kBlack);

levy->SetLineWidth(2);

levy->Draw("P");
    
    // Format pion histogram
    hRho_pions->SetLineColor(kRed);
    hRho_pions->SetMarkerColor(kRed);
    hRho_pions->SetMarkerStyle(20);
    hRho_pions->SetMarkerSize(0.8);
    hRho_pions->SetLineWidth(1);
    
    
    hRho_pions->Fit(levy,"R"); 
     cout << "\npion pairs Levy fit of AU-AU collision, levy fit:" << endl;
    cout << "alpha = " << levy->GetParameter(0)
     << " ± " << levy->GetParError(0) << endl;

     cout << "R = " << levy->GetParameter(1)
     << " ± " << levy->GetParError(1) << " fm" << endl;
    
    
    // Format all particles histogram
    hRho_all->SetLineColor(kBlue);
    hRho_all->SetMarkerColor(kBlue);
    hRho_all->SetMarkerStyle(21);
    hRho_all->SetMarkerSize(0.8);
    hRho_all->SetLineWidth(1);
    
    // Reset initial parameters
    levy->SetParameters(2.0, 4.3, 0.85); 
    levy->FixParameter(0,2.0);
    hRho_all->Fit(levy,"R"); 
    
    //fitting parameters
    cout << "\nAll particle pairs of AU-AU collision, levy fit:" << endl;
cout << "alpha = " << levy->GetParameter(0)
     << " ± " << levy->GetParError(0) << endl;
cout << "R     = " << levy->GetParameter(1)
     << " ± " << levy->GetParError(1) << " fm" << endl;
        
    
    hRho_pions->GetXaxis()->SetTitle("#rho [fm]");
hRho_pions->GetYaxis()->SetTitle("D(#rho)");

hRho_all->GetXaxis()->SetTitle("#rho [fm]");
hRho_all->GetYaxis()->SetTitle("D(#rho)");
    
    // Find suitable Y-axis range
    double min_val = min(hRho_pions->GetMinimum(1e-10), hRho_all->GetMinimum(1e-10));
    
    hRho_pions->SetMinimum(1e-8);
hRho_pions->SetMaximum(max_val*2);

hRho_pions->Draw("P");

    
    // Draw all particles on top
    hRho_all->Draw("P SAME");

    // Create legend
    TLegend* leg = new TLegend(0.15,0.70,0.45,0.88);
leg->SetBorderSize(0);
leg->SetFillStyle(0);

leg->AddEntry(hRho_pions,"Identical pion pairs","PE");
leg->AddEntry(hRho_all,"All particle pairs","PE");
leg->AddEntry(levy,"Gauss. fit","L");

leg->Draw();
    
    // Add title
    TLatex* title = new TLatex();
    title->SetNDC();
    title->SetTextFont(42);
    title->SetTextSize(0.045);
    title->DrawLatex(0.15, 0.92, "Hybird vHlle + SMASH");
    
    // Add system information
    TPaveText* infoBox = new TPaveText(0.60, 0.70, 0.95, 0.8, "NDC");
    infoBox->SetFillColor(0);
    infoBox->SetBorderSize(1);
    infoBox->SetTextFont(42);
    infoBox->SetTextSize(0.035);
    infoBox->SetTextAlign(12);
    infoBox->AddText(Form("%s , #sqrt{s_{NN}} = %.2f TeV", 
                          collision_system, sqrt_sNN));
    infoBox->AddText(Form("Centrality : %s", 
                          centrality));
    infoBox->AddText(Form("  %.2f < p_{T} [GeV/c] < %.2f", pT_min, pT_max));
    infoBox->AddText(Form("  |#eta| < %.1f", eta_cut));
    infoBox->AddText(Form("  %.2f < k_{T}[GeV/c] < %.2f", kT_min, kT_max));
    
    infoBox->Draw();
    // Save plot
    c1->SaveAs("hbt_comparison_pions_vs_all.png");
    
    
    // Save histograms
    hRho_pions->Write();
    hRho_all->Write();
    
    // Save canvases
    c1->Write();

    
    cout << "\n========================================" << endl;
    cout << "Output files created:" << endl;
    cout << "  - hbt_comparison_pions_vs_all.png" << endl;
    cout << "  - hbt_comparison_results.root" << endl;
    cout << "========================================" << endl;
    
}
