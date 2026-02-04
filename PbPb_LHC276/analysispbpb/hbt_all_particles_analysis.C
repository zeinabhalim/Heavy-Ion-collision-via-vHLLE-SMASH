// hbt_comparison_analysis.C
// Comparison of HBT analysis for pion pairs vs all particle pairs

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

/*// ===============================
// Lévy source function
// D(rho) = N * rho^2 * exp[-1/2 (rho/R)^alpha]
// ===============================
double LevySource(double *x, double *par) {
    double rho = x[0];
    double alpha = par[0];
    double R = par[1];
    double N = par[2];
    
    if (rho <= 0 || alpha <= 0 || R <= 0) return 0.0;
    
    double exponent = -0.5 * TMath::Power(rho / R, alpha);
    return N * rho * rho * TMath::Exp(exponent);
}*/

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
    const char* collision_system = "Pb Pb";
    const char* centrality = "20-30%";
    double sqrt_sNN = 2.76;
    
    // Kinematic cuts applied to all particles
    double pT_min = 0.15;     // GeV/c
    double pT_max = 2.0;      // GeV/c
    double eta_cut = 0.8;     // |η| < 1.0
    double kT_min = 0.3;    // GeV/c
    double kT_max = 0.4;    // GeV/c
    

    
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
    
    // Create histogram
    const int n_bins = 100;
    double log_min = log10(0.1);
    double log_max = log10(1000.0);
    double bins[n_bins + 1];
    
    for (int i = 0; i <= n_bins; i++) {
        bins[i] = pow(10, log_min + (log_max - log_min) * i / n_bins);
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
                
                // Pair cuts
                if (kT < kT_min || kT > kT_max) continue;
                
                
                pairs_count++;
                
                // Spatial separation
                double dx = p1.x - p2.x;
                double dy = p1.y - p2.y;
                double dz = p1.z - p2.z;
                double dt = p1.t - p2.t;
                
                // LCMS transformation
                double betaL = kz / (0.5 * (p1.E + p2.E));
                double gammaL = 1.0 / sqrt(1.0 - betaL * betaL);
                double dz_lcms = gammaL * (dz - betaL * dt);
                
                // Bertsch-Pratt coordinates
                double r_out  = (kx / kT) * dx + (ky / kT) * dy;
                double r_side = (-ky / kT) * dx + (kx / kT) * dy;
                double r_long = dz_lcms;
                
                double rho = sqrt(r_out * r_out + r_side * r_side + r_long * r_long);
                
                hRho->Fill(rho);
            }
        }
    }
    
    // Normalize histogram
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
void hbt_all_particles_analysis() {
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
    
  
    // =========================================
    // Create comparison plot
    // =========================================
    TCanvas* c1 = new TCanvas("c1", "HBT Comparison: Pions vs All Particles", 1000, 800);
    c1->SetLogx();
    c1->SetLogy();
    c1->SetGridx();
    c1->SetGridy();
    
    // Format pion histogram
    hRho_pions->SetLineColor(kRed);
    hRho_pions->SetMarkerColor(kRed);
    hRho_pions->SetMarkerStyle(20);
    hRho_pions->SetMarkerSize(0.8);
    hRho_pions->SetLineWidth(1);
    
    // Format all particles histogram
    hRho_all->SetLineColor(kBlue);
    hRho_all->SetMarkerColor(kBlue);
    hRho_all->SetMarkerStyle(21);
    hRho_all->SetMarkerSize(0.8);
    hRho_all->SetLineWidth(1);
    
    // Find suitable Y-axis range
    double max_val = max(hRho_pions->GetMaximum(), hRho_all->GetMaximum());
    double min_val = min(hRho_pions->GetMinimum(1e-10), hRho_all->GetMinimum(1e-10));
    
    // Draw pion pairs first
    hRho_pions->SetMaximum(max_val * 10);
    hRho_pions->SetMinimum(min_val * 0.1);
    hRho_pions->Draw("E1");
    
    // Draw all particles on top
    hRho_all->Draw("E1 SAME");
    
    // Create legend
    TLegend* leg = new TLegend(0.15, 0.70, 0.45, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    
    leg->AddEntry(hRho_pions, "Identical pion pairs", "PE");
   
    leg->AddEntry(hRho_all, "All particle pairs", "PE");
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
