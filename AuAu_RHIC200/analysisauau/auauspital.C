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
double rho_min = 0.5;
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

void Draw_R_vs_mT_3Pads(
    map<double, map<double, vector<double>>>& Rmap,
    map<double, map<double, vector<double>>>& Rerrmap,
    map<double, map<double, vector<double>>>& mTmap,
    vector<double>& rho_min_list,
    vector<double>& B_values)
{
    TCanvas* c = new TCanvas("c_R_mT_3pads","R vs mT",1200,400);
    c->Divide(3,1);

    vector<int> colors = {kRed, kBlue, kGreen+2};

    for(int i_rho = 0; i_rho < rho_min_list.size(); i_rho++)
    {
        double rho_min = rho_min_list[i_rho];

        c->cd(i_rho+1);
        gPad->SetGrid();

        TLegend* leg = new TLegend(0.50,0.65,0.88,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        bool first = true;

        for(int i_B = 0; i_B < B_values.size(); i_B++)
        {
            double B = B_values[i_B];

            auto& Rvals  = Rmap[rho_min][B];
            auto& Rerrs  = Rerrmap[rho_min][B];
            auto& mTvals = mTmap[rho_min][B];

            int n = Rvals.size();
            if(n == 0) continue;

            TGraphErrors* g = new TGraphErrors(
                n,
                &mTvals[0],
                &Rvals[0],
                nullptr,
                &Rerrs[0]
            );

            g->SetMarkerStyle(20 + i_B);
            g->SetMarkerColor(colors[i_B]);
            g->SetLineColor(colors[i_B]);
            g->SetLineWidth(2);

            if(first)
            {
                g->SetTitle(Form("#rho_{min} = %.1f; m_{T} [GeV]; R [fm]", rho_min));
                g->SetMinimum(3.2);
                g->SetMaximum(5.5);
                g->Draw("AP");
                first = false;
            }
            else
            {
                g->Draw("P SAME");
            }

            leg->AddEntry(g, Form("B = %.0f", B), "P");
        }

        leg->Draw();
        // ===== rho_min label box =====
TPaveText* rhoBox = new TPaveText(0.15, 0.75, 0.40, 0.88, "NDC");

rhoBox->SetFillStyle(0);
rhoBox->SetBorderSize(0);
rhoBox->SetTextFont(43);
rhoBox->SetTextSize(18);
rhoBox->SetTextAlign(12);

rhoBox->AddText(Form("#rho_{min} = %.1f fm", rho_min));

rhoBox->Draw();
    }

    c->SaveAs("R_vs_mT_3pads.png");
}

void Draw_R_vs_mT_ByB(
    map<double, map<double, vector<double>>>& Rmap,
    map<double, map<double, vector<double>>>& Rerrmap,
    map<double, map<double, vector<double>>>& mTmap,
    vector<double>& rho_min_list,
    vector<double>& B_values)
{
    int nPads = B_values.size();
    TCanvas* c = new TCanvas("c_R_mT_ByB","R vs mT by B",1200,400);
    c->Divide(nPads,1); // One pad per B value

    vector<int> colors = {kRed, kBlue, kGreen+2}; // for rho_min curves

    for(int i_B = 0; i_B < B_values.size(); i_B++)
    {
        double B = B_values[i_B];

        c->cd(i_B+1);
        gPad->SetGrid();

        TLegend* leg = new TLegend(0.55,0.65,0.88,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        bool first = true;

        for(int i_rho = 0; i_rho < rho_min_list.size(); i_rho++)
        {
            double rho_min = rho_min_list[i_rho];

            auto& Rvals  = Rmap[rho_min][B];
            auto& Rerrs  = Rerrmap[rho_min][B];
            auto& mTvals = mTmap[rho_min][B];

            int n = Rvals.size();
            if(n == 0) continue;

            TGraphErrors* g = new TGraphErrors(
                n,
                &mTvals[0],
                &Rvals[0],
                nullptr,
                &Rerrs[0]
            );

            g->SetMarkerStyle(20 + i_rho);
            g->SetMarkerColor(colors[i_rho]);
            g->SetLineColor(colors[i_rho]);
            g->SetLineWidth(2);

            if(first)
            {
                g->SetTitle(Form("B = %.0f; m_{T} [GeV]; R [fm]", B));
                g->SetMinimum(3.2);
                g->SetMaximum(5.5);
                g->Draw("AP");
                first = false;
            }
            else
            {
                g->Draw("P SAME");
            }

            leg->AddEntry(g, Form("#rho_{min} = %.1f", rho_min), "P");
        }

        leg->Draw();

        // Optional: B value label box
        TPaveText* BBox = new TPaveText(0.15, 0.75, 0.40, 0.88, "NDC");
        BBox->SetFillStyle(0);
        BBox->SetBorderSize(0);
        BBox->SetTextFont(43);
        BBox->SetTextSize(18);
        BBox->SetTextAlign(12);
        BBox->AddText(Form("B = %.0f", B));
        BBox->Draw();
    }

    c->SaveAs("R_vs_mT_byB.png");
}

// ===============================
// Main function
// ===============================
void auaubackref() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1::AddDirectory(kFALSE);
    
    // =========================================
    // Configuration
    // =========================================
    const char* oscar_file = "/home/zeinab/Documents/vhlle-smash/hybrid/AuAu_RHIC200/sampler.out/cent0_5/particle_lists_0.oscar";
    
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
    TFile* outFile = new TFile("hbt_comparison_results.root","RECREATE");
    ofstream fout("levy_fit_resultsWref30fm.txt");
    
    fout << "\n==================================================================================================\n";
fout << "                 Levy Fit Results using Levy reader file ,rho range( 1-30 fm)  \n";
fout << "=========================================================================================================\n";

     fout << setw(12) << "kT_range"
     << setw(17) << "alpha"
     << setw(17) << "R [fm]"
     << setw(17) << "lambda"
     << setw(15) << "chi2"
     << setw(12) << "NDF"
     << setw(14) << "C.L. percent"
      << setw(15) << "N_pairs"
     << "\n";

    // =========================================
    // Analyze pion pairs
    // =========================================
    vector<double> kT_bins;

double kT_min_val = 0.175;
double kT_max_val = 0.775;
double step = 0.05;

for(double k = kT_min_val; k <= kT_max_val + 1e-6; k += step)
{
    kT_bins.push_back(k);
}

int nBins = kT_bins.size()-1;
myLevy_reader = new Levy_reader("levy_proj3D_values.dat");

//storing the mt vs fitting paramters 
vector<double> mt_vals, mt_err;
vector<double> alpha_vals, alpha_err;
vector<double> R_vals, R_err;
vector<double> N_vals, N_err;
vector<double> B_values = {1600,2500, 3600};

//vector<double> R_scan, lambda_scan, rho_max_vals;
//vector<int> N_scan;

map<double, map<double, vector<double>>> R_map;
map<double, map<double, vector<double>>> R_map_err;
map<double, map<double, vector<double>>> rho_map;
map<double, map<double, vector<double>>> lambda_map;
map<double, map<double, vector<double>>> lambda_map_err;
map<double, vector<double>> chi2_map;

map<double, map<double, vector<double>>> R_vs_mT;
map<double, map<double, vector<double>>> Rerr_vs_mT;
map<double, map<double, vector<double>>> mT_map;
//map<double, vector<double>> rho_map;

for(int ibin=0; ibin<nBins; ibin++)
{
    double kT_min = kT_bins[ibin];
    double kT_max = kT_bins[ibin+1];

    cout<<"Processing kT: "<<kT_min<<" - "<<kT_max<<endl;

    int pion_pairs_count = 0;

TH1F* hRho_pions = AnalyzePairs(
    particles,
    kT_min,
    kT_max,
    true,
    Form("hRho_pions_%d",ibin),
    "Pion pairs",
    pion_pairs_count
);

   // =========================================
// Create canvas
// =========================================
TCanvas* c1 = new TCanvas("c1","HBT Comparison",1400,900);
c1->SetLogx();
c1->SetLogy();
c1->SetGridx();
c1->SetGridy();

// Histogram styling
hRho_pions->SetLineColor(kRed);
hRho_pions->SetMarkerColor(kRed);
hRho_pions->SetMarkerStyle(20);
hRho_pions->SetMarkerSize(0.8);

// Axis labels
hRho_pions->GetXaxis()->SetTitle("#rho [fm]");
hRho_pions->GetYaxis()->SetTitle("D(#rho)");

// Find range
double max_val = hRho_pions->GetMaximum();

// Draw histogram first
hRho_pions->SetMaximum(max_val*2);

// =========================================
// Perform levy fit
  // =========================================

// Initialize Levy reader object


TF1 *levy = new TF1("levy", LevySource3D, 0.7, 40.0, 3);

levy->SetParNames("alpha","R","N");


// Starting values
levy->SetParameters(1.704, 4.48, 1.013);

// Limits for fitting
levy->SetParLimits(0, 0.5, 2.0);   // alpha
levy->SetParLimits(1, 2.0, 15.0);  // R
levy->SetParLimits(2, 0.5, 2.0);   // N

double fit_min = 1.0;
double fit_max = 30.0;
hRho_pions->Fit(levy, "QRN", "", fit_min, fit_max);

// stat. check
double chi2   = levy->GetChisquare();
double ndf    = levy->GetNDF();
double chi2ndf = chi2 / ndf;
double CL     = TMath::Prob(chi2, ndf);

// Compute mean kT and mT
double kT_mean = 0.5 * (kT_min + kT_max);
double mpi = 0.13957;
double mT = sqrt(kT_mean*kT_mean + mpi*mpi);

// Store values
mt_vals.push_back(mT);
mt_err.push_back(0.0); 

alpha_vals.push_back(levy->GetParameter(0));
alpha_err.push_back(levy->GetParError(0));

R_vals.push_back(levy->GetParameter(1));
R_err.push_back(levy->GetParError(1));

N_vals.push_back(levy->GetParameter(2));
N_err.push_back(levy->GetParError(2));

// =============================================
// Stability scan for levy fitting parameters
// =============================================
vector<int> N_events_list = {10,30,100,300,1000,3000,10000};
vector<double> rho_min_list = {0.5, 0.7, 1.0};


vector<TF1*> scan_functions;

int idx = 0;

for(double rho_min_scan : rho_min_list)
{
    for(double B : B_values)
    {
        double rho_max_val = sqrt(B / mT);

            TF1* levy_scan = new TF1(
                Form("levy_scan_%d", idx++),
                LevySource3D,
                rho_min_scan,
                rho_max_val,
                3
            );

            levy_scan->SetParameters(1.7,4.5,1.0);

            // Fit in scan range
            hRho_pions->Fit(levy_scan, "QRN", "", rho_min_scan, rho_max_val);

           /* levy_scan->SetLineStyle(2); // dashed
            levy_scan->SetLineColor(kBlue);
            levy_scan->SetLineWidth(1);*/
          

         double chi2_scan = levy_scan->GetChisquare();
         double ndf_scan  = levy_scan->GetNDF();

        double chi2ndf_scan = (ndf_scan > 0) ? chi2_scan / ndf_scan : 0;
        
        

// store per rho_min

R_map[rho_min_scan][B].push_back(levy_scan->GetParameter(1));
R_map_err[rho_min_scan][B].push_back(levy_scan->GetParError(1));
rho_map[rho_min_scan][B].push_back(rho_max_val);
lambda_map[rho_min_scan][B].push_back(levy_scan->GetParameter(2));
lambda_map_err[rho_min_scan][B].push_back(levy_scan->GetParError(2));

// store per mt
R_vs_mT[rho_min_scan][B].push_back(levy_scan->GetParameter(1));
Rerr_vs_mT[rho_min_scan][B].push_back(levy_scan->GetParError(1));
mT_map[rho_min_scan][B].push_back(mT);

chi2_map[rho_min_scan].push_back(chi2ndf_scan);
//rho_map[rho_min_scan].push_back(rho_max_val);

            scan_functions.push_back(levy_scan);
        
    }
}
// ----- Create two drawing functions -----
// Global fit curves
TF1* f_full = new TF1("f_full", LevySource3D, 0.7, 40.0, 3);
f_full->SetParameters(levy->GetParameters());
f_full->SetLineColor(kRed);
f_full->SetLineStyle(2);
f_full->SetLineWidth(2);
f_full->SetNpx(5000);

TF1* f_fit = new TF1("f_fit", LevySource3D, 1.0, 30 , 3);
f_fit->SetParameters(levy->GetParameters());
f_fit->SetLineColor(kRed);
f_fit->SetLineStyle(1);
f_fit->SetLineWidth(3);
f_fit->SetNpx(1000);

hRho_pions->Draw("P");

// Global fit
f_full->Draw("SAME");
f_fit->Draw("SAME");

// Scan fits
/*for(auto f : scan_functions)
{
    f->Draw("SAME");
}*/
// Print results to console
cout << "\npion pairs Levy fit of AU-AU collision, levy fit:" << endl;
cout << "alpha = " << levy->GetParameter(0) << " ± " << levy->GetParError(0) << endl;
cout << "R     = " << levy->GetParameter(1) << " ± " << levy->GetParError(1) << " fm" << endl;
cout << "lambda= " << levy->GetParameter(2) << " ± " << levy->GetParError(2) << endl;
cout << "chi2/NDF = " << chi2 << "/" << ndf << "  C.L. = " << CL << endl;


// ----- Legend -----
TLegend* leg = new TLegend(0.60, 0.20, 0.90, 0.40);
leg->SetBorderSize(0);
leg->SetFillStyle(0);
leg->AddEntry(hRho_pions, " pion pairs", "PE");
leg->AddEntry(f_fit,      "Levy fit (1-30 fm)", "L");
leg->AddEntry(f_full,      "Levy fit (0.5-40 fm)", "L");
leg->Draw();

// ----- Title -----
TLatex* title = new TLatex();
title->SetNDC();
title->SetTextFont(42);
title->SetTextSize(0.045);

gPad->GetListOfPrimitives()->Print();
// ----- Information box -----
/*TPaveText* infoBox = new TPaveText(0.60, 0.65, 0.99, 0.85, "NDC");
infoBox->SetTextSize(0.035);
infoBox->SetFillStyle(0);   
infoBox->SetBorderSize(0); 
infoBox->SetTextFont(42);
infoBox->SetTextAlign(12);
infoBox->AddText(Form("%s , #sqrt{s_{NN}} = %.0f GeV", collision_system, sqrt_sNN));
infoBox->AddText(Form("Centrality : %s", centrality));
infoBox->AddText(Form("  %.2f < p_{T} [GeV/c] < %.2f", pT_min, pT_max));
infoBox->AddText(Form("  |#eta| < %.1f", eta_cut));
infoBox->AddText(Form("  %.2f < k_{T}[GeV/c] < %.2f", kT_min, kT_max));
infoBox->Draw();*/

// ----- Fit results box -----
TPaveText* fitBox = new TPaveText(0.60, 0.65, 0.99, 0.85, "NDC");
fitBox->SetFillStyle(0);   
fitBox->SetBorderSize(0); 
fitBox->SetTextFont(42);
fitBox->SetTextSize(0.035);
fitBox->SetTextAlign(12);
fitBox->AddText(Form("  %.2f < k_{T}[GeV/c] < %.2f", kT_min, kT_max));
fitBox->AddText(Form("#alpha = %.3f #pm %.3f", levy->GetParameter(0), levy->GetParError(0)));
fitBox->AddText(Form("R = %.3f #pm %.3f fm", levy->GetParameter(1), levy->GetParError(1)));
fitBox->AddText(Form("#lambda = %.3f #pm %.3f", levy->GetParameter(2), levy->GetParError(2)));
fitBox->AddText(Form("#chi^{2}/NDF = %.2f / %.0f", chi2, ndf));
fitBox->AddText(Form("C.L. = %.2f%%", CL * 100));
fitBox->Draw();

//loading fit results
fout << setw(12) << Form("%.3f-%.3f", kT_min, kT_max)
     << setw(20) << Form("%.5f±%.5f", levy->GetParameter(0), levy->GetParError(0))
     << setw(20) << Form("%.5f±%.5f", levy->GetParameter(1), levy->GetParError(1))
     << setw(20) << Form("%.5f±%.5f", levy->GetParameter(2), levy->GetParError(2))
     << setw(10) << Form("%.2f", chi2)
     << setw(10) << Form("%.0f", ndf)
     << setw(13) << Form("%.5f%%", CL*100)
          << setw(15) << pion_pairs_count   // NEW COLUMN
     << "\n";

       // Save canvas as PNG and also to ROOT file
        c1->SaveAs(Form("averAngular_Dist_kT_%.2f-%.2f.png",kT_min,kT_max));
        c1->Write();

    }
    fout.close();
    outFile->Close();
    
    // mt v.s. fitting paramters canvas 
    TCanvas* cSummary = new TCanvas("cSummary","Fit Parameters vs mT",1200,900);
cSummary->Divide(2,2);

// Convert to arrays (needed for TGraphErrors)
int n = mt_vals.size();

// -------- alpha --------
cSummary->cd(1);
TGraphErrors* gAlpha = new TGraphErrors(n,
    &mt_vals[0], &alpha_vals[0],
    &mt_err[0],  &alpha_err[0]);

gAlpha->SetTitle("#alpha vs m_{T}; m_{T} [GeV]; #alpha");
gAlpha->SetMarkerStyle(20);
gAlpha->Draw("AP");

// -------- R --------
cSummary->cd(2);
TGraphErrors* gR = new TGraphErrors(n,
    &mt_vals[0], &R_vals[0],
    &mt_err[0],  &R_err[0]);

gR->SetTitle("R vs m_{T}; m_{T} [GeV]; R [fm]");
gR->SetMarkerStyle(21);
gR->Draw("AP");

// -------- N (lambda) --------
cSummary->cd(3);

TGraphErrors* gN = new TGraphErrors(n,
    &mt_vals[0], &N_vals[0],
    &mt_err[0],  &N_err[0]);

gN->SetTitle("#lambda vs m_{T}; m_{T} [GeV]; #lambda");
gN->SetMarkerStyle(22);

gN->Draw("AP");   // Draw FIRST

// the 4th pad
cSummary->cd(4);

// Optional: clean look
gPad->SetFillColor(0);
gPad->SetFrameBorderMode(0);

// Create info box
TPaveText* infoBox = new TPaveText(0.10, 0.10, 0.90, 0.90, "NDC");

infoBox->SetTextSize(0.05);
infoBox->SetFillColor(0);
infoBox->SetBorderSize(1);
infoBox->SetTextFont(42);
infoBox->SetTextAlign(12);

// Add text
infoBox->AddText("System Information");
infoBox->AddText("---------------------------");

infoBox->AddText(Form("%s", collision_system));
infoBox->AddText(Form("#sqrt{s_{NN}} = %.0f GeV", sqrt_sNN));
infoBox->AddText(Form("Centrality: %s", centrality));

infoBox->AddText(" ");
infoBox->AddText("Kinematic cuts:");

infoBox->AddText(Form("%.2f < p_{T} < %.2f GeV/c", pT_min, pT_max));
infoBox->AddText(Form("|#eta| < %.1f", eta_cut));
double rhorangmax = 30;
double rhorangmin = 1;
infoBox->AddText(Form("%.1f < #rho < %.1f fm",rhorangmin, rhorangmax));

infoBox->Draw();

// stability check canvas 
// -------------------
auto DrawRhoMinGraphs = [&](const map<double, map<double, vector<double>>>& values_map,
                           const map<double, map<double, vector<double>>>& errors_map,
                           const map<double, map<double, vector<double>>>& rho_map,
                           const string& yTitle,
                           const string& fileName,
                           const vector<double>& rho_min_list,
                           const vector<double>& B_values)
{
    TCanvas* c = new TCanvas(Form("c_%s", fileName.c_str()),
                             Form("%s vs #rho_{max}", yTitle.c_str()), 1200, 900);
    c->Divide(2,2);

    vector<int> colors = {kRed, kBlue, kGreen+2};

    // ===== FIXED AXIS RANGES =====
    double xmin = 40;
    double xmax = 130;

    double ymin_R = 3.2;
    double ymax_R = 5.4;

    int pad = 1;

    for(double rho_min : rho_min_list)
    {
        c->cd(pad);

        TLegend* leg = new TLegend(0.60, 0.65, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        int color_idx = 0;

        for(double B : B_values)
        {
            const auto& vals = values_map.at(rho_min).at(B);
            const auto& errs = errors_map.at(rho_min).at(B);
            const auto& rhos = rho_map.at(rho_min).at(B);

            int n = vals.size();
            if(n == 0) continue;

            TGraphErrors* g = new TGraphErrors(n,
                                               &rhos[0], &vals[0],
                                               nullptr, &errs[0]);

            g->SetMarkerStyle(20 + color_idx);
            g->SetMarkerColor(colors[color_idx]);
            g->SetLineColor(colors[color_idx]);
            g->SetLineWidth(2);

            if(color_idx == 0)
            {
                // ===== Titles & axes =====
                g->SetTitle(Form("%s vs #rho_{max}", yTitle.c_str()));
                g->GetXaxis()->SetTitle("#rho_{max} [fm]");
                g->GetYaxis()->SetTitle(yTitle.c_str());

                // ===== FIXED AXES =====
                g->GetXaxis()->SetLimits(xmin, xmax);

                if(yTitle == "R [fm]")
                {
                    g->SetMinimum(ymin_R);
                    g->SetMaximum(ymax_R);
                }

                g->Draw("AP");
            }
            else
            {
                g->Draw("P SAME");
            }

            leg->AddEntry(g, Form("B = %.0f", B), "LP");

            color_idx++;
        }

        leg->Draw();

        // ===== rho_min box =====

TPaveText* rhoBox = new TPaveText(0.15, 0.75, 0.40, 0.88, "NDC");

rhoBox->SetFillStyle(0);   
rhoBox->SetBorderSize(0);  
rhoBox->SetTextFont(43);
rhoBox->SetTextSize(18);   
rhoBox->SetTextAlign(12);

rhoBox->AddText(Form("#rho_{min} = %.1f fm", rho_min));

rhoBox->Draw();

        pad++;
    }

if(pad <= 4) {
        c->cd(pad);
        gPad->SetFillColor(0);
        gPad->SetFrameBorderMode(0);
        TPaveText* infoBox = new TPaveText(0.1,0.1,0.9,0.9,"NDC");
     
infoBox->AddText("System Information");
infoBox->AddText("---------------------------");

infoBox->AddText(Form("%s", collision_system));
infoBox->AddText(Form("#sqrt{s_{NN}} = %.0f GeV", sqrt_sNN));
infoBox->AddText(Form("Centrality: %s", centrality));

infoBox->AddText(" ");
infoBox->AddText("Kinematic cuts:");

infoBox->AddText(Form("%.2f < p_{T} < %.2f GeV/c", pT_min, pT_max));
infoBox->AddText(Form("|#eta| < %.1f", eta_cut));
infoBox->Draw();
    }

    c->SaveAs(Form("%s.png", fileName.c_str()));
    c->Write();
};
std::vector<double> rho_min_list = {0.5, 0.7, 1.0};

DrawRhoMinGraphs(R_map, R_map_err, rho_map,
                 "R [fm]", "R_vs_rhomax",
                 rho_min_list, B_values);

DrawRhoMinGraphs(lambda_map, lambda_map_err, rho_map,
                 "#lambda", "lambda_vs_rhomax",
                 rho_min_list, B_values);
                 
Draw_R_vs_mT_3Pads(R_vs_mT, Rerr_vs_mT, mT_map,
                   rho_min_list, B_values);
                   
Draw_R_vs_mT_ByB(R_vs_mT, Rerr_vs_mT, mT_map, rho_min_list, B_values);

    cout << "\n========================================" << endl;
    cout << "Output files created:" << endl;
    cout << "  - hbt_comparison_pions30fm.png" << endl;
    cout << "  - hbt_comparison_results.root" << endl;
    cout << "========================================" << endl;
    
}

