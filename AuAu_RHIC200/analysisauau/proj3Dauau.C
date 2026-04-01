#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TStyle.h>
#include <RooFit.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <iomanip>
#include "Levy_proj_reader.h"

using namespace std;

Levy_reader* myLevy_reader;

// ===============================
// Particle structure
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
        if (p_abs <= fabs(pz)) return (pz>0?999.0:-999.0);
        return 0.5*TMath::Log((p_abs+pz)/(p_abs-pz));
    }
};

// ===============================
// Levy projection
// ===============================
double LevyProj1DFunc(const double *x, const double *par)
{
    double alpha = par[0];
    double R     = par[1];
    double N     = par[2];

    double Rcc = R * pow(2., 1./alpha);

    return (2.*N / Rcc) * myLevy_reader->getValue_1d(alpha, x[0]/Rcc);
}

// ===============================
// OSCAR file loader
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
    double pT_max = 1.0;      
    double eta_cut = 1.0;     



// ===============================
// HBTResults structure
// ===============================
struct HBTResults {
    TH1F* hOut;
    TH1F* hSide;
    TH1F* hLong;
    TH1F* hRho;
    int total_pairs; 
};
// ===============================
// Analyze pairs
// ===============================
HBTResults AnalyzePairs(const vector<Particle>& particles,
                        double kT_min,
                        double kT_max,
                        bool pions_only)
{
    map<int, vector<Particle>> events;
    for(auto& p : particles){
        if(pions_only && abs(p.pid)!=211) continue;
        events[p.event_id].push_back(p);
    }

   // =======================
// Histogram binning
// =======================
const int n_bins = 150; 
    double rho_min = 0.7;
    double rho_max = 40.0;
    double bins[n_bins + 1];
    
    double r = pow(rho_max / rho_min, 1.0 / n_bins);
    bins[0] = rho_min;
    for (int i = 1; i <= n_bins; i++) {
        bins[i] = bins[i-1] * r;
    }
  
TH1F* hRho = new TH1F("hRho", "D(#rho)", n_bins, bins);
hRho->Sumw2();

TH2F* hRhoVsOut  = new TH2F("hRhoVsOut",  "", n_bins, bins, n_bins, bins);
TH2F* hRhoVsSide = new TH2F("hRhoVsSide", "", n_bins, bins, n_bins, bins);
TH2F* hRhoVsLong = new TH2F("hRhoVsLong", "", n_bins, bins, n_bins, bins);

    hRhoVsOut->Sumw2();
    hRhoVsSide->Sumw2();
    hRhoVsLong->Sumw2();
    
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
                     hRhoVsOut->Fill(fabs(r_out), rho);
                     hRhoVsSide->Fill(fabs(r_side), rho);
                   hRhoVsLong->Fill(fabs(r_long), rho);
            }
        }
    }

auto normalize1D = [](TH1F* h) {
    double total = h->Integral(1, h->GetNbinsX()+1); 

    h->Scale(1.0 / total, "width"); 

    for(int i=1; i<=h->GetNbinsX(); i++) {
        double content = h->GetBinContent(i);
        double error   = h->GetBinError(i);
    }
};


    TH1F* hOut  = (TH1F*) hRhoVsOut->ProjectionX("hOut");
    TH1F* hSide = (TH1F*) hRhoVsSide->ProjectionX("hSide");
    TH1F* hLong = (TH1F*) hRhoVsLong->ProjectionX("hLong");


    HBTResults res;
    res.hOut = hOut; res.hSide = hSide; res.hLong = hLong; res.hRho = hRho;
    res.total_pairs = pairs_count;

    normalize1D(hOut);
    normalize1D(hSide);
    normalize1D(hLong);

    
    return res;

   
}

// ---------------------------
// Main 3D Levy Analysis Function
// ---------------------------
void auau3D()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1::AddDirectory(kFALSE);

    const char* oscar_file = "/home/zeinab/Documents/vhlle-smash/hybrid/sampler.out/test_run/particle_lists.oscar";

    // --- Load particles ---
    vector<Particle> particles = LoadOSCAR(oscar_file);
    if(particles.empty()){
        cout << "No particles loaded!" << endl;
        return;
    }

    // --- Initialize Levy reader ---
    myLevy_reader = new Levy_reader("levy_proj3D_values.dat");

    // --- kT bins ---
    vector<double> kT_bins;

double kT_min_val = 0.15;
double kT_max_val = 0.55;
double step = 0.05;

for(double k = kT_min_val; k <= kT_max_val + 1e-6; k += step)
{
    kT_bins.push_back(k);
}

int nBins = kT_bins.size()-1;

    const int NPAR = 7;

    // --- Loop over kT bins ---
    for(int ibin = 0; ibin < nBins; ibin++)
    {
        double kT_min = kT_bins[ibin];
        double kT_max = kT_bins[ibin+1];

        cout << "\n=== kT bin: " << kT_min << " - " << kT_max << " ===" << endl;

        // --- Analyze pairs ---
        HBTResults res = AnalyzePairs(particles, kT_min, kT_max, true);

        // --- Fit range ---
        const double fit_min = 0.7;
        const double fit_max = 40.0;

        // --- Chi2 functor ---
       struct Chi2Functor {
    HBTResults res;
    double fit_min, fit_max; 

    // constructor to take 3 arguments
    Chi2Functor(const HBTResults& r, double min, double max) 
        : res(r), fit_min(min), fit_max(max) {}

    double operator()(const double* p) const {
        double totalChi2 = 0;
        TH1F* hists[3] = {res.hOut, res.hSide, res.hLong};
        
        double max_fit_x[3] = {fit_max, 30.0, fit_max}; 

        for(int i=0; i<3; i++) {
            double alpha = p[0];
            double R     = p[i+1];
            double N     = p[i+4];
            double pars[3] = {alpha, R, N};

            for(int b=1; b <= hists[i]->GetNbinsX(); b++) {
                double x = hists[i]->GetBinCenter(b);
                
                // Use the specific cutoff
                if(x < fit_min || x > max_fit_x[i]) continue;
                
                double obs = hists[i]->GetBinContent(b);
                double err = hists[i]->GetBinError(b);
                if(err <= 0) continue;

                double exp = LevyProj1DFunc(&x, pars);
                totalChi2 += pow((obs - exp)/err, 2);
            }
        }
        return totalChi2;
    }
};


        Chi2Functor chi2func(res, fit_min, fit_max);
        ROOT::Math::Functor f(chi2func, NPAR);

        ROOT::Math::Minimizer* minimizer =
       ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

        minimizer->SetFunction(f);
        minimizer->SetMaxFunctionCalls(50000);
        minimizer->SetMaxIterations(50000);
        minimizer->SetTolerance(1e-6);

        // --- Initial parameters ---
minimizer->SetLimitedVariable(0, "alpha", 1.4, 0.05, 0.8, 2.0); 


minimizer->SetLimitedVariable(1, "Rout",  4.0, 0.1, 1.0, 20.0);
minimizer->SetLimitedVariable(2, "Rside", 10.0, 0.1, 1.0, 30.0);
minimizer->SetLimitedVariable(3, "Rlong", 6.0, 0.1, 1.0, 20.0);

// separate normalizations
minimizer->SetLimitedVariable(4, "N_out",  0.9, 0.01, 0.0, 5.0);
minimizer->SetLimitedVariable(5, "N_side", 0.5, 0.01, 0.0, 5.0);
minimizer->SetLimitedVariable(6, "N_long", 0.9, 0.01, 0.0, 5.0);

        minimizer->Minimize();

        const double* p = minimizer->X();
        const double* err = minimizer->Errors();
        double chi2 = minimizer->MinValue();
        int nPoints = res.hOut->GetNbinsX() + res.hSide->GetNbinsX() + res.hLong->GetNbinsX();
        int nPar = 7;
       int ndf = nPoints - nPar;
         cout << " chi2/NDF = " << chi2 << "/" << ndf << endl;

        double CL = (ndf>0) ? TMath::Prob(chi2, ndf) : 0;

        // --- Print results ---
        bool success = minimizer->Status() == 0;
        if (!success) cout << "WARNING: Fit did not converge!" << endl;

        cout << "Fit results:\n";
        cout << " alpha = " << p[0] << " ± " << err[0] << endl;
        cout << " Rout  = " << p[1] << " ± " << err[1] << endl;
        cout << " Rside = " << p[2] << " ± " << err[2] << endl;
        cout << " Rlong = " << p[3] << " ± " << err[3] << endl;
cout << " N_out  = " << p[4] << " ± " << err[4] << endl;
cout << " N_side = " << p[5] << " ± " << err[5] << endl;
cout << " N_long = " << p[6] << " ± " << err[6] << endl;
        cout << " chi²/NDF = " << chi2 << "/" << ndf << " (C.L. = " << CL*100 << "%)" << endl;
        
         //check points
      cout << "Pairs in Out: " << res.hOut->GetEntries() << endl;
      cout << "Pairs in Side: " << res.hSide->GetEntries() << endl;
      cout << "Pairs in Long: " << res.hLong->GetEntries() << endl;


        // --- Plot ---
        TCanvas* c = new TCanvas(Form("c_%d", ibin),
                                 Form("kT: %.2f-%.2f GeV/c, N_{pairs} = %d",
                                      kT_min, kT_max, res.total_pairs),
                                 1800, 600);
        c->Divide(3,1);

        TH1F* hdir[3] = {res.hOut, res.hSide, res.hLong};
        const char* names[3] = {"out", "side", "long"};

        double xmin[3] = {0.7, 0.7, 0.7};
        double xmax[3] = {35.0, 25.0, 35.0};
        double ymin[3] = {1e-5, 1e-5, 1e-5};
        double ymax[3] = {100, 100, 100};

        for (int idir = 0; idir < 3; ++idir) {
            c->cd(idir+1);
            gPad->SetLogx();
            gPad->SetLogy();

            TH1F* h = hdir[idir];
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.8);

            h->GetXaxis()->SetRangeUser(xmin[idir], xmax[idir]);

            // Compute y-axis max
            double maxy = 0;
            for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
                double x = h->GetXaxis()->GetBinCenter(ib);
                if (x < xmin[idir] || x > xmax[idir]) continue;
                double y = h->GetBinContent(ib);
                if (y > maxy) maxy = y;
            }
           
            h->Draw("E");
            
            res.hOut ->GetXaxis()->SetTitle("#rho_{out} [fm]");
            res.hSide->GetXaxis()->SetTitle("#rho_{side} [fm]");
             res.hLong->GetXaxis()->SetTitle("#rho_{long} [fm]");

            h ->GetYaxis()->SetTitle("D(#rho)");

            // Draw Levy fit
            double par[3] = {p[0], p[idir+1], p[idir+4]};
            TF1* flevy = new TF1(Form("f_%s_%d", names[idir], ibin),
                                 LevyProj1DFunc, xmin[idir], xmax[idir], 3);
            flevy->SetParameters(par);
            flevy->SetLineColor(kRed);
            flevy->SetLineWidth(3);
            flevy->Draw("SAME");

            // Parameter box
            TPaveText* box = new TPaveText(0.15,0.55,0.50,0.88,"NDC");
            box->SetFillColor(0);
            box->SetTextSize(0.035);
            box->AddText(Form("#alpha = %.3f #pm %.3f", p[0], err[0]));
            box->AddText(Form("R_{%s} = %.2f #pm %.2f fm", names[idir], p[idir+1], err[idir+1]));
            box->AddText(Form("#lambda_{%s} = %.2f #pm %.2f", names[idir], p[idir+1], err[idir+1]));
            box->AddText(Form("#chi^{2}/NDF = %.2f / %.d", chi2, ndf));
            box->AddText(Form("C.L. = %.2f%%", CL*100));
                       

            box->Draw();
        }

        c->cd(1);
        TPaveText* info = new TPaveText(0.15,0.85,0.45,0.95,"NDC");
        info->AddText(Form("kT: %.2f-%.2f GeV/c", kT_min, kT_max));
        info->AddText(Form("N_{pairs} = %d", res.total_pairs));
  
        info->Draw();

        c->SaveAs(Form("kT_%.2f-%.2f.png", kT_min, kT_max));

        delete minimizer;
    }
    
}
