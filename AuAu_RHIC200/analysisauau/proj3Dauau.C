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
#include <TRandom3.h>
#include "Levy_proj_reader.h"

using namespace std;

Levy_reader* myLevy_reader;

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
// Particle structure
// ===============================
struct Particle {
    Int_t event_id;
    Int_t pid;
    Int_t charge;
    Int_t ID;

    Double_t px, py, pz, E;
    Double_t x, y, z, t;

    Float_t pT() const { return TMath::Sqrt(px*px + py*py); }

    Float_t p() const { return TMath::Sqrt(px*px + py*py + pz*pz); }

    Float_t eta() const {
        Float_t p_abs = p();
        if (p_abs <= fabs(pz)) return (pz > 0 ? 999.0 : -999.0);
        return 0.5 * TMath::Log((p_abs + pz) / (p_abs - pz));
    }
};

//Filling from oscar file 

vector<Particle> LoadOSCAR(const char* fname)
{
    vector<Particle> particles;

    ifstream fin(fname);
    string line;
    int event_id = -1;

    while (getline(fin, line))
    {
        if (line.empty() || line[0]=='#') {
            if (line.find("# event") != string::npos) {
                string tmp;
                stringstream ss(line);
                ss >> tmp >> tmp >> event_id;
            }
            continue;
        }

        stringstream ss(line);

        Particle p;

        double mass, E;
        int pdg, ID, charge, ncoll;
        double form_time;
        double time_last_coll;
        double xsecfac;
        double proc_id_origin;
        double proc_type_origin;
        double pdg_mother1;
        double pdg_mother2;
        double baryon_number;
   

ss >> p.t
   >> p.x >> p.y >> p.z
   >> mass >> E
   >> p.px >> p.py >> p.pz
   >> pdg >> ID >> charge >> ncoll
   >> form_time
   >> xsecfac >> proc_id_origin >> proc_type_origin
   >> time_last_coll
   >> pdg_mother1 >> pdg_mother2 >> baryon_number;

        if (ss.fail()) continue;


        // on-shell SMASH energy correction
        double p2 = p.px*p.px + p.py*p.py + p.pz*p.pz;
        p.E = sqrt(p2 + mass*mass);
        
        // emission time correction
        TRandom3 *randGen = new TRandom3(0); 
       double tau_0 = (time_last_coll > 0) ? time_last_coll : time_last_coll - form_time ; // approx. last spread time 
       //double gamma = p.E / mass; 
       double delta_t = randGen->Exp(fabs(tau_0)); // Exponential decay mimics decoupling

       p.t = form_time + delta_t;

        p.pid = pdg;
        p.ID = ID;
        p.charge = charge;
        p.event_id = event_id;

        particles.push_back(p);
    }

    cout << "Loaded " << particles.size() << " particles\n";

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
const int n_bins = 300; 
    double rho_min = 0.1;
    double rho_max = 1e12;
    double bins[n_bins + 1];
    
    double r = pow(rho_max / rho_min, 1.0 / n_bins);
    bins[0] = rho_min;
    for (int i = 1; i <= n_bins; i++) {
        bins[i] = bins[i-1] * r;
    }
  
TH1F* hRho = new TH1F("hRho", "D(#rho)", n_bins, bins);
hRho->Sumw2();

TH1F* hOut  = new TH1F("hOut",  "", n_bins, bins);
TH1F* hSide = new TH1F("hSide", "", n_bins, bins);
TH1F* hLong = new TH1F("hLong", "", n_bins, bins);


hOut->Sumw2();
hSide->Sumw2();
hLong->Sumw2();
    
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
                 double c = 0.15; 
                 double qmCut = sqrt(c * mT);
           
                  // relative momentum  in the LCMS
                 double qx = p1.px - p2.px;
                 double qy = p1.py - p2.py;
                 double qzL = (4*pow((p1.pz * p2.E - p2.pz * p1.E),2))/(pow((p1.E+p2.E),2)-pow((p1.pz+p2.pz),2));

                 double_t qLCMS = TMath::Sqrt(qx*qx + qy*qy + qzL);

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
                
                /*/ 1. Pair kinematics in Lab Frame
double K0 = 0.5 * (p1.E + p2.E);
double beta = kz / K0;
double gamma = 1.0 / sqrt(1.0 - beta * beta);

// 3. Transform dz and dt to the LCMS frame
// This is the "Longitudinal" boost
double dz_L = gamma * (dz - beta * dt);
double dt_L = gamma * (dt - beta * dz);

// 4. Calculate Bertsch-Pratt coordinates
double phi = atan2(ky, kx);
double beta_T = kT / sqrt(K0*K0 - kz*kz); // Transverse velocity in LCMS

// r_side: purely geometric
double r_side = -TMath::Sin(phi) * dx + TMath::Cos(phi) * dy;

// r_out: includes the emission duration correction
double r_out_spatial = TMath::Cos(phi) * dx + TMath::Sin(phi) * dy;
double r_out = r_out_spatial - beta_T * dt_L;

// r_long: is simply the boosted dz
                double r_long = (K0*dz - kz*dt) / sqrt(K0*K0 - kz*kz);
//double r_long = dz_L;

double rho = sqrt(r_out * r_out + r_side * r_side + r_long * r_long);*/

                // LCMS transformation
                double_t K0 = 0.5 * (p1.E + p2.E); //pair energy
                double kp = K0*K0 - kz*kz;
                 double phi = atan2(ky, kx); // azimuthal angle
                
                // Bertsch-Pratt coordinates
                double r_out  = TMath::Cos(phi) * dx + TMath::Sin(phi) * dy - (kT/kp) * (K0*dt - kz*dz);
                double r_side = - TMath::Sin(phi) * dx + TMath::Cos(phi) * dy;
                double r_long = (K0*dz - kz*dt) / sqrt(kp);
                
                double rho = sqrt(r_out * r_out + r_side * r_side + r_long * r_long);
               
                   hRho->Fill(rho);
                   hOut->Fill(fabs(r_out));
                   hSide->Fill(fabs(r_side));
                   hLong->Fill(fabs(r_long));
            }
        }
    }

    HBTResults res;
    res.hOut = hOut; res.hSide = hSide; res.hLong = hLong; res.hRho = hRho;
    res.total_pairs = pairs_count;

    
    return res;
  
}

// ---------------------------
// Main Function
// ---------------------------
void auau3D()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TH1::AddDirectory(kFALSE);

    const char* oscar_file = "/home/zeinab/Documents/vhlle-smash/hybrid/AuAu_RHIC200/smash.out/cent0_5/particle_lists.oscar";
    

    // --- Load particles ---
    vector<Particle> particles = LoadOSCAR(oscar_file);
    if(particles.empty()){
        cout << "No particles loaded!" << endl;
        return;
    }
    
    //check particls at freeze out
   // 1. Declare histograms
TH1F* hEta_eta    = new TH1F("hEta_eta",   "#eta (eta meson)", 100, -5, 5);
TH1F* hEta_etap   = new TH1F("hEta_etap",  "#eta (eta')",      100, -5, 5);
TH1F* hEta_lambda = new TH1F("hEta_lambda","#eta (Lambda)",    100, -5, 5);
TH1F* hEta_k0s    = new TH1F("hEta_k0s",   "#eta (K0S)",       100, -5, 5);


// 2. Declare counters
int n_eta=0, n_etap=0, n_lambda=0, n_k0s=0;

// 3. Loop over particles and fill histograms
for(const auto& p : particles){
    if(p.pid==221) { n_eta++; hEta_eta->Fill(p.eta()); }
    if(p.pid==331) { n_etap++; hEta_etap->Fill(p.eta()); }
    if(p.pid==3122){ n_lambda++; hEta_lambda->Fill(p.eta()); }
    if(p.pid==310) { n_k0s++; hEta_k0s->Fill(p.eta()); }
}

// 4. Print counts
cout << "eta: " << n_eta << endl;
cout << "eta': " << n_etap << endl;
cout << "Lambda: " << n_lambda << endl;
cout << "K0S: " << n_k0s << endl;

// 5. Draw histograms in separate pads
// 5. Draw histograms in separate pads with labels and particle counts
TCanvas* cEta = new TCanvas("cEta","Particle #eta distributions",1000,800);
cEta->Divide(2,2); // 2x2 grid

// Function to create a TPaveText for particle counts
auto drawCountBox = [](int count) {
    TPaveText* box = new TPaveText(0.15,0.75,0.4,0.88,"NDC");
    box->SetFillColorAlpha(kWhite,0.0); // transparent
    box->SetBorderSize(0);
    box->SetTextAlign(12); // left/top
    box->AddText(Form("N = %d", count));
    return box;
};

// Pad 1: eta
cEta->cd(1);
hEta_eta->SetLineColor(kRed);
hEta_eta->GetXaxis()->SetTitle("#eta");
hEta_eta->GetYaxis()->SetTitle("Counts");
hEta_eta->Draw("HIST");
TPaveText* box1 = drawCountBox(n_eta);
box1->Draw();
TLegend* leg1 = new TLegend(0.65,0.75,0.88,0.88);
leg1->SetBorderSize(0);
leg1->SetFillStyle(0);
leg1->AddEntry(hEta_eta,"#eta","L");
leg1->Draw();

// Pad 2: eta'
cEta->cd(2);
hEta_etap->SetLineColor(kBlue);
hEta_etap->GetXaxis()->SetTitle("#eta");
hEta_etap->GetYaxis()->SetTitle("Counts");
hEta_etap->Draw("HIST");
TPaveText* box2 = drawCountBox(n_etap);
box2->Draw();
TLegend* leg2 = new TLegend(0.65,0.75,0.88,0.88);
leg2->SetBorderSize(0);
leg2->SetFillStyle(0);
leg2->AddEntry(hEta_etap,"#eta'","L");
leg2->Draw();

// Pad 3: Lambda
cEta->cd(3);
hEta_lambda->SetLineColor(kGreen+2);
hEta_lambda->GetXaxis()->SetTitle("#eta");
hEta_lambda->GetYaxis()->SetTitle("Counts");
hEta_lambda->Draw("HIST");
TPaveText* box3 = drawCountBox(n_lambda);
box3->Draw();
TLegend* leg3 = new TLegend(0.65,0.75,0.88,0.88);
leg3->SetBorderSize(0);
leg3->SetFillStyle(0);
leg3->AddEntry(hEta_lambda,"#Lambda","L");
leg3->Draw();

// Pad 4: K0S
cEta->cd(4);
hEta_k0s->SetLineColor(kMagenta);
hEta_k0s->GetXaxis()->SetTitle("#eta");
hEta_k0s->GetYaxis()->SetTitle("Counts");
hEta_k0s->Draw("HIST");
TPaveText* box4 = drawCountBox(n_k0s);
box4->Draw();
TLegend* leg4 = new TLegend(0.65,0.75,0.88,0.88);
leg4->SetBorderSize(0);
leg4->SetFillStyle(0);
leg4->AddEntry(hEta_k0s,"K^{0}_{S}","L");
leg4->Draw();

// Update canvas
cEta->Update();
// ==========================================
// Freeze-out time distributions (multi-species)
// ==========================================

// Pions
TH1F* hT_piplus  = new TH1F("hT_piplus",  "t #pi^{+}", 200, 0, 150);
TH1F* hT_piminus = new TH1F("hT_piminus", "t #pi^{-}", 200, 0, 150);
TH1F* hT_all     = new TH1F("hT_all",     "t all #pi", 200, 0, 150);

// Other species
TH1F* hT_eta     = new TH1F("hT_eta",     "t #eta",    200, 0, 150);
TH1F* hT_etap    = new TH1F("hT_etap",    "t #eta'",   200, 0, 150);
TH1F* hT_lambda  = new TH1F("hT_lambda",  "t #Lambda", 200, 0, 150);

// Fill histograms
for(const auto& p : particles)
{
    // --- Pions ---
    if(abs(p.pid) == 211){
        hT_all->Fill(p.t);
        if(p.pid == 211)  hT_piplus->Fill(p.t);
        if(p.pid == -211) hT_piminus->Fill(p.t);
    }

    // --- Other particles ---
    if(p.pid == 221)  hT_eta->Fill(p.t);
    if(p.pid == 331)  hT_etap->Fill(p.t);
    if(p.pid == 3122) hT_lambda->Fill(p.t);
}

// ==========================================
// Canvas
// ==========================================
TCanvas* cTime = new TCanvas("cTime", "Freeze-out time comparison", 1400, 1000);
cTime->Divide(2,2);

// ---  pi+ vs pi-
cTime->cd(1);
//gPad->SetLogy();
hT_piplus->GetYaxis()->SetTitle("Counts");
hT_piplus->SetLineColor(kRed);
hT_piminus->SetLineColor(kBlue);
hT_piminus->GetXaxis()->SetTitle("t [fm/c]");
hT_piplus->Draw("HIST");
hT_piminus->Draw("HIST SAME");

TLegend* legPi = new TLegend(0.6,0.7,0.88,0.88);
legPi->SetBorderSize(0);
legPi->SetFillStyle(0);
legPi->AddEntry(hT_piplus,"#pi^{+}","L");
legPi->AddEntry(hT_piminus,"#pi^{-}","L");
legPi->Draw();

// ---  eta
cTime->cd(2);
//gPad->SetLogy();
hT_eta->SetLineColor(kMagenta);
hT_eta->GetYaxis()->SetTitle("Counts");
hT_eta->GetXaxis()->SetTitle("t [fm/c]");
hT_eta->Draw("HIST");

TLegend* legP = new TLegend(0.6,0.7,0.88,0.88);
legP->SetBorderSize(0);
legP->SetFillStyle(0);
legP->AddEntry(hT_eta,"#eta","L");
legP->Draw();

// ---  eta meson
cTime->cd(3);
hT_etap->SetLineColor(kGreen+2);
hT_etap->GetYaxis()->SetTitle("Counts");
hT_etap->GetXaxis()->SetTitle("t [fm/c]");
hT_etap->Draw("HIST");

TLegend* legp = new TLegend(0.6,0.7,0.88,0.88);
legp->SetBorderSize(0);
legp->SetFillStyle(0);
legp->AddEntry(hT_etap,"#eta'","L");
legp->Draw();

// ---  Lambda
cTime->cd(4);
hT_lambda->GetYaxis()->SetTitle("Counts");
hT_lambda->SetLineColor(kOrange+7);
hT_lambda->GetXaxis()->SetTitle("t [fm/c]");
hT_lambda->Draw("HIST");

TLegend* legl = new TLegend(0.6,0.7,0.88,0.88);
legl->SetBorderSize(0);
legl->SetFillStyle(0);
legl->AddEntry(hT_lambda,"#Lambda","L");
legl->Draw();

cTime->SaveAs("Freezeout_time_all_species.png");

// ==========================================
// Momentum distributions for IDENTICAL pions
// ==========================================

// pi+
TH1F* hPx_piplus = new TH1F("hPx_piplus","p_{x} #pi^{+}",200,-2,2);
TH1F* hPy_piplus = new TH1F("hPy_piplus","p_{y} #pi^{+}",200,-2,2);
TH1F* hPz_piplus = new TH1F("hPz_piplus","p_{z} #pi^{+}",200,-10,10);

// pi-
TH1F* hPx_piminus = new TH1F("hPx_piminus","p_{x} #pi^{-}",200,-2,2);
TH1F* hPy_piminus = new TH1F("hPy_piminus","p_{y} #pi^{-}",200,-2,2);
TH1F* hPz_piminus = new TH1F("hPz_piminus","p_{z} #pi^{-}",200,-10,10);

for(const auto& p : particles)
{
    if(p.pid == 211){ // pi+
        hPx_piplus->Fill(p.px);
        hPy_piplus->Fill(p.py);
        hPz_piplus->Fill(p.pz);
    }

    if(p.pid == -211){ // pi-
        hPx_piminus->Fill(p.px);
        hPy_piminus->Fill(p.py);
        hPz_piminus->Fill(p.pz);
    }
}

TCanvas* cMomPi = new TCanvas("cMomPi","Momentum components of pions",1200,400);
cMomPi->Divide(3,1);
//p_x
cMomPi->cd(1);

hPx_piplus->SetLineColor(kRed);
hPx_piminus->SetLineColor(kBlue);

hPx_piplus->GetXaxis()->SetTitle("p_{x} [GeV/c]");
hPx_piplus->GetYaxis()->SetTitle("Counts");

hPx_piplus->Draw("HIST");
hPx_piminus->Draw("HIST SAME");

TLegend* legPx = new TLegend(0.65,0.75,0.88,0.88);
legPx->SetBorderSize(0);
legPx->SetFillStyle(0);
legPx->AddEntry(hPx_piplus,"#pi^{+}","L");
legPx->AddEntry(hPx_piminus,"#pi^{-}","L");
legPx->Draw();

//p_y
cMomPi->cd(2);

hPy_piplus->SetLineColor(kRed);
hPy_piminus->SetLineColor(kBlue);

hPy_piplus->GetXaxis()->SetTitle("p_{y} [GeV/c]");

hPy_piplus->Draw("HIST");
hPy_piminus->Draw("HIST SAME");

TLegend* legPy = new TLegend(0.65,0.75,0.88,0.88);
legPy->SetBorderSize(0);
legPy->SetFillStyle(0);
legPy->AddEntry(hPy_piplus,"#pi^{+}","L");
legPy->AddEntry(hPy_piminus,"#pi^{-}","L");
legPy->Draw();

//p_z
cMomPi->cd(3);

hPz_piplus->SetLineColor(kRed);
hPz_piminus->SetLineColor(kBlue);

hPz_piplus->GetXaxis()->SetTitle("p_{z} [GeV/c]");

hPz_piplus->Draw("HIST");
hPz_piminus->Draw("HIST SAME");

TLegend* legPz = new TLegend(0.65,0.75,0.88,0.88);
legPz->SetBorderSize(0);
legPz->SetFillStyle(0);
legPz->AddEntry(hPz_piplus,"#pi^{+}","L");
legPz->AddEntry(hPz_piminus,"#pi^{-}","L");
legPz->Draw();

cMomPi->SaveAs("Momentum_components_pions.png");


//========================================================================
//ANALYSIS PAIRS section
//=======================================================================
// --- Initialize Levy reader ---
    myLevy_reader = new Levy_reader("levy_proj3D_values.dat");

    // --- kT bins ---
    vector<double> kT_bins;

double kT_min_val = 0.15;
double kT_max_val = 0.85;
double step = 0.05;

for(double k = kT_min_val; k <= kT_max_val + 1e-6; k += step)
{
    kT_bins.push_back(k);
}

int nBins = kT_bins.size()-1;

vector<double> mt_vals, mt_err;
vector<double> alpha_vals, alpha_err;
vector<double> Rout_vals, Rout_err;
vector<double> Rside_vals, Rside_err;
vector<double> Rlong_vals, Rlong_err;
vector<double> lambda_vals, lambda_err;  


    // --- Loop over kT bins ---
    for(int ibin = 0; ibin < nBins; ibin++)
    {
        double kT_min = kT_bins[ibin];
        double kT_max = kT_bins[ibin+1];

        cout << "\n=== kT bin: " << kT_min << " - " << kT_max << " ===" << endl;

        // --- Analyze pairs ---
        HBTResults res = AnalyzePairs(particles, kT_min, kT_max, true);

        // --- Fit range ---
        const double fit_min = 0.1;
        double fit_max_x[3] = {100.0, 100.0, 100.0}; // out, side, long
        const int NPAR = 5;

struct LogLikelihood {
    HBTResults res;
    double fit_min, fit_max_x[3];
    int* pBinsUsed; 

    LogLikelihood(const HBTResults& r, double min, const double max_vals[3], int* binPtr)
        : res(r), fit_min(min), pBinsUsed(binPtr) {
        for(int i=0; i<3; i++) fit_max_x[i] = max_vals[i];
    }
  
    double operator()(const double* p) const {
        double logL = 0.0;
        int currentBins = 0; 
        
        TH1F* hists[3] = {res.hOut, res.hSide, res.hLong};

        for(int i=0; i<3; i++) {
            double pars[3] = {p[0], p[i+1], p[4]};
            double integral = hists[i]->Integral(0, hists[i]->GetNbinsX()+1); //NORMALIZATION 

            for(int b=1; b <= hists[i]->GetNbinsX(); b++) {
                double x = hists[i]->GetBinCenter(b);
                if(x < fit_min || x > fit_max_x[i]) continue;

                double expected = LevyProj1DFunc(&x, pars) * hists[i]->GetBinWidth(b) * integral;
                if(expected <= 1e-12) continue;

                double observed = hists[i]->GetBinContent(b);
                logL += (observed > 0) ? (expected - observed + observed * log(observed/expected)) : expected;
                currentBins++;
            }
        }
        
        if (pBinsUsed) *pBinsUsed = currentBins; 
        return logL;
    }
};

//counting used bins 
int actualBins = 0; 
LogLikelihood loglikfunc(res, fit_min, fit_max_x, &actualBins);
ROOT::Math::Functor f(loglikfunc, NPAR);

       ROOT::Math::Minimizer* minimizer =
       ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

        minimizer->SetFunction(f);
        minimizer->SetMaxFunctionCalls(50000);
        minimizer->SetMaxIterations(50000);
        minimizer->SetTolerance(1e-6);

        // --- fitting parameters ---
minimizer->SetLimitedVariable(0, "alpha", 1.5, 0.01, 0.5, 2.0);
minimizer->SetLimitedVariable(1, "Rout",  4.0, 0.1, 1.0, 30.0);
minimizer->SetLimitedVariable(2, "Rside", 2.0, 0.1, 1.0, 30.0);
minimizer->SetLimitedVariable(3, "Rlong", 4.0, 0.1, 1.0, 30.0);
minimizer->SetLimitedVariable(4, "N",1.0, 0.01, 0.0, 2.0);

        minimizer->Minimize();

        const double* p = minimizer->X();
        const double* err = minimizer->Errors();
       double chi2 =  minimizer->MinValue(); 
       int ndf = actualBins - NPAR; 
       double CL = (ndf > 0) ? TMath::Prob(chi2, ndf) : 0;
       
       // --- Compute mean kT and mT ---
double kT_mean = 0.5 * (kT_min + kT_max);
const double m_pi = 0.13957; 
double mT = sqrt(kT_mean * kT_mean + m_pi * m_pi);

// Store
mt_vals.push_back(mT);
mt_err.push_back(0.0);   // no error on mT (you can set to bin half‑width if needed)

alpha_vals.push_back(p[0]);
alpha_err.push_back(err[0]);

Rout_vals.push_back(p[1]);
Rout_err.push_back(err[1]);

Rside_vals.push_back(p[2]);
Rside_err.push_back(err[2]);

Rlong_vals.push_back(p[3]);
Rlong_err.push_back(err[3]);

lambda_vals.push_back(p[4]);
lambda_err.push_back(err[4]);


        // --- Print fit results ---
        bool success = minimizer->Status() == 0;
        if (!success) cout << "WARNING: Fit did not converge!" << endl;

        cout << "Fit results:\n";
        cout << " alpha = " << p[0] << " ± " << err[0] << endl;
        cout << " Rout  = " << p[1] << " ± " << err[1] << endl;
        cout << " Rside = " << p[2] << " ± " << err[2] << endl;
        cout << " Rlong = " << p[3] << " ± " << err[3] << endl;
        cout << " N  = " << p[4] << " ± " << err[4] << endl;
       cout << "-logL/NDF: " << chi2 << "/" << ndf << " (CL: " << CL * 100 << ")" << endl;
        
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

        for (int idir = 0; idir < 3; ++idir) {
               c->cd(idir+1);
               gPad->SetLogx();
               gPad->SetLogy();

               TH1F* h = hdir[idir];
               
          // scalling 
         h->Scale(1.0 / h->Integral(), "width");
          
           // Set the Range and Axis limits
           double xmin[3] = {0.1 , 0.1, 0.1};
           double xmax[3] = {1e12 , 1e12, 1e12};
            h->GetXaxis()->SetRangeUser(xmin[idir], xmax[idir]);
            h->SetMinimum(1e-10);  
            h->SetMaximum(1.0);  
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.8);
            
            h ->GetYaxis()->SetTitle("D(#rho)");
            res.hOut ->GetXaxis()->SetTitle("#rho_{out} [fm]");
            res.hSide->GetXaxis()->SetTitle("#rho_{side} [fm]");
            res.hLong->GetXaxis()->SetTitle("#rho_{long} [fm]");
            
            h->Draw("E");
  
         //Drawing fit function
          double R_par[3] = {p[1], p[2], p[3]}; // Rout, Rside, Rlong
          double N = p[4]; 
          double par[3] = {p[0], R_par[idir], N };

            TF1* flevy = new TF1(Form("f_%s_%d", names[idir], ibin),
                                 LevyProj1DFunc, xmin[idir], xmax[idir], 3);
            flevy->SetParameters(par);
            flevy->SetLineColor(kRed);
            flevy->SetLineWidth(3);
            flevy->Draw("SAME");

            // fit Parameter box
            TPaveText* box = new TPaveText(0.35, 0.15, 0.65, 0.40, "NDC");

box->SetFillStyle(0);     // transparent
box->SetBorderSize(0);    // no border
box->SetTextSize(0.035);
box->SetTextAlign(22);    // center text

box->AddText(Form("#alpha = %.3f #pm %.3f", p[0], err[0]));
box->AddText(Form("R_{%s} = %.2f #pm %.2f fm", names[idir], p[idir+1], err[idir+1]));
box->AddText(Form("#lambda = %.2f #pm %.2f", p[4], err[4]));
box->AddText(Form("#chi^{2}/NDF = %.2f / %d", chi2, ndf));
box->AddText(Form("C.L. = %.2f%%", CL*100));

box->Draw();
        }

        c->cd(1);
        TPaveText* info = new TPaveText(0.65, 0.80, 0.95, 0.95, "NDC");

info->SetFillStyle(0);     // transparent
info->SetBorderSize(0);    // no border
info->SetTextSize(0.04);
info->SetTextAlign(32);    // right-aligned text

info->AddText(Form("kT: %.2f-%.2f GeV/c", kT_min, kT_max));
info->AddText(Form("N_{pairs} = %d", res.total_pairs));

info->Draw();
        c->SaveAs(Form("kT_%.2f-%.2f.png", kT_min, kT_max));

        delete minimizer;
    }
    // --- New canvas: mT dependence of all fit parameters ---
TCanvas *c_mT = new TCanvas("c_mT", "Fit Parameters vs m_T", 1600, 1200);
c_mT->Divide(2, 2);

int nPoints = mt_vals.size();

// ========== Pad 1 : alpha ==========
c_mT->cd(1);
TGraphErrors *gAlpha = new TGraphErrors(nPoints,
    &mt_vals[0], &alpha_vals[0],
    &mt_err[0],  &alpha_err[0]);
gAlpha->SetTitle("#alpha vs m_{T}; m_{T} [GeV]; #alpha");
gAlpha->SetMarkerStyle(20);
gAlpha->SetMarkerColor(kBlack);
gAlpha->GetHistogram()->GetYaxis()->SetRangeUser(1.5, 2.5);
gAlpha->Draw("AP");


// ========== Pad 2 : Rout, Rside, Rlong together ==========
c_mT->cd(2);

TGraphErrors *gRout = new TGraphErrors(nPoints,
    &mt_vals[0], &Rout_vals[0],
    &mt_err[0],  &Rout_err[0]);
gRout->SetMarkerStyle(21);
gRout->SetMarkerColor(kRed);
gRout->SetLineColor(kRed);
gRout->SetTitle("HBT radii vs m_{T}; m_{T} [GeV]; R [fm]");

TGraphErrors *gRside = new TGraphErrors(nPoints,
    &mt_vals[0], &Rside_vals[0],
    &mt_err[0],  &Rside_err[0]);
gRside->SetMarkerStyle(22);
gRside->SetMarkerColor(kBlue);
gRside->SetLineColor(kBlue);


TGraphErrors *gRlong = new TGraphErrors(nPoints,
    &mt_vals[0], &Rlong_vals[0],
    &mt_err[0],  &Rlong_err[0]);
gRlong->SetMarkerStyle(23);
gRlong->SetMarkerColor(kGreen+2);
gRlong->SetLineColor(kGreen+2);



gPad->Update();

// Set range
gRout->GetHistogram()->GetYaxis()->SetRangeUser(4.0, 35.0);
// Draw first graph to set axes, then overlay the others
gRout->Draw("AP");
gRside->Draw("P SAME");
gRlong->Draw("P SAME");

// Legend
TLegend *legRadii = new TLegend(0.65, 0.70, 0.88, 0.88);
legRadii->SetBorderSize(0);
legRadii->SetFillStyle(0);
legRadii->AddEntry(gRout, "R_{out}", "P");
legRadii->AddEntry(gRside, "R_{side}", "P");
legRadii->AddEntry(gRlong, "R_{long}", "P");
legRadii->Draw();

// ========== Pad 3 : lambda (N) ==========
c_mT->cd(3);
TGraphErrors *gLambda = new TGraphErrors(nPoints,
    &mt_vals[0], &lambda_vals[0],
    &mt_err[0],  &lambda_err[0]);
gLambda->SetTitle("#lambda vs m_{T}; m_{T} [GeV]; #lambda");
gLambda->SetMarkerStyle(20);
gLambda->SetMarkerColor(kBlack);
gLambda->Draw("AP");

// ========== Pad 4 : system information ==========
c_mT->cd(4);
gPad->SetFillColor(0);
gPad->SetFrameBorderMode(0);

TPaveText *sysInfo = new TPaveText(0.10, 0.10, 0.90, 0.90, "NDC");
sysInfo->SetTextSize(0.045);
sysInfo->SetFillColor(0);
sysInfo->SetBorderSize(1);
sysInfo->SetTextFont(42);
sysInfo->SetTextAlign(12);

sysInfo->AddText("System Information");
sysInfo->AddText("---------------------------");
sysInfo->AddText(Form("%s", collision_system));
sysInfo->AddText(Form("#sqrt{s_{NN}} = %.0f GeV", sqrt_sNN));
sysInfo->AddText(Form("Centrality: %s", centrality));
sysInfo->AddText(" ");
sysInfo->AddText("Kinematic cuts:");
sysInfo->AddText(Form("%.2f < p_{T} < %.2f GeV/c", pT_min, pT_max));
sysInfo->AddText(Form("|#eta| < %.1f", eta_cut));
sysInfo->AddText(Form("k_{T} range: %.2f - %.2f GeV/c",
                kT_min_val, kT_max_val));
sysInfo->Draw();

// Optional: save the canvas
c_mT->SaveAs("Proj3D_FitParameters_vs_mT.png");
c_mT->Write();   // if you have a ROOT file open
    
}
