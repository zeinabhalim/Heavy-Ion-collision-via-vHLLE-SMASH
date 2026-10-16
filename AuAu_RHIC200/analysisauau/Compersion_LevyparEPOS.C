// ============================================================
// compare_SMASH_EPOS_all.C
// Single 3x2 canvas comparing the user's own wide-format
// analysis CSV with the paper's long-format supplementary CSV.
//   (a) Rout    (b) Rside   (c) Rlong
//   (d) alpha   (e) lambda  (f) system info
// Run:  root -l -b -q compare_SMASH_EPOS_all.C
// ============================================================

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TH1.h>
#include <TAxis.h>
#include <TPad.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

// ------------------------------------------------------------
// Record
// ------------------------------------------------------------
struct CSVPoint { double x, y, ey_lo, ey_hi; };

// ------------------------------------------------------------
// Splitter: comma, tab, space, CR, LF, semicolon are separators
// ------------------------------------------------------------
vector<string> SplitAny(const string& line)
{
    vector<string> out;
    string cur;
    for (char c : line) {
        if (c == ',' || c == '\t' || c == ' ' ||
            c == '\r' || c == '\n' || c == ';') {
            if (!cur.empty()) { out.push_back(cur); cur.clear(); }
        } else cur.push_back(c);
    }
    if (!cur.empty()) out.push_back(cur);
    return out;
}

// ------------------------------------------------------------
// Long-format loader (paper supplementary CSV)
//   plot_type,centrality,source,x,y,y_err_low,y_err_high
// ------------------------------------------------------------
vector<CSVPoint> LoadLong(const string& filename,
                          const string& targetPlot,
                          const string& targetCent,
                          const string& targetSource)
{
    vector<CSVPoint> pts;
    ifstream f(filename);
    if (!f.is_open()) { cerr << "Cannot open " << filename << endl; return pts; }

    string line;
    getline(f, line);   // skip header

    while (getline(f, line)) {
        if (line.empty()) continue;
        vector<string> F = SplitAny(line);
        if (F.size() < 7) continue;

        if (F[0] != targetPlot)   continue;
        if (F[1] != targetCent)   continue;
        if (F[2] != targetSource) continue;

        try {
            CSVPoint p;
            p.x     = stod(F[3]);
            p.y     = stod(F[4]);
            p.ey_lo = fabs(stod(F[5]));
            p.ey_hi = fabs(stod(F[6]));
            pts.push_back(p);
        } catch (...) {}
    }
    return pts;
}

vector<CSVPoint> LoadLongAnyName(const string& filename,
                                 const vector<string>& candidates,
                                 const string& targetCent,
                                 const string& targetSource)
{
    for (const auto& name : candidates) {
        vector<CSVPoint> pts = LoadLong(filename, name, targetCent, targetSource);
        if (!pts.empty()) return pts;
    }
    return {};
}

// ------------------------------------------------------------
// Wide-format loader (user's own analysis CSV)
//   columns read BY NAME: mT plus the requested quantity
// ------------------------------------------------------------
vector<CSVPoint> LoadWide(const string& filename,
                          const string& quantity,
                          const string& lambdaCutoff = "1mm")
{
    vector<CSVPoint> pts;
    ifstream f(filename);
    if (!f.is_open()) { cerr << "Cannot open " << filename << endl; return pts; }

    string headerLine;
    if (!getline(f, headerLine)) return pts;

    vector<string> cols = SplitAny(headerLine);
    if (!cols.empty() && cols[0].size() >= 3 &&
        (unsigned char)cols[0][0] == 0xEF &&
        (unsigned char)cols[0][1] == 0xBB &&
        (unsigned char)cols[0][2] == 0xBF)
        cols[0] = cols[0].substr(3);

    auto findCol = [&](const string& name) -> int {
        for (size_t i = 0; i < cols.size(); ++i)
            if (cols[i] == name) return (int)i;
        return -1;
    };

    string valCol, errCol;
    if      (quantity == "alpha")  { valCol = "alpha";  errCol = "alpha_err"; }
    else if (quantity == "Rout")   { valCol = "Rout";   errCol = "Rout_err";  }
    else if (quantity == "Rside")  { valCol = "Rside";  errCol = "Rside_err"; }
    else if (quantity == "Rlong")  { valCol = "Rlong";  errCol = "Rlong_err"; }
    else if (quantity == "lambda") { valCol = lambdaCutoff + "_lambda";
                                     errCol = lambdaCutoff + "_lambda_err"; }
    else { cerr << "Unknown quantity " << quantity << endl; return pts; }

    int iMT  = findCol("mT");
    int iVal = findCol(valCol);
    int iErr = findCol(errCol);
    if (iMT < 0 || iVal < 0 || iErr < 0) {
        cerr << "[" << filename << "] missing column(s): "
             << "mT=" << iMT << " " << valCol << "=" << iVal
             << " " << errCol << "=" << iErr << endl;
        return pts;
    }

    int maxIdx = std::max(std::max(iMT, iVal), iErr);

    string line;
    while (getline(f, line)) {
        if (line.empty()) continue;
        vector<string> F = SplitAny(line);
        if ((int)F.size() <= maxIdx) continue;
        try {
            CSVPoint p;
            p.x     = stod(F[iMT]);
            p.y     = stod(F[iVal]);
            p.ey_lo = p.ey_hi = fabs(stod(F[iErr]));
            pts.push_back(p);
        } catch (...) {}
    }
    return pts;
}

// ------------------------------------------------------------
// Build TGraphErrors from CSVPoint
// ------------------------------------------------------------
TGraphErrors* MakeGraph(const vector<CSVPoint>& pts,
                        Color_t color, Style_t marker, Size_t msize)
{
    const int n = (int)pts.size();
    if (n == 0) return nullptr;

    vector<double> vx(n), vy(n), vey(n);
    for (int i = 0; i < n; ++i) {
        vx[i]  = pts[i].x;
        vy[i]  = pts[i].y;
        vey[i] = 0.5 * (pts[i].ey_lo + pts[i].ey_hi);
    }

    TGraphErrors* g = new TGraphErrors(n, vx.data(), vy.data(),
                                          nullptr,   vey.data());
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerSize(msize);
    g->SetLineWidth(2);
    return g;
}

// ------------------------------------------------------------
// One panel
//   panelLabel     : "(a)", "(b)", ... top-right, non-bold
//   showLegend     : only on the first panel
//   wideLabel      : legend text for the user's wide CSV
//   paperLabel     : legend text for the paper CSV
//   lambdaCutoff   : "1mm", "5mm", "5cm", "50cm"
//   paperSource    : "EPOS3" for radii/alpha,
//                    "EPOS3_rho1mm" for lambda at 1 mm
// ------------------------------------------------------------
void DrawOnePad(const string& wideFile,
                const string& paperFile,
                const string& quantity,
                const vector<string>& paperPlotCandidates,
                const char*   ylabel,
                const char*   panelLabel,
                double        yMin, double yMax,
                const char*   wideLabel,
                const char*   paperLabel,
                bool          showLegend      = false,
                const string& lambdaCutoff    = "1mm",
                const string& paperSource     = "EPOS3",
                const string& paperCentrality = "0-10")
{
    gPad->SetLeftMargin  (0.17);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin   (0.10);
    gPad->SetRightMargin (0.04);

    // ---- Load both series ----
    vector<CSVPoint> wide  = LoadWide(wideFile, quantity, lambdaCutoff);
    vector<CSVPoint> paper = LoadLongAnyName(paperFile, paperPlotCandidates,
                                             paperCentrality, paperSource);

    cout << "  [" << quantity << "] wide=" << wide.size()
         << " paper=" << paper.size() << "\n";

    TGraphErrors* gWide  = MakeGraph(wide,  kBlue, 21, 1.0);
    TGraphErrors* gPaper = MakeGraph(paper, kRed,  24, 1.0);

    // ---- Frame ----
    TH1* frame = gPad->DrawFrame(0.15, yMin, 0.9, yMax);
    frame->GetXaxis()->SetTitle("m_{T} [GeV/c^{2}]");
    frame->GetYaxis()->SetTitle(ylabel);
    frame->GetXaxis()->SetTitleSize  (0.060);
    frame->GetXaxis()->SetLabelSize  (0.050);
    frame->GetYaxis()->SetTitleSize  (0.060);
    frame->GetYaxis()->SetLabelSize  (0.050);
    frame->GetYaxis()->SetTitleOffset(1.30);
    frame->GetXaxis()->SetTitleOffset(1.00);

    // ---- Draw ----
    if (gWide)  gWide ->Draw("PE SAME");
    if (gPaper) gPaper->Draw("PE SAME");

    // ---- Panel label: top-right, non-bold ----
    TLatex lab;
    lab.SetNDC();
    lab.SetTextFont(42);
    lab.SetTextSize(0.060);
    lab.SetTextAlign(31);
    lab.DrawLatex(0.95, 0.92, panelLabel);

    // ---- Legend: only on the first pad ----
    if (showLegend) {
        TLegend* leg = new TLegend(0.45, 0.72, 0.93, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.045);
        if (gWide)  leg->AddEntry(gWide,  wideLabel,  "P");
        if (gPaper) leg->AddEntry(gPaper, paperLabel, "P");
        leg->Draw();
    }
}

// ------------------------------------------------------------
// Info panel (f)
// ------------------------------------------------------------
void DrawInfoPad()
{
    gPad->SetLeftMargin  (0.17);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin   (0.10);
    gPad->SetRightMargin (0.04);

    TH1* frame = gPad->DrawFrame(0.0, 0.0, 1.0, 1.0);
    frame->GetXaxis()->SetNdivisions(0);
    frame->GetYaxis()->SetNdivisions(0);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0);
    frame->GetXaxis()->SetTitleSize(0);
    frame->GetYaxis()->SetTitleSize(0);

    const double x1 = 0.17, y1 = 0.15, x2 = 0.96, y2 = 0.90;

    TPaveText* box = new TPaveText(x1, y1, x2, y2, "NDC");
    box->SetFillStyle(0);
    box->SetBorderSize(0);
    box->SetTextFont(42);
    box->SetTextAlign(12);
    box->SetTextSize(0.055);
    box->AddText("Collision System");
    box->AddText("---------------------------");
    box->AddText("Au+Au,  #sqrt{s_{NN}} = 200 GeV");
    box->AddText("centrality 0-10%");
    box->AddText(" ");
    box->AddText("#pi^{#pm}#pi^{#pm} pairs");
    box->AddText(" ");
    box->AddText("0.15 < p_{T} < 1.00 GeV/c");
    box->AddText("|#eta| < 1.0");
    box->AddText(" ");
    box->AddText("#rho_{max}^{#lambda} = 1 mm");
    box->Draw();
}

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
void compare_SMASH_EPOS_all()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    // Wide-format file: your own analysis output.
    const string wideFile  = "EPOS_Levy_fit_parameters_cutoff.csv";
    // Long-format file: paper supplementary CSV.
    const string paperFile = "10052_2026_15526_MOESM1_ESM.csv";

    cout << "\n=== wide  : " << wideFile  << " ===\n";
    cout << "=== paper : " << paperFile << " ===\n";

    TCanvas* cAll = new TCanvas("cAll",
                                "EPOS analysis vs paper: radii and shape parameters",
                                2000, 1200);
    cAll->Divide(3, 2);

    // ---- Top row: radii ----
    cAll->cd(1);
    DrawOnePad(wideFile, paperFile, "Rout", {"R_out"},
               "R_{out} [fm]", "(a)", 2.0, 10.0,
               "EPOS (0-5%)", "EPOS3 (0-10%)",
               true, "1mm", "EPOS3");

    cAll->cd(2);
    DrawOnePad(wideFile, paperFile, "Rside", {"R_side"},
               "R_{side} [fm]", "(b)", 2.0, 10.0,
               "EPOS (0-5%)", "EPOS3 (0-10%)",
               false, "1mm", "EPOS3");

    cAll->cd(3);
    DrawOnePad(wideFile, paperFile, "Rlong", {"R_long"},
               "R_{long} [fm]", "(c)", 3.0, 15.0,
               "EPOS (0-5%)", "EPOS3 (0-10%)",
               false, "1mm", "EPOS3");

    // ---- Bottom row: shape parameters ----
    cAll->cd(4);
    DrawOnePad(wideFile, paperFile, "alpha", {"alpha"},
               "#alpha", "(d)", 1.40, 1.75,
               "EPOS (0-5%)", " EPOS3 (0-10%)",
               false, "1mm", "EPOS3");

    cAll->cd(5);
    DrawOnePad(wideFile, paperFile, "lambda", {"lambda_rhotest", "lambda"},
               "#lambda_{1mm}", "(e)", 0.80, 0.98,
               "EPOS (0-5%, 1 mm)", "EPOS3 (0-10%, 1 mm)",
               false, "1mm", "EPOS3_rho1mm");

    // ---- (f) info pad ----
    cAll->cd(6);
    DrawInfoPad();

    cAll->SaveAs("Compare_EPOS05_analysis_vs_paperEpos3.png");
    cAll->SaveAs("Compare_EPOS05_analysis_vs_paperEpos3.pdf");

    cout << "\nSaved Compare_EPOS_analysis_vs_paper_combined.{png,pdf}\n";
}
