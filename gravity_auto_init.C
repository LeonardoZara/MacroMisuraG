// gravity_auto_init.C
// Macro ROOT: fit automatico di una sinusoide smorzata con drift lineare
//
// Modello:
//   y(t) = A sin(omega t + phi) exp(-(t-x0)/tau) + q + m t
//
// La macro stima automaticamente i parametri iniziali:
//   q,m     : fit lineare preliminare della baseline
//   A       : metà del range dei dati corretti
//   omega   : distanza media tra massimi locali
//   phi     : dal primo punto corretto
//   tau     : decadimento approssimato dell'inviluppo dei picchi
//
// Uso:
//   .x gravity_auto_init.C
// oppure:
//   gravity_auto_init("gravity.dat", 0.0, 1500.0, false);

#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TAxis.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

struct InitPars {
  double A;
  double omega;
  double phi;
  double tau;
  double q;
  double m;
};

static double Clamp(double v, double lo, double hi) {
  if (v < lo) return lo;
  if (v > hi) return hi;
  return v;
}

static InitPars EstimateInitialParameters(TGraph *g, double xmin, double xmax) {
  const int n = g->GetN();
  double *x = g->GetX();
  double *y = g->GetY();

  // 1) Fit lineare preliminare per stimare drift q + m x
  TF1 *fLin = new TF1("fLin_auto_init","pol1", xmin, xmax);
  g->Fit(fLin, "RQ0"); // R: range, Q: quiet, 0: non disegna

  double q = fLin->GetParameter(0);
  double m = fLin->GetParameter(1);

  // 2) Dati corretti dalla baseline
  vector<double> xc, yc;
  xc.reserve(n);
  yc.reserve(n);

  for (int i=0; i<n; ++i) {
    if (x[i] < xmin || x[i] > xmax) continue;
    xc.push_back(x[i]);
    yc.push_back(y[i] - (q + m*x[i]));
  }

  int nc = (int)xc.size();
  if (nc < 5) {
    cerr << "Pochi punti nel range: uso parametri di fallback." << endl;
    return {1.0, 0.02, 0.0, 500.0, q, m};
  }

  // 3) Ampiezza iniziale dal range dei dati corretti
  auto mm = minmax_element(yc.begin(), yc.end());
  double ymin = *mm.first;
  double ymax = *mm.second;
  double A = 0.5*(ymax - ymin);
  if (A <= 0 || !isfinite(A)) A = 1.0;

  // 4) Trova massimi locali positivi per stimare il periodo
  vector<double> peakX;
  vector<double> peakY;

  // soglia minima per evitare piccoli picchi dovuti al rumore
  double threshold = 0.20 * A;

  for (int i=1; i<nc-1; ++i) {
    if (yc[i] > yc[i-1] && yc[i] > yc[i+1] && yc[i] > threshold) {
      peakX.push_back(xc[i]);
      peakY.push_back(yc[i]);
    }
  }

  double omega = 0.02; // fallback
  if (peakX.size() >= 2) {
    vector<double> periods;
    for (size_t k=1; k<peakX.size(); ++k) {
      double T = peakX[k] - peakX[k-1];
      if (T > 0) periods.push_back(T);
    }
    if (!periods.empty()) {
      sort(periods.begin(), periods.end());
      double Tmed = periods[periods.size()/2];
      omega = 2.0*TMath::Pi()/Tmed;
    }
  }

  // 5) Fase iniziale dal primo punto corretto
  double s = Clamp(yc[0]/A, -1.0, 1.0);
  double phi = asin(s) - omega*xc[0];

  // 6) Stima grossolana di tau dall'inviluppo dei picchi
  double tau = 0.5*(xmax - xmin); // fallback
  if (peakX.size() >= 2) {
    double x1 = peakX.front();
    double y1 = fabs(peakY.front());
    double x2 = peakX.back();
    double y2 = fabs(peakY.back());

    if (y1 > 0 && y2 > 0 && y2 < y1 && x2 > x1) {
      tau = -(x2 - x1)/log(y2/y1);
    }
  }

  // Limiti di sicurezza
  if (!isfinite(tau) || tau <= 0) tau = 0.5*(xmax-xmin);
  if (!isfinite(omega) || omega <= 0) omega = 0.02;
  if (!isfinite(phi)) phi = 0.0;

  delete fLin;
  return {A, omega, phi, tau, q, m};
}

void gravity_auto_init(const char *filename="gravity.dat",
                       double xmin=0.0,
                       double xmax=1210.0,
                       bool useErrors=false) {

  gStyle->SetOptFit(1111);

  TGraph *graph = nullptr;

  if (useErrors) {
    graph = new TGraphErrors(filename);
  } else {
    graph = new TGraph(filename);
  }

  if (!graph || graph->GetN() <= 0) {
    cerr << "Errore: impossibile leggere il file " << filename << endl;
    return;
  }

  graph->SetTitle("Oscillazione smorzata;tempo t;posizione y");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.6);
  graph->SetLineColor(kBlack);

  InitPars p = EstimateInitialParameters(graph, xmin, xmax);

  cout << "\nParametri iniziali stimati automaticamente:" << endl;
  cout << "A      = " << p.A << endl;
  cout << "omega  = " << p.omega << " rad/s" << endl;
  cout << "phi    = " << p.phi << endl;
  cout << "tau    = " << p.tau << " s" << endl;
  cout << "q      = " << p.q << endl;
  cout << "m      = " << p.m << endl;

  TF1 *f = new TF1("f",
    "([0]*sin(x*[1]+[2]))*exp(-(x-[6])/[3])+[4]+[5]*x",
    xmin, xmax);

  f->SetParNames("A","omega","phi","tau","q","m","x0");
  f->SetParameters(p.A, p.omega, p.phi, p.tau, p.q, p.m, xmin);

  // x0 viene fissato al bordo sinistro per rendere tau più stabile
  f->FixParameter(6, xmin);

  // Limiti ragionevoli per aiutare la convergenza
  f->SetParLimits(0, 0.05*p.A, 5.0*p.A);
  f->SetParLimits(1, 0.2*p.omega, 5.0*p.omega);
  f->SetParLimits(3, 0.05*(xmax-xmin), 10.0*(xmax-xmin));

  graph->Fit(f, "RME"); 
  // R: range
  // M: miglior minimizzatore quando possibile
  // E: stima errori migliore

  // Baseline finale
  TF1 *baseline = new TF1("baseline","[0]+[1]*x", xmin, xmax);
  baseline->SetParameters(f->GetParameter(4), f->GetParameter(5));
  baseline->SetLineColor(kGreen+2);
  baseline->SetLineStyle(2);
  baseline->SetLineWidth(2);

  // Costruisce grafico corretto dalla baseline
  int n = graph->GetN();
  double *x = graph->GetX();
  double *y = graph->GetY();

  vector<double> xcorr, ycorr;
  xcorr.reserve(n);
  ycorr.reserve(n);

  for (int i=0; i<n; ++i) {
    if (x[i] < xmin || x[i] > xmax) continue;
    xcorr.push_back(x[i]);

    // Correzione del drift:
    // sottraggo solo la pendenza m*x, NON l'intercetta q.
    // In questo modo la baseline nel secondo plot diventa orizzontale,
    // ma resta al livello q invece di andare a zero.
    double qFit = f->GetParameter(4);
    double mFit = f->GetParameter(5);
    ycorr.push_back(y[i] - mFit*x[i]);
  }

  TGraph *gCorr = new TGraph((int)xcorr.size(), xcorr.data(), ycorr.data());
  gCorr->SetTitle("Dati corretti dalla pendenza;tempo t;y - m t");
  gCorr->SetMarkerStyle(20);
  gCorr->SetMarkerSize(0.6);

  TF1 *fCorr = new TF1("fCorr",
    "[0]*sin(x*[1]+[2])*exp(-(x-[5])/[3]) + [4]",
    xmin, xmax);

  fCorr->SetParNames("A","omega","phi","tau","q","x0");
  fCorr->SetParameters(f->GetParameter(0),
                       f->GetParameter(1),
                       f->GetParameter(2),
                       f->GetParameter(3),
                       f->GetParameter(4),
                       xmin);
  fCorr->FixParameter(5, xmin);
  fCorr->SetLineColor(kRed);
  fCorr->SetLineWidth(2);

  gCorr->Fit(fCorr, "RME");

  double omega = fCorr->GetParameter(1);
  double omegaErr = fCorr->GetParError(1);
  double period = 2.0*TMath::Pi()/omega;
  double periodErr = period * omegaErr/omega;

  cout << "\nRisultati finali dopo correzione baseline:" << endl;
  cout << "A      = " << fCorr->GetParameter(0)
       << " +/- " << fCorr->GetParError(0) << endl;
  cout << "omega  = " << omega
       << " +/- " << omegaErr << " rad/s" << endl;
  cout << "T      = " << period
       << " +/- " << periodErr << " s" << endl;
  cout << "tau    = " << fCorr->GetParameter(3)
       << " +/- " << fCorr->GetParError(3) << " s" << endl;
  cout << "q      = " << f->GetParameter(4)
       << " +/- " << f->GetParError(4) << endl;
  cout << "m      = " << f->GetParameter(5)
       << " +/- " << f->GetParError(5) << endl;
  cout << "chi2/ndf original = " << f->GetChisquare()/f->GetNDF() << endl;
  cout << "chi2/ndf corretto = " << fCorr->GetChisquare()/fCorr->GetNDF() << endl;

  TCanvas *c = new TCanvas("c_gravity_auto","gravity auto init",1200,550);
  c->Divide(2,1);

  c->cd(1);
  graph->Draw("AP");
  f->SetLineColor(kRed);
  f->SetLineWidth(2);
  f->Draw("SAME");
  baseline->Draw("SAME");

  TLegend *leg1 = new TLegend(0.12,0.72,0.52,0.88);
  leg1->AddEntry(graph,"dati","p");
  leg1->AddEntry(f,"fit completo","l");
  leg1->AddEntry(baseline,"baseline lineare","l");
  leg1->Draw();

  c->cd(2);
  gCorr->Draw("AP");

  // Baseline corretta: dopo aver sottratto solo m*x,
  // resta una retta orizzontale al livello q.
  TF1 *baselineFlat = new TF1("baselineFlat","[0]", xmin, xmax);
  baselineFlat->SetParameter(0, f->GetParameter(4));
  baselineFlat->SetLineColor(kGreen+2);
  baselineFlat->SetLineStyle(2);
  baselineFlat->SetLineWidth(2);
  baselineFlat->Draw("SAME");

  fCorr->Draw("SAME");

  TLegend *leg2 = new TLegend(0.12,0.72,0.58,0.90);
  leg2->AddEntry(gCorr,"dati corretti: y - m x","p");
  leg2->AddEntry(baselineFlat,"baseline orizzontale q","l");
  leg2->AddEntry(fCorr,"fit smorzato + q","l");
  leg2->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.035);
  latex.DrawLatex(0.45,0.24,Form("#omega = %.5g #pm %.2g rad/s", omega, omegaErr));
  latex.DrawLatex(0.45,0.18,Form("T = %.5g #pm %.2g s", period, periodErr));
  latex.DrawLatex(0.45,0.12,Form("#tau = %.5g #pm %.2g s",
                                  fCorr->GetParameter(3),
                                  fCorr->GetParError(3)));

  c->SaveAs("gravity_auto_init.pdf");
  c->SaveAs("gravity_auto_init.png");
}

void gravity() {
  gravity_auto_init();
}
