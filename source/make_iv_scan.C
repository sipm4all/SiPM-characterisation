#pragma once

#include "makeiv.C"

bool iv_fit_oneshot = true;

double iv_fit_function(double *x, double *p)
{
  double x2 = (x[0] - p[3]) * (x[0] - p[3]);
  double x4 = x2 * x2;
  double x6 = x2 * x2 * x2;
  double below = p[0] * std::exp(p[2] * x[0]) + p[1] * std::exp(p[2] * x[0] * x[0]);
  double transition = 0.5 * (std::erf( (x[0] - p[3]) / p[4]) + 1.);
  double above = p[5] * p[5] * x2 + p[6] * p[6] * x4 + p[7] * p[7] * x6;
  return below + transition * above;
}

TF1 *get_iv_fit_function(double vbreak)
{
  auto f = new TF1("f_iv_fit_function", iv_fit_function, 0., 100., 7);
  f->SetParameter(0, 1.e-10);
  f->SetParameter(1, 1.e-10);
  f->SetParameter(2, -1.e-10);
  f->SetParameter(3, vbreak);
  f->SetParameter(4, 0.1);
  f->SetParameter(5, 1.e-10);
  f->SetParameter(6, 1.e-10);
  f->SetParameter(7, 1.e-10);
  return f;
}

TF1 *
iv_fit_better(TGraphErrors *g, float vstart, int ntry = 100)
{
  TVirtualFitter::SetMaxIterations(kMaxInt);
  auto fall = get_iv_fit_function(vstart);

  for (int itry = 0; itry < ntry; ++itry) {
  
    auto fbelow = get_iv_fit_function(vstart);
    for (int i = 3; i < fbelow->GetNpar(); ++i) fbelow->FixParameter(i, 0.);
    g->Fit(fbelow, "0q", "", 0., vstart - 1.);
    
    auto fabove = get_iv_fit_function(vstart);
    for (int i = 0; i < 3; ++i) fabove->FixParameter(i, 0.);
    fabove->FixParameter(3, vstart);
    fabove->FixParameter(4, 0.1);
    g->Fit(fabove, "0q", "", vstart + 1., vstart + 3.);
    
    for (int i = 0; i < 3; ++i) fall->SetParameter(i, fbelow->GetParameter(i));
    for (int i = 3; i < fbelow->GetNpar(); ++i) fall->SetParameter(i, fabove->GetParameter(i));
    fall->FixParameter(3, vstart);
    fall->FixParameter(4, 0.1);
    g->Fit(fall, "0q", "", 0., vstart + 3.);
    fall->ReleaseParameter(3);
    fall->ReleaseParameter(4);
    g->Fit(fall, "0q", "", 0., vstart + 3.);
    g->Fit(fall, "0q", "", 0., vstart + 3.);
    g->Fit(fall, "0q", "", 0., vstart + 3.);
    
    delete fbelow;
    delete fabove;

    std::cout << fall->GetParameter(3) << " +- " << fall->GetParError(3) << std::endl;
    
    if ( std::fabs(fall->GetParameter(3) - vstart) < 0.005 ) break;
    
    vstart = fall->GetParameter(3);
  }
  
  return fall;


}

TF1 *
iv_fit(TGraphErrors *g, float vstart, bool oneshot = true)
{
  TVirtualFitter::SetMaxIterations(kMaxUInt);


  std::cout << " --- IV FIT : " << g << std::endl;
    //    g->Print("all");
  
  //  auto f = new TF1("f", "expo(0) + pol0(2) + (x > [3]) * ( [4] * pow(abs(x - [3]), [5]) )", 0., 100.);



  auto f = new TF1("fSURF", "expo(0) + pol0(2) + (x > [3]) * ( [4] * pow(abs(x - [3]), [5]) + [6] * pow(abs(x - [3]), [7]) )", 0., 100.);
    
  for (int itry = 0; itry < 100; ++itry) {
  
  // fit background
  for (int i = 0; i < f->GetNpar(); ++i) f->FixParameter(i, 0.);
  f->ReleaseParameter(0);
  f->ReleaseParameter(1);
  f->SetParLimits(1, -1.e6, 0.);
  g->Fit(f, "0q", "", vstart - 2., vstart - 1.);

  // fit sensor
  for (int i = 0; i < f->GetNpar(); ++i) f->ReleaseParameter(i);
  f->FixParameter(0, f->GetParameter(0));
  f->FixParameter(1, f->GetParameter(1));
  f->SetParameter(2, 1.e-12);
  f->SetParLimits(2, 0., 1.e-10);
  f->SetParameter(3, vstart);
  f->SetParLimits(3, vstart - 1., vstart + 1.);
  f->SetParameter(4, 1.e-12);
  f->SetParLimits(4, 0., 1.e-6);
  f->FixParameter(5, 2.);
  //  f->SetParLimits(5, 0., 10.);
  f->SetParameter(6, 1.e-12);
  //  f->SetParLimits(6, 0., 1.e-6);
  f->FixParameter(7, 4.);
  //  f->SetParLimits(7, 0., 10.);
  g->Fit(f, "0q", "", vstart + 1., vstart + 2.);
  
  // fit all
  for (int i = 0; i < f->GetNpar(); ++i) f->ReleaseParameter(i);
  f->SetParLimits(1, -1.e3, 0.);
  f->SetParLimits(2, 0., 1.e-10);
  f->SetParLimits(3, vstart - 3., vstart + 3.);
  f->SetParLimits(4, 0., 1.e-6);
  f->FixParameter(5, 2.);
  //  f->SetParLimits(5, 0., 10.);
  // R+++++  f->SetParLimits(6, 0., 1.e-6);
  f->FixParameter(7, 4.);
  //  f->SetParLimits(7, 0., 10.);
  g->Fit(f, "0", "", 0., vstart + 3.);
  f->ReleaseParameter(5);
  f->ReleaseParameter(7);
  int res = g->Fit(f, "0M", "", 0., vstart + 3.);
  std::cout << " --->>> " << res << std::endl;

  f->SetRange(0., vstart + 5.);

  if (iv_fit_oneshot || fabs(f->GetParameter(3) - vstart) < 0.01) break;

  vstart = f->GetParameter(3);
  
  }

  //  f->SetParameter(4, 0.);
  //  f->DrawClone("same");
  std::cout << "vdb = " << f->GetParameter(3) << " +- " << f->GetParError(3) << std::endl;
  
  return f;
}

TGraphErrors *
make_iv_scan(std::string filename, std::string filenameOPEN = "", std::string filenameZERO = "", float vbreak = -1.)
{
  auto g = makeiv(filename, filenameZERO, true, true, false);
  if (!g) return nullptr;
  if (filenameOPEN.empty()) {
    //    make_csv(g, filename + ".processed.csv");
    return g;
  }

  auto gOPEN = makeiv(filenameOPEN, "", true, true, false);
  if (!gOPEN) return nullptr;
  auto fOPEN = (TF1 *)gROOT->GetFunction("pol1");    
  gOPEN->Fit(fOPEN);
  for (int i = 0; i < g->GetN(); ++i) {
    double offset = fOPEN->Eval(g->GetX()[i]);
    g->GetY()[i] -= offset;
  }

  //  make_csv(g, filename + ".processed.opensubtracted.csv");
  
  return g;
}

TCanvas *
draw_iv_scan(std::string dirname,
             std::string tagname,
             std::vector<std::string> rows,
             std::vector<std::string> cols,
             std::vector<std::string> outliers,
             TCanvas *c = nullptr,
             int marker = 20,
             int color = kAzure-3,
             float vmin = 0.,
             float vmax = 100.,
             float vcut = 100.,
             bool subtractopen = true,
             float subtractsurface = -1.)
{

  if (!c) {
    c = new TCanvas("c", "c", 800, 800);
    c->DrawFrame(vmin, 1.e-13, vmax, 1.e-3, ";bias voltage (V);current (A)");
    c->SetLogy();
  }

  //  auto cmap = new TCanvas("cmap", "cmap", 1600, 800);
  //  cmap->Divide(8, 4);

  auto hbreak = new TH1F("hbreak", "", 1000., 0., 100.);
  
  int ipad = 1;
  for (auto row : rows) {
    for (auto col : cols) {
      //      cmap->cd(ipad)->DrawFrame(vmin, 1.e-11, vmax, 1.e-3, ";bias voltage (V);current (A)");
      //      cmap->cd(ipad)->SetLogy();
      ipad++;
      std::string chname = row + col;
      if (std::find(outliers.begin(), outliers.end(), chname) != outliers.end()) continue;
      std::string filename = dirname + "/" + tagname + "_243K_" + chname + ".ivscan.csv";
      std::string filenameOPEN = dirname + "/" + tagname + "_243K_OPEN-" + chname + ".ivscan.csv";

      auto g = make_iv_scan(filename, subtractopen ? filenameOPEN : "", "", subtractsurface);
      
      if (subtractsurface > 0.) {
        auto fSURF = iv_fit(g, subtractsurface);
        hbreak->Fill(fSURF->GetParameter(3.));
        fSURF->SetParameter(4, 0.);
        fSURF->SetRange(0., 100.);
        for (int i = 0; i < g->GetN(); ++i) {
          double surfcurr = fSURF->Eval(g->GetX()[i]);
          g->GetY()[i] -= surfcurr;
        }    
      }

      auto gg = new TGraphErrors;
      gg->SetName(g->GetName());
      gg->SetMarkerStyle(marker);
      gg->SetMarkerColor(color);
      gg->SetLineColor(color);
      for (int i = 0; i < g->GetN(); ++i) {
        if (g->GetX()[i] > vcut) continue;
        auto n = gg->GetN();
        gg->SetPoint(n, g->GetX()[i], g->GetY()[i]);
        gg->SetPointError(n, g->GetEX()[i], g->GetEY()[i]);
      }

      gg->Draw("samelp");

    }
  }

  new TCanvas;
  hbreak->Draw();
  
  return c;
}

