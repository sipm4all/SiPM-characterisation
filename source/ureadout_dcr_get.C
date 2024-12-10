#pragma once

TGraphErrors *
ureadout_dcr_get(const std::string filename, const std::string whatx, const std::string whaty = "dead_rate")
{
  auto fin = TFile::Open(filename.c_str());
  if (!fin || !fin->IsOpen()) return nullptr;
  auto tin = (TTree *)fin->Get("ureadout_dcr_scan");
  if (!tin) return nullptr;
  std::cout << tin->GetEntries() << std::endl;

  
  
  int base_threshold, threshold, bias_dac;
  float bias_voltage, raw_rate, raw_ratee, dead_rate, dead_ratee, fit_rate, fit_ratee;
  tin->SetBranchAddress("bias_dac", &bias_dac);
  tin->SetBranchAddress("bias_voltage", &bias_voltage);
  tin->SetBranchAddress("base_threshold", &base_threshold);
  tin->SetBranchAddress("threshold", &threshold);
  tin->SetBranchAddress("raw_rate", &raw_rate);
  tin->SetBranchAddress("raw_ratee", &raw_ratee);
  tin->SetBranchAddress("dead_rate", &dead_rate);
  tin->SetBranchAddress("dead_ratee", &dead_ratee);
  tin->SetBranchAddress("fit_rate", &fit_rate);
  tin->SetBranchAddress("fit_ratee", &fit_ratee);

  auto g = new TGraphErrors;

  std::map<std::string, float> x, y, ey;
  
  for (int iev = 0; iev < tin->GetEntries(); ++iev) {
    tin->GetEntry(iev);

    x["bias_voltage"] = bias_voltage;
    x["threshold"] = threshold;

    y["raw_rate"] = raw_rate;
    y["dead_rate"] = dead_rate;
    y["fit_rate"] = fit_rate;
    
    ey["raw_rate"] = raw_ratee;
    ey["dead_rate"] = dead_ratee;
    ey["fit_rate"] = fit_ratee;
    
    if (!x.count(whatx)) {
      std::cout << "unknown whatx: " << whatx << std::endl;
      continue;
    }

    if (!y.count(whaty) || !ey.count(whaty)) {
      std::cout << "unknown whaty: " << whaty << std::endl;
      continue;
    }
    g->SetPoint(iev, x[whatx], y[whaty]);
    g->SetPointError(iev, 0., ey[whaty]);
  }

  g->Sort();
  return g;
}
