#pragma once

#include "make_iv_scan.C"
#include "ureadout_dcr_get.C"
#include "../utils/graphutils.C"

namespace database {

std::string fields[8] = {
  "run", "quality", "step", "setup", "iv-mux-1", "iv-mux-2", "dcr-chip-2", "dcr-chip-3" 
};

std::vector<std::string> channels = {
  "A1", "A2", "A3", "A4",
  "B1", "B2", "B3", "B4",
  "C1", "C2", "C3", "C4"
}; 

std::string basedir = ".";
std::map<std::string, std::map<std::string, std::map<std::string, std::string>>> boards;

void read_database(std::string fname);
void dump_filenames(std::string fname);
std::map<std::string, std::string> get_filename(std::string board, std::string channel, std::string step);
TGraphErrors *get_dcr_vbias_scan(std::string board, std::string channel, std::string step, int marker = 1, int color = 1);
TGraphErrors *get_iv_scan(std::string board, std::string channel, std::string step, int marker = 1, int color = 1);
TGraphErrors *get_gain(std::string board, std::string channel, std::string step, int marker = 1, int color = 1);

TGraphErrors *
get_gain(std::string board, std::string channel, std::string step, int marker, int color)
{
  /** get IV **/
  auto giv = get_iv_scan(board, channel, step, marker, color);
  if (!giv) return nullptr;
  /** subtract surface current **/
  float ave = 0.;
  for (int i = 0; i < 5; ++i) ave += giv->GetY()[i];
  graphutils::y_shift(giv, ave / 5.);

  /** get DCR **/
  auto gdcr = get_dcr_vbias_scan(board, channel, step);
  if (!gdcr) return nullptr;

  /** compute gain **/
  auto ggain = graphutils::ratio(giv, gdcr);
  graphutils::y_scale(ggain, 1. / TMath::Qe());

  return ggain;
}

void
dump_filenames(std::string fname)
{
  ofstream fout(fname);
  for (auto &val1 : boards) {
    auto board = val1.first;
    if (board == ".") continue;
    for (auto &val2 : val1.second) {
      auto step = val2.first;
      for (auto &channel : channels) {
        auto fnames = get_filename(board, channel, step);
        for (auto fname : fnames)
          fout << fname.second << std::endl;
      }
    }
  }
  fout.close();
}

TGraphErrors *
get_dcr_vbias_scan(std::string board, std::string channel, std::string step, int marker, int color)
{
  auto fname = get_filename(board, channel, step);
  if (fname["dcr-vbias"].empty()) return nullptr;
  auto g = ureadout_dcr_get(fname["dcr-vbias"], "bias_voltage", "dead_rate");
  graphutils::set_style(g, marker, color);
  return g;
}

TGraphErrors *
get_iv_scan(std::string board, std::string channel, std::string step, int marker, int color)
{
  auto fname = get_filename(board, channel, step);
  if (fname["iv"].empty() || fname["iv-open"].empty()) return nullptr;
  auto g = make_iv_scan(fname["iv"], fname["iv-open"]);
  graphutils::set_style(g, marker, color);
  return g;
}

std::map<std::string, std::string>
get_filename(std::string board, std::string channel, std::string step)
{
  if (!boards.count(board)) return {};
  if (!boards[board].count(step)) return {};
  
  auto ivrun = boards[board][step]["iv-run"];
  auto ivsetup = boards[board][step]["iv-setup"];
  auto ivmux = boards[board][step]["iv-mux"];
  std::string ivfname = basedir + "/" + ivrun + "/" + ivsetup + "/iv/HAMA3_sn0_mux" + ivmux + "/HAMA3_sn0_243K_" + channel + ".ivscan.csv";
  std::string ivofname = basedir + "/" + ivrun + "/" + ivsetup + "/iv/HAMA3_sn0_mux" + ivmux + "/HAMA3_sn0_243K_OPEN-" + channel + ".ivscan.csv";

  auto dcrrun = boards[board][step]["dcr-run"];
  auto dcrsetup = boards[board][step]["dcr-setup"];
  auto dcrchip = boards[board][step]["dcr-chip"];

  std::string dcrthrfname = basedir + "/" + dcrrun + "/" + dcrsetup + "/dcr/HAMA3-chip" + dcrchip + "/rate/threshold_scan/chip" + dcrchip + "-" + channel + ".ureadout_dcr_scan.tree.root";
  std::string dcrbiasfname = basedir + "/" + dcrrun + "/" + dcrsetup + "/dcr/HAMA3-chip" + dcrchip + "/rate/vbias_scan/chip" + dcrchip + "-" + channel + ".ureadout_dcr_scan.tree.root";

  return { {"iv", ivfname} , {"iv-open", ivofname} , {"dcr-threshold", dcrthrfname} , {"dcr-vbias", dcrbiasfname} };
}

void
read_database(std::string fname)
{
  std::ifstream ifs(fname);
  std::string line;
  std::string board;
  while (std::getline(ifs, line)) {
    if (line[0] == '#' || line[0] == ' ') continue;
    std::stringstream ss(line);
    std::map<std::string, std::string> linedata;
    std::string data;
    for (int i = 0; i < 8; ++i) {
      ss >> data;
      linedata[fields[i]] = data;
    }
    
    auto run = linedata["run"];
    auto step = linedata["step"];
    auto setup = linedata["setup"]; 

    board = linedata["iv-mux-1"];
    boards[board][step]["iv-run"] = run;
    boards[board][step]["iv-setup"] = setup;
    boards[board][step]["iv-mux"] = "1";

    board = linedata["iv-mux-2"];
    boards[board][step]["iv-run"] = run;
    boards[board][step]["iv-setup"] = setup;
    boards[board][step]["iv-mux"] = "2";

    board = linedata["dcr-chip-2"];
    boards[board][step]["dcr-run"] = run;
    boards[board][step]["dcr-setup"] = setup;
    boards[board][step]["dcr-chip"] = "2";

    board = linedata["dcr-chip-3"];
    boards[board][step]["dcr-run"] = run;
    boards[board][step]["dcr-setup"] = setup;
    boards[board][step]["dcr-chip"] = "3";
    
  }

}

}
