#pragma once

#include "database.C"
#include "../utils/general_utility.h"
#include "../utils/tree_database.h"

namespace database
{
  namespace laser
  {
    //  --- Constants ---
    //  Physics constants
    const double coarse_to_ns = 3.125; // ns
    //  --- Utility data structures ---
    //  --- Data structure in database reference
    std::vector<std::string> fields = {"run", "status", "step", "board", "channels"};
    std::vector<std::string> sub_run_fields = {"channel", "bcrconfig", "opmode", "deltathreshold", "vbias", "directory", "notes"};
    //  --- Tree batabase
    tree_database database;
    std::unordered_map<int, std::vector<std::string>> info_database;
    //  --- Possible positions
    std::vector<std::string> laser_positions = {"empty", "center", "top-left", "top-right", "bottom-left", "bottom-right"};
    //  --- Global variables
    std::string basedir = "./Data/";

    //  --- Declarations ---
    //  --- I/O
    //  --- --- database
    void read_database(std::string database_input_file);
    void read_sub_run_database(std::string database_input_file, std::string type = "target");
    void read_all_sub_runs(std::string target_run_tag);
    void read_all_sub_runs();
    //  --- --- data
    void download_run(std::string campaign, std::string target_run_tag = ".");
    void unzip_sub_run(std::string target_file, std::string target_dir = ".", std::string campaign = "2024-laser-window");
    void unzip_all_sub_runs(std::string target_run_tag);
    void update_local_data_repository(std::string current_campaign = "2024-laser-window");
    //  --- Getters
    std::unordered_map<std::string, std::vector<std::string>> get_run_infos(std::string run);
    std::unordered_map<std::string, std::vector<std::string>> get_sub_run_infos(std::string sub_run);
    std::unordered_map<std::string, std::vector<std::string>> get_sub_run_infos(int sub_run_id);
    std::unordered_map<std::string, std::vector<std::string>> get_sub_run_infos(std::string board, std::string step, std::string channel, std::string vbias, std::string type, std::string position);
    std::vector<std::pair<int, std::string>> get_all_laser_runs() { return database.get_children(-1, 4, {{3, "laser"}}); }
    std::vector<std::pair<int, std::string>> get_all_laser_runs(std::string quality);
    std::vector<std::pair<int, std::string>> get_all_laser_runs(std::string board, std::string step) { return database.get_children(-1, 4, {{3, "laser"}, {2, step}, {1, board}}); }
    std::pair<int, std::string> get_latest_laser_run(std::string board, std::string step);
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> get_sub_runs_infos(std::string run, std::string channel, std::string vbias, std::string type, std::string position);
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> get_sub_runs_infos(std::string board, std::string step, std::string channel, std::string vbias, std::string type, std::string position) { return get_sub_runs_infos(get_latest_laser_run(board, step).second, channel, vbias, type, position); }
    //  --- General utilities
    std::map<std::string, TGraphErrors *> get_graph_sensor_from_channel(std::string run, std::string sensor, std::string type, std::string position, TGraphErrors *(*graph_getter)(std::string, std::string, std::string, std::string));
    //  --- Measurements
    std::array<std::array<float, 2>, 2> get_sub_run_value(std::string sub_run);
    std::array<float, 2> get_sub_run_sig(std::string sub_run) { return get_sub_run_value(sub_run)[0]; }
    std::array<float, 2> get_sub_run_bkg(std::string sub_run) { return get_sub_run_value(sub_run)[1]; }
    TGraphErrors *get_pPDE_vs_vbias(std::string run, std::string channel, std::string type, std::string position);
    TGraphErrors *get_bkg_vs_vbias(std::string run, std::string channel, std::string type, std::string position);
    TGraphErrors *get_pPDE_vs_bkg(std::string run, std::string channel, std::string type, std::string position);
    std::map<std::string, TGraphErrors *> get_sensor_pPDE_vs_vbias(std::string run, std::string sensor, std::string type, std::string position) { return get_graph_sensor_from_channel(run, sensor, type, position, get_pPDE_vs_vbias); };
    std::map<std::string, TGraphErrors *> get_sensor_bkg_vs_vbias(std::string run, std::string sensor, std::string type, std::string position) { return get_graph_sensor_from_channel(run, sensor, type, position, get_bkg_vs_vbias); };
    std::map<std::string, TGraphErrors *> get_sensor_pPDE_vs_bkg(std::string run, std::string sensor, std::string type, std::string position) { return get_graph_sensor_from_channel(run, sensor, type, position, get_pPDE_vs_bkg); };
    std::map<std::string, TGraphErrors *> get_stability_check(std::vector<std::string> run_list);
    std::map<std::string, TGraphErrors *> get_stability_check(std::vector<std::pair<int, std::string>> run_list);
    void show_stability_check(std::vector<std::string> run_list);
    void show_stability_check(std::vector<std::pair<int, std::string>> run_list);
    //  Implementations
    //  ---  I/O
    //  --- --- database
    void read_database(std::string database_input_file)
    {
      //  Start reading the file
      std::ifstream data_stream(database_input_file);
      std::string current_line;
      while (std::getline(data_stream, current_line))
      {
        //  Skip comment characters
        if (current_line[0] == '#' || current_line[0] == ' ')
          continue;
        //  Read database
        std::stringstream string_in_stream(current_line);
        std::unordered_map<std::string, std::string> data_by_field;
        std::string current_data;
        for (auto current_field : fields)
        {
          string_in_stream >> current_data;
          data_by_field[current_field] = current_data;
        }
        //  Record quality
        auto current_entry = database.find_or_create_node(-1, data_by_field["board"], data_by_field["step"], static_cast<std::string>("laser"), data_by_field["run"]);
        database.find_or_create_node(current_entry, data_by_field["status"]);

        //  Record channels
        std::vector<std::string> channels;
        for (size_t i = 0; i < data_by_field["channels"].length(); i += 3)
        {
          // Extract the two-character group
          char letter = data_by_field["channels"][i];
          char symbol = data_by_field["channels"][i + 1];
          if (symbol == '*')
          {
            // Expand when encountering '*'
            channels.push_back(std::string(1, letter) + '1');
            channels.push_back(std::string(1, letter) + '2');
            channels.push_back(std::string(1, letter) + '3');
            channels.push_back(std::string(1, letter) + '4');
          }
          else
          {
            // Add the two-character group as is
            channels.push_back(std::string(1, letter) + symbol);
          }
        }
        for (auto current_channel : channels)
          info_database[database.find_or_create_node(current_entry, static_cast<std::string>("channels"))].push_back(current_channel);
      }
    }
    void read_sub_run_database(std::string database_input_file, std::string type = "target")
    {
      //  Start reading the file
      std::ifstream data_stream(database_input_file);
      if (!data_stream)
        return;
      std::string current_line;
      while (std::getline(data_stream, current_line))
      {
        //  Skip comment characters
        if (current_line[0] == '#' || current_line[0] == ' ')
          continue;
        //  Read database
        std::stringstream string_in_stream(current_line);
        std::unordered_map<std::string, std::string> data_by_field;
        std::string current_data;
        for (auto current_field : sub_run_fields)
        {
          string_in_stream >> current_data;
          data_by_field[current_field] = current_data;
        }
        //  Recover global run and sub_run number
        std::string global_run_tag = data_by_field["directory"].substr(40, 15);
        auto run_tags = database.get_name_ids(global_run_tag);
        if (run_tags.size() == 0)
        {
          cerr << "[WARNING][datbase::laser::read_sub_run_database] no run entry in general database, abort loading sub_run" << endl;
          return;
        }
        if (run_tags.size() > 1)
        {
          cerr << "[WARNING][datbase::laser::read_sub_run_database] too many run entries in general database, abort loading sub_run" << endl;
          return;
        }
        auto current_entry_father = run_tags[0];
        std::string sub_run_tag = data_by_field["directory"].substr(64, 15);
        std::string channel = data_by_field["channel"];
        auto current_entry = database.find_or_create_node(current_entry_father, sub_run_tag);
        database.find_or_create_node(current_entry, channel);
        for (auto current_field : sub_run_fields)
          info_database[database.find_or_create_node(current_entry, current_field)].push_back(data_by_field[current_field]);
        info_database[database.find_or_create_node(current_entry, static_cast<std::string>("type"))].push_back(type);
        info_database[database.find_or_create_node(current_entry, static_cast<std::string>("global_run_tag"))].push_back(global_run_tag);
      }
    }
    void read_all_sub_runs(std::string target_run_tag)
    {
      //  Load the reference sensor
      read_sub_run_database(basedir + "/" + target_run_tag + "/database.reference.A1.txt", "reference");
      //  Load the target sensors
      auto run_tags = database.get_name_ids(target_run_tag);
      if (run_tags.size() == 0)
      {
        cerr << "[WARNING][datbase::laser::read_sub_run_database] no run entry in general database, abort loading sub_run" << endl;
        return;
      }
      if (run_tags.size() > 1)
      {
        cerr << "[WARNING][datbase::laser::read_sub_run_database] too many run entries in general database, abort loading sub_run" << endl;
        return;
      }
      auto channels_id = database.find_or_create_node(run_tags[0], static_cast<std::string>("channels"));
      for (auto current_channel : info_database[channels_id])
        read_sub_run_database(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".txt", "target");
      for (auto current_channel : info_database[channels_id])
        read_sub_run_database(basedir + "/" + target_run_tag + "/database.proto." + current_channel + ".txt", "proto");
    }
    void read_all_sub_runs()
    {
      auto all_laser_runs = get_all_laser_runs();
      for (auto current_run_entry : all_laser_runs)
        read_all_sub_runs(current_run_entry.second);
    }
    //  --- --- data
    void download_run(std::string campaign, std::string target_run_tag = ".")
    {
      std::system(("mkdir -p " + basedir + "/" + target_run_tag + "/").c_str());
      //  Download the reference sensor
      std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.reference." + "A1" + ".signal.tgz " + basedir + "/" + target_run_tag + "/").c_str());
      std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.reference." + "A1" + ".txt " + basedir + "/" + target_run_tag + "/").c_str());
      unzip_sub_run(basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".signal.tgz", basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".signal/");
      read_sub_run_database(basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".txt");
      auto current_run = get_run_infos(target_run_tag);
      for (auto current_channel : current_run["channels"])
      {
        std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz " + basedir + "/" + target_run_tag + "/").c_str());
        std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.target." + current_channel + ".txt " + basedir + "/" + target_run_tag + "/").c_str());
        unzip_sub_run(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz", basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal/");
        read_sub_run_database(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".txt");
      }
    }
    void unzip_sub_run(std::string target_file, std::string target_dir = ".", std::string campaign = "2024-laser-window")
    {
      std::system(("mkdir -p " + target_dir).c_str());
      std::system(("tar -xf " + target_file + " -C " + target_dir).c_str());
      std::system(("mv " + target_dir + "/home/eic/DATA/" + campaign + "/actual/*/subruns/* " + target_dir).c_str());
      std::system(("rm -r " + target_dir + "/home").c_str());
      std::system(("rm " + target_file).c_str());
    }
    void unzip_all_sub_runs(std::string target_run_tag)
    {
      //  Load the reference sensor
      unzip_sub_run(basedir + "/" + target_run_tag + "/database.reference." + "A1.signal" + ".tgz", basedir + "/" + target_run_tag + "/database.reference." + "A1.signal/");
      //  Load the target sensors
      auto current_run = get_run_infos(target_run_tag);
      for (auto current_channel : current_run["channels"])
        unzip_sub_run(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz", basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal/");
    }
    void update_local_data_repository(std::string current_campaign = "2024-laser-window")
    {
      auto all_downloaded_runs = utility::get_folders(basedir);
      auto all_recorded_laser_runs = get_all_laser_runs();
      for (auto current_run : all_recorded_laser_runs)
      {
        if (std::find(all_downloaded_runs.begin(), all_downloaded_runs.end(), current_run.second) != all_downloaded_runs.end())
          continue;
        auto current_run_infos = get_run_infos(current_run.second);
        if (!current_run_infos.count("good"))
          continue;
        download_run(current_campaign, current_run.second);
        unzip_all_sub_runs(current_run.second);
      }
    }
    //  --- Getters
    std::unordered_map<std::string, std::vector<std::string>> get_run_infos(std::string run)
    {
      std::unordered_map<std::string, std::vector<std::string>> result;
      auto list_of_infos = database.get_children(database.get_name_ids(run)[0]);
      for (auto current_info : list_of_infos)
        result[current_info.second] = info_database[current_info.first];
      return result;
    }
    std::unordered_map<std::string, std::vector<std::string>> get_sub_run_infos(std::string sub_run)
    {
      std::unordered_map<std::string, std::vector<std::string>> result;
      auto list_of_infos = database.get_children(-1, 6, {{5, sub_run}});
      for (auto current_info : list_of_infos)
        result[current_info.second] = info_database[current_info.first];
      return result;
    }
    std::unordered_map<std::string, std::vector<std::string>> get_sub_run_infos(int sub_run_id)
    {
      std::unordered_map<std::string, std::vector<std::string>> result;
      auto list_of_infos = database.get_children(sub_run_id);
      for (auto current_info : list_of_infos)
        result[current_info.second] = info_database[current_info.first];
      return result;
    }
    std::vector<std::pair<int, std::string>> get_all_laser_runs(std::string quality)
    {
      auto all_recorded_laser_runs = get_all_laser_runs();
      std::vector<std::pair<int, std::string>> filtered_run_list;
      for (auto current_run : all_recorded_laser_runs)
      {
        auto current_run_infos = get_run_infos(current_run.second);
        if (!current_run_infos.count(quality))
          continue;
        filtered_run_list.push_back(current_run);
      }
      return filtered_run_list;
    }
    std::pair<int, std::string> get_latest_laser_run(std::string board, std::string step)
    {
      auto all_laser_runs = get_all_laser_runs(board, step);
      std::sort(all_laser_runs.begin(), all_laser_runs.end(),
                [](const std::pair<int, std::string> &a, const std::pair<int, std::string> &b)
                { return a.second > b.second; });
      return all_laser_runs[0];
    }
    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> get_sub_runs_infos(std::string run, std::string channel, std::string vbias, std::string type, std::string position)
    {
      std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> result;
      auto sub_run_list = database.get_children(database.get_name_ids(run)[0]);
      for (auto current_sub_run : sub_run_list)
      {
        auto current_sub_run_infos = get_sub_run_infos(current_sub_run.first);
        if (current_sub_run_infos.empty())
          continue;
        if (channel != current_sub_run_infos["channel"][0])
          continue;
        if (vbias != current_sub_run_infos["vbias"][0])
          continue;
        if (position != current_sub_run_infos["notes"][0])
          continue;
        if (type != current_sub_run_infos["type"][0])
          continue;
        result[current_sub_run.second] = current_sub_run_infos;
      }
      return result;
    }
    //  --- General utilities
    std::map<std::string, TGraphErrors *> get_graph_sensor_from_channel(std::string run, std::string sensor, std::string type, std::string position, TGraphErrors *(*graph_getter)(std::string, std::string, std::string, std::string))
    {
      std::vector<TGraphErrors *> graphs_buffer;
      for (auto current_channel : sensor_to_channels[sensor])
        graphs_buffer.push_back(graph_getter(run, current_channel, type, position));
      return graphutils::average_graphs(graphs_buffer);
    }
    //  --- Measurements
    std::array<std::array<float, 2>, 2> get_sub_run_value(std::string sub_run)
    {
      std::array<std::array<float, 2>, 2> result; // 2x2 array to store signal and background values and errors

      // Construct the file path based on parameters and open the ROOT file
      auto sub_run_infos = get_sub_run_infos(sub_run);
      auto gneeral_run = get_sub_run_infos(sub_run);
      TFile *target_file = new TFile((basedir + "/" + sub_run_infos["global_run_tag"][0] + "/database." + sub_run_infos["type"][0] + "." + sub_run_infos["channel"][0] + ".signal/" + sub_run + "/decoded/signal.root").c_str());

      //  Check the file is there
      if (target_file->IsZombie())
        return {{{0., 0.}, {0., 0.}}};

      // Retrieve the histogram "hDelta" from the ROOT file
      auto current_hDelta = (TH1F *)(target_file->Get("hDelta"));

      // Initialize variables to hold errors for background and full signal calculations
      double bkg_err_1, bkg_err_2, full_err;

      // Compute background values and their errors for two separate regions within the histogram
      auto bkg_val_1 = current_hDelta->IntegralAndError(19, 35, bkg_err_1); // Background region 1
      auto bkg_val_2 = current_hDelta->IntegralAndError(67, 83, bkg_err_2); // Background region 2

      // Compute full signal value and error over the main signal region
      auto full_val = current_hDelta->IntegralAndError(35, 67, full_err);

      // Combine the two background values and calculate the total background error
      float bkg_val = (bkg_val_1 + bkg_val_2);                                    // Sum of background regions
      float bkg_err = TMath::Sqrt(bkg_err_1 * bkg_err_1 + bkg_err_2 * bkg_err_2); // Combined background error

      // Calculate the signal value and its error by subtracting background from full signal
      float sig_val = full_val - bkg_val;
      float sig_err = TMath::Sqrt(bkg_err * bkg_err + full_err * full_err); // Total error for signal

      // Convert background values and errors to nanoseconds using a scaling factor
      bkg_val *= 1.e9 / (100 * coarse_to_ns);
      bkg_err *= 1.e9 / (100 * coarse_to_ns);

      // Clean up by deleting the histogram and the ROOT file object to free memory
      delete current_hDelta;
      delete target_file;

      // Store signal and background values and errors in the result array
      result[0] = {sig_val, sig_err}; // Signal {value, error}
      result[1] = {bkg_val, bkg_err}; // Background {value, error}

      // Return the result array with signal and background data
      return result;
    }
    TGraphErrors *get_pPDE_vs_vbias(std::string run, std::string channel, std::string type, std::string position)
    {
      // Create a new TGraphErrors object to store the graph with error bars
      auto result = new TGraphErrors();

      // Set the x-axis and y-axis titles for the graph
      result->GetXaxis()->SetTitle("V_{bias} (V)");          // X-axis: bias voltage (V)
      result->GetYaxis()->SetTitle("Pseudo-efficiency (%)"); // Y-axis: pseudo-efficiency percentage

      // Retrieve all sub-run identifiers associated with the specified run
      auto sub_run_list = database.get_children(database.get_name_ids(run)[0]);

      // Loop through each sub-run's data, unpacking sub-run ID and name
      for (auto [current_sub_run_id, current_sub_run_name] : sub_run_list)
      {
        // Get information of the current sub-run
        auto current_sub_run_info = get_sub_run_infos(current_sub_run_id);

        // Skip this sub-run if its position, type, or channel does not match the specified parameters
        if (current_sub_run_info.empty())
          continue;
        if (current_sub_run_info["notes"][0] != position)
          continue;
        if (current_sub_run_info["type"][0] != type)
          continue;
        if (current_sub_run_info["channel"][0] != channel)
          continue;

        // Get the current point index in the TGraphErrors (used to add new points)
        auto iPnt = result->GetN();

        // Convert the "vbias" value from string to double for use as the x-axis value
        auto x_value = std::stod(current_sub_run_info["vbias"][0]);

        // Retrieve the y-value (pseudo-efficiency) and its associated error for the sub-run
        auto y_value = get_sub_run_sig(current_sub_run_name);

        // Set the x and y values for the current point in the graph
        result->SetPoint(iPnt, x_value, y_value[0]);

        // Set the y-axis error for the current point in the graph
        result->SetPointError(iPnt, 0, y_value[1]);
      }
      // Return the completed TGraphErrors object with the plotted data and errors
      return graphutils::average_same_x(result)["ave_err"];
    }
    TGraphErrors *get_bkg_vs_vbias(std::string run, std::string channel, std::string type, std::string position)
    {
      // Create a new TGraphErrors object to store the graph with error bars
      auto result = new TGraphErrors();

      // Set the x-axis and y-axis titles for the graph
      result->GetXaxis()->SetTitle("V_{bias} (V)");          // X-axis: bias voltage (V)
      result->GetYaxis()->SetTitle("Pseudo-efficiency (%)"); // Y-axis: pseudo-efficiency percentage

      // Retrieve all sub-run identifiers associated with the specified run
      auto sub_run_list = database.get_children(database.get_name_ids(run)[0]);

      // Loop through each sub-run's data, unpacking sub-run ID and name
      for (auto [current_sub_run_id, current_sub_run_name] : sub_run_list)
      {
        // Get information of the current sub-run
        auto current_sub_run_info = get_sub_run_infos(current_sub_run_id);

        // Skip this sub-run if its position, type, or channel does not match the specified parameters
        if (current_sub_run_info.empty())
          continue;
        if (current_sub_run_info["notes"][0] != position)
          continue;
        if (current_sub_run_info["type"][0] != type)
          continue;
        if (current_sub_run_info["channel"][0] != channel)
          continue;

        // Get the current point index in the TGraphErrors (used to add new points)
        auto iPnt = result->GetN();

        // Convert the "vbias" value from string to double for use as the x-axis value
        auto x_value = std::stod(current_sub_run_info["vbias"][0]);

        // Retrieve the y-value (pseudo-efficiency) and its associated error for the sub-run
        auto y_value = get_sub_run_bkg(current_sub_run_name);

        // Set the x and y values for the current point in the graph
        result->SetPoint(iPnt, x_value, y_value[0]);

        // Set the y-axis error for the current point in the graph
        result->SetPointError(iPnt, 0, y_value[1]);
      }
      // Return the completed TGraphErrors object with the plotted data and errors
      return graphutils::average_same_x(result)["ave_err"];
    }
    TGraphErrors *get_pPDE_vs_bkg(std::string run, std::string channel, std::string type, std::string position)
    {
      // Create a new TGraphErrors object to store the graph with error bars
      auto result = new TGraphErrors();

      // Set the x-axis and y-axis titles for the graph
      result->GetXaxis()->SetTitle("V_{bias} (V)");          // X-axis: bias voltage (V)
      result->GetYaxis()->SetTitle("Pseudo-efficiency (%)"); // Y-axis: pseudo-efficiency percentage

      // Retrieve all sub-run identifiers associated with the specified run
      auto sub_run_list = database.get_children(database.get_name_ids(run)[0]);

      // Loop through each sub-run's data, unpacking sub-run ID and name
      for (auto [current_sub_run_id, current_sub_run_name] : sub_run_list)
      {
        // Get information of the current sub-run
        auto current_sub_run_info = get_sub_run_infos(current_sub_run_id);

        // Skip this sub-run if its position, type, or channel does not match the specified parameters
        if (current_sub_run_info.empty())
          continue;
        if (current_sub_run_info["notes"][0] != position)
          continue;
        if (current_sub_run_info["type"][0] != type)
          continue;
        if (current_sub_run_info["channel"][0] != channel)
          continue;

        // Get the current point index in the TGraphErrors (used to add new points)
        auto iPnt = result->GetN();

        // Convert the "vbias" value from string to double for use as the x-axis value
        auto x_value = get_sub_run_bkg(current_sub_run_name);

        // Retrieve the y-value (pseudo-efficiency) and its associated error for the sub-run
        auto y_value = get_sub_run_sig(current_sub_run_name);

        // Set the x and y values for the current point in the graph
        result->SetPoint(iPnt, x_value[0], y_value[0]);

        // Set the y-axis error for the current point in the graph
        result->SetPointError(iPnt, x_value[1], y_value[1]);
      }
      // Return the completed TGraphErrors object with the plotted data and errors
      return result; // graphutils::average(result, 3);
    }
    std::map<std::string, TGraphErrors *> get_stability_check(std::vector<std::string> run_list)
    {
      std::map<std::string, TGraphErrors *> result;
      for (auto current_position : laser_positions)
      {
        result["pPDE::" + current_position] = new TGraphErrors();
        result["bkg::" + current_position] = new TGraphErrors();
        auto i_pnt = -1;
        for (auto current_run : run_list)
        {
          i_pnt++;
          auto current_timestamp = std::stoll(current_run);
          auto current_pPDE = get_pPDE_vs_vbias(current_run, "A1", "reference", current_position);
          auto current_bkg = get_bkg_vs_vbias(current_run, "A1", "reference", current_position);
          result["pPDE::" + current_position]->SetPoint(i_pnt, current_timestamp, current_pPDE->GetPointY(2));
          result["pPDE::" + current_position]->SetPointError(i_pnt, 0, current_pPDE->GetErrorY(2));
          result["bkg::" + current_position]->SetPoint(i_pnt, current_timestamp, current_pPDE->GetPointY(2));
          result["bkg::" + current_position]->SetPointError(i_pnt, 0, current_pPDE->GetErrorY(2));
        }
      }
      return result;
    }
    std::map<std::string, TGraphErrors *> get_stability_check(std::vector<std::pair<int, std::string>> run_list)
    {
      std::vector<std::string> run_list_string(run_list.size());
      std::transform(
          run_list.begin(),
          run_list.end(),
          run_list_string.begin(),
          [](const std::pair<int, std::string> &pair)
          { return pair.second; });
      return get_stability_check(run_list_string);
    }
    void show_stability_check(std::vector<std::string> run_list)
    {
      auto all_graphs = get_stability_check(run_list);
      TCanvas *plot_all = new TCanvas();
      plot_all->Divide(1, 2);
      plot_all->cd(1);
      for (auto current_position : laser_positions)
        all_graphs["pPDE::" + current_position]->Draw("ALP");
      //  result["bkg::" + current_position];
      plot_all->cd(2);
    }
    void show_stability_check(std::vector<std::pair<int, std::string>> run_list)
    {
      std::vector<std::string> run_list_string(run_list.size());
      std::transform(
          run_list.begin(),
          run_list.end(),
          run_list_string.begin(),
          [](const std::pair<int, std::string> &pair)
          { return pair.second; });
      return show_stability_check(run_list_string);
    }
  }
  /*
    namespace laser
    {
      //  --- Constants ---
      //  Physics constants
      const double coarse_to_ns = 3.125; // ns
      //  --- Utility data structures ---
      //  --- Data structure in database reference
      std::vector<std::string> fields = {
          "run", "status", "step", "board", "channels"};
      std::vector<std::string> criteria = {
          "run", "status", "channels"};
      std::vector<std::string> sub_run_fields = {
          "channel", "bcrconfig", "opmode", "deltathreshold", "vbias", "directory", "notes"};
      std::vector<std::string> sub_run_criteria = {
          "channel", "bcrconfig", "opmode", "deltathreshold", "vbias", "directory", "notes"};
      //  --- Global variables
      std::string basedir = "./Data/";
      //  --- Board, Status, Run, Information
      std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>>> database_memory;
      //  --- Run, Type, Channel, sub_run, Information
      std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>>>>> sub_runs_database_memory;


      //  --- Declarations ---
      //  --- I/O
      bool look_for_run(std::string global_run_tag);
      std::unordered_map<std::string, std::vector<std::string>> get_run(std::string global_run_tag);
      void download_run(std::string campaign, std::string target_run_tag = ".");
      void unzip_sub_run(std::string target_file, std::string target_dir = ".");
      void unzip_all_sub_runs(std::string target_run_tag);
      //  --- Quality assurance
      void show_database();
      void show_sub_run_database();
      //  --- Getters
      std::vector<std::string> get_channels(std::string board, std::string status, std::string run);
      std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> get_sub_runs(std::string run, std::string type, std::string channel);

      //  Implementations
      //  ---  I/O
      void read_database(std::string database_input_file)
      {
        //  Start reading the file
        std::ifstream data_stream(database_input_file);
        std::string current_line;
        while (std::getline(data_stream, current_line))
        {
          //  Skip comment characters
          if (current_line[0] == '#' || current_line[0] == ' ')
            continue;
          //  Read database
          std::stringstream string_in_stream(current_line);
          std::unordered_map<std::string, std::string> data_by_field;
          std::string current_data;
          for (auto current_field : fields)
          {
            string_in_stream >> current_data;
            data_by_field[current_field] = current_data;
          }
          //  Record quality
          database_memory[data_by_field["board"]][data_by_field["step"]][data_by_field["run"]]["status"].push_back(data_by_field["status"]);
          //  Record channels
          std::vector<std::string> channels;
          for (size_t i = 0; i < data_by_field["channels"].length(); i += 3)
          {
            // Extract the two-character group
            char letter = data_by_field["channels"][i];
            char symbol = data_by_field["channels"][i + 1];
            if (symbol == '*')
            {
              // Expand when encountering '*'
              channels.push_back(std::string(1, letter) + '1');
              channels.push_back(std::string(1, letter) + '2');
              channels.push_back(std::string(1, letter) + '3');
              channels.push_back(std::string(1, letter) + '4');
            }
            else
            {
              // Add the two-character group as is
              channels.push_back(std::string(1, letter) + symbol);
            }
          }
          for (auto current_channel : channels)
            database_memory[data_by_field["board"]][data_by_field["step"]][data_by_field["run"]]["channels"].push_back(current_channel);
        }
      }
      bool look_for_run(std::string global_run_tag)
      {
        for (auto [current_board, all_statuses] : database_memory)
          for (auto [current_status, all_runs] : all_statuses)
          {
            auto current_run = all_runs.find(global_run_tag);
            if (current_run != all_runs.end())
              return true;
          }
        return false;
      }
      std::unordered_map<std::string, std::vector<std::string>> get_run(std::string global_run_tag)
      {
        for (auto [current_board, all_statuses] : database_memory)
          for (auto [current_status, all_runs] : all_statuses)
          {
            auto current_run = all_runs.find(global_run_tag);
            if (current_run != all_runs.end())
              return current_run->second;
          }
        return {};
      }
      void read_sub_run_database(std::string database_input_file, std::string type = "target")
      {
        //  Check the run is there in the global database

        //  Start reading the file
        std::ifstream data_stream(database_input_file);
        if (!data_stream)
          return;
        std::string current_line;
        while (std::getline(data_stream, current_line))
        {
          //  Skip comment characters
          if (current_line[0] == '#' || current_line[0] == ' ')
            continue;
          //  Read database
          std::stringstream string_in_stream(current_line);
          std::unordered_map<std::string, std::string> data_by_field;
          std::string current_data;
          for (auto current_field : sub_run_fields)
          {
            string_in_stream >> current_data;
            data_by_field[current_field] = current_data;
          }
          //  Recover global run and sub_run number
          std::string global_run_tag = data_by_field["directory"].substr(40, 15);
          if (!look_for_run(global_run_tag))
            cerr << "[WARNING][datbase::laser::read_sub_run_database] no run entry in general databse, loading sub_run data as orphan" << endl;
          std::string sub_run_tag = data_by_field["directory"].substr(64, 15);
          std::string channel = data_by_field["channel"];
          for (auto current_field : sub_run_fields)
            sub_runs_database_memory[global_run_tag][type][channel][sub_run_tag][current_field].push_back(data_by_field[current_field]);
        }
      }
      void read_all_sub_runs(std::string target_run_tag)
      {
        if (!look_for_run(target_run_tag))
          cerr << "[WARNING][datbase::laser::read_all_sub_runs] no run entry in general databse, loading sub_run data as orphan" << endl;
        //  Load the reference sensor
        read_sub_run_database(basedir + "/" + target_run_tag + "/database.reference.A1.txt", "reference");
        //  Load the target sensors
        auto current_run = get_run(target_run_tag);
        for (auto current_channel : current_run["channels"])
          read_sub_run_database(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".txt", "target");
        for (auto current_channel : current_run["channels"])
          read_sub_run_database(basedir + "/" + target_run_tag + "/database.proto." + current_channel + ".txt", "target");
      }
      void read_all_sub_runs()
      {
        for (auto [current_board, all_statuses] : database_memory)
          for (auto [current_status, all_runs] : all_statuses)
            for (auto [current_run, all_infos] : all_runs)
              read_all_sub_runs(current_run);
      }
      void download_run(std::string campaign, std::string target_run_tag = ".")
      {
        if (!look_for_run(target_run_tag))
          cerr << "[WARNING][datbase::laser::download_run] no run entry in general databse, loading sub_run data as orphan" << endl;
        std::system(("mkdir -p " + basedir + "/" + target_run_tag + "/").c_str());
        //  Download the reference sensor
        std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.reference." + "A1" + ".signal.tgz " + basedir + "/" + target_run_tag + "/").c_str());
        std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.reference." + "A1" + ".txt " + basedir + "/" + target_run_tag + "/").c_str());
        unzip_sub_run(basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".signal.tgz", basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".signal/");
        read_sub_run_database(basedir + "/" + target_run_tag + "/database.reference." + "A1" + ".txt");
        auto current_run = get_run(target_run_tag);
        for (auto current_channel : current_run["channels"])
        {
          std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz " + basedir + "/" + target_run_tag + "/").c_str());
          std::system(("scp eic@eicdesk01:/home/eic/DATA/" + campaign + "/actual/" + target_run_tag + "/database.target." + current_channel + ".txt " + basedir + "/" + target_run_tag + "/").c_str());
          unzip_sub_run(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz", basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal/");
          read_sub_run_database(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".txt");
        }
      }
      void unzip_sub_run(std::string target_file, std::string target_dir = ".")
      {
        std::system(("mkdir -p " + target_dir).c_str());
        std::system(("tar -xf " + target_file + " -C " + target_dir).c_str());
        std::system(("mv " + target_dir + "/home/eic/DATA/2024-laser-window/actual/\*\/sub_runs/\* " + target_dir).c_str());
        std::system(("rm -r " + target_dir + "/home").c_str());
        std::system(("rm " + target_file).c_str());
      }
      void unzip_all_sub_runs(std::string target_run_tag)
      {
        if (!look_for_run(target_run_tag))
          cerr << "[WARNING][datbase::laser::unzip_all_sub_runs] no run entry in general databse, loading sub_run data as orphan" << endl;
        //  Load the reference sensor
        unzip_sub_run(basedir + "/" + target_run_tag + "/database.reference." + "A1.signal" + ".tgz", basedir + "/" + target_run_tag + "/database.reference." + "A1.signal/");
        //  Load the target sensors
        auto current_run = get_run(target_run_tag);
        for (auto current_channel : current_run["channels"])
          unzip_sub_run(basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal.tgz", basedir + "/" + target_run_tag + "/database.target." + current_channel + ".signal/");
      }
      //  --- Quality assurance
      void show_database()
      {
        cout << "[INFO][database::show_database] Starting print all database entries" << endl;
        for (auto [current_board, all_statuses] : database_memory)
          for (auto [current_status, all_runs] : all_statuses)
            for (auto [current_run, all_informations] : all_runs)
              for (auto [current_information, all_values] : all_informations)
                for (auto current_value : all_values)
                  cout << "[INFO][database::show_database] [Board] " << current_board << " - [Status] " << current_status << " - [run] " << current_run << " - [Info] " << current_information << " - [Value] " << current_value << endl;
      }
      void show_sub_run_database()
      {
        cout << "[INFO][database::show_sub_run_database] Starting print all database entries" << endl;
        for (auto [current_run, all_types] : sub_runs_database_memory)
          for (auto [current_type, all_channels] : all_types)
            for (auto [current_channel, all_sub_runs] : all_channels)
              for (auto [current_sub_run, all_informations] : all_sub_runs)
                for (auto [current_information, all_values] : all_informations)
                  for (auto current_value : all_values)
                    cout << "[INFO][database::show_sub_run_database] [Run] " << current_run << " - [Type] " << current_type << " - [Channel] " << current_channel << " - [sub_run] " << current_sub_run << " - [Info] " << current_information << " - [Value] " << current_value << endl;
      }
      //  --- Getters
      std::vector<std::string> get_channels(std::string board, std::string status, std::string run)
      {
        // TODO: Set limitations and protections for missing infos (no run, no board, no status)
        return database_memory[board][status][run]["channels"];
      }
      std::unordered_map<std::string, std::unordered_map<std::string, std::vector<std::string>>> get_sub_runs(std::string run, std::string type, std::string channel)
      {
        return sub_runs_database_memory[run][type][channel];
      }
      //  --- Measures
      //  Returns values for signal and background in large coincidence region
      std::array<std::array<float, 2>, 2> get_sub_run_value(std::string global_run, std::string sub_run, std::string channel, std::string type = "target")
      {
        std::array<std::array<float, 2>, 2> result;
        TFile *target_file = new TFile((basedir + "/" + global_run + "/database." + type + "." + channel + ".signal/" + sub_run + "/decoded/signal.root").c_str());
        auto current_hDelta = (TH1F *)(target_file->Get("hDelta"));
        double bkg_err_1, bkg_err_2, full_err;
        auto bkg_val_1 = current_hDelta->IntegralAndError(19, 35, bkg_err_1);
        auto bkg_val_2 = current_hDelta->IntegralAndError(67, 83, bkg_err_2);
        auto full_val = current_hDelta->IntegralAndError(35, 67, full_err);
        float bkg_val = (bkg_val_1 + bkg_val_2);
        float bkg_err = TMath::Sqrt(bkg_err_1 * bkg_err_1 + bkg_err_2 * bkg_err_2);
        float sig_val = full_val - bkg_val;
        float sig_err = TMath::Sqrt(bkg_err * bkg_err + full_err * full_err);
        bkg_val *= 1.e9 / (100 * coarse_to_ns);
        bkg_err *= 1.e9 / (100 * coarse_to_ns);
        delete current_hDelta;
        delete target_file;
        result[0] = {sig_val, sig_err};
        result[1] = {bkg_val, bkg_err};
        return result;
      }
      std::array<float, 2> get_sub_run_sig(std::string global_run, std::string sub_run, std::string channel, std::string type = "target")
      {
        return get_sub_run_value(global_run, sub_run, channel, type)[0];
      }
      std::array<float, 2> get_sub_run_bkg(std::string global_run, std::string sub_run, std::string channel, std::string type = "target")
      {
        return get_sub_run_value(global_run, sub_run, channel, type)[1];
      }
      //  --- --- General functions
      //  --- --- --- Getters of curves
      TGraphErrors *get_pPDE_vs_vbias(std::string run, std::string channel, std::string type, std::string position = "center")
      {
        auto result = new TGraphErrors();
        result->GetXaxis()->SetTitle("V_{bias} (V)");
        result->GetYaxis()->SetTitle("Pseudo-efficiency (%)");
        auto all_sub_runs_data = get_sub_runs(run, type, channel);
        for (auto [current_sub_run, all_informations] : all_sub_runs_data)
        {
          if (all_informations["notes"][0] != position)
            continue;
          auto iPnt = result->GetN();
          auto x_value = std::stod(all_informations["vbias"][0]);
          auto y_value = get_sub_run_sig(run, current_sub_run, channel, type);
          result->SetPoint(iPnt, x_value, y_value[0]);
          result->SetPointError(iPnt, 0, y_value[1]);
        }
        return result;
      }
      TGraphErrors *get_bkg_vs_vbias(std::string run, std::string channel, std::string type, std::string position = "center")
      {
        auto result = new TGraphErrors();
        result->GetXaxis()->SetTitle("V_{bias} (V)");
        result->GetYaxis()->SetTitle("DCR (Hz)");
        auto all_sub_runs_data = get_sub_runs(run, type, channel);
        for (auto [current_sub_run, all_informations] : all_sub_runs_data)
        {
          if (all_informations["notes"][0] != position)
            continue;
          auto iPnt = result->GetN();
          auto x_value = std::stod(all_informations["vbias"][0]);
          auto y_value = get_sub_run_bkg(run, current_sub_run, channel, type);
          result->SetPoint(iPnt, x_value, y_value[0]);
          result->SetPointError(iPnt, 0, y_value[1]);
        }
        return result;
      }
      TGraphErrors *get_pPDE_vs_bkg(std::string run, std::string channel, std::string type, std::string position = "center")
      {
        auto result = new TGraphErrors();
        result->GetXaxis()->SetTitle("DCR (Hz)");
        result->GetYaxis()->SetTitle("Pseudo-efficiency (%)");
        auto all_sub_runs_data = get_sub_runs(run, type, channel);
        for (auto [current_sub_run, all_informations] : all_sub_runs_data)
        {
          if (all_informations["notes"][0] != position)
            continue;
          auto iPnt = result->GetN();
          auto x_value = get_sub_run_bkg(run, current_sub_run, channel, type);
          auto y_value = get_sub_run_sig(run, current_sub_run, channel, type);
          result->SetPoint(iPnt, x_value[0], y_value[0]);
          result->SetPointError(iPnt, x_value[1], y_value[1]);
        }
        return result;
      }
      TGraphErrors *get_pPDE_merit_vs_vbias(std::string run, std::string channel, std::string type, std::string position = "center")
      {
        auto result = new TGraphErrors();
        result->GetXaxis()->SetTitle("V_{bias} (V)");
        result->GetYaxis()->SetTitle("Pseudo-efficiency figure of merit (%)");
        auto all_sub_runs_data = get_sub_runs(run, type, channel);
        for (auto [current_sub_run, all_informations] : all_sub_runs_data)
        {
          if (all_informations["notes"][0] != position)
            continue;
          auto iPnt = result->GetN();
          auto x_value = std::stod(all_informations["vbias"][0]);
          auto dcr_value = get_sub_run_bkg(run, current_sub_run, channel, type);
          auto pde_value = get_sub_run_sig(run, current_sub_run, channel, type);
          auto y_value_val = pde_value[0] / dcr_value[0];
          auto y_value_err = y_value_val * std::sqrt((pde_value[1] / pde_value[0]) * (pde_value[1] / pde_value[0]) + (dcr_value[1] / dcr_value[0]) * (dcr_value[1] / dcr_value[0]));
          std::array<float, 2> y_value = {y_value_val, y_value_err};
          result->SetPoint(iPnt, x_value, y_value[0]);
          result->SetPointError(iPnt, 0, y_value[1]);
        }
        return result;
      }

      std::array<float, 2>
      get_pPDE_at_overvoltage(std::string run, std::string type, std::string channel, float overvoltage)
      {
        return database::get_value_at_overvoltage<database::laser::get_pPDE_vs_vbias>(run, type, channel, overvoltage);
      }
      std::array<float, 2>
      get_bkg_at_overvoltage(std::string run, std::string type, std::string channel, float overvoltage)
      {
        return database::get_value_at_overvoltage<database::laser::get_bkg_vs_vbias>(run, type, channel, overvoltage);
      }

      std::vector<TGraphErrors *>
      get_general_TGraphs_pPDE(std::string sensor, std::vector<std::pair<std::string, std::vector<std::string>>> target_list_w_runs, float overvoltage, int marker = kFullCircle, int color = kBlue)
      {
        return get_general_TGraphs<float, get_pPDE_at_overvoltage>(sensor, target_list_w_runs, overvoltage, marker, color);
      }
      std::vector<TGraphErrors *>
      get_general_TGraphs_bkg(std::string sensor, std::vector<std::pair<std::string, std::vector<std::string>>> target_list_w_runs, float overvoltage, int marker = kFullCircle, int color = kBlue)
      {
        return get_general_TGraphs<float, get_bkg_at_overvoltage>(sensor, target_list_w_runs, overvoltage, marker, color);
      }
    }
    */
}
