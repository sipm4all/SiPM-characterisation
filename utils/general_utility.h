//  --- --- ---
//  General utility functions
//
//  Auth: Nicola Rubini
//  Mail: nicola.rubini@bo.infn.it
//
#pragma once

namespace utility
{
  //  Swap two values
  template <typename arg_type>
  inline void swap_values(arg_type &first_element, arg_type &second_element);
  //  Square sum
  template <typename arg_type>
  inline arg_type sq_sum(arg_type value);
  template <typename arg_type, typename... Args>
  inline arg_type sq_sum(arg_type first, Args... args);
  //  Average
  // template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
  // std::map<std::string, std::array<T, 2>> average(std::vector<std::array<T, 2>> list_of_measurements, bool skip_unfit_skim = false);
  //  Set precision of the number
  template <typename T>
  T round_digits(T value, int digits);
  //  Concatenate vectors
  template <typename T>
  std::vector<T> merge(const std::vector<T> &vec1, const std::vector<T> &vec2);
  //  Get folder list of a location
  std::vector<std::string> get_folders(const std::string &dir_path);

  //  Read a txt file with some numbers
  void readFileToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap);
  void readTxtToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap);
  void readCsvToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap);

  //  Calculate intercept between lines
  std::array<std::array<float, 2>, 2> get_intercept(float m1, float em1, float q1, float eq1, float m2, float em2, float q2, float eq2);
  std::array<std::array<float, 2>, 2> get_intercept(std::array<float, 2> m1, std::array<float, 2> q1, std::array<float, 2> m2, std::array<float, 2> q2) { return get_intercept(m1[0], m1[1], q1[0], q1[1], m2[0], m2[1], q2[0], q2[1]); };
  std::array<std::array<float, 2>, 2> get_intercept_x(float m, float em, float q, float eq) { return utility::get_intercept(m, em, q, eq, 0, 0, 0, 0); };
  std::array<std::array<float, 2>, 2> get_intercept_y(float m, float em, float q, float eq) { return {q, eq}; };
  std::array<std::array<float, 2>, 2> get_intercept_pol1_pol2(float x0_0, float ex0_0, float x1_0, float ex1_0, float x0_1, float ex0_1, float x1_1, float ex1_1, float x2_1, float ex2_1);
}

template <typename arg_type>
inline void utility::swap_values(arg_type &first_element, arg_type &second_element)
{
  auto memory = second_element;
  second_element = first_element;
  first_element = memory;
}

template <typename arg_type>
inline arg_type utility::sq_sum(arg_type value)
{
  return value * value; // Base case: square the single argument
}

template <typename arg_type, typename... Args>
inline arg_type utility::sq_sum(arg_type first, Args... args)
{
  return first * first + sq_sum(args...); // Recursive call
}

template <typename T>
T utility::round_digits(T value, int digits)
{
  //  Check the passed argument is a number
  static_assert(std::is_arithmetic<T>::value, "Type must be numeric");
  T factor = std::pow(10, digits);
  return std::ceil(value * factor) / factor;
}

template <typename T>
std::vector<T> utility::merge(const std::vector<T> &vec1, const std::vector<T> &vec2)
{
  std::vector<T> result = vec1;                          // Start with the first vector
  result.insert(result.end(), vec2.begin(), vec2.end()); // Concatenate the second vector
  return result;
}

std::vector<std::string> utility::get_folders(const std::string &dir_path)
{
  std::vector<std::string> folders;
  for (const auto &entry : std::filesystem::directory_iterator(dir_path))
    if (entry.is_directory())
      folders.push_back(entry.path().filename().string());
  return folders;
}

namespace utility
{
  //  Check again
  template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
  std::map<std::string, std::array<T, 2>> average(std::vector<std::array<T, 2>> list_of_measurements, bool skip_unfit_skim = false)
  {
    //  Final result
    std::map<std::string, std::array<T, 2>> result;
    std::vector<std::array<T, 2>> skimmed_data;

    //  Skim dataset from dangerous values
    skimmed_data = list_of_measurements;
    if (!skip_unfit_skim)
    {
      skimmed_data.clear();
      for (auto current_measurement : list_of_measurements)
      {
        if (std::isnan(current_measurement[0]) || std::isnan(current_measurement[1]))
          continue;
        if (std::isinf(current_measurement[0]) || std::isinf(current_measurement[1]))
          continue;
        skimmed_data.push_back(current_measurement);
      }
    }

    //  Measure average and error on the average
    for (auto current_measurement : skimmed_data)
    {
      result["ave_err"][0] += current_measurement[0];
      result["ave_rms"][0] += current_measurement[0];
      result["ave_sqe"][0] += current_measurement[0];
      result["sqa_err"][0] += current_measurement[0] * current_measurement[0];
      result["sqa_rms"][0] += current_measurement[0] * current_measurement[0];
      result["sqa_sqe"][0] += current_measurement[0] * current_measurement[0];
      result["ave_err"][1] += current_measurement[1];
      // result["ave_rms"][1] += current_measurement[1];
      result["ave_sqe"][1] += current_measurement[1] * current_measurement[1];
      result["sqa_err"][1] += current_measurement[1];
      // result["sqa_rms"][1] += current_measurement[1];
      result["sqa_sqe"][1] += current_measurement[1] * current_measurement[1];
    }
    result["ave_err"][0] /= skimmed_data.size();
    result["ave_rms"][0] /= skimmed_data.size();
    result["ave_sqe"][0] /= skimmed_data.size();
    result["sqa_err"][0] /= skimmed_data.size();
    result["sqa_rms"][0] /= skimmed_data.size();
    result["sqa_sqe"][0] /= skimmed_data.size();
    result["sqa_err"][0] = sqrt(result["sqa_err"][0]);
    result["sqa_rms"][0] = sqrt(result["sqa_rms"][0]);
    result["sqa_sqe"][0] = sqrt(result["sqa_sqe"][0]);
    result["ave_err"][1] /= skimmed_data.size();
    // result["ave_rms"][1] += current_measurement[1];
    result["ave_sqe"][1] = sqrt(result["ave_sqe"][1]);
    result["sqa_err"][1] /= skimmed_data.size();
    // result["sqa_rms"][1] += current_measurement[1];
    result["sqa_sqe"][1] = sqrt(result["sqa_sqe"][1]);

    //  Measure RMS

    return result;
  }
}

void utility::readCsvToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap)
{
  if (!std::filesystem::exists(filename))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File does not exist: " + filename);
  }

  std::ifstream file(filename);
  if (!file)
  {
    throw std::runtime_error("[ERROR][readFileToMap] Could not open file " + filename);
  }

  // Read header line
  std::string headerLine;
  if (!std::getline(file, headerLine))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File " + filename + " is empty or missing header row.");
  }

  // Ensure the string has at least 3 characters to check
  if (filename.size() < 3)
  {
    std::cerr << "[ERROR][utility::readFileToMap] Filename too short to determine file type: " << filename << std::endl;
    return;
  }

  // Parse column titles
  std::vector<std::string> columnTitles;
  std::istringstream headerStream(headerLine);
  std::string title;

  // Split header line by commas
  size_t start = 0, end = 0;
  while ((end = headerLine.find(',', start)) != std::string::npos)
  {
    columnTitles.push_back(std::regex_replace(headerLine.substr(start, end - start), std::regex("^\\s+|\\s+$"), ""));
    start = end + 1;
  }
  columnTitles.push_back(std::regex_replace(headerLine.substr(start), std::regex("^\\s+|\\s+$"), ""));

  for (const auto &col : columnTitles)
  {
    if (dataMap.find(col) == dataMap.end())
    {
      dataMap[col] = {}; // Initialize vector if not already present
    }
  }

  // Define the allowed characters regex pattern
  const std::regex validPattern(R"(^\s*[+-]?[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?\s*$)");

  // Read remaining lines
  std::string line;
  while (std::getline(file, line))
  {
    std::vector<std::string> rowValues;
    start = 0;

    // Split the line by commas
    while ((end = line.find(',', start)) != std::string::npos)
    {
      rowValues.push_back(line.substr(start, end - start));
      start = end + 1;
    }
    rowValues.push_back(line.substr(start));

    if (rowValues.size() != columnTitles.size())
    {
      throw std::runtime_error("[ERROR][readFileToMap] Mismatched data and header columns in file " + filename);
    }

    bool isValidLine = true;

    // Check each value for allowed characters
    for (const auto &value : rowValues)
    {
      if (!std::regex_match(value, validPattern))
      {
        isValidLine = false;
        break;
      }
    }

    // If the line is valid, add its values to the respective columns
    if (isValidLine)
    {
      for (size_t i = 0; i < columnTitles.size(); ++i)
      {
        dataMap[columnTitles[i]].push_back(rowValues[i]);
      }
    }
    else
    {
      std::cout << "[WARNING] Skipping line for character incompatibility: " << line << std::endl;
    }
  }

  file.close();
}

void utility::readTxtToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap)
{
  if (!std::filesystem::exists(filename))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File does not exist: " + filename);
  }

  std::ifstream file(filename);
  if (!file)
  {
    throw std::runtime_error("[ERROR][readFileToMap] Could not open file " + filename);
  }

  // Read header line
  std::string headerLine;
  if (!std::getline(file, headerLine))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File " + filename + " is empty or missing header row.");
  }

  // Parse column titles
  std::istringstream headerStream(headerLine);
  std::vector<std::string> columnTitles;
  std::string title;
  while (headerStream >> title)
  {
    columnTitles.push_back(title);
    if (dataMap.find(title) == dataMap.end())
    {
      dataMap[title] = {}; // Initialize vector if not already present
    }
  }

  // Define the allowed characters regex pattern
  const std::regex validPattern(R"(^\s*[+-]?[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)?\s*$)");

  // Read remaining lines
  std::string line;
  while (std::getline(file, line))
  {
    std::istringstream lineStream(line);
    std::vector<std::string> rowValues;
    std::string value;
    bool isValidLine = true;

    // Process each column in the line
    for (size_t i = 0; i < columnTitles.size(); ++i)
    {
      if (!(lineStream >> value))
      {
        throw std::runtime_error("[ERROR][readFileToMap] Mismatched data and header columns.");
      }

      // Check if the value contains only allowed characters
      if (!std::regex_match(value, validPattern))
      {
        isValidLine = false;
        break; // Skip the rest of the line
      }

      rowValues.push_back(value);
    }

    // If the line is valid, add its values to the respective columns
    if (isValidLine)
    {
      for (size_t i = 0; i < rowValues.size(); ++i)
      {
        dataMap[columnTitles[i]].push_back(rowValues[i]);
      }
    }
    else
      cout << "[WARNING] Skipping line for character incompatibility" << endl;
  }

  file.close();
}

void utility::readFileToMap(const std::string &filename, std::map<std::string, std::vector<std::string>> &dataMap)
{
  if (!std::filesystem::exists(filename))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File does not exist: " + filename);
  }

  std::ifstream file(filename);
  if (!file)
  {
    throw std::runtime_error("[ERROR][readFileToMap] Could not open file " + filename);
  }

  // Read header line
  std::string headerLine;
  if (!std::getline(file, headerLine))
  {
    throw std::runtime_error("[ERROR][readFileToMap] File " + filename + " is empty or missing header row.");
  }

  // Ensure the string has at least 3 characters to check
  if (filename.size() < 3)
  {
    std::cerr << "[ERROR][utility::readFileToMap] Filename too short to determine file type: " << filename << std::endl;
    return;
  }

  // Get the last three characters of the string
  std::string fileExtension = filename.substr(filename.size() - 3);

  // Convert to lowercase for case-insensitive comparison
  for (char &c : fileExtension)
    c = std::tolower(c);

  // Perform actions based on the file extension
  if (fileExtension == "txt")
  {
    utility::readTxtToMap(filename, dataMap);
  }
  else if (fileExtension == "csv")
  {
    utility::readCsvToMap(filename, dataMap);
  }
  else
  {
    std::cerr << "[WARNING][utility::readFileToMap] Unknown file type: " << filename << ", treating as txt" << std::endl;
  }
}

std::array<std::array<float, 2>, 2> utility::get_intercept(float m1, float em1, float q1, float eq1, float m2, float em2, float q2, float eq2)
{
  std::array<std::array<float, 2>, 2> result;
  result[0][0] = (q2 - q1) / (m1 - m2);
  result[1][0] = (m1 * q2 - m2 * q1) / (m1 - m2);
  auto denominator_err_abs = sqrt(em1 * em1 + em2 * em2);
  auto denominator_err_rel = fabs(denominator_err_abs / (m1 - m2));
  auto numerator1_err_abs = sqrt(eq1 * eq1 + eq2 * eq2);
  auto numerator1_err_rel = fabs(numerator1_err_abs / (q2 - q1));
  auto numerator2_1_err_rel = sqrt((em1 / m1) * (em1 / m1) + (eq2 / q2) * (eq2 / q2));
  auto numerator2_1_err_abs = m1 * q2 * numerator2_1_err_rel;
  auto numerator2_2_err_rel = sqrt((em2 / m2) * (em2 / m2) + (eq1 / q1) * (eq1 / q1));
  auto numerator2_2_err_abs = m2 * q1 * numerator2_1_err_rel;
  auto numerator2_err_abs = sqrt(numerator2_1_err_abs * numerator2_1_err_abs + numerator2_2_err_abs * numerator2_2_err_abs);
  auto numerator2_err_rel = fabs(numerator2_err_abs / (m1 * q2 - m2 * q1));
  result[0][1] = result[0][0] * sqrt(denominator_err_rel * denominator_err_rel + numerator1_err_rel * numerator1_err_rel);
  result[1][1] = result[1][0] * sqrt(denominator_err_rel * denominator_err_rel + numerator2_err_rel * numerator2_err_rel);
  return result;
}

std::array<std::array<float, 2>, 2> utility::get_intercept_pol1_pol2(float x0_0, float ex0_0, float x1_0, float ex1_0, float x0_1, float ex0_1, float x1_1, float ex1_1, float x2_1, float ex2_1)
{
  std::array<std::array<float, 2>, 2> result = {0};
  return result;
  /*
  std::array<std::array<float, 2>, 2> result;
  auto delta = sqrt((x1_1 - x1_0) * (x1_1 - x1_0) - 4 * x2_1 * (x0_1 - x0_0));
  result[0][0] = (-(x1_1 - x1_0) + delta) / (2 * x2_1);
  result[1][0] = -1;

  auto intercept_before = pol1_before->GetParameter(0);
  auto angularc_before = pol1_before->GetParameter(1);
  auto c0_after = pol1_after->GetParameter(0);
  auto c1_after = pol1_after->GetParameter(1);
  auto c2_after = pol1_after->GetParameter(2);
  auto delta = sqrt((c1_after - angularc_before) * (c1_after - angularc_before) - 4 * c2_after * (c0_after - intercept_before));
  auto intercept_before_econtrib = pol1_before->GetParError(0) / (delta);
  auto angularc_before_econtrib = pol1_before->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
  auto c0_after_econtrib = pol1_after->GetParError(0) / (delta);
  auto c1_after_econtrib = -pol1_after->GetParError(1) * ((c1_after - angularc_before) / (delta) + 1) / (2 * c2_after);
  auto c2_after_econtrib = pol1_after->GetParError(2) * ((c0_after - intercept_before) / (c2_after * delta) - (-delta - c1_after + angularc_before) / (2 * c2_after * c2_after));
  auto intercept_val = result[0] = (-(c1_after - angularc_before) + delta) / (2 * c2_after);
  auto intercept_err = result[1] = sqrt(intercept_before_econtrib * intercept_before_econtrib + angularc_before_econtrib * angularc_before_econtrib + c0_after_econtrib * c0_after_econtrib + c1_after_econtrib * c1_after_econtrib + c2_after_econtrib * c2_after_econtrib);

  auto denominator_err_abs = sqrt(em1 * em1 + em2 * em2);
  auto denominator_err_rel = fabs(denominator_err_abs / (m1 - m2));
  auto numerator1_err_abs = sqrt(eq1 * eq1 + eq2 * eq2);
  auto numerator1_err_rel = fabs(numerator1_err_abs / (q2 - q1));
  auto numerator2_1_err_rel = sqrt((em1 / m1) * (em1 / m1) + (eq2 / q2) * (eq2 / q2));
  auto numerator2_1_err_abs = m1 * q2 * numerator2_1_err_rel;
  auto numerator2_2_err_rel = sqrt((em2 / m2) * (em2 / m2) + (eq1 / q1) * (eq1 / q1));
  auto numerator2_2_err_abs = m2 * q1 * numerator2_1_err_rel;
  auto numerator2_err_abs = sqrt(numerator2_1_err_abs * numerator2_1_err_abs + numerator2_2_err_abs * numerator2_2_err_abs);
  auto numerator2_err_rel = fabs(numerator2_err_abs / (m1 * q2 - m2 * q1));
  result[0][1] = result[0][0] * sqrt(denominator_err_rel * denominator_err_rel + numerator1_err_rel * numerator1_err_rel);
  result[1][1] = result[1][0] * sqrt(denominator_err_rel * denominator_err_rel + numerator2_err_rel * numerator2_err_rel);
  return result;
  */
}

std::array<std::array<float, 2>, 2> get_intercept_tf1(TF1 *first_target, TF2 *second_target)
{
  std::array<std::array<float, 2>, 2> result = {0};
  return result;
}