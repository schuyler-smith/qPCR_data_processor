/*
 *
 * Author:      Schuyler D. Smith
 * Function:    smart_chip_analyzer
 * Purpose:     process outputs from the SmartChip qPCR
 *
 */


#include "sca.cpp"
#include "version.hpp"
#include <lyra/lyra.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <dirent.h>
#ifndef _WIN32
  #include <sys/stat.h>
#endif

void file_check(std::string file) {
  if (!file.empty()) {
    std::ifstream check(file);
    if (!check.is_open()) {
      std::cerr << "Error: file " << file << " does not exist" << std::endl;
    }
  }
}

void make_dir(const std::string& directoryPath) {
  std::stringstream ss(directoryPath);
  std::string word;
  std::string mkdir_path;

  while (std::getline(ss, word, '/')) {
    mkdir_path += word + '/';
    #ifdef _WIN32
      _mkdir(mkdir_path.c_str());
    #else
      mkdir(mkdir_path.c_str(), 0777);
    #endif

  }
}


int main(int argc, char *argv[]) {
	std::string input;
	std::string output;
  std::string value = "Ct";
  std::string negative_control = "NEG";
  std::string standard_id = "STD";
  std::string non_template_control = "NTC";
  bool        headers = true;
  float       efficiency_min = 1.70;
  float       efficiency_max = 2.20;
  float       r_sqared_threshold = 0.85;
  std::string replacement_stds;
  std::string gene_magnitudes;
  // help flag
  bool show_help    = false;
  bool show_version = false;
  // declare options
  auto cli
    = lyra::help(show_help)
    | lyra::opt([&](bool){ show_version += 1; }) 
      ["-v"]["--version"]
      ("Version information.")
    | lyra::opt( input, "").required()
      ["-i"]["--input"]
      ("Input CSV file. Output from the smartchip qPCR. (required)")
    | lyra::opt( output, "").optional()
      ["-o"]["--output"]
      ("Output directory path (default uses input path), and/or prefix for output files (default uses input filename).")
    | lyra::opt( value, "Ct")
      ["-c"]["--Ct"]
      ("Column name in input file to use for the cycle thresholds.")
    | lyra::opt( headers, "true" )
      ["-j"]["--headers"]
      ("Whether the input file has a header row.")
    | lyra::opt( negative_control, "NEG" )
      ["-n"]["--negcontrol"]
      ("Sample identifiers for the negative controls.")
    | lyra::opt( standard_id, "STD")
      ["-s"]["--standard"]
      ("Sample identifiers for the standards.")
    | lyra::opt( non_template_control, "NTC" )
      ["-t"]["--nontemplate"]
      ("Sample identifiers for the non-template controls.")
    | lyra::opt( efficiency_min, "1.70" )
      ["-e"]["--effmin"]
      ("Minimum value for qPCR efficiency for PASS.")
    | lyra::opt( efficiency_max, "2.20" )
      ["-E"]["--effmax"]
      ("Maximum value for qPCR efficiency for PASS.")
    | lyra::opt( r_sqared_threshold, "0.85" )
      ["-R"]["--rsquare"]
      ("Minimum value for r-squared of Ct value and log-abundances to receive PASS.")
    | lyra::opt( replacement_stds, "" ).optional()
      ["-r"]["--replacements"]
      ("Replacement values for when standards FAIL quality check.")
    | lyra::opt( gene_magnitudes, "" ).optional()
      ["-m"]["--magnitudes"]
      ("File with maps (gene:value) for gene abundance coefficients.")
  ;

  // Check that the arguments where valid:
  auto result = cli.parse({ argc, argv });
  if (!result) {
    std::cerr << "Error in command line: " << result.message() << std::endl;
    return 1;
  }
  if (show_help) {
    std::cout
      << PROGRAM_VERSION
      << cli << "\n";
    return 0;
  }
  if (show_version) {
    std::cout << PROGRAM_VERSION;
    return 0;
  }

  // create input file array:
  std::vector<std::string> inputs;
  DIR* directory = opendir(input.c_str());
  if (directory != nullptr) {
    dirent* entry;
    while ((entry = readdir(directory)) != nullptr) {
      if (entry->d_type == DT_REG) {
        inputs.push_back(input + '/' + entry->d_name);
      }
    }
    closedir(directory);
  } else {
    inputs.push_back(input);
  }

  // Check all input files exist
  for (std::string file : inputs) {
    file_check(file);
  }
  std::string files[] = { replacement_stds, gene_magnitudes };
  for (std::string file : files) {
    file_check(file);
  }

  if (output.empty()) {
    output = input.substr(0, input.find_last_of("/\\"));
  }
  if (!(output.back() == '/' || output.back() == '\\')) {
    output += "/";
  }
  output += "sma_output/";
  make_dir(output);

  for (std::string input_file : inputs) {
    std::string output_file;
    std::string filename;
    filename = input_file.substr(input_file.find_last_of("/\\") + 1);
    filename = filename.substr(0, filename.find_last_of("."));
    output_file = output + filename;

    std::cout << input_file << "\n";
    std::cout << output_file << "\n";
    
    smart_chip_analyzer(
      input_file, 
      output_file, 
      value, 
      negative_control, 
      standard_id, 
      non_template_control, 
      headers, 
      efficiency_min, 
      efficiency_max, 
      r_sqared_threshold,
      replacement_stds,
      gene_magnitudes
    );
  }

	return(0);
}