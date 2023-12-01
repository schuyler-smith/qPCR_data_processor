/*
 *
 * Author:      Schuyler D. Smith
 * Function:    smart_chip_analyzer
 * Purpose:     process outputs from the SmartChip qPCR
 *
 */


#include "scaclass.cpp"
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

int main(int argc, char *argv[]) {
	std::string input;
	std::string output;
  std::string value = "Ct";
  std::string assay = "Assay";
  std::string sample = "Sample";
  std::string efficiency = "Efficiency";
  std::string negative_control = "NEG";
  std::string standard_id = "STD";
  std::string non_template_control = "NTC";
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
    | lyra::opt( sample, "Sample")
      ["-s"]["--sample"]
      ("Column name in input file to use for the Sample IDs.")
    | lyra::opt( assay, "Assay")
      ["-a"]["--assay"]
      ("Column name in input file to use for the Assay names.")
    | lyra::opt( efficiency, "Efficiency")
      ["-f"]["--efficiency"]
      ("Column name in input file to use for the Efficiency scores.")
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
  inputs = SmartchipInfra::create_file_array(input);

  for (std::string input_file : inputs) {
    SmartchipParameters sma(input_file);
    sma.set_replacement_stds(replacement_stds);
    sma.set_gene_magnitudes(gene_magnitudes);
    if (!output.empty()) {
      sma.set_output_dir(output);
    }
    sma.set_assay_colname(assay);
    sma.set_sample_colname(sample);
    sma.set_qPCR_ct_colname(value);
    sma.set_efficiency_colname(efficiency);
    sma.set_negative_control(negative_control);
    sma.set_standard_id(standard_id);
    sma.set_non_template_control(non_template_control);
    sma.set_efficiency_min(efficiency_min);
    sma.set_efficiency_max(efficiency_max);
    sma.set_r_sqared_threshold(r_sqared_threshold);
    SmartchipAnalyzer sma_report(sma);
    sma_report.build_reports();
  }

	return(0);
}