// class for processes in include/sca.cpp

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#ifdef _WIN32
  #include <dirent.h>
#else
  #include <sys/stat.h>
#endif

#include "maths_sds.hpp"
#include "maps_sds.hpp"
#include "defs_sds.hpp"
#include "outputs.hpp"


class SmartchipIngest {

  public:
    bool        headers;
    std::string output_path;
    std::string val_var;
    void file_check(const std::string file);
  
  protected:

  private:
    vstring     id  = {"Assay", "Sample"};
    void make_dir(const std::string& directoryPath);
    std::vector<std::string> create_input_array(std::string input);
};

void SmartchipIngest::file_check(const std::string file) {
  if (!file.empty()) {
    std::ifstream check(file);
    if (!check.is_open()) {
      std::cerr << "Error: file '" << file << "' does not exist." << std::endl;
    }
  }
}

void SmartchipIngest::make_dir(const std::string& directoryPath) {
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

std::vector<std::string> SmartchipIngest::create_input_array(std::string input) {
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
  return inputs;
}



// class SmartchipQC: public SmartchipIngest {

//   public:
//     bool        headers;
//     std::string qpcr_output;
//     std::string val_var;
  
//   protected:
//     um_str_vdbl group_Ct  = map_variable_vec_numeric(qpcr_output, id, {val_var}, headers);
//     um_str_dbl  Ct_means  = um_mean(group_Ct);
//     um_str_dbl  NTC_means;
//     um_str_dbl  STD_means;
//     um_str_str  QC_NTC;
//     um_str_dbl  NEG_means;
//     um_str_str  QC_NEG;
//     um_str_str  std_QC;
//     um_str_dbl  rsqr_map;
//     um_str_dbl  std_efficiency_map;
//     um_str_pair_dbl_dbl regression_map;

//   private:
//     vstring     id  = {"Assay", "Sample"};
// };



// // Create reports for assays, samples, and all
// class SmartchipReport {

//   public:
//     void create_reports()

// };

// SmartchipReport.create_reports(output, assays, groupID, group_QC, group_assay, group_sample, group_copyN, Ct_perc_below,
//   regression_map, group_efficiency, std_efficiency_map, rsqr_map, std_QC, NEG_means, QC_NEG, NTC_means, STD_means, QC_NTC);