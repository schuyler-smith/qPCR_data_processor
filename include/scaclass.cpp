/*
 *
 * Author:      Schuyler D. Smith
 * Function:    test
 * Purpose:     testing
 *
 */

#ifndef SCACLASS
#define SCACLASS


// class for processes in include/sca.cpp
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <dirent.h>
#ifndef _WIN32
  #include <sys/stat.h>
#endif
#include <stdexcept>

#include "defs_sds.hpp"
#include "maths_sds.hpp"
#include "maps_sds.hpp"
#include "defs_sds.hpp"
#include "outputs.hpp"

namespace SmartchipInfra {
  void file_check(const std::string& file) {
    if (!file.empty()) {
      std::ifstream check(file);
      if (!check.is_open()) {
        throw std::invalid_argument("Error: file '" + file + "' not found.");
      }
    }
  }

  std::string basename(const std::string& filepath) {
    std::string basename;
    basename = filepath.substr(filepath.find_last_of("/\\") + 1);
    basename = basename.substr(0, basename.find_last_of("."));
    return(basename);
  }

  std::vector<std::string> create_file_array(const std::string& path) {
    std::vector<std::string> files;
    DIR* directory = opendir(path.c_str());
    if (directory != nullptr) {
      dirent* entry;
      while ((entry = readdir(directory)) != nullptr) {
        if (entry->d_type == DT_REG) {
          files.push_back(path + '/' + entry->d_name);
        }
      }
      closedir(directory);
    } else {
      files.push_back(path);
    }
    return(files);
  }

  std::vector<std::string> read_csv_headers(const std::string& file_name) {
    std::string csv_file_name = file_name;
    std::ifstream csv_file(csv_file_name);
    if (!csv_file.is_open()) {
      std::cerr << "Failed to open the CSV file." << std::endl;
    }
    std::string header_line;
    std::vector<std::string> columns;
    std::string column;
    if (std::getline(csv_file, header_line)) {
      std::istringstream header_stream(header_line);
      while (std::getline(header_stream, column, ',')) {
        columns.push_back(column);
      }
    } else {
      std::cerr << "Failed to read the header line from the CSV file." << std::endl;
    }
    csv_file.close();
    return(columns);
  }

  bool contained_in_vector(const std::vector<std::string>& vector, const std::string& string) {
    if (std::find(vector.begin(), vector.end(), string) == vector.end()) {
      return false;
    }
    return true;
  }

  void check_headers(const std::string& csv, const std::vector<std::string>& required_headers) {
    std::vector<std::string> file_headers;
    file_headers = SmartchipInfra::read_csv_headers(csv);
    for (const std::string& header : required_headers) {
      if (!SmartchipInfra::contained_in_vector(file_headers, header)) {
        throw std::invalid_argument("'" + csv + "' missing required field '" + header + "'.");
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
}

//
// SmartchipIngest
//

class SmartchipIngest {
  public:
    std::string data;
    std::string input_path;
    std::string output_dir;
    std::string output_file;
    std::string replacement_stds_path;
    std::string gene_magnitudes_path;
    
    SmartchipIngest(
      const std::string&  qPCR_data_path
    ) : input_path(qPCR_data_path) {
      construct_ingest();
    }

    void construct_ingest();
    void set_output_dir(const std::string&);
    void set_input(const std::string&);
    void set_replacement_stds(const std::string&);
    void set_gene_magnitudes(const std::string&);
  
  protected:
    um_str_flo  gene_magnitudes;
};

void SmartchipIngest::construct_ingest() {
  set_input(input_path);
  output_dir = input_path.substr(0, input_path.find_last_of("/\\"));
  output_dir += "/sca_output/";
  output_file = output_dir;
  output_file += SmartchipInfra::basename(input_path);
}

void SmartchipIngest::set_output_dir(const std::string& path) {output_dir = path;}
void SmartchipIngest::set_input(const std::string& qPCR_data_path) {
  input_path = qPCR_data_path;
  try {
    SmartchipInfra::file_check(input_path);
  }
  catch (const std::exception& e) {
    std::cerr << "Error in input_path: " << e.what() << std::endl;
  }
  data = input_path;
}

void SmartchipIngest::set_replacement_stds(const std::string& path) {
  try {
    SmartchipInfra::file_check(path);
    replacement_stds_path = path;
  }
  catch (const std::exception& e) {
    std::cerr << "Error in set_replacement_stds(): " << e.what() << std::endl;
  }
}

void SmartchipIngest::set_gene_magnitudes(const std::string& path) {
  try {
    SmartchipInfra::file_check(path);
    gene_magnitudes = file_map(path);
  }
  catch (const std::exception& e) {
    std::cerr << "Error in set_gene_magnitudes(): " << e.what() << std::endl;
  }
}

//
// Class SmartchipParameters
//

class SmartchipParameters : public SmartchipIngest {
  public:
    std::string assay_colname;
    std::string sample_colname;
    std::string ct_colname;
    std::string negative_control_id;
    std::string standard_id;
    std::string non_template_id;
    std::string efficiency_colname;
    float       efficiency_min;
    float       efficiency_max;
    float       r_sqared_threshold;

    SmartchipParameters(
      const std::string& qPCR_data_path
    ) : SmartchipIngest(qPCR_data_path) {
      construct_params();
    }
    SmartchipParameters(
      const SmartchipIngest ingest
    ) : SmartchipParameters(ingest.input_path) {}

    // void check_Smartchip_headers();

    void construct_params();
    void set_assay_colname(const std::string&);
    void set_sample_colname(const std::string&);
    void set_qPCR_ct_colname(const std::string&);
    void set_negative_control(const std::string&);
    void set_standard_id(const std::string&);
    void set_non_template_control(const std::string&);
    void set_efficiency_colname(const std::string&);
    void set_efficiency_min(const float&);
    void set_efficiency_max(const float&);
    void set_r_sqared_threshold(const float&);
};

void SmartchipParameters::construct_params() {
  assay_colname         = "Assay";
  sample_colname        = "Sample";
  ct_colname            = "Ct";
  negative_control_id   = "NEG";
  standard_id           = "STD";
  non_template_id       = "NTC";
  non_template_id       = "Efficiency";
  efficiency_min        = 1.70;
  efficiency_max        = 2.20;
  r_sqared_threshold    = 0.85;
}

void SmartchipParameters::set_assay_colname(const std::string& x)         {assay_colname = x; SmartchipInfra::check_headers(data, {x});}
void SmartchipParameters::set_sample_colname(const std::string& x)        {sample_colname = x; SmartchipInfra::check_headers(data, {x});}
void SmartchipParameters::set_qPCR_ct_colname(const std::string& x)       {ct_colname = x; SmartchipInfra::check_headers(data, {x});}
void SmartchipParameters::set_efficiency_colname(const std::string& x)    {efficiency_colname = x; SmartchipInfra::check_headers(data, {x});}
void SmartchipParameters::set_negative_control(const std::string& x)      {negative_control_id = x;}
void SmartchipParameters::set_standard_id(const std::string& x)           {standard_id = x;}
void SmartchipParameters::set_non_template_control(const std::string& x)  {non_template_id = x;}
void SmartchipParameters::set_efficiency_min(const float& x)              {efficiency_min = x;}
void SmartchipParameters::set_efficiency_max(const float& x)              {efficiency_max = x;}
void SmartchipParameters::set_r_sqared_threshold(const float& x)          {r_sqared_threshold = x;}

// void SmartchipParameters::check_Smartchip_headers() {
//   try {
//     SmartchipInfra::check_headers(data, {sample_colname, assay_colname, efficiency_colname, ct_colname});
//   }
//   catch (const std::exception& e) {
//     std::cerr << "Error in input_files: " << e.what() << std::endl;
//   }
// }

//
// Class SmartchipExtract
//

class SmartchipExtract : public SmartchipParameters {
  public:
    SmartchipExtract(
      const std::string& qPCR_data_path
    ) : SmartchipParameters(qPCR_data_path) {
      construct_extract();
    }
    SmartchipExtract(
      const SmartchipIngest ingest
    ) : SmartchipExtract(ingest.input_path) {}
    SmartchipExtract(
      const SmartchipParameters parameters
    ) : SmartchipExtract(parameters.input_path) {}
    
    vstring       assays;
  
  // protected:
    um_str_vstr   assay_group;
    um_str_str    group_assay;
    um_str_vdbl   group_Ct;
    um_str_dbl    Ct_means;
    um_str_vdbl   array_Ct;
    um_str_dbl    Ct_perc_below;
    um_str_dbl    Ct_sd;
    um_str_str    group_sample;
    um_str_dbl    group_efficiency;
    um_str_vstr   replacement_assay_group;
    um_str_vdbl   replacement_group_Ct;
    vstring       groups;

  private:
    void construct_extract();
    std::vector<std::string> id = {assay_colname, sample_colname};
    void extract_assay();
    void extract_groups();
};

void SmartchipExtract::construct_extract() {
  // check_Smartchip_headers();
  assay_group       = map_variable_vec(data, {assay_colname}, id);
  group_assay       = map_variable(data, id, assay_colname);
  group_Ct          = map_variable_vec_numeric(data, id, {ct_colname});
  Ct_means          = um_mean(group_Ct);
  array_Ct          = map_variable_vec_numeric(data, {assay_colname}, {ct_colname});
  Ct_perc_below     = um_percent_below_threshold(array_Ct, 33);
  Ct_sd             = um_sd(group_Ct);
  group_sample      = map_variable(data, id, sample_colname);
  group_efficiency  = um_mean(map_variable_vec_numeric(data, id, {efficiency_colname}));
  if (!replacement_stds_path.empty()) {
    replacement_assay_group = map_variable_vec(replacement_stds_path, {assay_colname}, id);
    replacement_group_Ct    = map_variable_vec_numeric(replacement_stds_path, id, {ct_colname});
  }
  extract_assay();
  extract_groups();
}

void SmartchipExtract::extract_assay() {
  for (auto it : assay_group) {assays.push_back(it.first);}
  std::sort(assays.begin(), assays.end(), [](const std::string& lhs, const std::string& rhs) {
    const auto result = std::mismatch(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), 
      [](const unsigned char lhs, const unsigned char rhs) {return std::tolower(lhs) == std::tolower(rhs);});
    return result.second != rhs.cend() && (result.first == lhs.cend() || std::tolower(*result.first) < std::tolower(*result.second));
  });
}

void SmartchipExtract::extract_groups() {
  for (auto const& it : group_assay) {groups.push_back(it.first);}
  std::sort(groups.begin(), groups.end(), [](const std::string& lhs, const std::string& rhs) {
    const auto result = std::mismatch(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), 
      [](const unsigned char lhs, const unsigned char rhs) {return std::tolower(lhs) == std::tolower(rhs);});
    return result.second != rhs.cend() && (result.first == lhs.cend() || std::tolower(*result.first) < std::tolower(*result.second));
  });
}

//
// Class SmartchipTransform
//

class SmartchipTransform : public SmartchipExtract {
  public:
    SmartchipTransform(
      const std::string& qPCR_data_path
    ) : SmartchipExtract(qPCR_data_path) {
      construct_transform();
    }
    SmartchipTransform(
      const SmartchipIngest ingest
    ) : SmartchipTransform(ingest.input_path) {}
    SmartchipTransform(
      const SmartchipParameters parameters
    ) : SmartchipTransform(parameters.input_path) {}
    SmartchipTransform(
      const SmartchipExtract extracted
    ) : SmartchipTransform(extracted.input_path) {}

  // protected:
    um_str_dbl  NTC_means;
    um_str_dbl  STD_means;
    um_str_str  QC_NTC;
    um_str_dbl  NEG_means;
    um_str_str  QC_NEG;
    um_str_str  std_QC;
    um_str_dbl  rsqr_map;
    um_str_dbl  std_efficiency_map;
    um_str_pair_dbl_dbl regression_map;  
    um_str_str  group_QC;
    um_str_vdbl group_copyN;

  private:
    void construct_transform();
    void quality_check_NTC(const std::string&);
    void quality_check_NEG(const std::string&);
    void quality_check_STD(const std::string&);
    void quality_check_EFF(const std::string&);
    void extract_log_value_and_Ct(um_str_vstr, const std::string&, vdouble&, vdouble&);
    void regression_analysis(const std::string&, const vdouble&, const vdouble&);
    void calculate_copyN(const std::string&);
};
// 
void SmartchipTransform::construct_transform() {
  for (auto assay : assays) {
    vdouble  log_abundances;
    vdouble  Ct_values;
    quality_check_NTC(assay);
    quality_check_NEG(assay);
    extract_log_value_and_Ct(assay_group, assay, log_abundances, Ct_values);
    regression_analysis(assay, log_abundances, Ct_values);
    quality_check_STD(assay);
    if (std_QC[assay] == "FAIL") {
      if (replacement_assay_group.find(assay) != replacement_assay_group.end()) {
        log_abundances.clear();
        Ct_values.clear();
        extract_log_value_and_Ct(replacement_assay_group, assay, log_abundances, Ct_values);
        regression_analysis(assay, log_abundances, Ct_values);
      }
    }
  }
  for (auto group : groups) {
    quality_check_EFF(group);
    calculate_copyN(group);
  }
}

void SmartchipTransform::quality_check_NTC(const std::string& assay) {
  NTC_means[assay] = Ct_means[search_vsrting(assay_group[assay], negative_control_id)];
  STD_means[assay] = Ct_means[search_vsrting(assay_group[assay], standard_id + "1")];
  if ((NTC_means[assay] - STD_means[assay]) < 3) {
    QC_NTC[assay] = "FAIL";
  } else {
    QC_NTC[assay] = "PASS";
  }
}

void SmartchipTransform::quality_check_NEG(const std::string& assay) {
  if (negative_control_id == "none") {
    NEG_means[assay] = std::numeric_limits<double>::quiet_NaN();
    QC_NEG[assay] = "NONE";
  } else {
    NEG_means[assay] = Ct_means[search_vsrting(assay_group[assay], negative_control_id)];
    if (NEG_means[assay] < 35) {
      QC_NEG[assay] = "FAIL";
    } else {
      QC_NEG[assay] = "PASS";
    }
  }
}

void SmartchipTransform::quality_check_STD(const std::string &assay) {
  if (std_efficiency_map[assay] >= efficiency_min 
    && std_efficiency_map[assay] <= efficiency_max 
    && rsqr_map[assay] >= r_sqared_threshold) {
    std_QC[assay] = "PASS";
  } else {
    std_QC[assay] = "FAIL";
  }
}

void SmartchipTransform::extract_log_value_and_Ct(
  um_str_vstr assay_group,
  const std::string& assay, 
  vdouble& log_abundances, 
  vdouble& Ct_values
) {
  for (auto group : assay_group[assay]) {
    if (group.find(standard_id) != std::string::npos) {
      char log_abundance_value = group.back();
      vdouble group_Cts = group_Ct[group];
      Ct_values.insert(std::end(Ct_values), std::begin(group_Cts), std::end(group_Cts));
      for (int i=0; i < group_Cts.size(); i++) {
        log_abundances.push_back(log_abundance_value - 48);
      }
    }
  }
}

void SmartchipTransform::regression_analysis(const std::string& assay, const vdouble& log_abundances, const vdouble& Ct_values) {
  if (std::none_of(Ct_values.begin(), Ct_values.end(), [](double ct) {return !std::isnan(ct);})) {
    regression_map[assay]     = {0,0};
    rsqr_map[assay]           = 0;
    std_efficiency_map[assay] = 0;
  } else {
    regression_map[assay]     = lm(log_abundances, Ct_values);
    rsqr_map[assay]           = coef_determination(log_abundances, Ct_values);
    std_efficiency_map[assay] = std::pow(10, -1/regression_map[assay].second);
  }
}

void SmartchipTransform::quality_check_EFF(const std::string& group) {
  if (group_efficiency[group] >= efficiency_min) {
    group_QC[group] = "PASS";
  } else {
    group_QC[group] = "FAIL";
  }
}

void SmartchipTransform::calculate_copyN(const std::string& group) {
  float gene_coefficient;
  std::string a = group_assay[group];
  vdouble copy_N;
  for (auto Ct : group_Ct[group]) {
    double N = std::pow(10, (Ct - regression_map[a].first)/regression_map[a].second);
    copy_N.push_back(N);
  }
  if (gene_magnitudes.find(a) != gene_magnitudes.end()) {
    gene_coefficient = gene_magnitudes[a];
  } else { gene_coefficient = 1; }
  copy_N = magnify(copy_N, gene_coefficient);
  group_copyN[group] = copy_N;
}

//
// Class SmartchipAnalyzer
//

class SmartchipAnalyzer : public SmartchipTransform {
  public:
    SmartchipAnalyzer(
      const std::string& qPCR_data_path
    ) : SmartchipTransform(qPCR_data_path) {
      construct_load();
    }
    SmartchipAnalyzer(
      const SmartchipIngest ingest
    ) : SmartchipTransform(ingest.input_path) {}
    SmartchipAnalyzer(
      const SmartchipParameters parameters
    ) : SmartchipTransform(parameters.input_path) {}
    SmartchipAnalyzer(
      const SmartchipExtract extracted
    ) : SmartchipAnalyzer(extracted.input_path) {}
    SmartchipAnalyzer(
      const SmartchipTransform transformed
    ) : SmartchipAnalyzer(transformed.input_path) {}

    void build_reports();

  private:
    void construct_load();
};

void SmartchipAnalyzer::construct_load() {
}

void SmartchipAnalyzer::build_reports() {
  SmartchipInfra::make_dir(output_dir);
  create_reports(output_file, assays, 
    groups, group_QC, group_assay, group_sample, group_copyN, Ct_perc_below,
    regression_map, group_efficiency, std_efficiency_map, rsqr_map, 
    std_QC, NEG_means, QC_NEG, NTC_means, STD_means, QC_NTC);
}

#endif
