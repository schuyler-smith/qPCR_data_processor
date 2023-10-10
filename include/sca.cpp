/*
 *
 * Author:  Schuyler D. Smith
 * Function:  smart_chip_analyzer
 * Purpose: process outputs from the SmartChip qPCR
 *
 */

#ifndef SCA
#define SCA

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <cctype>
#include <map>
#include <unordered_map>
#include <regex>
#include <numeric>
#include <cmath>

#include "maths_sds.hpp"
#include "maps_sds.hpp"
#include "defs_sds.hpp"
#include "outputs.hpp"

void smart_chip_analyzer(
  std::string   qpcr_data, 
  std::string   output_path, 
  std::string   val_var, 
  std::string   negative_control,  
  std::string   standard, 
  std::string   non_template_control, 
  bool          headers,
  float         efficiency_min,
  float         efficiency_max,
  float         r_sqared_threshold,
  std::string   replacement_stds,
  std::string   gene_magnitudes_file
) {
  vstring       id                  = {"Assay", "Sample"};
  um_str_vstr   assay_group         = map_variable_vec(qpcr_data, {"Assay"}, id, headers);
  um_str_str    group_assay         = map_variable(qpcr_data, id, "Assay", headers);
  um_str_vdbl   group_Ct            = map_variable_vec_numeric(qpcr_data, id, {val_var}, headers);
  um_str_dbl    Ct_means            = um_mean(group_Ct);
  um_str_vdbl   array_Ct            = map_variable_vec_numeric(qpcr_data, {"Assay"}, {val_var}, headers);
  um_str_dbl    Ct_perc_below       = um_percent_below_threshold(array_Ct, 33);
  um_str_dbl    Ct_sd               = um_sd(group_Ct);
  um_str_str    group_sample        = map_variable(qpcr_data, id, "Sample", headers);
  um_str_dbl    group_efficiency    = um_mean(map_variable_vec_numeric(qpcr_data, id, {"Efficiency"}, headers));

  um_str_vstr   replacement_assay_group;
  um_str_vdbl   replacement_group_Ct;
  if (!replacement_stds.empty()) {
    replacement_assay_group = map_variable_vec(replacement_stds, {"Assay"}, id, headers);
    replacement_group_Ct    = map_variable_vec_numeric(replacement_stds, id, {val_var}, headers);
  }

  vstring assays;
  for (auto it : assay_group) {assays.push_back(it.first);}
  std::sort(assays.begin(), assays.end(), [](const std::string& lhs, const std::string& rhs) {
    const auto result = std::mismatch(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), 
      [](const unsigned char lhs, const unsigned char rhs) {return std::tolower(lhs) == std::tolower(rhs);});
    return result.second != rhs.cend() && (result.first == lhs.cend() || std::tolower(*result.first) < std::tolower(*result.second));
  });

// QC - Ct of NTC needs to be 3+ Cq greater than Control1 for each assay
  um_str_dbl  NTC_means;
  um_str_dbl  STD_means;
  um_str_str  QC_NTC;
  um_str_dbl  NEG_means;
  um_str_str  QC_NEG;
  um_str_str  std_QC;
  um_str_dbl  rsqr_map;
  um_str_dbl  std_efficiency_map;
  um_str_pair_dbl_dbl regression_map;
  for (auto assay : assays) {
    NTC_means[assay] = Ct_means[search_vsrting(assay_group[assay], non_template_control)];
    STD_means[assay] = Ct_means[search_vsrting(assay_group[assay], standard + "1")];
    if ((NTC_means[assay]-STD_means[assay]) < 3) {QC_NTC[assay] = "FAIL";} else {QC_NTC[assay] = "PASS";}

    if (negative_control == "none") {
      NEG_means[assay] = std::numeric_limits<double>::quiet_NaN();
      QC_NEG[assay] = "NONE";
    } else {
      NEG_means[assay] = Ct_means[search_vsrting(assay_group[assay], negative_control)];
      if (NEG_means[assay] < 35) {
        QC_NEG[assay] = "FAIL";
      } else {
        QC_NEG[assay] = "PASS";
      }
    }

    vdouble  log_abundances;
    vdouble  Ct_values;
    for (auto group : assay_group[assay]) {
      if (group.find(standard) != std::string::npos) {
        char log_abundance_value = group.back();
        vdouble group_Cts = group_Ct[group];
        Ct_values.insert(std::end(Ct_values), std::begin(group_Cts), std::end(group_Cts));
        for (int i=0; i < group_Cts.size(); i++) {
          log_abundances.push_back(log_abundance_value - 48);
        }
      }
    }

    if (std::none_of(Ct_values.begin(), Ct_values.end(), [](double ct) {return !std::isnan(ct);})) {
      regression_map[assay]     = {0,0};
      rsqr_map[assay]           = 0;
      std_efficiency_map[assay] = 0;
    } else {
      regression_map[assay]     = lm(log_abundances, Ct_values);
      rsqr_map[assay]           = coef_determination(log_abundances, Ct_values);
      std_efficiency_map[assay] = std::pow(10, -1/regression_map[assay].second);
    }
    if (std_efficiency_map[assay] >= efficiency_min 
     && std_efficiency_map[assay] <= efficiency_max 
     && rsqr_map[assay] >= r_sqared_threshold) {
      std_QC[assay] = "PASS";
    } else {
      std_QC[assay] = "FAIL";
    }
    if (std_QC[assay] == "FAIL") {
      if (replacement_assay_group.find(assay) != replacement_assay_group.end()) {
        std::cout << assay << "\n";
        for (auto group : replacement_assay_group[assay]) {
          if (group.find(standard) != std::string::npos) {
            char log_abundance_value = group.back();
            vdouble group_Cts = replacement_group_Ct[group];
            Ct_values.insert(std::end(Ct_values), std::begin(group_Cts), std::end(group_Cts));
            for (int i=0; i < group_Cts.size(); i++) {
              log_abundances.push_back(log_abundance_value-48);
            }
          }
        }
        regression_map[assay]     = lm(log_abundances, Ct_values);
        rsqr_map[assay]           = coef_determination(log_abundances, Ct_values);
        std_efficiency_map[assay] = std::pow(10, -1/regression_map[assay].second);
        std_QC[assay]             = "FAIL_REPLACED";
      }
    }
  }

  um_str_str  group_QC;
  um_str_vdbl group_copyN;
  float       gene_coefficient;
  um_str_flo  gene_magnitudes = file_map(gene_magnitudes_file);
  // create group IDs
  vstring     groupID;
  for (auto const& it : group_assay) {groupID.push_back(it.first);}
  std::sort(groupID.begin(), groupID.end(), [](const std::string& lhs, const std::string& rhs) {
    const auto result = std::mismatch(lhs.cbegin(), lhs.cend(), rhs.cbegin(), rhs.cend(), 
      [](const unsigned char lhs, const unsigned char rhs) {return std::tolower(lhs) == std::tolower(rhs);});
    return result.second != rhs.cend() && (result.first == lhs.cend() || std::tolower(*result.first) < std::tolower(*result.second));
  });
  for (auto group : groupID) {
    std::string assay = group_assay[group];
    if (group_efficiency[group] >= efficiency_min) {
      group_QC[group] = "PASS";
    } else {
      group_QC[group] = "FAIL";
    }

    // adjust ct values for start of curve
    vdouble copy_N;
    for (auto Ct : group_Ct[group]) {
      double N = std::pow(10, (Ct - regression_map[assay].first)/regression_map[assay].second);
      copy_N.push_back(N);
    }
    if (gene_magnitudes.find(assay) != gene_magnitudes.end()) {
      gene_coefficient = gene_magnitudes[assay];
    } else { gene_coefficient = 1; }
    copy_N = magnify(copy_N, gene_coefficient);
    group_copyN[group] = copy_N;
  }


// Create reports for assays, samples, and all
  create_reports(output_path, assays, 
    groupID, group_QC, group_assay, group_sample, group_copyN, Ct_perc_below,
    regression_map, group_efficiency, std_efficiency_map, rsqr_map, 
    std_QC, NEG_means, QC_NEG, NTC_means, STD_means, QC_NTC);
}

#endif