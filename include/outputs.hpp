/*
 *
 * Author:  Schuyler D. Smith
 * Function:  smart_chip_analyzer
 * Purpose: process outputs from the SmartChip qPCR
 *
 */

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

#ifndef OUTPUTS
#define OUTPUTS


void create_assay_report(
  std::string output,
  vstring     assays,
  um_str_dbl  std_efficiency_map,
  um_str_pair_dbl_dbl regression_map,
  um_str_dbl  rsqr_map,
  um_str_str  std_QC,
  um_str_dbl  NEG_means,
  um_str_str  QC_NEG,
  um_str_dbl  NTC_means,
  um_str_dbl  STD_means,
  um_str_str  QC_NTC,
  um_str_dbl  Ct_perc_below
) {
  std::ofstream assay_report_file(output + "_assay_QC_report.csv");
  assay_report_file 
    << "Assay" << ","
    << "STD_Efficiency" << ","
    << "Slope" << ","
    << "Intercept" << ","
    << "Rsqr" << ","
    << "QC_StdCurve" << ","
    << "NEG_Ct" << ","
    << "QC_NEG" << ","
    << "NTC_diff" << ","
    << "QC_NTC" << ","
    << "Percent_Positive_Samples" << "\n";
  for (auto assay : assays) {
    assay_report_file 
      << assay << ","
      << std_efficiency_map[assay] << ","
      << regression_map[assay].second << ","
      << regression_map[assay].first << ","
      << rsqr_map[assay] << ","
      << std_QC[assay] << ","
      << NEG_means[assay] << ","
      << QC_NEG[assay] << ","
      << NTC_means[assay] - STD_means[assay] << ","
      << QC_NTC[assay] << ","
      << std::round(Ct_perc_below[assay]*100) << "\n";
  }
}

void create_sample_report(
  std::string output,
  vstring     groupID,
  um_str_str  group_QC,
  um_str_str  group_assay,
  um_str_str  group_sample,
  um_str_vdbl group_copyN,
  um_str_dbl  group_efficiency
) {   
  std::ofstream sample_report_file(output + "_sample_QC_report.csv");
  sample_report_file 
    << "Assay," 
    << "Sample," 
    << "Mean_Copy_N," 
    << "Sd_Copy_N," 
    << "Mean_Efficiency," 
    << "QCSample," 
    << "\n";
  for (auto group : groupID) {
    sample_report_file 
      << group_assay[group] << ","
      << group_sample[group] << ","
      << mean(group_copyN[group]) << ","
      << sd(group_copyN[group]) << ","
      << group_efficiency[group] << ","
      << group_QC[group] 
      << "\n";
  }
}

void create_LIMS_report(
  std::string output,
  vstring     groupID,
  um_str_str  group_QC,
  um_str_str  group_assay,
  um_str_str  group_sample,
  um_str_vdbl group_copyN,
  um_str_dbl  group_efficiency
) {   
  std::ofstream LIMS_report_file(output + "_LIMS_report.csv");
  int row = 0;
  LIMS_report_file 
    << ","
    << "Number,"
    << "Assay,"
    << "Cycle,"
    << "FunctionalGroup,"
    << "GeneClass,"
    << "Measure," 
    << "JIC,"
    << "Sample," 
    << "meanCopyN," 
    << "stderr_CopyN," 
    << "Mean_Efficiency," 
    << "QCSample," 
    << "\n";
  for (auto group : groupID) {
    ++row;
    LIMS_report_file
      << row
      << ","
      << ","
      << ","
      << ","
      << ","
      << group_assay[group] << ","
      << ","
      << group_sample[group] << ","
      << mean(group_copyN[group]) << ","
      << sd(group_copyN[group]) << ","
      << group_efficiency[group] << ","
      << group_QC[group] 
      << "\n";
  }
}

void create_full_report(
  std::string output,
  vstring     groupID,
  um_str_str  group_QC,
  um_str_str  group_assay,
  um_str_str  group_sample,
  um_str_vdbl group_copyN,
  um_str_dbl  group_efficiency,
  um_str_dbl  std_efficiency_map,
  um_str_dbl  rsqr_map,
  um_str_str  std_QC,
  um_str_dbl  NEG_means,
  um_str_str  QC_NEG,
  um_str_dbl  NTC_means,
  um_str_dbl  STD_means,
  um_str_str  QC_NTC
) {   
  std::ofstream all_report_file(output + "_sample_qpcr_output_with_assay_info_qc.csv");
  all_report_file 
    << "Assay," 
    << "Sample," 
    << "Mean_Copy_N," 
    << "stderr," 
    << "meanEffi," 
    << "QCSample," 
    << "STD_Efficiency," 
    << "Rsqr," 
    << "QC_StdCurve," 
    << "NEG_Ct," 
    << "QC_NEG," 
    << "NTC_diff," 
    << "QC_NTC," 
    << "\n";
  for (auto group : groupID) {
    std::string assay = group_assay[group];
    all_report_file 
      << assay << ","
      << group_sample[group] << ","
      << mean(group_copyN[group]) << ","
      << sd(group_copyN[group]) << ","
      << group_efficiency[group] << ","
      << group_QC[group] << ","
      << std_efficiency_map[assay] << ","
      << rsqr_map[assay] << ","
      << std_QC[assay] << ","
      << NEG_means[assay] << ","
      << QC_NEG[assay] << ","
      << NTC_means[assay] - STD_means[assay] << ","
      << QC_NTC[assay] 
      << "\n";
  }
}

void create_reports(
  std::string output,
  vstring     assays,
  vstring     groupID,
  um_str_str  group_QC,
  um_str_str  group_assay,
  um_str_str  group_sample,
  um_str_vdbl group_copyN,
  um_str_dbl  Ct_perc_below,
  um_str_pair_dbl_dbl regression_map,
  um_str_dbl  group_efficiency,
  um_str_dbl  std_efficiency_map,
  um_str_dbl  rsqr_map,
  um_str_str  std_QC,
  um_str_dbl  NEG_means,
  um_str_str  QC_NEG,
  um_str_dbl  NTC_means,
  um_str_dbl  STD_means,
  um_str_str  QC_NTC
) {
  create_assay_report(output, assays, std_efficiency_map, regression_map, rsqr_map, 
    std_QC, NEG_means, QC_NEG, NTC_means, STD_means, QC_NTC, Ct_perc_below);
  // create_sample_report(output, groupID, group_QC, group_assay, group_sample,
  //     group_copyN, group_efficiency);
  create_LIMS_report(output, groupID, group_QC, group_assay, group_sample,
    group_copyN, group_efficiency);
  create_full_report(output, groupID, group_QC, group_assay, group_sample, 
    group_copyN, group_efficiency, std_efficiency_map, rsqr_map, 
    std_QC, NEG_means, QC_NEG, NTC_means, STD_means, QC_NTC);
}

#endif