/*
 *
 * Author:      Schuyler D. Smith
 * Function:    smart_chip_analyzer
 * Purpose:     process outputs from the SmartChip qPCR
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

#include "sca.hpp"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "gopt.h"




typedef std::unordered_map<std::string, std::vector<std::string> > um_str_vstr;

auto map_headers(std::string input_csv_file) {
    std::ifstream csv_file(input_csv_file);
    std::string line;
    std::string csv_line;
    std::vector<std::string> row_data;
    std::unordered_map<std::string, int> headers;
    while(std::getline(csv_file, line)) {
        std::stringstream line_(line);
        std::getline(line_, csv_line, '\t');
        row_data = string_split(csv_line,',');
        for (int i = 0; i < row_data.size(); ++i) {
            headers[row_data[i]] = i;
        }
        break;
    }
    return headers;
}

auto map_variable_vec(std::string input_csv_file, std::vector<std::string> variables, std::vector<std::string> values, bool headers=true) {
    std::ifstream csv_file(input_csv_file);
    std::string line;
    std::vector<std::string> row_data;
    um_str_vstr variable_map;
    std::unordered_map<std::string, int> colnames;
    colnames = map_headers(input_csv_file);
    while(std::getline(csv_file, line)) {
        if (headers) {headers=false; continue;}
        row_data = string_split(line,',');
        std::string key, value;
        for (int i = 0; i < variables.size(); ++i) {
            key = key.append(row_data[colnames[variables[i]]]);
        }
        for (int i = 0; i < values.size(); ++i) {
            value = value.append(row_data[colnames[values[i]]]);
        }
        if (value.empty()) {value = "NAN";}
        if (variable_map.count(key)) {
            if (std::find(variable_map[key].begin(), variable_map[key].end(), value) != variable_map[key].end()) {continue;
            } else {variable_map[key].push_back(value);}
        } else{variable_map.insert(std::pair<std::string, std::vector<std::string> >(key, {value}));}
    }
    return variable_map;
}

auto map_variable(std::string input_csv_file, std::vector<std::string> variables, std::string value, bool headers=true) {
    std::ifstream csv_file(input_csv_file);
    std::string line;
    std::vector<std::string> row_data;
    std::unordered_map<std::string, std::string> variable_map;
    std::unordered_map<std::string, int> colnames;
    colnames = map_headers(input_csv_file);
    while(std::getline(csv_file, line)) {
        if (headers) {headers=false; continue;}
        row_data = string_split(line,',');
        std::string key;
        for (int i = 0; i < variables.size(); ++i) {
            key.append(row_data[colnames[variables[i]]]);
        }
        variable_map[key] = row_data[colnames[value]];
    }
    return variable_map;
}

auto map_csv(std::string input_csv_file, std::string variable, bool headers=true) {
    std::ifstream csv_file(input_csv_file);
    std::string line;
    std::vector<std::string> row_data;
    std::unordered_map<std::string, std::string> csv_map;
    std::unordered_map<std::string, int> colnames;
    colnames = map_headers(input_csv_file);
    while(std::getline(csv_file, line)) {
        if (headers) {headers=false; continue;}
        row_data = string_split(line,',');
        std::string key, value;
        for (int i = 0; i < colnames.size(); ++i) {
            if (i == colnames[variable]) {
                key = row_data[i];
            } else {
                if (value.empty()) {
                    value = row_data[i];
                } else {
                    value.append("," + row_data[i]);
                }
            }
        }
        csv_map[key] = value;
    }
    return csv_map;
}

auto um_mean(um_str_vstr val_map) {
    std::unordered_map<std::string, float> means;
    for (auto it : val_map) {
        std::vector<float> values(it.second.size());
        std::transform(it.second.begin(), it.second.end(), values.begin(), [](const std::string& value) {return std::stof(value);});
        means[it.first] = mean(values);
    }
    return means;
}

auto um_sd(um_str_vstr val_map) {
    std::unordered_map<std::string, float> means;
    std::unordered_map<std::string, float> sd;
    means = um_mean(val_map);
    for (auto it : val_map) {
        std::vector<float> values(it.second.size());
        std::transform(it.second.begin(), it.second.end(), values.begin(), [](const std::string& value) {return std::stof(value);});
        float variance = 0;
        for (auto v : values) {
            variance += std::pow((v - means[it.first]), 2);
        }
        sd[it.first] = std::sqrt(variance);
    }
    return means;
}


auto smart_chip_analyzer(std::string input_csv_file, std::string assays_file, std::string output, std::string val_var, std::vector<std::string> cols = {"Assay", "Sample"}, std::string neg_control = "H12", std::string standard = "STD", std::string non_template = "NTC", bool headers=true) {
    um_str_vstr assay_cols = map_variable_vec(input_csv_file, {"Assay"}, cols, headers);
    std::unordered_map<std::string, std::string> cols_assay = map_variable(input_csv_file, cols, "Assay", headers);
    std::unordered_map<std::string, std::string> cols_sample = map_variable(input_csv_file, cols, "Sample", headers);
    um_str_vstr cols_Ct = map_variable_vec(input_csv_file, cols, {val_var}, headers);
    std::unordered_map<std::string, float> Ct_means = um_mean(cols_Ct);
    std::unordered_map<std::string, float> Ct_sd = um_sd(cols_Ct);
    std::unordered_map<std::string, float> cols_efficiency = um_mean(map_variable_vec(input_csv_file, cols, {"Efficiency"}, headers));

    std::vector<std::string> keys;
    for (auto const& it : cols_assay) {keys.push_back(it.first);}
    std::sort(keys.begin(), keys.end(), numeric_string_compare);

    std::vector<std::string> assay;
    for (auto it : assay_cols) {assay.push_back(it.first);}
    sort_numeric_strings(assay);

// QC - Ct of NTC needs to be 3+ Cq greater than Control1 for each assay
    std::unordered_map<std::string, float> NTC_means;
    std::unordered_map<std::string, float> STD_means;
    std::unordered_map<std::string, std::string> QC_NTC;
    for (auto it : assay) {
        NTC_means[it] = Ct_means[search_vsrting(assay_cols[it], non_template)];
        STD_means[it] = Ct_means[search_vsrting(assay_cols[it], standard + "1")];
        if ((NTC_means[it]-STD_means[it]) < 3) {QC_NTC[it] = "FAIL";} else {QC_NTC[it] = "PASS";}
    }
    
// QC - Ct of negative control needs to be >35 or not amplified for each assay
    std::unordered_map<std::string, float> NEG_means;
    std::unordered_map<std::string, std::string> QC_NEG;
    for (auto it : assay) {
        NEG_means[it] = Ct_means[search_vsrting(assay_cols[it], neg_control)];
        if (NEG_means[it] < 35) {QC_NEG[it] = "FAIL";} else {QC_NEG[it] = "PASS";}
    }

// Calculate standard curves - pull out slope, intercept and R2 for efficiency calculation and QC
    std::unordered_map<std::string, int> log_abundance_map;
    for (int i=0; i <= 10; i++) {
        log_abundance_map[standard+std::to_string(i)] = i;
    }
    std::unordered_map<std::string, std::pair<float, float> > regression_map;
    std::unordered_map<std::string, float> rsqr_map;
    std::unordered_map<std::string, float> std_efficiency_map;
    std::unordered_map<std::string, std::string> std_QC;
    for (auto it : assay) {
        std::vector<float> assay_means;
        std::vector<float> log_abundance;
        for (auto i : assay_cols[it]) {
            if (i.find(standard) != std::string::npos) {
                assay_means.push_back(Ct_means[i]);
                log_abundance.push_back(log_abundance_map[cols_sample[i]]);
            }
        }
        regression_map[it] = lm(log_abundance, assay_means);
        rsqr_map[it] = coef_determination(log_abundance, assay_means);
        std_efficiency_map[it] = std::pow(10, -1/regression_map[it].second);
        // QC - PCR efficiency must be between 1.75-2.1 and R2 > .975     
        if (std_efficiency_map[it] > 1.75 && std_efficiency_map[it] < 2.1 && rsqr_map[it] > 0.975) {std_QC[it] = "PASS";} else {std_QC[it]="FAIL";}
    }

    std::ofstream assay_report_file(output + "_assay_QC_report.csv");
    assay_report_file << "Assay," << "STD_Efficiency," << "Rsqr," << "QC_StdCurve," << "NEG_Ct," << "QC_NEG," << "NTC_diff," << "QC_NTC," << "\n";
    for (auto it : assay) {
        assay_report_file << it << ",";
        assay_report_file << std_efficiency_map[it] << ",";
        assay_report_file << rsqr_map[it] << ",";
        assay_report_file << std_QC[it] << ",";
        assay_report_file << NEG_means[it] << ",";
        assay_report_file << QC_NEG[it] << ",";
        assay_report_file << NTC_means[it] - STD_means[it] << ",";
        assay_report_file << QC_NTC[it] << ",";
        assay_report_file << std::endl;
    }

    std::unordered_map<std::string, std::string> assays_map = map_csv(assays_file, "Assay");
    std::ofstream sample_report_file(output + "_sample_QC_report.csv");
    sample_report_file << "Assay," << "Sample," << "Mean_Copy_N," << "Sd_Copy_N," << "Mean_Efficiency," << "QCSample,";
    sample_report_file  << "Cycle," << "FunctionalGroup," << "GeneClass," << "Measure," << "ActualPrimerSet," << "\n";
    for (auto it : keys) {
        std::string key_assay = cols_assay[it];
        sample_report_file << key_assay << ",";
        sample_report_file << cols_sample[it] << ",";
            std::vector<float> copy_N;
            for (auto i : cols_Ct[it]) {
                copy_N.push_back(std::pow(10, (std::stof(i) - regression_map[key_assay].first)/regression_map[key_assay].second));
            }
        sample_report_file << mean(copy_N) << ",";
        sample_report_file << sd(copy_N) << ",";
        sample_report_file << cols_efficiency[it] << ",";
            std::string sample_QC = "FAIL";
            if (cols_efficiency[it] >= 1.75) {sample_QC="PASS";}
        sample_report_file << sample_QC << ",";
        sample_report_file << assays_map[key_assay] << "\n";
    }

    
    std::ofstream all_report_file(output + "_sample_qpcr_output_with_assay_info_qc.csv");
    all_report_file << "Assay," << "Sample," << "meanCt," << "stderr," << "meanEffi," << "QCSample,";
    all_report_file << "STD_Efficiency," << "Rsqr," << "QC_StdCurve," << "NEG_Ct," << "QC_NEG," << "NTC_diff," << "QC_NTC,";
    all_report_file  << "Cycle," << "FunctionalGroup," << "GeneClass," << "Measure," << "ActualPrimerSet," << "\n";
    for (auto it : keys) {
        std::string key_assay = cols_assay[it];
        sample_report_file << key_assay << ",";
        all_report_file << cols_sample[it] << ",";
            std::vector<float> copy_N;
            for (auto i : cols_Ct[it]) {
                copy_N.push_back(std::pow(10, (std::stof(i) - regression_map[key_assay].first)/regression_map[key_assay].second));
            }
        all_report_file << mean(copy_N) << ",";
        all_report_file << sd(copy_N) << ",";
        all_report_file << cols_efficiency[it] << ",";
            std::string sample_QC = "FAIL";
            if (cols_efficiency[it] >= 1.75) {sample_QC="PASS";}
        all_report_file << sample_QC << ",";
        all_report_file << std_efficiency_map[key_assay] << ",";
        all_report_file << rsqr_map[key_assay] << ",";
        all_report_file << std_QC[key_assay] << ",";
        all_report_file << NEG_means[key_assay] << ",";
        all_report_file << QC_NEG[key_assay] << ",";
        all_report_file << NTC_means[key_assay] - STD_means[key_assay] << ",";
        all_report_file << QC_NTC[key_assay] << ",";
        all_report_file << assays_map[key_assay] << "\n";
    }

    return 0;
}





int main(int argc, char **argv)
{
	std::string input;
	std::string assays;
	std::string output;
    std::string value;
    std::string neg_control;
    std::string standard;
    std::string non_template;
    bool headers;
    // declare options
    struct option options[9];

    options[0].long_name  = "help";
    options[0].short_name = 'h';
    options[0].flags      = GOPT_ARGUMENT_FORBIDDEN;

    options[1].long_name  = "input";
    options[1].short_name = 'i';
    options[1].flags      = GOPT_ARGUMENT_REQUIRED;

    options[2].long_name  = "assays";
    options[2].short_name = 'a';
    options[2].flags      = GOPT_ARGUMENT_REQUIRED;

    options[3].long_name  = "output";
    options[3].short_name = 'o';
    options[3].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[4].long_name  = "value-var";
    options[4].short_name = 'v';
    options[4].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[5].long_name  = "neg_control";
    options[5].short_name = 'n';
    options[5].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[6].long_name  = "standard";
    options[6].short_name = 's';
    options[6].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[7].long_name  = "ntc";
    options[7].short_name = 'c';
    options[7].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[8].long_name  = "headers";
    options[8].short_name = 'r';
    options[8].flags      = GOPT_ARGUMENT_OPTIONAL;

    options[9].flags      = GOPT_LAST;

    argc = gopt(argv, options);
    gopt_errors(argv[0], options);
    if (options[0].count) {
        fprintf (stdout, "see the manual\n");
        exit (EXIT_SUCCESS);
    }
    if (options[1].count) {
        input = options[1].argument;
    } else {exit(EXIT_FAILURE);}
    if (options[2].count) {
        assays = options[2].argument;
    } else {
        exit(EXIT_FAILURE);}
    if (options[3].count) {
        output = options[3].argument;
    } else {
        output = input.substr(0, input.find_last_of(".")); 
    }
    if (options[4].count) {
        value = options[4].argument;
    } else {value = "Ct";}
    if (options[5].count) {
        neg_control = options[5].argument;
    } else {neg_control = "NEG";}
    if (options[6].count) {
        standard = options[6].argument;
    } else {standard = "STD";}
    if (options[7].count) {
        non_template = options[7].argument;
    } else {non_template = "NTC";}
    if (options[8].count) {
        headers = options[8].argument;
    } else {headers = true;}

	smart_chip_analyzer(input, assays, output, value, {"Assay", "Sample"}, neg_control, standard, non_template, headers);
	return(0);
}