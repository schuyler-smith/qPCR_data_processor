/*
 *
 * Author:  Schuyler D. Smith
 *
 */

#ifndef MAP_SDS
#define MAP_SDS

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


typedef std::unordered_map<std::string, std::vector<std::string> >  um_str_vstr;
typedef std::multimap<std::string, std::string>                     mm_str_str;
typedef std::unordered_map<std::string, std::vector<float> >        um_str_vflo;
typedef std::unordered_map<std::string, std::vector<double> >       um_str_vdbl;

auto extract_keys(std::unordered_map<std::string, float> const& input_map)
{
	std::vector<std::string> keys;
	for (auto const& element : input_map) {keys.push_back(element.first);}
  return keys;
}

auto extract_values(std::unordered_map<std::string, float> const& input_map)
{
	std::vector<float> values;
	for (auto const& element : input_map) {values.push_back(element.second);}
	return values;
}

auto map_headers(std::string input_csv_file)
{
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

auto map_variable_vec(std::string input_csv_file, std::vector<std::string> variables, std::vector<std::string> values, bool headers=true)
{
  std::ifstream csv_file(input_csv_file);
  std::string line;
  std::vector<std::string> row_data;
  std::unordered_map<std::string, int> colnames;
  colnames = map_headers(input_csv_file);
  um_str_vstr variable_map;
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

auto map_variable_vec_numeric(std::string input_csv_file, std::vector<std::string> variables, std::vector<std::string> values, bool headers=true)
{
  std::ifstream csv_file(input_csv_file);
  std::string line;
  std::vector<std::string> row_data;
  std::unordered_map<std::string, int> colnames;
  colnames = map_headers(input_csv_file);
  um_str_vdbl variable_map;
  double value;
  while(std::getline(csv_file, line)) {
    if (headers) {headers=false; continue;}
    row_data = string_split(line,',');
    std::string key, val;
    for (int i = 0; i < variables.size(); ++i) {
      key = key.append(row_data[colnames[variables[i]]]);
    }
    for (int i = 0; i < values.size(); ++i) {
      val = row_data[colnames[values[i]]];
      if (val.empty()) {val = "NAN";}
      value = std::stof(val);
    }
    if (variable_map.count(key)) {
      if (std::find(variable_map[key].begin(), variable_map[key].end(), value) != variable_map[key].end()) {continue;
      } else {variable_map[key].push_back(value);}
    } else{variable_map.insert(std::pair<std::string, std::vector<double> >(key, {value}));}
  }
  return variable_map;
}

auto um_percent_below_threshold(um_str_vdbl val_map, double threshold)
{
  std::unordered_map<std::string, double> perc;
  for (auto it : val_map) {
    double i = 0;
    std::vector<double> values(it.second);
    for (auto val : values) {
      if (val <= threshold) {
        ++i;
      }
    }
    perc[it.first] = i/values.size();
  }
  return perc;
}

auto um_mean(um_str_vdbl val_map)
{
  std::unordered_map<std::string, double> means;
  for (auto it : val_map) {
    std::vector<double> values(it.second);
    means[it.first] = mean(values);
  }
  return means;
}

auto um_sd(um_str_vdbl val_map)
{
  std::unordered_map<std::string, double> means;
  std::unordered_map<std::string, double> sd;
  means = um_mean(val_map);
  for (auto it : val_map) {
    std::vector<double> values(it.second);
    double variance = 0;
    for (auto v : values) {
      variance += std::pow((v - means[it.first]), 2);
    }
    sd[it.first] = std::sqrt(variance);
  }
  return sd;
}

auto map_variable(std::string input_csv_file, std::vector<std::string> variables, std::string value, bool headers=true)
{
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

auto map_value_as_numeric(std::string input_csv_file, std::vector<std::string> variables, std::string value, bool headers=true)
{
  std::unordered_map<std::string, double> numeric_map;
  std::unordered_map<std::string, std::string> val_map = map_variable(input_csv_file, variables, value, headers);
  for (auto it : val_map) {
    if (it.second.empty()) {it.second = "NAN";}
    numeric_map[it.first] = std::stof(it.second);
  }
  return numeric_map;
}

auto map_csv(std::string input_csv_file, std::string variable, bool headers=true)
{
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


void print_mm(mm_str_str& myContainer)
{
  for (auto pr : myContainer) {std::cout << pr.first << ", " << pr.second << '\n';}
}

void print_um(um_str_vstr const &m)
{
  for (auto const &pair: m) {
    std::string sec;
    for (int i = 0; i < pair.second.size(); ++i) {
      sec.append(", ");
      sec.append(pair.second[i]);
    }
    std::cout << "{" << pair.first 
    << ": " 
    << sec 
    << "}\n";
  }
}


auto file_map(std::string file_name)
{
  std::unordered_map<std::string, float> out_map;
  std::ifstream file(file_name);
  std::string line;
  while(std::getline(file, line))
  {
    std::vector<std::string> words;
    words = string_split(line,':');
    out_map[words[0]] = std::stof(words[1]);
  }
  return out_map;
}


#endif

