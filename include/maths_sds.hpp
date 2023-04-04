/*
 *
 * Author:  Schuyler D. Smith
 *
 */

#ifndef MATHS_SDS
#define MATHS_SDS


#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <math.h>
#include <tuple>


auto which_nan(std::vector<double> values) {
  std::vector<int> nan_vals;
  for (int i = 0; i < values.size(); ++i) {
    if (values[i] != values[i]) {
      nan_vals.push_back(i);
    }
  }
  return nan_vals;
}

auto rm_nan(std::vector<double> values, std::vector<int> nan_indices = {}) {
  if (nan_indices.size() == 0) {nan_indices = which_nan(values);}
  if (nan_indices.size() > 0) {
    std::sort(nan_indices.begin(), nan_indices.end(), std::greater<int>()); 
    for (auto nan_i : nan_indices) {
      values.erase(values.begin() + nan_i);
    }
  }
  if (values.size() == 0) {values.push_back(NAN);}
  return values;
}
bool is_numeric(const std::string& s) {
  return !s.empty() && std::find_if (s.begin(), s.end(), [](unsigned char c) { return !std::isdigit(c); }) == s.end();
}

auto mean(std::vector<double> values, bool rm_na = true) {
  if (rm_na) {values = rm_nan(values);}
  if (values.empty()) {return std::numeric_limits<double>::quiet_NaN();}
  double mean;
  auto const n = static_cast<int>(values.size());
  mean = std::accumulate(values.begin(), values.end(), 0.0) / n;
  return mean;
}

auto variance(std::vector<double> x, bool rm_na = true) {
  if (rm_na) {x = rm_nan(x);}
  if (x.empty()) {return std::numeric_limits<double>::quiet_NaN();}
  double X = mean(x);
  auto const n = static_cast<int>(x.size());
  double variance = 0;
  for (auto v : x) {
    variance += std::pow((v - X), 2);
  }
  variance = variance/(n-1);
  return variance;
}

auto sd(std::vector<double> values, bool rm_na = true) {
  return std::sqrt(variance(values, rm_na));
}

auto magnify(std::vector<double> values, float scalar, bool rm_na = true) {
  if (rm_na) {values = rm_nan(values);}
  for (int i = 0; i < values.size(); i++) {
    values[i] *= scalar;
  }
  return values;
}

auto rm_nan_pairs(std::vector<double> x, std::vector<double> y) {
  std::vector<int>    na_vals     = which_nan(x);
  std::vector<int>    y_na_vals   = which_nan(y);
  na_vals.insert(std::end(na_vals), std::begin(y_na_vals), std::end(y_na_vals));
  std::vector<double>  X   = rm_nan(x, na_vals);
  std::vector<double>  Y   = rm_nan(y, na_vals);
  return std::make_tuple(X, Y);
}

auto covariance(std::vector<double> x, std::vector<double> y, bool rm_na = true) {
  if (rm_na) {tie(x, y) = rm_nan_pairs(x, y);}
  double X = mean(x);
  double Y = mean(y);
  const auto N = x.size();
  double covar = 0;
  for (int i = 0; i < N; ++i) {
    covar += (x[i]-X) * (y[i] - Y);
  }
  return covar/(N-1);
}

auto lm(std::vector<double> X, std::vector<double> Y) {
  tie(X, Y) = rm_nan_pairs(X, Y);
  const auto n   = X.size();
  // for (int i = 0; i < n; ++i) {Y[i] -= 1;}
  const auto sX  = std::accumulate(X.begin(), X.end(), 0.0);
  const auto sY  = std::accumulate(Y.begin(), Y.end(), 0.0);
  const auto sXX = std::inner_product(X.begin(), X.end(), X.begin(), 0.0);
  const auto sXY = std::inner_product(X.begin(), X.end(), Y.begin(), 0.0);
  const auto b   = ((n * sXY) - (sX * sY)) / ((n * sXX) - (sX * sX));
  const auto a   = (sY - (b * sX)) / n;
  return std::make_pair(a,b);
}

auto coef_determination(std::vector<double> X, std::vector<double> Y) {
  tie(X, Y) = rm_nan_pairs(X, Y);
  const auto n  = X.size();
  double Y_hat = mean(Y);
  std::pair<double, double> ab = lm(X, Y);
  std::vector<double> res = Y;
  for (auto& it : res) {
    it -= Y_hat;
    it = it * it;
  }
  const auto SSE = std::accumulate(res.begin(), res.end(), 0.0);
  std::vector<double> pred = X;
  for (auto& it : pred) {
    it = it * ab.second + ab.first;
    it -= Y_hat;
    it = it * it;
  }
  const auto SSR  = std::accumulate(pred.begin(), pred.end(), 0.0);
  return SSR/SSE;
}

void sort_numeric_strings(std::vector<std::string> &v) {
  std::vector<int> temp;
  for (auto it : v) {temp.push_back(std::stoi(it));}
  std::sort(temp.begin(), temp.end()); 
  for (int i = 0; i < temp.size(); ++i) {
    v[i] = std::to_string(temp[i]);
  }
}

auto string_split(std::string str, const char& delim) {
  std::string word;
  std::vector<std::string> result_vector;
  str.erase(std::remove(str.begin(), str.end(), '\n'), str.cend());
  for (std::string::const_iterator character = str.begin(); character != str.end(); character++) {
    if (*character == delim) {
      // if (!word.empty()) {
        result_vector.push_back(word);
        word.clear();
      // }
    } else {word += *character;}
  }
  if (!word.empty()) {result_vector.push_back(word);}
  return result_vector;
}

auto search_vsrting(std::vector<std::string> v, std::string pattern) {
  std::string match;
  for (int i = 0; i < v.size(); ++i) {
    if (v[i].find(pattern) != std::string::npos) {
      match = v[i];
      break;
    } else{match = "NA";}
  }
  return match;
}

bool is_not_digit(char c) {
  return !std::isdigit(c);
}

bool numeric_string_compare(const std::string& s1, const std::string& s2) {
  // handle empty strings...
  std::string::const_iterator it1 = s1.begin(), it2 = s2.begin();
  if (std::isdigit(s1[0]) && std::isdigit(s2[0])) {
    int n1, n2;
    std::stringstream ss(s1);
    ss >> n1;
    ss.clear();
    ss.str(s2);
    ss >> n2;
    if (n1 != n2) return n1 < n2;
    it1 = std::find_if (s1.begin(), s1.end(), is_not_digit);
    it2 = std::find_if (s2.begin(), s2.end(), is_not_digit);
  }
  return std::lexicographical_compare(it1, s1.end(), it2, s2.end());
}


#endif

