/*
 *
 * Author:  Schuyler D. Smith
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

#ifndef DEFS_SDS
#define DEFS_SDS

typedef std::vector<std::string>                                    vstring;
typedef std::vector<int>                                            vint;
typedef std::vector<float>                                          vfloat;
typedef std::vector<double>                                         vdouble;
typedef std::unordered_map<std::string, std::string>                um_str_str;
typedef std::unordered_map<std::string, int>                        um_str_int;
typedef std::unordered_map<std::string, float>                      um_str_flo;
typedef std::unordered_map<std::string, double>                     um_str_dbl;
typedef std::unordered_map<std::string, std::vector<std::string> >  um_str_vstr;
typedef std::unordered_map<std::string, std::vector<float> >        um_str_vflo;
typedef std::unordered_map<std::string, std::vector<double> >       um_str_vdbl;
typedef std::unordered_map<std::string, std::pair<float, float> >   um_str_pair_flo_flo;
typedef std::unordered_map<std::string, std::pair<double, double> > um_str_pair_dbl_dbl;

#endif