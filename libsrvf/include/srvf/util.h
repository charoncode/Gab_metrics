/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#ifndef SRVF_UTIL_H
#define SRVF_UTIL_H 1

#include "pointset.h"
#include "matrix.h"

#include <cstddef>
#include <vector>
#include <map>


namespace srvf
{

namespace util
{

std::vector<double> 
linspace(double a, double b, size_t n);

std::vector<double> 
random_vector(size_t len, double first=0.0, int dir=0);

std::vector<double> 
unique (std::vector<double> v1, std::vector<double> v2, double thresh=1e-6);

Pointset diff(const Pointset &X);
Pointset diff(const Pointset &X, const std::vector<double> &tv);

std::map<size_t,size_t> build_lookup_map (
  const std::vector<double> &tv1, const std::vector<double> &tv2 );

} // namespace srvf::util

} // namespace srvf

#endif // SRVF_UTIL_H
