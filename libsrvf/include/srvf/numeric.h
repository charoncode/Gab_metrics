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
#ifndef SRVF_NUMERIC_H
#define SRVF_NUMERIC_H 1

#include <cmath>


namespace srvf
{

namespace numeric
{

/**
 * Floating point comparison test.
 *
 * Reference:
 * http://www.cygnus-software.com/papers/comparingfloats/comparingfloats.htm
 */
inline bool almost_equal(double a, double b, 
                         double abs_thresh=1e-6, double rel_thresh=1e-3)
{
  double abs_diff = fabs(a - b);
  if (abs_diff <= abs_thresh) return true;

  double amag = fabs(a);
  double bmag = fabs(b);
  double largest = (amag > bmag ? amag : bmag);

  return (abs_diff <= rel_thresh * largest);
}

} // namespace numeric

} // namespace srvf

#endif // SRVF_NUMERIC_H
