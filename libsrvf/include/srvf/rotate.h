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
#ifndef SRVF_ROTATE_H
#define SRVF_ROTATE_H 1

#include "srvf.h"
#include "matrix.h"


namespace srvf
{

Matrix optimal_rotation (const Srvf &Q1, const Srvf &Q2);

Matrix optimal_rotation(const Srvf &Q1, const Srvf &Q2, 
  double a, double b, double c, double d, 
  size_t Q1_start_idx=(size_t)(-1), size_t Q2_start_idx=(size_t)(-1) );

} // namespace srvf

#endif // SRVF_ROTATE_H
