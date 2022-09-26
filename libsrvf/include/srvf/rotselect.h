/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  FSU Statistical Shape Analysis and Modeling Group
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
#ifndef ROTSELECT_H
#define ROTSELECT_H 1

#include "paretoset.h"
#include <srvf/matrix.h>
#include <srvf/srvf.h>

#include <vector>


namespace srvf
{
namespace pmatch
{


/**
 * Computes optimal rotations for all partial matches in the given grid.
 *
 * The rotations are computed relative to the given parametrizations 
 * of \a Q1 and \a Q2.  If the curves were parametrized differently, then 
 * this function would generally yield different results.
 */
std::vector<Matrix> select_rotations (
  const Srvf &Q1, const Srvf &Q2, 
  const std::vector<double> &tv1, const std::vector<double> &tv2 );


} // namespace pmatch
} // namespace srvf

#endif // ROTSELECT_H
