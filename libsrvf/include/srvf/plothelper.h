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
#ifndef PLOTHELPER_H
#define PLOTHELPER_H 1

#include <srvf/plf.h>
#include <srvf/srvf.h>

#include <FL/Fl.h>


namespace srvf
{
namespace plot
{

/**
 * High-level convenience routine for plotting a vector of 1-D Plf's.
 */
void plot_1d_plfs(const std::vector<Plf> &v, 
  const std::vector<Color> &colors,
  size_t x=0, size_t y=0, size_t w=800, size_t h=400, 
  const char *title="-");

/**
 * High-level convenience routine for plotting a vector of 1-D Srvf's.
 */
void plot_1d_srvfs(const std::vector<Srvf> &v, 
  const std::vector<Color> &colors,
  size_t x=0, size_t y=0, size_t w=800, size_t h=400, 
  const char *title="-");


} // namespace plot
} // namespace srvf

#endif // PLOTHELPER_H
