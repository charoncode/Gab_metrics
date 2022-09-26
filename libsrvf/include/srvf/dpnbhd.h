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
#ifndef DPNBHD_H
#define DPNBHD_H 1

#include <cstddef>

namespace srvf
{

#ifndef DP_NBHD_DIM
#define DP_NBHD_DIM 7 
#endif

#if DP_NBHD_DIM == 17
#define DP_NBHD_SIZE 191
extern size_t dp_nbhd[DP_NBHD_SIZE][2];
#elif DP_NBHD_DIM == 12
#define DP_NBHD_SIZE 91
extern size_t dp_nbhd[DP_NBHD_SIZE][2];
#elif DP_NBHD_DIM == 10
#define DP_NBHD_SIZE 63
extern size_t dp_nbhd[DP_NBHD_SIZE][2];
#elif DP_NBHD_DIM == 7
#define DP_NBHD_SIZE 35
extern size_t dp_nbhd[DP_NBHD_SIZE][2];
#else // DP_NBHD_DIM = 6 (default)
#define DP_NBHD_SIZE 23
extern size_t dp_nbhd[DP_NBHD_SIZE][2];
#endif

} // namespace srvf

#endif // DPNBHD_H
