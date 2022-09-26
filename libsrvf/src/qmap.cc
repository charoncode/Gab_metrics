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
#include <srvf/qmap.h>
#include <srvf/pointset.h>
#include <srvf/util.h>

#include <cstddef>
#include <cmath>


namespace srvf
{

/**
 * Returns a new \c Srvf representing the square-root velocity function of 
 * the given \c Plf.
 *
 * \param F a \c Plf
 * \return a new \c Srvf representing the SRVF of \a F
 */
Srvf plf_to_srvf(const Plf &F)
{
  Pointset dF=srvf::util::diff(F.samps(),F.params());
  for (size_t i=0; i<dF.npts(); ++i)
  {
    double rnqi=sqrt(dF.norm(i));
    if (rnqi>1e-6)
    {
      dF.scale(i,1.0/rnqi);
    }
    else
    {
      for (size_t j=0; j<dF.dim(); ++j)
      {
        dF[i][j]=0.0;
      }
    }
  }

  return Srvf(dF,F.params());
}


Plf srvf_to_plf(const Srvf &Q)
{
  Pointset samps(Q.dim(), Q.ncp());

  for (size_t i=0; i<Q.dim(); ++i)
  {
    samps[0][i] = 0.0;
  }
  
  for (size_t i=1; i<samps.npts(); ++i)
  {
    double nvi = Q.samps().norm(i-1);
    double dt = Q.params()[i] - Q.params()[i-1];
    weighted_sum(samps, Q.samps(), i-1, i-1, 1.0, nvi*dt, samps, i);
  }

  return Plf(samps, Q.params());
}

} // namespace srvf
