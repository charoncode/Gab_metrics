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
#include <srvf/opencurves.h>
#include <srvf/plf.h>
#include <srvf/rotate.h>
#include <srvf/reparam.h>

#include <cmath>
#include <iostream>


namespace srvf
{

namespace opencurves
{

double shape_distance (const Srvf &Q1, const Srvf &Q2, 
                       bool optimize_rots, 
                       bool optimize_scale, 
                       bool optimize_reparams, 
                       size_t nrounds, 
                       size_t dp_grid_width, 
                       size_t dp_grid_height)
{
  Srvf Q1l(Q1), Q2l(Q2);

  if (optimize_scale)
  {
    Q1l.scale_to_unit_norm();
    Q2l.scale_to_unit_norm();
  }

  if (optimize_rots && optimize_reparams)
  {
    std::vector<double> tv1, tv2;
    if (dp_grid_width > 0) 
      tv1 = srvf::util::linspace(Q1.domain_lb(),Q1.domain_ub(),dp_grid_width);
    if (dp_grid_height > 0) 
      tv2=srvf::util::linspace(Q2.domain_lb(),Q2.domain_ub(),dp_grid_height);
    
    const std::vector<double> &tv1r = (dp_grid_width > 0 ? tv1 : Q1.params());
    const std::vector<double> &tv2r = (dp_grid_height > 0 ? tv2 : Q2.params());

    Matrix R = optimal_rotation(Q1l, Q2l);
    Q2l.rotate(R);

    for (size_t i=0; i<nrounds; ++i)
    {
      Plf G = optimal_reparam(Q1l, Q2l, tv1r, tv2r);
      Q2l = gamma_action(Q2l, G);
      R = optimal_rotation(Q1l, Q2l);
      Q2l.rotate(R);
    }

    if (optimize_scale)
      return sphere_distance(Q1l, Q2l);
    else
      return l2_distance(Q1l, Q2l);
  }
  else if (optimize_rots)
  {
    Matrix R = optimal_rotation(Q1l, Q2l);
    Q2l.rotate(R);

    if (optimize_scale)
      return sphere_distance(Q1l, Q2l);
    else
      return l2_distance(Q1l, Q2l);
  }
  else if (optimize_reparams)
  {
    std::vector<double> tv1, tv2;
    if (dp_grid_width > 0) 
      tv1 = srvf::util::linspace(Q1.domain_lb(),Q1.domain_ub(),dp_grid_width);
    if (dp_grid_height > 0) 
      tv2=srvf::util::linspace(Q2.domain_lb(),Q2.domain_ub(),dp_grid_height);
    
    const std::vector<double> &tv1r = (dp_grid_width > 0 ? tv1 : Q1.params());
    const std::vector<double> &tv2r = (dp_grid_height > 0 ? tv2 : Q2.params());

    Plf G = optimal_reparam(Q1l, Q2l, tv1r, tv2r);
    Q2l = gamma_action(Q2l, G);

    if (optimize_scale)
      return sphere_distance(Q1l, Q2l);
    else
      return l2_distance(Q1l, Q2l);
  }
  else
  {
    if (optimize_scale)
      return sphere_distance(Q1l, Q2l);
    else
      return l2_distance(Q1l, Q2l);
  }
}

/**
 * Computes the Karcher mean of a collection of points on the unit 
 * sphere in \f$ L^2(I,R^n) \f$.
 */
Srvf karcher_mean(const std::vector<Srvf> &Qs, double tol, size_t max_iters, 
                  bool optimize_rots, bool optimize_reparams)
{
  // Corner case:  if Qs is empty, return an empty Srvf
  if (Qs.empty())
  {
    Srvf empty;
    return empty;
  }

  double radius = srvf::l2_norm(Qs[0]);
  double stepsize = 0.25;

  // Initialize mean to one of the Qs.  This is arbitrary, and not necessarily 
  // the best solution.  But (hopefully) good enough.
  Srvf Mu(Qs[0]);
  
  for (size_t iter=1; (max_iters==0 || iter<=max_iters); ++iter)
  {
    // The update direction is the average of the shooting directions.
    Srvf update_dir(Mu.domain_lb(), Mu.domain_ub(), 
                    std::vector<double>(Mu.dim(),0.0));

#pragma omp parallel
{
    #pragma omp for schedule (static) nowait
    for (size_t i=0; i<Qs.size(); ++i)
    {
      Srvf Svi;

      if (optimize_rots && optimize_reparams)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        // Rotate
        Matrix Ri = optimal_rotation(Mu, Qsi);
        Qsi.rotate(Ri);

        // Reparametrize
        Plf Gi = optimal_reparam(Mu, Qsi);
        Srvf Qsir = gamma_action(Qsi, Gi);

        // And rotate again
        Ri = optimal_rotation(Mu, Qsir);
        Qsir.rotate(Ri);

        // Compute current shooting vector
        Svi = shooting_vector(Mu, Qsir);
      }
      else if (optimize_rots)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        Matrix Ri = optimal_rotation(Mu, Qsi);
        Qsi.rotate(Ri);
        Svi = shooting_vector(Mu, Qsi);
      } 
      else if (optimize_reparams)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        Plf Gi = optimal_reparam(Mu, Qsi);
        Srvf Qsir = gamma_action(Qsi, Gi);
        Svi = shooting_vector(Mu, Qsir);
      }
      else
      {
        Svi = shooting_vector(Mu, Qs[i]);
      }

      #pragma omp critical
      update_dir = srvf::linear_combination(update_dir, Svi, 1.0, 1.0);
    }
}

    update_dir.scale(1.0 / (double)Qs.size());
    double udr = l2_norm(update_dir);
    std::cout << iter << ") karcher_mean(): gradient norm = " << udr 
              << std::endl;
    if (udr < tol)
    {
      // Norm of gradient is small enough to stop here
      break;
    }
    else
    {
      // Update mean estimate and keep going
      Mu = srvf::linear_combination(Mu, update_dir, 1.0, stepsize);
      Mu.scale(radius / l2_norm(Mu));
    }
  }

  return Mu;
}


/**
 * Returns the shooting vector from \a Q1 to \a Q2.
 *
 * \a Q1 and \a Q2 must have the same \f$ L^2 \f$ norm.
 * 
 * The shooting vector from \a Q1 to \a Q2 is the tangent vector at \a Q1 
 * of a spherical geodesic from \a Q1 to \a Q2.  Its length is equal to the 
 * geodesic distance between the two points.
 */
Srvf shooting_vector(const Srvf &Q1, const Srvf &Q2)
{
  double radius = srvf::l2_norm(Q1);
  double radius2 = radius * radius;

  // Check for case where Q1 and Q2 are close to zero
  if (radius2 < 1e-8)
  {
    return Srvf(Q1);
  }

  // Angle between Q1 and Q2
  double cos_theta = srvf::l2_product(Q1, Q2) / radius2;
  if (cos_theta >  1.0) cos_theta =  1.0;
  if (cos_theta < -1.0) cos_theta = -1.0;
  double theta = acos(cos_theta);

  // Shooting vector
  Srvf sv = srvf::linear_combination(Q1, Q2, -cos_theta, 1.0);

  // Scale shooting vector to have length equal to great-circle distance
  double norm_sv = srvf::l2_norm(sv);
  if (norm_sv > 1e-6)
  {
    sv.scale(radius * theta / norm_sv);
  }

  return sv;
}

} //namespace opencurves
} // namespace srvf
