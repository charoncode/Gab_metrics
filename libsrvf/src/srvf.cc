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
#include <srvf/srvf.h>
#include <srvf/plf.h>
#include <srvf/matrix.h>
#include <srvf/interp.h>

#include <cmath>


namespace srvf
{

/**
 * Evaluate this \c Srvf at a single point.
 *
 * \param t the parameter value at which the function will be evaluated
 * \return a \c Pointset containing the function value at \a t
 */
Point Srvf::evaluate(double t) const
{
  Pointset res = 
    srvf::interp::interp_const(samps(),params(),std::vector<double>(1,t));
  return res[0];
}

/**
 * Evaluate this \c Srvf at one or more points.
 *
 * The \a tv must be sorted in non-decreasing order.
 *
 * The result will be stored in \a result, which must have the correct 
 * dimension and at least \c tv.size() points.
 *
 * \param tv the parameter values at which the function will be evaluated
 * \param result receives result
 */
Pointset Srvf::evaluate(const std::vector<double> &tv) const
{
  return srvf::interp::interp_const(samps(), params(), tv);
}

/**
 * Apply a rotation to this \c Srvf.
 *
 * \param R a \c DxD rotation matrix, where \c D=this->dim().
 */
void Srvf::rotate(const Matrix &R)
{
  samps_.rotate(R);
}

/**
 * Apply a uniform scaling to this \c Srvf.
 *
 * \param sf the scale factor
 */
void Srvf::scale(double sf)
{
  samps_.scale(sf);
}

/**
 * Radial projection to the unit sphere in Hilbert space.
 *
 * If the current norm is zero, then this \c Srvf represents the zero 
 * function and no action is taken.
 */
void Srvf::scale_to_unit_norm()
{
  double nrm = l2_norm(*this);
  if (nrm > 1e-6)
  {
    scale(1.0 / nrm);
  }
}


////////////////////////////////////////////////////////////////////////////
///////////////////////// friend functions /////////////////////////////////
////////////////////////////////////////////////////////////////////////////


/**
 * Computes the L^2 norm of the given \c Srvf.
 *
 * The \f$L^2\f$ norm of a function \f$f:[a,b]\to R^n\f$ is defined as
 * \f[ |f| = \left( \int_a^b |f(t)| dt \right)^{1/2} \f]
 *
 * \param Q the \c Srvf
 * \return the L^2 norm of \a Q
 */
double l2_norm(const Srvf &Q)
{
  double res=0.0;
  for (size_t i=0; i<Q.ncp()-1; ++i)
  {
    double nqi=0.0;
    for (size_t j=0; j<Q.dim(); ++j)
    {
      double x=Q.samps()[i][j];
      nqi += x*x;
    }
    double dt = Q.params()[i+1]-Q.params()[i];
    res += nqi * dt;
  }
  return sqrt(res);
}

/**
 * Computes the L^2 inner product of two \c Srvf's.
 *
 * The \f$L^2\f$ inner product of two functions 
 * \f$f_1,f_2:[a,b]\to R^n\f$ is defined as
 *
 * \f[ \left\langle f_1,f_2 \right\rangle = 
 *     \int_a^b \left\langle f_1(t),f_2(t) \right\rangle dt \f]
 *
 * If either \a Q1 or \a Q2 has fewer than 2 sample points (i.e. is empty, 
 * or is defined on an interval of length 0), then the result is 0.
 *
 * \param Q1 the first \c Srvf
 * \param Q2 the second \c Srvf.  Must be defined on the same interval 
 *        as \a Q1.
 * \return the L^2 inner product of \a Q1 and \a Q2
 */
double l2_product(const Srvf &Q1, const Srvf &Q2)
{
  if (Q1.ncp()<2 || Q2.ncp()<2) return 0.0;

  double Q1a=Q1.params()[0], Q1b=Q1.params()[Q1.ncp()-1];
  double Q2a=Q2.params()[0], Q2b=Q2.params()[Q2.ncp()-1];

  if (Q1.dim()!=Q2.dim())
    throw std::invalid_argument("Q1 and Q2 must have the same dimension");
  if (fabs(Q1a-Q2a)>1e-4 || fabs(Q1b-Q2b)>1e-4)
    throw std::invalid_argument("Q1 and Q2 must have the same domain");

  size_t dim=Q1.dim();
  double tlast=(Q1a<Q2a ? Q1a : Q2a);
  size_t i1=0, i2=0;
  double res=0.0;

  while ((i1<Q1.ncp()-1) && (i2<Q2.ncp()-1))
  {
    double ipi=0.0;
    for (size_t j=0; j<dim; ++j)
    {
      ipi += Q1.samps()[i1][j]*Q2.samps()[i2][j];
    }

    double tnext1=Q1.params()[i1+1];
    double tnext2=Q2.params()[i2+1];
    double dt;
    if (tnext1<tnext2)
    {
      dt=tnext1-tlast;
      tlast=tnext1;
      ++i1;
      if (fabs(tnext1-tnext2)<1e-6) ++i2;
    }
    else
    {
      dt=tnext2-tlast;
      tlast=tnext2;
      ++i2;
      if (fabs(tnext1-tnext2)<1e-6) ++i1;
    }

    res += ipi*dt;
  }

  return res;
}

/**
 * Compute the L^2 distance between \a Q1 and \a Q2.
 *
 * The \f$L^2\f$ distance between two functions 
 * \f$f_1,f_2:[a,b]\to R^n\f$ is defined as
 *
 * \f[ |f_1-f_2| = \left(\int_a^b |f_1(t)-f_2(t)|^2 dt \right)^{1/2}\f]
 *
 * \param Q1 the first \c Srvf
 * \param Q2 the second \c Srvf.  Must have the same dimension and domain 
 *           as \a Q1.
 * \return the L^2 distance between \a Q1 and \a Q2
 */
double l2_distance(const Srvf &Q1, const Srvf &Q2)
{
  // TODO: refactor (this is basically the same thing as l2_product())
  if (Q1.ncp()<2 || Q2.ncp()<2) return 0.0;

  double Q1a=Q1.params()[0], Q1b=Q1.params()[Q1.ncp()-1];
  double Q2a=Q2.params()[0], Q2b=Q2.params()[Q2.ncp()-1];

  if (Q1.dim()!=Q2.dim())
    throw std::invalid_argument("Q1 and Q2 must have the same dimension");
  if (fabs(Q1a-Q2a)>1e-4 || fabs(Q1b-Q2b)>1e-4)
    throw std::invalid_argument("Q1 and Q2 must have the same domain");

  size_t dim=Q1.dim();
  double tlast=(Q1a<Q2a ? Q1a : Q2a);
  size_t i1=0, i2=0;
  double res=0.0;

  while ((i1<Q1.ncp()-1) && (i2<Q2.ncp()-1))
  {
    double ipi=0.0;
    for (size_t j=0; j<dim; ++j)
    {
      double dqi = Q1.samps()[i1][j]-Q2.samps()[i2][j];
      ipi += dqi*dqi;
    }

    double tnext1=Q1.params()[i1+1];
    double tnext2=Q2.params()[i2+1];
    double dt;
    if (tnext1<tnext2)
    {
      dt=tnext1-tlast;
      tlast=tnext1;
      ++i1;
      if (fabs(tnext1-tnext2)<1e-6) ++i2;
    }
    else
    {
      dt=tnext2-tlast;
      tlast=tnext2;
      ++i2;
      if (fabs(tnext1-tnext2)<1e-6) ++i1;
    }

    res += ipi*dt;
  }

  return sqrt(res);
}

/**
 * Computes the great circle distance between two SRVFs.
 *
 * The great circle distance between two functions 
 * \f$ q_1,q_2\in L^2([a,b],R^n) \f$
 *
 * having the same norm is given by 
 *
 * \f[ \theta=\cos^{-1}\left(\frac{\left\langle q_1, q_2 \right\rangle}
 *                                {|q_1||q_2|} \right) \f]
 *
 * \param Q1 the first \c Srvf
 * \param Q2 the second \c Srvf.  Must have the same dimension, domain, 
 *           and L^2 norm as \a Q1.
 * \return the great circle distance between \a Q1 and \a Q2.
 */
double sphere_distance(const Srvf &Q1, const Srvf &Q2)
{
  double nrm1=l2_norm(Q1);
  double nrm2=l2_norm(Q2);

  // TODO: is this a good tolerance?
  if (fabs(nrm1-nrm2)>5e-3)
    throw std::invalid_argument("Difference between norms exceeds tolerance.");

  if (fabs(nrm1*nrm2)>1e-6)
  {
    double ip=l2_product(Q1,Q2);
    double ct=ip/(nrm1*nrm2);
    if (ct<-1.0) ct=-1.0;
    if (ct>1.0)  ct=1.0;
    return acos(ct);
  }
  else
  {
    return 0.0;
  }
}

/**
 * Linear combination of \a Q1 and \a Q2 with the specified weights.
 *
 * Returns a new \c Srvf representing a linear combination of \a Q1 and 
 * \a Q2 with the specified weights \a w1 and \a w2.
 *
 * If \a Q1 is empty, then the result will be \c w2*Q2.  Similarly, if \a Q2 
 * is empty, then the result will be \c w1*Q1.
 *
 * \a Q1 and \a Q2 must be defined on the same interval.
 *
 * \param Q1 the first \c Srvf
 * \param Q2 the second \c Srvf
 * \param w1 the first weight
 * \param w2 the second weight
 * \return a new \c Srvf representing \f$ w_1Q_1 + w_2Q_2 \f$.
 */
Srvf linear_combination(const Srvf &Q1, const Srvf &Q2, 
                double w1, double w2)
{
  if ((fabs(Q1.domain_lb()-Q2.domain_lb()) > 1e-6) ||
      (fabs(Q1.domain_ub()-Q2.domain_ub()) > 1e-6) )
  { throw std::invalid_argument("Q1 and Q2 must have the same domain."); }
  if (Q1.dim() != Q2.dim())
  { throw std::invalid_argument("Q1 and Q2 must have the same dimension."); }


  if (Q1.is_empty())
  {
    Srvf res(Q2);
    res.samps().scale(w2);
    return res;
  }
  else if (Q2.is_empty())
  {
    Srvf res(Q1);
    res.samps().scale(w1);
    return res;
  }
  else
  {
    // Get new changepoints
    std::vector<double> params=srvf::util::unique(Q1.params(), Q2.params());
    size_t dim=Q1.dim();
    size_t ncp=params.size();

    // Evaluate Q1 and Q2 at the midpoint of each interval
    std::vector<double> tv(ncp-1);
    for (size_t i=0; i<ncp-1; ++i)
    {
      tv[i] = 0.5*(params[i] + params[i+1]);
    }
    Pointset samps1 = Q1.evaluate(tv);
    Pointset samps2 = Q2.evaluate(tv);

    // New samps is linear combination of samps1 and samps2
    Pointset samps(dim, ncp-1);
    for (size_t i=0; i<ncp-1; ++i)
    {
      weighted_sum(samps1,samps2,i,i,w1,w2,samps,i);
    }

    return Srvf(samps, params);
  }
}

/**
 * Returns a refinement of \a Q whose changepoint parameters include \a tv.
 *
 * \a tv is a 1-row matrix containing additional changepoints.
 *
 * The resulting \c Srvf will have the union of \c Q.params() and \a tv as 
 * its changepoint parameters.  On each of the resulting subintervals, the 
 * value of the new \c Srvf will be the same as the value of \a Q.
 *
 * \param Q reference to an existing \c Srvf
 * \param tv a \c Matrix containing the new changepoint parameters
 * \return a new \c Srvf representing the specified refinement of \a Q
 */
Srvf refinement(const Srvf &Q, const std::vector<double> &new_params)
{
  std::vector<double> params=srvf::util::unique(Q.params(), new_params);
  std::vector<double> tv(params.size()-1);
  for (size_t i=0; i<tv.size(); ++i)
  {
    tv[i] = 0.5*(params[i]+params[i+1]);
  }
  Pointset samps = Q.evaluate(tv);

  return Srvf(samps,params);
}

/**
 * Returns a new \c Srvf representing the action of \a gamma on \a Q.
 *
 * If \a Q is defined on the interval \f$ [a,b] \f$, then \a gamma must 
 * represent a non-decreasing map \f$ \gamma:[a,b]\to[a,b] \f$ sending 
 * \f$ a \mapsto a \f$ and \f$ b \mapsto b \f$.
 *
 * The action of \f$ \gamma \f$ on \a Q is defined as
 * \f[ Q*\gamma = \sqrt{\dot{\gamma}}(Q\circ\gamma) \f]
 *
 * \param Q an existing \c Srvf
 * \param gamma a \c Plf representing a reparametrization function
 */
Srvf gamma_action(const Srvf &Q, const Plf &gamma)
{
  if (gamma.dim() != 1)
  { throw std::invalid_argument("gamma must be 1-dimensional"); }

  if ((gamma.samps()[0][0] < Q.domain_lb()-1e-6) ||
      (gamma.samps()[gamma.ncp()-1][0] > Q.domain_ub()+1e-4))
  { throw std::invalid_argument(
      "Range of gamma must be contained in domain of Q"); }

  std::vector<double> xparams = gamma.preimages(Q.params());
  std::vector<double> params = srvf::util::unique(gamma.params(),xparams);
  Srvf Qr = refinement(Q,gamma.samps().to_vector());
  Pointset samps(Qr.samps());

  size_t dim=Q.dim();
  size_t ncp=params.size();
  for (size_t i=0; i<ncp-1; ++i)
  {
    double gdi = (Qr.params()[i+1]-Qr.params()[i]) / (params[i+1]-params[i]);
    double rgdi = sqrt(gdi);
    
    for (size_t j=0; j<dim; ++j)
    {
      samps[i][j] *= rgdi;
    }
  }

  return Srvf(samps,params);
}

/**
 * Returns a constant-speed parametrization of \a Q.
 */
Srvf constant_speed_param(const Srvf &Q)
{
  double Qnorm = l2_norm(Q);
  double int_width = Q.domain_ub() - Q.domain_lb();
  
  // Check for the zero function
  if (Qnorm == 0.0 || int_width == 0.0)
  {
    return Srvf(Q.domain_lb(), Q.domain_ub(), std::vector<double>(Q.dim(),0.0));
  }

  double v = Qnorm / int_width;

  // Initialize result to empty SRVF
  Srvf Qres;
  std::vector<double> partial_integrals;
  size_t last_nonzero=0; bool have_nonzero=false;

  // Build sample points
  for (size_t i=0; i<Q.samps().npts(); ++i)
  {
    double vi = Q.samps().norm(i);
    double dt = Q.params()[i+1] - Q.params()[i];
    double cur_pi = dt * vi * vi;

    if (fabs(cur_pi) < 1e-5) continue;

    if (have_nonzero && Q.samps().on_same_ray(i,last_nonzero))
    {
      partial_integrals[partial_integrals.size()-1] += cur_pi;
    }
    else
    {
      // Add new sample point and scale it to the correct norm
      Qres.samps().push_back(Q.samps()[i]);
      Qres.samps().scale(Qres.samps().npts()-1, v / vi);
      partial_integrals.push_back(cur_pi);
    }

    last_nonzero = i;
    have_nonzero = true;
  }
  
  // Build paramter vector
  Qres.params().push_back(Q.domain_lb());
  for (size_t i=0; i<partial_integrals.size(); ++i)
  {
    double tlast = Qres.params()[Qres.params().size()-1];
    double dt = partial_integrals[i] / (v*v);
    Qres.params().push_back(tlast + dt);
  }

  return Qres;
}

} // namespace srvf
