/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012 FSU Statistical Shape Analysis and Modeling Group
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
#ifndef SRVF_PLOT_H
#define SRVF_PLOT_H 1

#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/matrix.h>
#include "render.h"

namespace srvf
{
namespace plot
{

/**
 * Abstract base for plot classes.
 *
 * A \c Plot object is responsible for managing a collection of \c Plf 
 * and \c Srvf instances, and for defining the geometry of all of 
 * these functions.  \c Plot objects do not handle viewing.
 */
class Plot
{ 
public:
  
  Plot()
   : plot_radius_(1.0)
  { }
  
  virtual void 
  render(Renderer &r) = 0;

  virtual void 
  pan_view(double dx, double dy) { }

  virtual void 
  scale_view(double dx, double dy) { }

  virtual void
  rotate_view(double dx, double dy) { }

  virtual void 
  insert(const Plf &F, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    plfs_.push_back(F);
    plf_visibilities_.push_back(true);
    plf_colors_.push_back(c);
    plf_thicknesses_.push_back(thickness);
    plf_modes_.push_back(mode);
    plf_intervals_.push_back (
      std::pair<double,double>(F.domain_lb(), F.domain_ub()) );
    bounding_boxes_.push_back(F.bounding_box());
    for (size_t i=0; i<F.samps().npts(); ++i)
    {
      plot_radius_ = std::max(plot_radius_, F.samps()[i].norm());
    }
  }

  virtual void 
  insert(const Srvf &Q, Color c, double thickness=1.0, DrawingMode mode=LINES)
  {
    srvfs_.push_back(Q);
    srvf_colors_.push_back(c);
    srvf_thicknesses_.push_back(thickness);
    srvf_modes_.push_back(mode);
  }

  Plf &get_plf(size_t idx)
  { return plfs_[idx]; }

  virtual void set_plf_visibility(size_t idx, bool v)
  { plf_visibilities_[idx] = v; }
  virtual bool get_plf_visibility(size_t idx)
  { return plf_visibilities_[idx]; }

  virtual void set_plf_color(size_t idx, Color c)
  { plf_colors_[idx] = c; }
  virtual Color get_plf_color(size_t idx)
  { return plf_colors_[idx]; }

  virtual void set_plf_thickness(size_t idx, double v)
  { plf_thicknesses_[idx] = v; }
  virtual double get_plf_thickness(size_t idx)
  { return plf_thicknesses_[idx]; }

  virtual void set_plf_mode(size_t idx, DrawingMode m)
  { plf_modes_[idx] = m; }
  virtual DrawingMode get_plf_mode(size_t idx)
  { return plf_modes_[idx]; }

  virtual void set_plf_interval(size_t idx, std::pair<double,double> ivl)
  { plf_intervals_[idx] = ivl; }
  virtual std::pair<double,double> get_plf_interval(size_t idx)
  { return plf_intervals_[idx]; }

protected:
  
  std::vector<Plf> plfs_;
  std::vector<Srvf> srvfs_;
  std::vector<bool> plf_visibilities_;
  std::vector<Color> plf_colors_;
  std::vector<Color> srvf_colors_;
  std::vector<double> plf_thicknesses_;
  std::vector<double> srvf_thicknesses_;
  std::vector<DrawingMode> plf_modes_;
  std::vector<DrawingMode> srvf_modes_;
  std::vector<std::pair<double,double> > plf_intervals_;
  std::vector<std::vector<Point> > bounding_boxes_;
  double plot_radius_;
};


class Plot2D : public Plot
{ 
public:

  Plot2D() 
   : trans_x_(0.0), trans_y_(0.0), rot_theta_(0.0), scale_(1.0)
  { }

  virtual void render(Renderer &r);

  virtual void 
  pan_view(double dx, double dy);

  virtual void 
  scale_view(double dx, double dy);

  virtual void
  rotate_view(double dx, double dy);

private:
  
  double trans_x_;
  double trans_y_;
  double rot_theta_;
  double scale_;
};


class Plot3D : public Plot
{ 
public:

  Plot3D() 
   : hrot_(0.0), vrot_(0.0), scale_(1.0)
  { }

  virtual void render(Renderer &r);

  virtual void 
  scale_view(double dx, double dy);

  virtual void
  rotate_view(double dx, double dy);

private:
  
  double hrot_;
  double vrot_;
  double scale_;
};


/**
 * A specialized plot for 1-D functions.
 */
class FunctionPlot : public Plot
{ 
public:

  virtual void render(Renderer &r);
};


} // namespace plot
} // namespace srvf

#endif // SRVF_PLOT_H
