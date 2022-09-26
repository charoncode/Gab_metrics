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
#ifndef MATCHVIEW_H
#define MATCHVIEW_H 1

#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/plot.h>
#include <srvf/rotate.h>
#include <srvf/paretoset.h>

#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <GL/gl.h>
#include <GL/glu.h>

#include <limits>


class MatchView : public Fl_Gl_Window
{
public:

  MatchView(int x, int y, int w, int h)
   : Fl_Gl_Window(x, y, w, h),
     prev_x_(0), prev_y_(0), button_state_(0),
     do_rotations_(false)
  { }

  MatchView(int x, int y, int w, int h, const char *l)
   : Fl_Gl_Window(x, y, w, h, l),
     prev_x_(0), prev_y_(0), button_state_(0),
     do_rotations_(false)
  { }

  void set_plot(srvf::plot::Plot *p)
  { 
    plot_ = p; 
    redraw();
  }

  void set_matches(const srvf::pmatch::ParetoSet &S)
  {
    matches_ = S;
  }

  void set_interval(size_t idx, double a, double b)
  {
    plot_->set_plf_interval(idx, std::pair<double,double>(a,b));

    // Changing the domain also changes the bounding box and centroid, 
    // so the curve needs to be re-centered
    srvf::Plf &F = plot_->get_plf(idx);
    srvf::Point X(F.dim(), 0.0);
    double npts=0.0;
    for (size_t i=0; i<F.ncp(); ++i)
    {
      if (F.params()[i] >= a && F.params()[i] <= b)
      {
        X += F.samps()[i];
        ++npts;
      }
    }
    if (npts > 0.0) X *= (-1.0 / npts);
    F.translate( X );
    redraw();
  }

  void set_match(size_t bucket_idx, size_t match_idx)
  {
    if ( (bucket_idx < matches_.nbuckets()) && 
         (match_idx < matches_[bucket_idx].size()) )
    {
      double a = matches_[bucket_idx][match_idx].a;
      double b = matches_[bucket_idx][match_idx].b;
      double c = matches_[bucket_idx][match_idx].c;
      double d = matches_[bucket_idx][match_idx].d;

      set_interval (0, a, b);
      set_interval (1, c, d);

      if (do_rotations_)
      {
        srvf::Srvf Q1 = srvf::plf_to_srvf(plot_->get_plf(0));
        srvf::Srvf Q2 = srvf::plf_to_srvf(plot_->get_plf(1));
        srvf::Matrix R = srvf::optimal_rotation(Q1, Q2, a, b, c, d);
        plot_->get_plf(1).rotate(R);
      }
      redraw();
    }
  }

  void show_f1(bool v)
  { plot_->set_plf_visibility(0,v); redraw(); }

  void show_f2(bool v)
  { plot_->set_plf_visibility(1,v); redraw(); }

  void do_rotations(bool v)
  { do_rotations_ = v; }

  size_t nbuckets()
  { return matches_.nbuckets(); }

  size_t bucket_size(size_t bucket_idx)
  { return matches_[bucket_idx].size(); }

  double match_dist(size_t bucket_idx, size_t match_idx)
  { 
    if (bucket_idx < matches_.nbuckets() && 
        match_idx < matches_[bucket_idx].size())
      return matches_[bucket_idx][match_idx].dist; 
    else
      return std::numeric_limits<double>::infinity();
  }

  double match_length(size_t bucket_idx, size_t match_idx)
  { 
    if (bucket_idx < matches_.nbuckets() && 
        match_idx < matches_[bucket_idx].size())
      return matches_[bucket_idx][match_idx].length(); 
    else
      return -std::numeric_limits<double>::infinity();
  }

  virtual void draw()
  {
    renderer_.device_width(w());
    renderer_.device_height(h());
    renderer_.viewport(0, 0, w(), h());
    
    renderer_.clear_color(srvf::plot::Color(1.0,1.0,1.0));
    renderer_.clear();

    glPushMatrix();
    plot_->render(renderer_);
    glPopMatrix();
  }


  int handle(int event)
  {
    double dx, dy;

    switch(event) {
      case FL_PUSH:
        prev_x_ = Fl::event_x();
        prev_y_ = Fl::event_y();
        if      ( Fl::event_button() == 1 ) button_state_ |= 1; // left down
        else if ( Fl::event_button() == 2 ) button_state_ |= 2; // middle down
        else if ( Fl::event_button() == 3 ) button_state_ |= 4; // right down
        return 1;
      case FL_DRAG:
        dx = (double)(Fl::event_x() - prev_x_) / (double)w();
        dy = (double)(Fl::event_y() - prev_y_) / (double)h();
        prev_x_ = Fl::event_x();
        prev_y_ = Fl::event_y();

        switch(button_state_)
        {
          case 1: // left button down only
            plot_->pan_view(dx, dy);
            break;
          case 2: // middle button down only
            break;
          case 4: // right button down only
            plot_->rotate_view(dx, dy);
            break;
          case 5: // left and right buttons down
            plot_->scale_view(dx, dy);
            break;
        }
        redraw();
        return 1;
      case FL_RELEASE:   
        if      ( Fl::event_button() == 1 ) button_state_ &= 6; // left up
        else if ( Fl::event_button() == 2 ) button_state_ &= 5; // middle up
        else if ( Fl::event_button() == 3 ) button_state_ &= 3; // right up
        return 1;
      case FL_FOCUS :
      case FL_UNFOCUS :
        return 1;  // we want keyboard events
      case FL_KEYDOWN:
        return 0;
      case FL_SHORTCUT:
        return 0;
      default:
        return Fl_Gl_Window::handle(event);
    }
  }

private:  
  srvf::plot::Plot *plot_;
  srvf::plot::OpenGlRenderer renderer_;
  srvf::pmatch::ParetoSet matches_;

  // Mouse state
  int prev_x_;
  int prev_y_;
  int button_state_;

  bool do_rotations_;
};

#endif // MATCHVIEW_H
