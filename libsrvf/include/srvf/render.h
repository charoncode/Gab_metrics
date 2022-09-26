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
#ifndef SRVF_RENDER_H
#define SRVF_RENDER_H 1

#include <srvf/matrix.h>
#include <cstdlib>

namespace srvf
{
namespace plot
{


struct Color
{
  Color(float r, float g, float b)
   : red(r), green(g), blue(b)
  { }

  float red;
  float green;
  float blue;
};


enum DrawingMode
{
  POINTS=0,
  LINES=1,
  TUBES=2
};
  

/**
 * Abstract base class for renderers.
 */
class Renderer
{
public:

  virtual size_t device_width() = 0;
  virtual void device_width(size_t w) = 0;

  virtual size_t device_height() = 0;
  virtual void device_height(size_t h) = 0;

  virtual void viewport(double x, double y, double w, double h) = 0;
  virtual void ortho(double left, double right, 
                       double bottom, double top, 
                       double near, double far) = 0;
  virtual void frustum(double left, double right, 
                       double bottom, double top, 
                       double near, double far) = 0;

  virtual void clear_color(Color c) = 0;
  virtual void clear() = 0;
  virtual void begin(DrawingMode mode) = 0;
  virtual void vertex(double x, double y) = 0;
  virtual void vertex(double x, double y, double z) = 0;
  virtual void end() = 0;
  
  virtual void set_color(Color c) = 0;
  virtual void set_thickness(double t) = 0;

  virtual void translate(double x, double y) = 0;
  virtual void translate(double x, double y, double z) = 0;
  virtual void scale(double sfx, double sfy) = 0;
  virtual void scale(double sfx, double sfy, double sfz) = 0;
  virtual void rotate(double angle) = 0;
  virtual void rotate(double angle, double ax, double ay, double az) = 0;

protected:
  
  Renderer() { }
};

/**
 * OpenGL renderer.
 */
class OpenGlRenderer : public Renderer
{
public:
  
  OpenGlRenderer();

  virtual size_t device_width(){ return dev_width_; }
  virtual void device_width(size_t w){ dev_width_ = w; }

  virtual size_t device_height(){ return dev_height_; }
  virtual void device_height(size_t h){ dev_height_ = h; }

  virtual void viewport(double x, double y, double w, double h);
  virtual void ortho(double left, double right, 
                       double bottom, double top, 
                       double near, double far);
  virtual void frustum(double left, double right, 
                       double bottom, double top, 
                       double near, double far);

  virtual void clear_color(Color c);
  virtual void clear();
  virtual void begin(DrawingMode mode);
  virtual void vertex(double x, double y);
  virtual void vertex(double x, double y, double z);
  virtual void end();
  
  virtual void set_color(Color c);
  virtual void set_thickness(double t);

  virtual void translate(double x, double y);
  virtual void translate(double x, double y, double z);
  virtual void scale(double sfx, double sfy);
  virtual void scale(double sfx, double sfy, double sfz);
  virtual void rotate(double angle);
  virtual void rotate(double angle, double ax, double ay, double az);

private:
  DrawingMode mode_;
  double thickness_;
  size_t dev_width_;
  size_t dev_height_;
};


} // namespace plot
} // namespace srvf

#endif // SRVF_RENDER_H
