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
#include "matchview.h"
#include "ui.h"

#include <srvf/partialmatch.h>
#include <srvf/paretoset.h>

#include <srvf/matrix.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/fileio.h>

#include <FL/Fl.H>

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

#define DEFAULT_GRID_WIDTH 0
#define DEFAULT_GRID_HEIGHT 0


static void do_usage(const char *progname)
{
  std::cout << "USAGE: " << progname << " [OPTIONS] file1 file2" << std::endl;
  std::cout << "where" << std::endl;
  std::cout << "  file1 contains samples for the first curve" << std::endl;
  std::cout << "  file2 contains samples for the second curve" << std::endl;
  std::cout << "Recognized options are:" << std::endl;
  std::cout << "  -W n\t\tuse grid width n" << std::endl;
  std::cout << "  -H n\t\tuse grid height n" << std::endl;
  std::cout << "  -r\t\toptimize over rotations" << std::endl;
  std::cout << "  -p\t\tsample point matrices in file1 and file2 " 
               "are in point-per-column ordering" << std::endl;
  std::cout << "  -o outfile\twrite matches to outfile" << std::endl;
  std::cout << "  -h\t\tshow this message" << std::endl;
}

int main( int argc, char **argv ){
  size_t grid_width = DEFAULT_GRID_WIDTH;
  size_t grid_height = DEFAULT_GRID_HEIGHT;
  bool do_rotations = false;
  double salukwadze_weight = 1.0;
  const char *output_filename = "matches.mat";
  srvf::Pointset::PackingMethod packing = srvf::Pointset::POINT_PER_ROW;
  int opt;

  while( (opt=getopt(argc, argv, "W:H:rw:po:h")) != -1 ){
    switch( opt ){
      case 'W':
        grid_width = atoi(optarg);
        break;
      case 'H':
        grid_height = atoi(optarg);
        break;
      case 'r':
        do_rotations = true;
        break;
      case 'w':  // shape distance weight for Salukwadze distance
        salukwadze_weight = atof(optarg);
        break;
      case 'p':  // sample point matrix = point per column
        packing = srvf::Pointset::POINT_PER_COLUMN;
        break;
      case 'o':
        output_filename = optarg;
        break;
      case 'h':
        do_usage(argv[0]);
        return 0;
      case '?':
        do_usage(argv[0]);
        return -1;
    }
  }
  if (argc - optind != 2) {
    do_usage(argv[0]);
    return -1;
  }

  std::ifstream ifs1(argv[optind]);
  std::ifstream ifs2(argv[optind+1]);

  std::vector<srvf::Matrix> F1data = srvf::io::load_csv ( ifs1, ' ', '\n' );
  std::vector<srvf::Matrix> F2data = srvf::io::load_csv ( ifs2, ' ', '\n' );

  ifs1.close();
  ifs2.close();

  if (F1data.size() == 0)
  {
    std::cerr << "Failed to load matrix from " << argv[optind]
              << "; exiting." << std::endl;
    return -1;
  }
  if (F2data.size() == 0)
  {
    std::cerr << "Failed to load matrix from " << argv[optind+1]
              << "; exiting." << std::endl;
    return -1;
  }

  srvf::Pointset F1samps(F1data[0], packing);
  srvf::Pointset F2samps(F2data[0], packing);

  srvf::Plf F1(F1samps);
  srvf::Plf F2(F2samps);

  std::cout << "F1 has " << F1.ncp() << " changepoints" << std::endl;
  std::cout << "F2 has " << F2.ncp() << " changepoints" << std::endl;

  if (F1.dim() != 2 && F1.dim() != 3)
  {
    std::cerr << "Invalid dimension: curves must have dimension 2 or 3." 
              << std::endl;
    return -1;
  }
  if (F1.dim() != F2.dim())
  {
    std::cerr << "F1 and F2 have different dimensions." << std::endl;
    return -1;
  }

  double L1 = F1.arc_length();
  double L2 = F2.arc_length();
  double L = std::min(L1, L2);
  F1.scale( 1.0 / L );
  F2.scale( 1.0 / L );

  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);

  srvf::pmatch::ParetoSet S = srvf::pmatch::find_matches (
    Q1, Q2, do_rotations, grid_width, grid_height);

  std::vector<double> reg_errs = S.regression_errors();
  for (size_t i=0; i<reg_errs.size(); ++i)
  {
    if (S[i].empty()) continue;

    std::cout << S[i][0].length() << " "
              << reg_errs[i] << std::endl;
  }

  std::ofstream ofs(output_filename);
  ofs << "# srvf_pmatch output" << std::endl;
  ofs << "# File 1: " << argv[optind] << std::endl;
  ofs << "# File 2: " << argv[optind+1] << std::endl;
  ofs << "# Salukwadze dist: " 
      << S.salukwadze_dist(salukwadze_weight) 
      << std::endl << std::endl;

  ofs << "# name: salukwadze_dist" << std::endl;
  ofs << "# type: matrix" << std::endl;
  ofs << "# rows: 1" << std::endl;
  ofs << "# columns: 1" << std::endl;
  ofs << S.salukwadze_dist(salukwadze_weight) 
      << std::endl << std::endl;

  ofs << "# name: pareto_set" << std::endl;
  ofs << "# type: matrix" << std::endl;
  ofs << "# rows: " << S.total_size() << std::endl;
  ofs << "# columns: 5" << std::endl;
  for (size_t i=0; i<S.nbuckets(); ++i)
  {
    for (size_t j=0; j<S[i].size(); ++j)
    {
      ofs << S[i][j].a << " " << S[i][j].b << " " 
          << S[i][j].c << " " << S[i][j].d << " "
          << S[i][j].dist << std::endl;
    }
  }
  ofs.close();

  srvf::plot::Plot *plot = NULL;
  if (F1.dim() == 2)
    plot = new srvf::plot::Plot2D();
  else // dim = 3
    plot = new srvf::plot::Plot3D();

  plot->insert(F1, srvf::plot::Color(0.0, 0.0, 1.0), 2.0);
  plot->insert(F2, srvf::plot::Color(1.0, 0.0, 0.0), 2.0);

  UserInterface ui;

  Fl_Window *win = ui.make_window();
  ui.match_view->do_rotations(do_rotations);
  ui.match_view->set_plot(plot);
  ui.match_view->set_matches(S);
  ui.slider_match_length->range(0.0, (double)(S.nbuckets()-1));
  ui.slider_match_length->step(1.0);
  ui.slider_match_length->value((double)(S.nbuckets()-1));
  ui.set_match();

  win->show();
  Fl::run();

  delete plot;
  return 0;
}

