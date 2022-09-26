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
#include <srvf/partialmatch.h>
#include <srvf/paretoset.h>

#include <srvf/matrix.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/interp.h>
#include <srvf/fileio.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>

#define DEFAULT_GRID_WIDTH 31
#define DEFAULT_GRID_HEIGHT 31


static void do_usage(const char *progname)
{
  std::cout << "USAGE: " << progname << " [OPTIONS] file1 file2" << std::endl;
  std::cout << "where" << std::endl;
  std::cout << "  file1 contains samples for the first curve" << std::endl;
  std::cout << "  file2 contains samples for the second curve" << std::endl;
  std::cout << "Recognized options are:" << std::endl;
  std::cout << "  -A r\t\tautomatically calculate grid size so that " 
               "there is a breakpoint for every rth sample point"
            << std::endl;
  std::cout << "  -W n\t\tuse grid width n" << std::endl;
  std::cout << "  -H n\t\tuse grid height n" << std::endl;
  std::cout << "  -r\t\toptimize over rotations" << std::endl;
  std::cout << "  -w w\t\tuse w as ratio of partiality weight "
               "to shape distance weight in Salukwadze distance " 
               "(default is 1.0)" << std::endl;
  std::cout << "  -k t\t\tuse t as regression error threshold when finding "
               "the knee of the Pareto curve" << std::endl;
  std::cout << "  -p sample point matrices in file1 and file2 " 
               "are in point-per-column ordering" << std::endl;
  std::cout << "  -o file\twrite matches to file" << std::endl;
  std::cout << "  -h\t\tshow this message" << std::endl;
}


// Returns the index of the element of v which is closest to t.
static size_t find_nearest_(double t, const std::vector<double> &v)
{
  int idx = srvf::interp::lookup(v, t);
  if (idx < 0) return 0;

  size_t res = (size_t)idx;
  if (res+1 < v.size() && (fabs(t - v[res]) > fabs(v[res+1] - t)))
  {
    ++res;
  }
  return res;
}


// Program entry point
int main( int argc, char **argv ){
  size_t grid_width = DEFAULT_GRID_WIDTH;
  size_t grid_height = DEFAULT_GRID_HEIGHT;
  bool do_rotations = false;
  double sample_to_grid_ratio = -1.0;
  double salukwadze_ratio = 1.0;
  double regression_thresh = 0.003;
  const char *output_filename = "matches.csv";
  srvf::Pointset::PackingMethod packing = srvf::Pointset::POINT_PER_ROW;
  int opt;

  // Get the command line arguments
  while( (opt=getopt(argc, argv, "A:W:H:rw:k:po:h")) != -1 ){
    switch( opt ){
      case 'A':  // set sample to grid ratio
        sample_to_grid_ratio = atof(optarg);
        break;
      case 'W':  // set matching grid width
        grid_width = atoi(optarg);
        break;
      case 'H':  // set matching grid height
        grid_height = atoi(optarg);
        break;
      case 'r':  // optimize over rotations
        do_rotations = true;
        break;
      case 'w':  // partiality : distance weight ratio for Salukwadze distance
        salukwadze_ratio = atof(optarg);
        break;
      case 'k':  // regression error threshold
        regression_thresh = atof(optarg);
        break;
      case 'p':  // sample point matrix = point per column
        packing = srvf::Pointset::POINT_PER_COLUMN;
        break;
      case 'o':  // set output filename
        output_filename = optarg;
        break;
      case 'h':  // print the help message
        do_usage(argv[0]);
        return 0;
      case '?':  // unrecognized option
        do_usage(argv[0]);
        return -1;
    }
  }
  if (argc - optind != 2) {
    do_usage(argv[0]);
    return -1;
  }

  // Load sample points from file into Matrix instances
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

  
  // Create Pointset's from the Matrix instances
  srvf::Pointset F1samps(F1data[0], packing);
  srvf::Pointset F2samps(F2data[0], packing);

  // Create the curves F1 and F2
  srvf::Plf F1(F1samps);
  srvf::Plf F2(F2samps);


  // Scale F1 and F2 down by a common factor
  double L1 = F1.arc_length();
  double L2 = F2.arc_length();
  double L = std::min(L1, L2);
  F1.scale( 1.0 / L );
  F2.scale( 1.0 / L );


  // Compute the SRVFs of F1 and F2
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);

  // Calculate matching grid dimensions
  if (sample_to_grid_ratio > 0.0)
  {
    grid_width = F1.ncp() / sample_to_grid_ratio;
    grid_height = F2.ncp() / sample_to_grid_ratio;
  }
  std::cout << "Using grid dimensions " << grid_width << " x " 
            << grid_height << std::endl;

  // Find the Pareto-optimal partial matches
  srvf::pmatch::ParetoSet S = srvf::pmatch::find_matches (
    Q1, Q2, do_rotations, grid_width, grid_height);

  // Write the results to file
  std::ofstream ofs(output_filename);
  if (!ofs)
  {
    std::cerr << "Failed to open " << output_filename << " for writing."
              << std::endl;
    return -1;
  }

  ofs << "# srvf_pmatch output" << std::endl;
  ofs << "# File 1: " << argv[optind] << std::endl;
  ofs << "# File 2: " << argv[optind+1] << std::endl;
  ofs << "# Grid dimensions: " << grid_width << " x " 
      << grid_height << std::endl;
  ofs << "# Salukwadze dist: " 
      << S.salukwadze_dist(salukwadze_ratio) 
      << std::endl << std::endl;

  size_t knee_idx = S.find_knee(regression_thresh);
  ofs << "# Best match at knee: "
      << find_nearest_(S[knee_idx][0].a, F1.params()) << " "
      << find_nearest_(S[knee_idx][0].b, F1.params()) << " "
      << find_nearest_(S[knee_idx][0].c, F1.params()) << " "
      << find_nearest_(S[knee_idx][0].d, F1.params()) << " "
      << S[knee_idx][0].dist << std::endl;
  ofs << "# Partiality of knee: "
      << (2.0 - S[knee_idx][0].length()) << std::endl << std::endl;

  ofs << "# name: salukwadze_dist" << std::endl;
  ofs << "# type: matrix" << std::endl;
  ofs << "# rows: 1" << std::endl;
  ofs << "# columns: 1" << std::endl;
  ofs << S.salukwadze_dist(salukwadze_ratio) 
      << std::endl << std::endl;

  ofs << "# name: pareto_set" << std::endl;
  ofs << "# type: matrix" << std::endl;
  ofs << "# rows: " << S.total_size() << std::endl;
  ofs << "# columns: 5" << std::endl;
  for (size_t i=0; i<S.nbuckets(); ++i)
  {
    for (size_t j=0; j<S[i].size(); ++j)
    {
      size_t ai = find_nearest_(S[i][j].a, F1.params());
      size_t bi = find_nearest_(S[i][j].b, F1.params());
      size_t ci = find_nearest_(S[i][j].c, F2.params());
      size_t di = find_nearest_(S[i][j].d, F2.params());

      ofs << ai << " " << bi << " " << ci << " " << di << " "
          << S[i][j].dist << std::endl;
    }
  }
  ofs.close();

  return 0;
}

