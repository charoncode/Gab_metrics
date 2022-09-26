/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  Daniel Robinson
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
#include <srvf/matrix.h>
#include <srvf/plf.h>
#include <srvf/srvf.h>
#include <srvf/qmap.h>
#include <srvf/opencurves.h>
#include <srvf/fileio.h>

#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>


static void do_usage(const char *progname)
{
  std::cout << "USAGE: " << progname << " [OPTIONS] file1 file2" << std::endl;
  std::cout << "where" << std::endl;
  std::cout << "  file1 contains samples for the first curve" << std::endl;
  std::cout << "  file2 contains samples for the second curve" << std::endl;
  std::cout << "Recognized options are:" << std::endl;
  std::cout << "  -W n\tuse n for DP grid width" << std::endl;
  std::cout << "  -H n\tuse n for DP grid height" << std::endl;
  std::cout << "  -r\t\toptimize over rotations" << std::endl;
  std::cout << "  -s\t\trescale SRVFs to unit norm" << std::endl;
  std::cout << "  -e\t\telastic distance (optimize over reparametrizations)" 
            << std::endl;
  std::cout << "  -p sample point matrices in file1 and file2 " 
               "are in point-per-column ordering" << std::endl;
  std::cout << "  -h\t\tshow this message" << std::endl;
}

// Program entry point
int main( int argc, char **argv ){
  bool do_rotations = false;
  bool do_reparams = false;
  bool do_rescale = false;
  size_t grid_width = 0;
  size_t grid_height = 0;
  srvf::Pointset::PackingMethod packing = srvf::Pointset::POINT_PER_ROW;
  int opt;

  // Get the command line arguments
  while( (opt=getopt(argc, argv, "W:H:rseph")) != -1 ){
    switch( opt ){
      case 'W':  // DP grid width
        grid_width = (size_t)atoi(optarg);
        break;
      case 'H':  // DP grid height
        grid_height = (size_t)atoi(optarg);
        break;
      case 'r':  // optimize over rotations
        do_rotations = true;
        break;
      case 's':  // rescale to unit norm
        do_rescale = true;
        break;
      case 'e':  // optimize over reparametrizations
        do_reparams = true;
        break;
      case 'p':  // sample point matrix = point per column
        packing = srvf::Pointset::POINT_PER_COLUMN;
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
  if (F1samps.dim() != F2samps.dim())
  {
    std::cerr << "Sample point sets in " << argv[optind] 
              << " and " << argv[optind+1]
              << " have different dimensions; exiting."
              << std::endl;
    return -1;
  }

  // Create the curves F1 and F2
  srvf::Plf F1(F1samps);
  srvf::Plf F2(F2samps);

  // Compute the SRVFs
  srvf::Srvf Q1 = srvf::plf_to_srvf(F1);
  srvf::Srvf Q2 = srvf::plf_to_srvf(F2);

  double shape_dist = srvf::opencurves::shape_distance (
    Q1, Q2, do_rotations, do_rescale, do_reparams, 1, grid_width, grid_height);
  std::cout << "Shape distance = " << shape_dist << std::endl;

  return 0;
}
