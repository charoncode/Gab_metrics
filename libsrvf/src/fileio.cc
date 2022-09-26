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
#include <srvf/fileio.h>
#include <srvf/matrix.h>
#include <srvf/exceptions.h>

#include <cstddef>
#include <cctype>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>


// Returns the first non-space character in the line.  
// If there are no non-space characters in the line, returns 0.
static char first_nonspace_char_(const std::string &s)
{
  for (size_t i=0; i<s.size(); ++i)
    if (!isspace(s[i])) return s[i];

  return 0;
}


namespace srvf
{
namespace io
{

/**
 * Load a collection of sample sets from a CSV file.
 *
 * Any line beginning with a non-numeric character is ignored.
 */
std::vector<Matrix> load_csv (std::istream &is, char fieldsep, char linesep)
{
  std::vector<Matrix> res;
  std::string line;
  std::vector<double> entries;

  size_t this_matrix_rows=0;
  size_t this_matrix_cols=0;
  while (std::getline(is, line, linesep))
  {
    char testchar = first_nonspace_char_(line);
    if (testchar == '.' || testchar == '-' || isdigit(testchar))
    {
      size_t this_row_cols=0;
      double vi;
      std::istringstream iss(line);
      while (iss >> vi)
      {
        entries.push_back(vi);
        ++this_row_cols;
        iss.ignore(1,fieldsep);
      }

      if (this_matrix_rows > 0)
      {
        if (this_row_cols != this_matrix_cols)
        {
          throw BadFormat("inconsistent row length");
        }
      }
      else
      {
        this_matrix_cols=this_row_cols;
      }
      ++this_matrix_rows;
      line.clear();
    }
    else
    {
      if (this_matrix_rows>0)
      {
        res.push_back(Matrix(this_matrix_rows,this_matrix_cols,entries));
        entries.clear();
        this_matrix_rows=0;
        this_matrix_cols=0;
      }
    }
  }

  if (this_matrix_rows>0)
  {
    res.push_back(Matrix(this_matrix_rows,this_matrix_cols,entries));
  }
  return res;
}

/**
 * Save a collection of sample sets to a CSV file.
 */
void save_csv (std::ostream &os, 
               const std::vector<Matrix> &data,
               char fieldsep, char linesep)
{
  for (size_t i=0; i<data.size(); ++i)
  {
    for (size_t r=0; r<data[i].rows(); ++r)
    {
      for (size_t c=0; c<data[i].cols(); ++c)
      {
        os << data[i](r,c);
        if (c<data[i].cols()-1)
        {
          os << fieldsep;
        }
      }
      os << linesep;
    }
    os << linesep;
  }
}

} // namespace srvf::io
} // namespace srvf
