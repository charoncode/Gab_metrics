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
#ifndef SRVF_EXCEPTIONS_H
#define SRVF_EXCEPTIONS_H 1

#include <stdexcept>

namespace srvf
{

class BadFormat : public std::runtime_error 
{ 
public:
  explicit BadFormat(const std::string &msg)
   : runtime_error(msg)
  { }
};

class AlgorithmFailure : public std::runtime_error 
{ 
public:
  explicit AlgorithmFailure(const std::string &msg)
   : runtime_error(msg)
  { }
};

class UnsupportedOperation : public std::runtime_error 
{ 
public:
  explicit UnsupportedOperation(const std::string &msg)
   : runtime_error(msg)
  { }
};

}

#endif // SRVF_EXCEPTIONS_H
