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
#ifndef MINMAXHEAP_H
#define MINMAXHEAP_H 1

#include <cstddef>
#include <vector>
#include <iostream>
#include <stdexcept>

namespace srvf{

/**
 * A min-max heap.
 * 
 * An implementation of a min-max heap, as described by Atkinson et al. 
 * in "Min-Max Heaps and Generalized Priority Queues", Communications of 
 * the ACM, Vol. 29 No. 10, October 1986.
 *
 * The min-max heap supports the following operations:
 * <ul>
 *   <li> Create (\f$ O(n) \f$ operations)
 *   <li> Insert (\f$ O(\log(n)) \f$ operations)
 *   <li> Delete Min (\f$ O(\log(n)) \f$ operations)
 *   <li> Delete Max (\f$ O(\log(n)) \f$ operations)
 * </ul>
 *
 * This is implemented as a class template.  The template parameter \c T 
 * must be a type for which \c operator< is defined.
 */
template <class T>
class MinMaxHeap
{
public:

  /**
   * Creates an empty heap.
   */
  MinMaxHeap()
   : data_()
  { }

  /**
   * Creates a heap containing the given elements.
   */
  MinMaxHeap(const std::vector<T> &data)
   : data_(data)
  { minmax_heapify_(); }

  /**
   * Returns the number of elements in the heap.
   */
  inline size_t size() const
  {
    return data_.size();
  }

  /**
   * Indicates whether or not the heap is empty.
   */
  inline bool empty() const
  {
    return data_.empty();
  }

  /**
   * Returns a reference to the maximum element in the heap.
   */
  inline T &max()
  {
    if (empty()) throw std::logic_error("max() called on an empty heap");

    return data_[get_max_index_()];
  }

  /**
   * Returns a \c const reference to the maximum element in the heap.
   */
  inline const T &max() const
  {
    if (empty()) throw std::logic_error("max() called on an empty heap");

    return data_[get_max_index_()];
  }

  /**
   * Returns a reference to the minimum element in the heap.
   */
  inline T &min()
  {
    if (empty()) throw std::logic_error("min() called on an empty heap");
    
    return data_[0];
  }

  /**
   * Returns a \c const reference to the minimum element in the heap.
   */
  inline const T &min() const
  {
    if (empty()) throw std::logic_error("min() called on an empty heap");
    
    return data_[0];
  }

  /**
   * Inserts the given element into the heap.
   */
  void insert(const T &X)
  {
    data_.push_back(X);
    bubble_up_(data_.size()-1);
  }

  /**
   * Removes all elements from the heap.
   */
  void clear()
  {
    data_.clear();
  }

  /**
   * Removes the maximum element from the heap.
   */
  void remove_max()
  {
    if (empty()) throw std::logic_error("remove_max() called on an empty heap");
    
    if (data_.size() > 1)
    {
      size_t max_idx = get_max_index_();
      if (max_idx + 1 < data_.size()) 
        data_[max_idx] = data_.back();
      data_.pop_back();
      trickle_down_(max_idx);
    }
    else 
    {
      data_.clear();
    }
  }

  /**
   * Removes the minimum element from the heap.
   */
  void remove_min()
  {
    if (empty()) throw std::logic_error("remove_min() called on an empty heap");
    
    if (data_.size() > 1)
    {
      data_[0] = data_.back();
      data_.pop_back();
      trickle_down_(0);
    }
    else 
    {
      data_.clear();
    }
  }

  /**
   * Returns a \c const reference to the \c vector containing the elements 
   * of the heap.
   */
  inline const std::vector<T> &data()
  {
    return data_;
  }

  /**
   * Returns a reference to an element of the heap by its index.
   */
  inline T &operator[] (size_t idx)
  {
    if (idx >= data_.size()) 
      throw std::out_of_range("Index idx is out of bounds.");

    return data_[idx];
  }

  /**
   * Returns a \c const reference to an element of the heap by its index.
   */
  inline const T &operator[] (size_t idx) const
  {
    if (idx >= data_.size()) 
      throw std::out_of_range("Index idx is out of bounds.");

    return data_[idx];
  }

  /**
   * Print a text representation of the heap to stdout.
   */
  void print()
  {
    for (size_t i=0; i<data_.size(); ++i)
    {
      std::cout << data_[i];
      if (i+1 >= data_.size()) break;

      if (level_(i+1) == level_(i))
        std::cout << " ";
      else
        std::cout << std::endl;
    }
    std::cout << std::endl;
  }

private:
  std::vector<T> data_;

  /**
   * Indicates whether the given index is the index of an 
   * existing node in the heap.
   */
  inline bool exists_(size_t i)
  { return (i < data_.size()); }

  /**
   * Returns an invalid node index.
   */
  inline size_t nonexistant_node_()
  { return data_.size(); }

  /**
   * Return the level of the node with index \a i.
   */
  inline size_t level_(size_t i)
  { 
    size_t res=0;
    if (i>0)
    {
      ++i;
      while (i){
        i>>=1;
        ++res;
      }
      return res-1;
    }
    return res;
  }

  /**
   * Returns the index of the largest element in the heap.
   */
  inline size_t get_max_index_()
  {
    int idx = (data_.size() == 1 ? 0 : 1);
    if (data_.size() > 2 && data_[idx] < data_[2]) idx = 2;
    return idx;
  }

  /**
   * Indicates whether the given node is on a min level.
   */
  inline bool is_on_min_level_(size_t i)
  { return (level_(i) % 2 == 0); }

  /**
   * Indicates whether the given node is on a max level.
   */
  inline bool is_on_max_level_(size_t i)
  { return !is_on_min_level_(i); }

  /**
   * Return the index of the parent of the node with index \a i.
   */
  inline size_t parent_(size_t i)
  { return (i>0 ? (i-1)/2 : 0); }

  /**
   * Return the index of the left child of the node with index \a i.
   */
  inline size_t left_child_(size_t i)
  { return 2*i + 1; }

  /**
   * Return the index of the right child of the node with index \a i.
   */
  inline size_t right_child_(size_t i)
  { return 2*i + 2; }

  /**
   * Return the index of the smallest child of the given node.
   */
  size_t find_smallest_child_(size_t i)
  {
    size_t res = nonexistant_node_();
    size_t c;

    if (exists_(c=left_child_(i)))
      if (!exists_(res) || data_[c] < data_[res]) res = c;
    if (exists_(c=right_child_(i)))
      if (!exists_(res) || data_[c] < data_[res]) res = c;

    return res;
  }

  /**
   * Return the index of the largest child of the given node.
   */
  size_t find_largest_child_(size_t i)
  {
    size_t res = nonexistant_node_();
    size_t c;

    if (exists_(c=left_child_(i)))
      if (!exists_(res) || data_[res] < data_[c]) res = c;
    if (exists_(c=right_child_(i)))
      if (!exists_(res) || data_[res] < data_[c]) res = c;

    return res;
  }

  /**
   * Return the index of the smallest of the given nodes children and 
   * grandchildren.
   *
   * If the given node has no children, returns an invalid index.
   */
  size_t find_smallest_child_or_grandchild_(size_t i)
  {
    size_t res = nonexistant_node_();
    size_t c, cc;

    if ( exists_(c=left_child_(i)) )
    {
      if ( !exists_(res) || data_[c] < data_[res]) res = c;
      if ( exists_(cc=left_child_(c)) )
        if ( !exists_(res) || data_[cc] < data_[res] ) res = cc;
      if ( exists_(cc=right_child_(c)) )
        if ( !exists_(res) || data_[cc] < data_[res] ) res = cc;
    }
    if ( exists_(c=right_child_(i)) )
    {
      if ( !exists_(res) || data_[c] < data_[res]) res = c;
      if ( exists_(cc=left_child_(c)) )
        if ( !exists_(res) || data_[cc] < data_[res] ) res = cc;
      if ( exists_(cc=right_child_(c)) )
        if ( !exists_(res) || data_[cc] < data_[res] ) res = cc;
    }

    return res;
  }

  /**
   * Return the index of the largest of the given nodes children and 
   * grandchildren.
   *
   * If the given node has no children, returns an invalid index.
   */
  size_t find_largest_child_or_grandchild_(size_t i)
  {
    size_t res = nonexistant_node_();
    size_t c, cc;

    if ( exists_(c=left_child_(i)) )
    {
      if ( !exists_(res) || data_[res] < data_[c]) res = c;
      if ( exists_(cc=left_child_(c)) )
        if ( !exists_(res) || data_[res] < data_[cc] ) res = cc;
      if ( exists_(cc=right_child_(c)) )
        if ( !exists_(res) || data_[res] < data_[cc] ) res = cc;
    }
    if ( exists_(c=right_child_(i)) )
    {
      if ( !exists_(res) || data_[res] < data_[c]) res = c;
      if ( exists_(cc=left_child_(c)) )
        if ( !exists_(res) || data_[res] < data_[cc] ) res = cc;
      if ( exists_(cc=right_child_(c)) )
        if ( !exists_(res) || data_[res] < data_[cc] ) res = cc;
    }

    return res;
  }

  inline void swap_(size_t i, size_t j)
  {
    //std::cout << "swap(" << i << "," << j << ")" << std::endl;
    T tmp = data_[i];
    data_[i] = data_[j];
    data_[j] = tmp;
  }

  void trickle_down_(size_t i)
  {
    if (is_on_min_level_(i))
      trickle_down_min_(i);
    else
      trickle_down_max_(i);
  }

  void trickle_down_min_(size_t i)
  {
    size_t m = find_smallest_child_or_grandchild_(i);
    if ( exists_(m) )
    {
      if (m >= 4*i+3)
      {
        // m is a grandchild of i
        if (data_[m] < data_[i])
        {
          swap_(m, i);
          size_t p = parent_(m);
          if (data_[p] < data_[m]) swap_(m,p);
          trickle_down_min_(m);
        }
      }
      else
      {
        // m is a child of i
        if (data_[m] < data_[i]) swap_(m,i);
      }
    }
  }

  void trickle_down_max_(size_t i)
  {
    size_t m = find_largest_child_or_grandchild_(i);
    if ( exists_(m) )
    {
      if (m >= 4*i+3)
      {
        // m is a grandchild of i
        if (data_[i] < data_[m])
        {
          swap_(m, i);
          size_t p = parent_(m);
          if (data_[m] < data_[p]) swap_(m,p);
          trickle_down_max_(m);
        }
      }
      else
      {
        // m is a child of i
        if (data_[i] < data_[m]) swap_(m,i);
      }
    }
  }

  void bubble_up_(size_t i)
  {
    if (is_on_min_level_(i))
    {
      size_t p = parent_(i);
      if (exists_(p) && data_[p] < data_[i])
      {
        swap_(i,p);
        bubble_up_max_(p);
      }
      else bubble_up_min_(i);
    }
    else
    {
      size_t p = parent_(i);
      if (exists_(p) && data_[i] < data_[p])
      {
        swap_(i,p);
        bubble_up_min_(p);
      }
      else bubble_up_max_(i);
    }
  }

  void bubble_up_min_(size_t i)
  {
    size_t p = parent_(i);
    size_t gp = parent_(p);

    if (gp < p && data_[i] < data_[gp])
    {
      swap_(i,gp);
      bubble_up_min_(gp);
    }
  }

  void bubble_up_max_(size_t i)
  {
    size_t p = parent_(i);
    size_t gp = parent_(p);

    if (gp < p && data_[gp] < data_[i])
    {
      swap_(i,gp);
      bubble_up_max_(gp);
    }
  }

  /**
   * Bottom-up min-max heap construction.
   */
  void minmax_heapify_()
  {
    for (size_t i=1; i<=data_.size(); ++i)
      trickle_down_(data_.size()-i);
  }
};

} // namespace srvf

#endif // MINMAXHEAP_H
