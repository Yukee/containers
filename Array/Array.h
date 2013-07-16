#ifndef ARRAY_H
#define ARRAY_H

#include <stdexcept>
#include <ostream>
#include <cstdarg> // variable number of arguments for the constructors
#include <cstdio> // memcpy

#include "../Vector.h"

namespace ij
{

typedef unsigned int INT;

template <class T=int>
class Array
{
public:

    // initialize a 2D array containing 1 data point
    Array();

    // usage: number of directions d followed by unsigned integers
    // indicating the dimension of the arrray in each direction
    // WARNING: this constructor can't check if d agrees with the number
    // of following arguments you pass, or even if the type of these
    // arguments (unsigned int) is correct!
    // so use it carefully (or even better: don't use it)
    Array(const unsigned int d, ...);

    // initialize an Array of dimensions D
    Array(const Vector<INT> & D);

    // initialize an Array of dimensions D filled with value v
    Array(const Vector<INT> &D, const T & v);

    ~Array();

    // getters
    Vector<INT> get_dim();
    unsigned int get_size();

    // check if arrays have the same dimensions
    bool are_compatible(const Array<T> & rhs) const;

    // check if arrays have the same dimensions and values
    inline friend bool operator==(const Array<T> & a1, const Array<T> & a2)
    {
      // same dimensions
      bool are_eq = a1.are_compatible(a2);

      // same values
      for(INT it=0;it<a1.N_ && are_eq;++it) are_eq = (a1.data_[it] == a2.data_[it]);

      return are_eq;
    }

    inline Array<T> & operator=(const Array<T> & rhs)
    {
      if(!(rhs == *this))
        {
            D_ = rhs.D_;
            N_ = rhs.N_;

            if(data_) delete[] data_;
            data_ = NULL;
            data_ = new T[N_];
            for(INT it=0;it<N_;++it) data_[it] = rhs.data_[it];
        }
      return *this;
    }

    inline T & operator[](const INT i) const
    {
        if(i >= N_) throw std::out_of_range("In Array<T>::operator[](unsigned int i)");
        return data_[i];
    }

    inline T & operator[](const Vector<INT> & r)
    {
        INT d = D_.size();

        if(r.size() != d) throw std::invalid_argument("In Array<T>::operator[](Vector<unsigned int> r) size of r differs from the number of dimensions of the Array");
        for(INT i=0;i<d;i++) if(D_[i] <= r[i]) throw std::out_of_range("In Array<T>::operator[](Vector<unsigned int> r)");

        int element = r[d-1];
            for(unsigned int i=0; i<d-1; i++)
            {
                element*=D_[d-2-i];
                element+=r[d-2-i];
            }
        return data_[element];
    }

  //void resize(const unsigned int d, ...);

private :

  // pointer to data stored in the array
  T *data_;

  // dimensions of the array
  Vector<INT> D_;

  // number of elements
  unsigned int N_;
};

template <class T>
Array<T>::Array()
{
    D_.resize(2);
    D_[0] = 1;
    D_[1] = 1;

    N_ = D_[0]*D_[1];

    data_ = new T[N_];
}

template <class T>
Array<T>::Array(const unsigned int d, ...)
{
    D_.resize(d);

    va_list vl;
    va_start(vl, d);

    N_ = 1;
    for(unsigned int i=0;i<d;i++)
    {
        D_[i] = va_arg(vl, unsigned int);
        N_ *= D_[i];
    }

    data_ = new T[N_];
}

template <class T>
Array<T>::Array(const Vector<INT> &D)
{
    D_ = D;
    N_ = 1;
    for(INT i=0;i<D_.size();i++) N_ *= D_[i];

    data_ = new T[N_];
}

template <class T>
Array<T>::Array(const Vector<INT> &D, const T &v)
{
    D_ = D;
    N_ = 1;
    for(INT i=0;i<D_.size();i++) N_ *= D_[i];

    data_ = new T[N_];

    for(INT it=0;it<N_;++it) data_[it] = v;
}

template <class T>
Vector<INT> Array<T>::get_dim()
{
    return D_;
}

template <class T>
unsigned int Array<T>::get_size()
{
    return N_;
}

template <class T>
bool Array<T>::are_compatible(const Array<T> &rhs) const
{
    // same number of elements
  bool are_comp = (N_ == rhs.N_);

  // same number of directions
  are_comp = (D_.size() == rhs.D_.size());

  // same dimensions
  for(unsigned int i=0;i<D_.size() && are_comp;i++) are_comp = (D_[i] == rhs.D_[i]);

  return are_comp;
}

template <class T>
Array<T>::~Array()
{
  if(data_) delete[] data_;
  data_ = NULL;
}

}

#endif // ARRAY_H
