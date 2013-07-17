#ifndef ARRAY_H
#define ARRAY_H

#include <stdexcept>
#include <ostream>
#include <cstdarg> // variable number of arguments for the constructors
#include <vector>

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
    Array(const Vector<INT> &D);

    // initialize an Array of dimensions D filled with value v
    Array(const Vector<INT> &D, const T & v);

    // getters
    Vector<INT> get_dim() const;
    INT get_size() const;

    // check if arrays have the same dimensions and values
    inline friend bool operator==(const Array<T> & a1, const Array<T> & a2)
    {
      // same dimensions
        bool are_eq = (a1.D_ == a2.D_);

      // same values
      for(INT it=0;it<a1.D_.size() && are_eq;++it) are_eq = (a1.data_[it] == a2.data_[it]);

      return are_eq;
    }

    inline Array<T> & operator=(const Array<T> & rhs)
    {
      if(!(rhs == *this))
        {
            D_ = rhs.D_;

            data_ = rhs.data_;
        }
      return *this;
    }

    // access operators
    T & operator[](const INT i) const;
    T & operator[](const INT i);

    T & operator[](const Vector<INT> & r) const;
    T & operator[](const Vector<INT> & r);

    T & at(const INT i) const;
    T & at(const INT i);

    T & at(const Vector<INT> & r) const;
    T & at(const Vector<INT> & r);

    // usual operations
    template <class U> friend Array<U> & operator+(const Array<U> & a1, const Array<U> & a2);
    template <class U> friend Array<U> & operator*(const Array<U> & a1, const Array<U> & a2);

private :

  // data stored in the array
  std::vector<T> data_;

  // dimensions of the array
  Vector<INT> D_;
};

template <class T>
Array<T>::Array()
{
    D_.resize(2, 1);
    data_.resize(2);
}

template <class T>
Array<T>::Array(const unsigned int d, ...)
{
    D_.resize(d);

    va_list vl;
    va_start(vl, d);

    INT n = 1;
    for(unsigned int i=0;i<d;i++)
    {
        D_[i] = va_arg(vl, unsigned int);
        n *= D_[i];
    }

    data_.resize(n);
}

template <class T>
Array<T>::Array(const Vector<INT> &D)
{
    D_ = D;
    INT n = 1;
    for(INT i=0;i<D_.size();i++) n *= D_[i];
    data_.resize(n);
}

template <class T>
Array<T>::Array(const Vector<INT> &D, const T &v)
{
    D_ = D;
    INT n = 1;
    for(INT i=0;i<D_.size();i++) n *= D_[i];
    data_.resize(n, v);
}

template <class T>
Vector<INT> Array<T>::get_dim() const
{
    return D_;
}

template <class T>
INT Array<T>::get_size() const
{
    return data_.size();
}

template <class T>
T & Array<T>::operator[](const INT i) const
{
    return data_[i];
}

template <class T>
T & Array<T>::operator[](const INT i)
{
    return data_[i];
}

template <class T>
T & Array<T>::at(const INT i) const
{
    return data_.at(i);
}

template <class T>
T & Array<T>::at(const INT i)
{
    return data_.at(i);
}

template <class T>
T & Array<T>::operator[](const Vector<INT> & r) const
{
    INT n = D_.size();
    INT elem = r[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= D_[n-i-2];
        elem += r[n-i-2];
    }

    return data_[elem];
}

template <class T>
T & Array<T>::operator[](const Vector<INT> & r)
{
    INT n = D_.size();
    INT elem = r[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= D_[n-i-2];
        elem += r[n-i-2];
    }

    return data_[elem];
}

template <class T>
T & Array<T>::at(const Vector<INT> & r) const
{
    if(r.size() != D_.size()) throw std::invalid_argument("In ij::Array<T>::at");

    for(INT i=0;i<r.size();i++) if(r[i] >= D_[i]) throw std::out_of_range("In ij::Array<T>::at");

    INT n = D_.size();
    INT elem = r[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= D_[n-i-2];
        elem += r[n-i-2];
    }

    return data_[elem];
}

template <class T>
T & Array<T>::at(const Vector<INT> & r)
{
    if(r.size() != D_.size()) throw std::invalid_argument("In ij::Array<T>::at");

    for(INT i=0;i<r.size();i++) if(r[i] >= D_[i]) throw std::out_of_range("In ij::Array<T>::at");

    INT n = D_.size();
    INT elem = r[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= D_[n-i-2];
        elem += r[n-i-2];
    }

    return data_[elem];
}

template <class T>
Array<T> & operator+(const Array<T> & a1, const Array<T> & a2)
{
    INT N = a1.D_.size();
    if(N != a2.D_.size()) throw std::invalid_argument("In ij::Array<T>::operator+");

    Array<T> temp (N);
    for(int it=0;it<N;++it) temp[it] = a1[it] + a2[it];

    return temp;
}

template <class T>
Array<T> & operator*(const Array<T> & a1, const Array<T> & a2)
{
    INT N = a1.D_.size();
    if(N != a2.D_.size()) throw std::invalid_argument("In ij::Array<T>::operator*");

    Array<T> temp (N);
    for(int it=0;it<N;++it) temp[it] = a1[it] * a2[it];

    return temp;
}

}

#endif // ARRAY_H
