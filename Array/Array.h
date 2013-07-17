#ifndef ARRAY_H
#define ARRAY_H

#include <stdexcept>
#include <ostream>
#include <cstdarg> // variable number of arguments for the constructors
#include <vector>

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
    Array(const std::vector<INT> &D);

    // initialize an Array of dimensions D filled with value v
    Array(const std::vector<INT> &D, const T & v);

    // getters
    std::vector<INT> get_dim();
    INT get_size();

    // check if arrays have the same dimensions and values
    inline friend bool operator==(const Array<T> & a1, const Array<T> & a2)
    {
      // same dimensions
        bool are_eq = (a1.D_ == a2.D_);

      // same values
      for(INT it=0;it<a1.N_ && are_eq;++it) are_eq = (a1.data_[it] == a2.data_[it]);

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
    T & operator[](const std::vector<INT> & r) const;
    T & at(const INT i) const;
    T & at(const std::vector<INT> & r) const;


private :

  // data stored in the array
  std::vector<T> data_;

  // dimensions of the array
  std::vector<INT> D_;
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
Array<T>::Array(const std::vector<INT> &D)
{
    D_ = D;
    INT n = 1;
    for(INT i=0;i<D_.size();i++) n *= D_[i];
    data_.resize(n);
}

template <class T>
Array<T>::Array(const std::vector<INT> &D, const T &v)
{
    D_ = D;
    INT n = 1;
    for(INT i=0;i<D_.size();i++) n *= D_[i];
    data_.resize(n, v);
}

template <class T>
std::vector<INT> Array<T>::get_dim()
{
    return D_;
}

template <class T>
INT Array<T>::get_size()
{
    return data_.size();
}

template <class T>
T Array<T>::operator[](const INT i) const
{
    return data_[i];
}

template <class T>
T Array<T>::at(const INT i) const
{
    return data_.at(i);
}

template <class T>
T Array<T>::operator[](const std::vector<INT> & r) const
{
    INT n = D_.size();
    INT elem = D_[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= r[n-i];
        elem += D_[n-i];
    }

    return data_[elem];
}

template <class T>
T Array<T>::at(const std::vector<INT> & r) const
{
    if(r.size() != D_.size()) throw std::invalid_argument("In Array<T>::at");

    for(INT i=0;i<r.size();i++) if(r[i] >= D_[i]) throw std::out_of_range("In Array<T>::at");

    INT n = D_.size();
    INT elem = D_[n-1];
    for(INT i=0;i<n-1;i++)
    {
        elem *= r[n-i];
        elem += D_[n-i];
    }

    return data_[elem];
}

}

#endif // ARRAY_H
