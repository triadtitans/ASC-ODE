#ifndef FILE_AUTODIFF
#define FILE_AUTODIFF

#include<math.h>
#include<iostream>
/**************************************************************************/
/* File:   autodiff.h                                                     */
/* Author: Linus Knoll                                                    */
/* Date:   26. March. 24                                                  */
/**************************************************************************/


// Automatic differentiation datatype


namespace ASC_ode
{

template <size_t D, typename SCAL = double>
class AutoDiffVec
{
  SCAL val;
  SCAL dval[D?D:1];
public:

  typedef AutoDiffVec<D,SCAL> TELEM;
  typedef SCAL TSCAL;

  AutoDiffVec() = default;

  AutoDiffVec  (SCAL aval) throw()
  {
    val = aval;
    for (size_t i = 0; i < D; i++)
      dval[i] = 0;
  }

  /// init object with (val, e_diffindex)
  AutoDiffVec  (SCAL aval, size_t diff_ind)  throw()
  {
    val = aval;
    for (size_t i = 0; i < D; i++)
      dval[i] = 0;
    
    dval[diff_ind] = 1;
  }

  

  AutoDiffVec & operator= (const AutoDiffVec & ad2) = default;

  /// assign constant value
  AutoDiffVec & operator= (SCAL aval) throw()
  {
    val = aval;
    for (size_t i = 0; i < D; i++)
      dval[i] = 0;
    return *this;
  }
  
  /// returns value
  inline SCAL Val () const throw() { return val; }
  
  /// returns partial derivative
  inline SCAL DVal (size_t i) const throw() { return dval[i]; }

  void LoadGradient (const SCAL * p) 
  {
    for (size_t i = 0; i < D; i++)
      dval[i] = p[i];
  }

  inline SCAL & Val () throw() { return val; }

  inline SCAL & DVal (size_t i) throw() { return dval[i];}
  
};

// Helper Functions

template<size_t D, typename SCAL>
std::ostream & operator<< (std::ostream & ost, const AutoDiffVec<D,SCAL> & x)
{
  ost << x.Val() << ", D = ";
  for (size_t i = 0; i < D; i++)
    ost << x.DVal(i) << " ";
  return ost;
}


// AutoDiff + AutoDiff
template<size_t D, typename SCAL>
AutoDiffVec<D, SCAL> operator+ (const AutoDiffVec<D, SCAL> & x, const AutoDiffVec<D, SCAL> & y)
{
    AutoDiffVec<D, SCAL> res;
    res.Val() = x.Val() + y.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x.DVal(i) + y.DVal(i);

    return res;
}

// minus AutoDiffVec
template<size_t D, typename SCAL>
AutoDiffVec<D,SCAL> operator- (const AutoDiffVec<D,SCAL> & x) throw()
{
  AutoDiffVec<D,SCAL> res;
  res.Val() = -x.Val();
  for (size_t i = 0; i < D; i++)
    res.DVal(i) = -x.DVal(i);
  return res;
}

// AutoDiff - AutoDiff
template<size_t D, typename SCAL>
AutoDiffVec<D, SCAL> operator- (const AutoDiffVec<D, SCAL> & x, const AutoDiffVec<D, SCAL> & y)
{
    AutoDiffVec<D, SCAL> res;
    res.Val() = x.Val() - y.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x.DVal(i) - y.DVal(i);

    return res;
}

// AutoDiff + SCALAR
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator+ (const AutoDiffVec<D, SCAL> & x, SCAL2 y)
{
    AutoDiffVec<D, SCAL> res;
    res.Val() = x.Val() + y;
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x.DVal(i);

    return res;
}

// AutoDiff - SCALAR
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator- (const AutoDiffVec<D, SCAL> & x, SCAL y)
{
    AutoDiffVec<D, SCAL> res;
    res.Value() = x.Val() - y;
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x.DVal(i);

    return res;
}

// SCALAR + AutoDiff
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator+ (SCAL2 x, const AutoDiffVec<D, SCAL> &y)
{
    AutoDiffVec<D, SCAL> res;
    res.Val() = x + y.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = y.DVal(i);

    return res;
}

// SCALAR - AutoDiff
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator- (SCAL x, const AutoDiffVec<D, SCAL> & y)
{
    AutoDiffVec<D, SCAL> res;
    res.Val() = y - x.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = -x.DVal(i);

    return res;
}

// SCALAR * AutoDiff
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator* (SCAL2 x, const AutoDiffVec<D, SCAL> & y)
{
  AutoDiffVec<D, SCAL> res;
    res.Val() = x * y.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x * y.DVal(i);

    return res;
}

// AutoDiff * SCALAR
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D, SCAL> operator* (const AutoDiffVec<D, SCAL> & y, SCAL2 x)
{
  AutoDiffVec<D, SCAL> res;
    res.Val() = y * x.Val();
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = y * x.DVal(i);

    return res;
}

// AutoDiff * AutoDiff
template<size_t D, typename SCAL>
AutoDiffVec<D, SCAL> operator* (const AutoDiffVec<D, SCAL> & x, const AutoDiffVec<D, SCAL> & y)
{
    AutoDiffVec<D, SCAL> res;
    SCAL x_val = x.Val();
    SCAL y_val = y.Val();

    res.Val() = x_val*y_val;
    for (size_t i = 0; i < D; i++)
        res.DVal(i) = x_val*y.DVal(i) + y_val*x.DVal(i);

    return res;
}

// Inverse of AutoDiffVec
template<size_t D, typename SCAL>
AutoDiffVec<D,SCAL> Inv (const AutoDiffVec<D,SCAL> & x)
{
  AutoDiffVec<D,SCAL> res(1.0 / x.Val());
  for (size_t i = 0; i < D; i++)
    res.DVal(i) = -x.DVal(i) / (x.Val() * x.Val());
  return res;
}

// AutoDiffVec / AutoDiffVec
template<size_t D, typename SCAL>
AutoDiffVec<D,SCAL> operator/ (const AutoDiffVec<D,SCAL> & x, const AutoDiffVec<D,SCAL> & y)
{
  return x * Inv (y);
}

// AutoDiffVec / double
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D,SCAL> operator/ (const AutoDiffVec<D,SCAL> & x, SCAL2 y)
{
  return (1.0/y) * x;
}

// SCALAR / AutoDiffVec
template<size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D,SCAL> operator/ (SCAL2 x, const AutoDiffVec<D,SCAL> & y)
{
  return x * Inv(y);
}

// AutoDiff += SCALAR
template <size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D,SCAL> & operator+= (AutoDiffVec<D,SCAL> & x, SCAL2 y) throw()
{
  x.Val() += y;
  return x;
}

// AutoDiff += AutoDiff
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> & operator+= (AutoDiffVec<D,SCAL> & x, AutoDiffVec<D,SCAL> y)
{
  x.Val() += y.Val();
  for (size_t i = 0; i < D; i++)
    x.DVal(i) += y.DVal(i);
  return x;
}

// AutoDiff -= AutoDiff
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> & operator-= (AutoDiffVec<D,SCAL> & x, AutoDiffVec<D,SCAL> y)
{
  x.Val() -= y.Val();
  for (size_t i = 0; i < D; i++)
    x.DVal(i) -= y.DVal(i);
  return x;

}


// AutoDiff -= SCALAR
template <size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D,SCAL> & operator-= (AutoDiffVec<D,SCAL> & x, SCAL2 y)
{
  x.Val() -= y;
  return x;
}

// AutoDiff *= AutoDiff
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> & operator*= (AutoDiffVec<D,SCAL> & x, AutoDiffVec<D,SCAL> y) 
{
  for (size_t i = 0; i < D; i++)
    x.DVal(i) = x.DVal(i)*y.Val() + x.Val() * y.DVal(i);
  x.Val() *= y.Val();
  return x;
}

// AutoDiff *= SCALAR
template <size_t D, typename SCAL, typename SCAL2>
AutoDiffVec<D,SCAL> & operator*= (AutoDiffVec<D,SCAL> & x, SCAL2 y) 
{
  x.Val() *= y;
  for (size_t i = 0; i < D; i++)
    x.DVal(i) *= y;
  return x;
}

// Autodiff /= SCALAR
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> & operator/= (AutoDiffVec<D,SCAL> & x, SCAL y) 
{
  SCAL iy = 1.0 / y;
  x.Val() *= iy;
  for (size_t i = 0; i < D; i++)
    x.DVal(i) *= iy;
  return x;
}

using std::sqrt;
template<size_t D, typename SCAL>
AutoDiffVec<D,SCAL> sqrt (const AutoDiffVec<D,SCAL> & x)
{
  AutoDiffVec<D,SCAL> res;
  res.Val() = sqrt(x.Val());
  for (size_t j = 0; j < D; j++)
    res.DVal(j) = 0.5 / res.Val() * x.DVal(j);
  return res;
}

using std::log;
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> log (AutoDiffVec<D,SCAL> x)
{
  AutoDiffVec<D,SCAL> res;
  res.Value() = log(x.Val());
  for (size_t k = 0; k < D; k++)
    res.DVal(k) = x.DVal(k) / x.Val();
  return res;
}

using std::exp;
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> exp (AutoDiffVec<D,SCAL> x)
{
  AutoDiffVec<D,SCAL> res;
  res.Val() = exp(x.Val());
  for (size_t k = 0; k < D; k++)
    res.DVal(k) = x.DVal(k) * res.Val();
  return res;
}

using std::pow;
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> pow (AutoDiffVec<D,SCAL> x, AutoDiffVec<D,SCAL> y )
{
  return exp(log(x)*y);
}

using std::sin;
template <size_t D, typename SCAL>
AutoDiffVec<D,SCAL> sin (AutoDiffVec<D,SCAL> x)
{
  AutoDiffVec<D,SCAL> res;
  res.Val() = sin(x.Val());
  SCAL c = cos(x.Val());
  for (size_t k = 0; k < D; k++)
    res.DVal(k) = x.DVal(k) * c;
  return res;
}


}

#endif