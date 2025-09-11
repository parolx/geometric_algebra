/*===========================================================================*\
*/ /**
  @file tensor.hpp

  General purpose C library.
  Tensors and matrices.

  @date Created: 02.05.2024
  @date Last Revision:

  @copyright Copyright 2024 E.Yanenko. All rights reserved.
*/

#ifndef __ETENSOR_HPP
#define __ETENSOR_HPP

#include "elib/etypes.hpp"
#include <array>
#include <stdexcept>
#include <vector>

namespace ey
{
namespace tens
{
/*===========================================================================*\
*/ /**
  @defgroup tens Tensors and matrices.
 */
///@{
/**
  @addtogroup spans Sub-tensors, access to tensor components.
  @ingroup tens
  @tparam D Amount of sub-tensor dimensions.
  @tparam R Tensor value type.
  @param s Pointer to array with sub-tensor dimensions.
  @param b Pointer to tensor data.
  @param n Component index.
  @param d Dimension index.
  @param bSigned Signed or unsigned result. Sign defines if upper (positive) or lower (negative) index.
  @see    Tensor<D,R>
 */
///@{
/**  Generic sub-tensor, dimension decreased by one.  */
template<nat D, typename R> struct SubTensor
{
  /**  Constructor.  */
  SubTensor(int const* s, R* b): m_dims(s), m_data(b) {;}
  /**  Access to sub-sub-tensor.  */
  SubTensor<D-1,R> operator[](nat n) const {return {m_dims+1, m_data+n*m_dims[2*D-2]};}
  /**  Direct access to tensor value.  */
  const R& operator()(nat n) const {return m_data[n];}
  /**  Size of certain dimension.  */
  int DSize(nat d, bool bSigned=false) const {return bSigned? m_dims[d] : std::abs(m_dims[d]);}
private:
  int const* m_dims;          ///< Pointer to array with sizes of sub-tensor dimensions.
  R* m_data;                  ///< Pointer to tensor data.
};
/** Vector-like sub-tensor, direct access to components.  */
template<typename R> struct SubTensor<1,R>
{
  /**  Constructor.  */
  SubTensor(int const* s, R* b): m_dims(s), m_data(b) {;}
  /**  Access to sub-sub-tensor.  */
  R& operator[](nat n) const {return m_data[n];}
  /**  Direct access to tensor value.  */
  const R& operator()(nat n) const {return m_data[n];}
  /**  Size of certain dimension.  */
  int DSize(nat d, bool bSigned=false) const {return bSigned? m_dims[d] : std::abs(m_dims[d]);}
private:
  int const* m_dims;          ///< Pointer to array with sizes of sub-tensor dimensions.
  R* m_data;                  ///< Pointer to tensor data.
};
/** Last sub-tensor, access to single scalar.  */
template<typename R> struct SubTensor<0,R>
{
  /**  Constructor.  */
  SubTensor(int const* s, R* b): m_data(b) {;}
  /**  Dimension 0 - access to single scalar.  */
  operator R&() const {return *m_data;}
  /**  Size of certain dimension.  */
  int DSize(nat d, bool bSigned=false) const {return 0;}
private:
  R* m_data;                  ///< Pointer to tensor data.
};
///@}
/*===========================================================================*\
*/ /**
  Tensor, multi-dimensional array with upper/lower indices.

  @tparam D Tensor dimension.
  @tparam R Tensor value type.
*/
template<nat D, typename R=real> class Tensor
{
public:
  /**
  Constructor.

  @param dims Array with sizes of tensor dimensions. Positive - upper index, negative - lower.
  */
  Tensor(const std::array<int,D>& dims) {Resize(dims);}
  /**
  Constructor.

  @param st Sub-tensor.
  */
  Tensor(const SubTensor<D,R>& st)
  {
    for (int i=0; i<D; ++i)
      m_dims[i]=st.DSize(i,true);
    DoResize();
    for (int i=0; i<Size(); ++i)
      m_data[i]=st(i);
  }
  /**
  Change tensor dimensions.

  @param dims Array with sizes of tensor dimensions. Positive - upper index, negative - lower.
  */
  void Resize(const std::array<int,D>& dims)
  {
    for (nat i=0; i<D; ++i)
      m_dims[i] = dims[i];
    DoResize();
  }
  /**
  Tensor / sub-tensor sizes.

  @param d Tensor or sub-tensor dimension.
  */
  nat Size(nat d=D) const {return D? m_dims[2*D-d] : 1;}
  /**
  Size of certain dimension.

  @param d Dimension index.
  @param bSigned Signed or unsigned result. Sign defines if upper (positive) or lower (negative) index.
  */
  int DSize(nat d, bool bSigned=false) const {return bSigned? m_dims[d] : std::abs(m_dims[d]);}
  /** Access to sub-tensors. */
  SubTensor<D-1,R> operator[](nat n) {return {std::addressof(m_dims[1]),std::addressof(m_data[n*m_dims[2*D-2]])};}
  /** Access to sub-tensors. */
  SubTensor<D-1,R const> operator[](nat n) const {return {std::addressof(m_dims[1]),std::addressof(m_data[n*m_dims[2*D-2]])};}
  /** Direct access to dimension sizes. */
  const std::array<int,2*D> Dims() const {return m_dims;}
  /** Direct access to tensor data. */
  R& operator()(nat n) {return m_data[n];}
  /** Direct access to tensor data. */
  const R& operator()(nat n) const {return m_data[n];}
  /** Tensor addition. */
  Tensor& operator +=(const Tensor& t)
  {
    if (t.Dims() != m_dims)
      throw std::domain_error("Attempt to add tensors of different size / dimension.");
    for (nat n=0; n<Size(); n++)
      m_data[n] += t(n);
    return *this;
  }
  /** Tensor subtraction. */
  Tensor& operator -=(const Tensor& t)
  {
    if (t.Dims() != m_dims)
      throw std::domain_error("Attempt to subtract tensors of different size / dimension.");
    for (nat n=0; n<Size(); n++)
      m_data[n] -= t(n);
    return *this;
  }
  /** Tensor scaling. */
  Tensor& operator *=(R r)
  {
    for (nat n=0; n<Size(); n++)
      m_data[n] *= r;
    return *this;
  }
  /** Tensor scaling. */
  Tensor& operator /=(R r)
  {
    if (r == 0.0)
      throw std::overflow_error("Tensor divide by zero.");
    r = (R)1.0/r;
    for (nat n=0; n<Size(); n++)
      m_data[n] *= r;
    return *this;
  }
protected:
  std::array<int,2*D> m_dims; ///< Array dimensions (till index D) and sizes of sub-arrays (till 2*D-1).
  std::vector<R> m_data;      ///< Array data.
  /** Perform resize. */
  void DoResize()
  {
    int s=1;
    if (D != 0)
    {
      for (int i=D-1; i>=0; --i)
      {
        s *= DSize(i);
        m_dims[2*D-1-i] = s;
      }
    }
    m_data.resize(s,(R)0);
  }
};
template<typename R=real> using Vector = Tensor<1,R>; ///< Vector.
template<typename R=real> using Matrix = Tensor<2,R>; ///< Matrix.
/*===========================================================================*\
*/ /**
  Returns unitary/scaled unitary matrix.

  @param v Value on main diagonal.
  @param nSize Matrix size.
  @return Returns unitary/scaled unitary matrix.
  @see    Matrix<nRows,nCols,R>
*/
template <typename R=real> Matrix<R> l(nat nSize, R v=1.0)
{
  Matrix<R> u({nSize,-(int)nSize});
  for (nat r=0; r<nSize; ++r)
      u[r][r] = v;
  return u;
}
/*===========================================================================*\
*/ /**
  Calculates transposed matrix.

  @return  Returns transposed matrix.
*/
template<typename R> Matrix<R> T(const Matrix<R>& m1)
{
  nat r,c;
  nat nRows=m1.DSize(0);
  nat nCols=m1.DSize(1);
  Matrix<R> m({m1.DSize(1,true),m1.DSize(0,true)});
  for (r=0; r<nRows; ++r)
  {
    for (c=0; c<nCols; ++c)
      m[r][c]=m1[c][r];
  }
  return m;
}
/**
  @addtogroup tensop1 Operators on tensors.
  @ingroup tens
  @tparam D Tensor dimension.
  @tparam R Tensor value type.
  @param t1 Tensor operand.
  @see    Tensor<D,R>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates negative tensor.

  @return  Returns negative tensor.
*/
template<nat D, typename R> Tensor<D,R> operator -(const Tensor<D,R>& t1)
{
  Tensor<D,R> t(t1);
  for (nat r=0; r<t.Size(); ++r)
    t(r)=-t1(r);
  return t;
}
/*===========================================================================*\
*/ /**
  Calculates scaled tensor.

  @return  Returns scaled tensor.
*/
template<nat D, typename R> inline Tensor<D,R> operator *(const Tensor<D,R>& t1, R r)
{
  Tensor<D,R> t(t1);
  t *= r;
  return t;
}
/*===========================================================================*\
*/ /**
  Calculates scaled tensor.

  @return  Returns scaled tensor.
*/
template<nat D, typename R> inline Tensor<D,R> operator *(R r, const Tensor<D,R>& t1)
{
  Tensor<D,R> t(t1);
  t *= r;
  return t;
}
/*===========================================================================*\
*/ /**
  Calculates scaled tensor.

  @return  Returns scaled tensor.
*/
template<nat D, typename R> inline Tensor<D,R> operator /(const Tensor<D,R>& t1, R r)
{
  Tensor<D,R> t(t1);
  t /= r;
  return t;
}
///@}
/**
  @addtogroup tensop2 Sums and products of tensors.
  @ingroup tens
  @tparam D Tensor dimension.
  @tparam R Tensor value type.
  @param t1,t2 Tensor operands.
  @see    Tensor<D,R>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates sum of tensors.

  @return  Returns sum.
*/
template<nat D, typename R> inline Tensor<D,R> operator +(const Tensor<D,R>& t1, const Tensor<D,R>& t2)
{
  Tensor<D,R> t(t1);
  t += t2;
  return t;
}
/*===========================================================================*\
*/ /**
  Calculates difference of tensors.

  @return  Returns difference.
*/
template<nat D, typename R> inline Tensor<D,R> operator -(const Tensor<D,R>& t1, const Tensor<D,R>& t2)
{
  Tensor<D,R> t(t1);
  t -= t2;
  return t;
}
/*===========================================================================*\
*/ /**
  Calculates product of tensors.

  @return  Returns product.
*/
template<nat nDim1, nat nDim2, typename R> Tensor<nDim1+nDim2-2,R> operator *(const Tensor<nDim1,R>& t1, const Tensor<nDim2,R>& t2)
{
  nat n,r,c,i,j;
  if (t1.DSize(nDim1-1,true)+t2.DSize(0,true) != 0)
    throw std::domain_error("Different tensor dimensions in multiplication.");

  nat s1=1;
  nat s2=1;
  std::array<int,nDim1+nDim2-2> dims;
  for (i=j=0; i<nDim1-1; ++i)
  {
    s1*=t1.DSize(i);
    dims[j++]=t1.DSize(i,true);
  }
  for (i=1; i<nDim2; ++i)
  {
    s2*=t2.DSize(i);
    dims[j++]=t2.DSize(i,true);
  }
  Tensor<nDim1+nDim2-2,R> t(dims);

  nat d=t2.DSize(0);
  for (n=r=0; r<s1; ++r)
  {
    for (c=0; c<s2; ++c,++n)
    {
      for (i=0; i<d; ++i)
        t(n) += t1(r*d+i)*t2(c+i*s2);
    }
  }
  return t;
}
///@}
/*===========================================================================*\
*/ /**
  Computes a solution of system of linear equations using Gauss elimination method.

  @param a Equation matrix, will be destroyed on output!
  @param x Right side of equations, will be replaced with solution on output.
          Supply unitary matrix to calculate inverse of a.
  @return Returns determinant of matrix. If determinant equals zero - system can't be solved.
  @note   Implements Gauss's elimination algorithm for square factor matrix nDim x nDim
          and rectangle result matrix nDim x nRes.
  @see    Tensor<D,R>
*/
template <typename R=real> R SolveGauss(Matrix<R>& a, Matrix<R>& x)
{
  nat r,r1,rb,c;
  R e,u;
  nat nDim=a.DSize(0);
  if ((a.DSize(1) != (int)nDim) || (x.DSize(0) != (int)nDim))
    throw std::domain_error("Wrong matrix dimensions in SolveGauss().");

  nat nRes=x.DSize(1);
  R rDet=1.0;

  for (r=0; r<nDim; r++)                // First step - initial matrix -> upper
  {                                     //   triangle matrix.
    u=a[r][r];
    rb=r;

    for (r1=r+1; r1<nDim; r1++)
    {
      if (std::abs(a[r1][r]) > std::abs(u))
      {
        u=a[r1][r];
        rb=r1;
      }
    }
    if (u == 0.0)
    {
      rDet=0.0;
      break;
    }
    if (rb != r)
    {                                   // Swap two rows in the matrix and vector.
      for (c=r; c<nDim; c++)            // Swap two rows in the matrix.
        ::std::swap(a[rb][c],a[r][c]);
      for (c=0; c<nRes; c++)            // Swap two rows in the vector.
        ::std::swap(x[rb][c],x[r][c]);
      rDet=-rDet;
    }
    // Diagonal element -> 1.
    rDet*=u;
    a[r][r]=1.0;
    for (c=r+1; c<nDim; c++)
      a[r][c] /= u;
    for (c=0; c<nRes; c++)
      x[r][c] /= u;
    for (r1=r+1; r1<nDim; r1++)         // All elements in the column
    {                                   //   under diagonal -> 0.
      e=a[r1][r];
      a[r1][r] = 0.0;
      for (c=r+1; c<nDim; c++)
        a[r1][c] -= e*a[r][c];
      for (c=0; c<nRes; c++)
        x[r1][c] -= e*x[r][c];
    }
  }
  if (r > 1)                            // Second step - upper triangle
  {
    for (r--; r>0; r--)                 // All elements in the column
    {
      for (r1=r-1; r1<CNT_END; r1--)    // All elements in the column
      {
        e=a[r1][r];
        a[r1][r] = 0.0;
        for (c=0; c<nRes; c++)
          x[r1][c] -= e*x[r][c];
      }
    }
  }
  return rDet;
}
///@}
};
};
#endif

