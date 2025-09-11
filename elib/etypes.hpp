/*===========================================================================*\
*/ /**
  @file etypes.hpp

  General purpose C library.
  Common definitions.

  @date Created: 02.05.2024
  @date Last Revision:

  @copyright Copyright 2024 E.Yanenko. All rights reserved.
*/

#ifndef __ETYPES_HPP
#define __ETYPES_HPP

#include <cstdlib>
#include <algorithm>
#include <limits>
#include <cmath>

namespace ey
{
/*===========================================================================*\
*/ /**
  @defgroup etypes Common definitions.
 */
///@{
typedef size_t nat;               ///< Natural numbers.
typedef long double real;         ///< Real values.

#define CNT_END ((nat)-1)         ///< Item count, which means total amount of items.

#define BIT_END (CNT_END & ~(CNT_END >> 1)) ///< Highest bit for natural numbers.

#define CNT_BIT (sizeof(nat)*8)   ///< Amount of bits for natural numbers.

#ifdef sign
#undef sign
#endif
/*===========================================================================*\
*/ /**
  Get sign.

  @tparam T argument type.
  @param x Argument.
  @return Sign of argument: -1 - negative, 1 - positive, 0 - zero.
*/
template <typename T> int sign(const T& x)
{
  return (x > (T)0)-(x < (T)0);
}
/*===========================================================================*\
*/ /**
  Compare real value with zero.

  @tparam T Type of arguments.
  @param a Values to compare.
  @param s Scale to compare.
  @param eps Relative error allowed.
  @param thr Absolute threshold.
  @return -1 - a < 0, 1 - a > 0, 0 - a approximately equal to 0.
*/
template<typename T> inline int CmpReal0(T a, T s=(T)1.0, T eps=std::numeric_limits<T>::epsilon(), T thr=std::numeric_limits<T>::min())
{
  s = std::max(thr,s*eps);
  if (a < -s)
    return -1;
  if (a > s)
    return 1;
  return 0;
}
/*===========================================================================*\
*/ /**
  Compare real values.

  @tparam T Type of arguments.
  @param a,b Values to compare.
  @param eps Relative error allowed.
  @param thr Absolute threshold.
  @return -1 - a < b, 1 - a > b, 0 - a approximately equal to b.
*/
template<typename T> inline int CmpReal(T a, T b, T eps=std::numeric_limits<T>::epsilon(), T thr=std::numeric_limits<T>::min())
{
  return CmpReal0(a-b,std::abs(a)+std::abs(b),eps,thr);
}
/*===========================================================================*\
*/ /**
  Computes sqrt(a^2 + b^2) without destructive underflow or overflow.

  @tparam T Type of arguments.
  @param a,b Input values to compare.
  @return sqrt(a^2+b^2).
*/
template<typename T> T pythag(T a, T b)
{
  T absa,absb;
  absa=std::abs(a);
  absb=std::abs(b);
  if (absa > absb)
  {
    absb/=absa;
    return absa*std::sqrt((T)1.0+absb*absb);
  }
  if (absb == 0.0)
    return 0.0;
  absa/=absb;
  return absb*std::sqrt((T)1.0+absa*absa);
}
///@}
};
#endif

