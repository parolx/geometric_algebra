/*===========================================================================*\
*/ /**
  @file geom.hpp

  General purpose C library.
  Geometrical algebra. Support for projection and conformal geometry.

  @date Created: 18.04.2024
  @date Last Revision:

  @copyright Copyright 2024 E.Yanenko. All rights reserved.
*/

#ifndef __EGEOM_HPP
#define __EGEOM_HPP

#include "elib/etypes.hpp"
#include "elib/tensor.hpp"
#include <ostream>
#include <map>
#include <array>
#include <stdexcept>

namespace ey
{
namespace ga
{
/*===========================================================================*\
*/ /**
  @defgroup geom Geometrical algebra.
 */
/**
  @addtogroup gspace Space of geometric algebra. Support for projection and conformal geometry.
  @ingroup geom
 */
///@{
template<typename Sp> class BElem;
template<typename Sp> class MVect;

/*===========================================================================*\
*/ /**
  Space of geometric algebra.

  @tparam _Pos Number of space dimensions with positive signature.
  @tparam _Neg Number of space dimensions with negative signature.
  @tparam _Proj Number of additional projection dimensions. 0 - usual vector space, 1 - projective geometry,
          2 - conformal geometry.
  @tparam _R Type of real values.
  @tparam _Cospace Boolean, true if it is co-space.
*/
template<nat _Pos=3,nat _Neg=0,nat _Proj=0,typename _R=real,bool _Cospace=false> class Space
{
public:
  using R = _R;               ///< Type of real values.
  using Cospace = Space<_Pos,_Neg,_Proj,_R,!_Cospace>;    ///< Co-space.
  /** Dimensions of sub-spaces. */
  enum dims: nat
  {
    Pos = _Pos,               ///< Amount of axes with positive signature.
    Neg = _Neg,               ///< Amount of axes with negative signature.
    Proj = _Proj,             ///< Amount of projection axes.
    Dim = _Pos+_Neg+_Proj,    ///< Total space dimension.
  };
  /** Grades. */
  enum grades: nat
  {
    scalar = 0,               ///< Grade of scalar.
    vect,                     ///< Grade of vector.
    bivect,                   ///< Grade of bivector.
    trivect,                  ///< Grade of trivector.
    pseudotrivect = Dim-3,    ///< Grade of pseudo-trivector.
    pseudobivect = Dim-2,     ///< Grade of pseudo-bivector.
    pseudovect = Dim-1,       ///< Grade of pseudo-vector.
    pseudoscalar = Dim        ///< Grade of pseudo-scalar.
  };
  /** Starting indices of sub-spaces. */
  enum indices: nat
  {
    weight = 0,               ///< Projective space.
    time = Proj,              ///< Minkovsky space.
    vect0 = Proj+(Neg>0? Pos : 0)   ///< Space vector.
  };
  /** Scalar and pseudo-scalar. */
  enum ones: nat
  {
    l = 0,                    ///< Scalar 1 basis space element bits. Can be
                              ///<   automatically converted to basis element.
    I = ((nat)1<<Dim)-1       ///< Pseudo-scalar I basis space element bis. Can be
                              ///<   automatically converted to basis element.
  };
  /**
  Get signature of certain coordinate.

  @param n Coordinate number.
  @return Coordinate signature.
  */
  static int Sign(nat n)
  {
    return n<Proj? (n? -1 : 1) : (n<Proj+Pos? 1 : -1);
  }
  /**
  Get grade according to bits.

  @param nBits Basis element bits.
  @return Grade of basis element.
  */
  static nat Grade(nat nBits)
  {
    nat nBit=1;
    nat nGrade=0;
    for (; nBit<=nBits; nBit<<=1)
    {
      if (nBits & nBit)
        nGrade++;
    }
    return nGrade;
  }
  /**
  Get element index according to element grade.

  @param nBits Basis element bits.
  @return Component index.
  */
  static nat Index(nat nBits);
  /**
  Get bits, according to given grade and index.

  @param nGrade Grade of basis element.
  @param idx Number of basis element of this grade.
  */
  static nat Bits(nat nGrade,nat idx);
  /**
  Get amount of components of given grade.

  @param nGrade Grade.
  @return Grade of basis element.
  */
  static nat Components(nat nGrade);
  /**
  Get vector.

  @param lst Initialization list.
  @return Vector.
  */
  static MVect<Space> Vector(const std::array<R,_Pos+_Neg>& lst = {})
  {
    R r;
    nat i=time;
    MVect<Space> mv;
    for (auto it=lst.cbegin(); (i < Dim) && (it != lst.cend()); ++it,++i)
    {
      r=*it;
      if (r != 0.0)
        mv[BElem<Space>(vect,i)]=r;
    }
    return mv;
  }
  /**
  Get point in projection space.

  @param lst Initialization list.
  @return Point.
  */
  static MVect<Space> inline Point(const std::array<R,_Pos+_Neg>& lst = {}) requires (Proj == 1)
  {
    MVect<Space> mv=Vector(lst);
    mv.w(1.0);
    return mv;
  }
  /**
  Get point in conformal space.

  @param lst Initialization list.
  @return Point.
  */
  static MVect<Space> Point(const std::array<R,_Pos+_Neg>& lst = {}) requires (Proj > 1)
  {
    MVect<Space> mv=Vector(lst);
    R s2=0.0;
    nat i=time;
    for (auto it=lst.cbegin(); (i < Dim) && (it != lst.cend()); ++it,++i)
      s2 += (*it)*(*it);
    mv.p(0.5-0.5*s2);
    mv.n(0.5+0.5*s2);
    return mv;
  }
  /**
  Get infinite point in conformal space.

  @return Point at infinity.
  */
  static constexpr MVect<Space> Infinity() requires (Proj > 1)
  {
    MVect<Space> mv(2);
    mv.p(-1.0);
    return mv;
  }
};
///@}
/**
  @addtogroup gspbelem Space basis elements support.
  @ingroup gspace
  @tparam _Pos Number of space dimensions with positive signature.
  @tparam _Neg Number of space dimensions with negative signature.
  @tparam _Proj Number of additional projection dimensions. 0 - usual vector space, 1 - projective geometry,
          2 - conformal geometry.
  @tparam _R Type of real values.
  @tparam _Cospace Boolean, true if it is co-space.
  @param nBits Basis element bits.
  @param nGrade Grade of element.
  @param idx Number of basis element of this grade.
  @see    Space<_Pos,_Neg,_Proj,_R>
 */
///@{
/*===========================================================================*\
*/ /**
  Get element index according to element grade.

  @return Component index.
*/
template<nat _Pos,nat _Neg,nat _Proj, typename _R, bool _Cospace> nat Space<_Pos,_Neg,_Proj,_R,_Cospace>::Index(nat nBits)
{
  nat idx,nGrade,nMask,nBit,nComb,nCombBit;

  idx=nGrade=0;
  for (nBit=nMask=nComb=1; nMask && (nMask<=nBits); nMask<<=1)
  // Check each bit and calculate amount of combinations nGrade from nBit.
  {
    if (nMask & nBits)
    {
      // Grade is equal to amount of set bits.
      nGrade++;
      nCombBit=nComb*(nBit-nGrade)/nGrade;
      // Increase index by combination count corresponding to set bit.
      idx+=nCombBit;
      nComb+=nCombBit;
    }
    else
      nComb=nComb*nBit/(nBit-nGrade);
    nBit++;
  }
  return idx;
}
/*===========================================================================*\
*/ /**
  Get bits, according to given grade and index.
*/
template<nat _Pos,nat _Neg,nat _Proj, typename _R, bool _Cospace> nat Space<_Pos,_Neg,_Proj,_R,_Cospace>::Bits(nat nGrade,nat idx)
{
  nat nBit,nComb,nCombBit;
  nat nBits=0;
  if ((nGrade > Dim) || (idx >= Components(nGrade)))
    return 0;
  for (; idx&&nGrade; nGrade--)
  {
    nCombBit=1;
    // Index is not zero, find how many bits are necessary to allow idx combinations.
    for (nBit=nComb=nGrade+1; nComb <= idx; )
    {
      nCombBit=nComb;
      nBit++;
      // Avoid overflow.
      nComb=nComb*nBit/(nBit-nGrade);
    }
    // This bit is necessarily set. This allows us to decrement the grade.
    nBits |= (nat)1<<(nBit-1);
    // Decrease index by combination count corresponding to set bit.
    idx-=nCombBit;
  }
  if (nGrade)
  // Set rest low bits, if index is zero and grade is still not zero.
  // This if() is necessary, because 1<<nGrade does not work at nGrade=32. Why???
    nBits |= ((nat)2<<(nGrade-1))-1;
  return nBits;
}
/*===========================================================================*\
*/ /**
  Get amount of components of given grade.

  @return Grade of basis element.
*/
template<nat _Pos,nat _Neg,nat _Proj, typename _R, bool _Cospace> nat Space<_Pos,_Neg,_Proj,_R,_Cospace>::Components(nat nGrade)
{
  nat nComb=Dim;
  // This is just amount of combinations nGrade from _Dim.
  if (nGrade > Dim)
    return 0;
  nat q=nGrade;
  if (2*q > Dim)
    q=Dim-q;
  if (!q)
    return 1;

  for (nat i=1; i<q; i++)
    nComb=nComb*(Dim-i)/(i+1);

  return nComb;
}
///@}
/**
  @addtogroup belem Basis elements.
  @ingroup geom
*/
///@{
/*===========================================================================*\
*/ /**
  Basis elements of geometric algebra.

  @tparam Sp Space of basis element.
*/
template<typename Sp> class BElem
{
public:
  /** Default constructor. */
  BElem():
    m_nRepr(0)
  {;}
  /**
  Constructor.

  @param nBits Bit representation of basis element.
  */
  BElem(nat nBits):
    m_nBits(nBits & Sp::I),
    m_nGrade(Sp::Grade(m_nBits))
  {;}
  /**
  Constructor.

  @param nGrade Grade of basis element.
  @param idx Number of basis element of this grade.
  */
  BElem(nat nGrade, nat idx):
    m_nBits(Sp::Bits(nGrade,idx)),
    m_nGrade(Sp::Grade(m_nBits))
  {;}
  /** Get bit representation of basis element. */
  nat Bits() const {return m_nBits;}
  /** Get grade of basis element. */
  nat Grade() const {return m_nGrade;}
  /**
  Conversion to natural number to be used for sorting.

  @return Combined m_nBits and m_nGrade.
  */
  operator const nat&() const {return m_nRepr;}

private:
  union
  {
    nat m_nRepr;                            ///< Full representation.
    struct
    {
      nat m_nBits: Sp::Dim;                 ///< Bit representation of basis element.
      nat m_nGrade: 8*sizeof(nat)-Sp::Dim;  ///< Grade of basis element.
    };
  };
};
///@}
/**
  @addtogroup mvect Multi-vectors.
  @ingroup geom
 */
///@{
/*===========================================================================*\
*/ /**
  Representation of multi-vector.

  @tparam Sp Space to which belongs multi-vector.
*/
template<typename Sp> class MVect: public std::map<BElem<Sp>,typename Sp::R>
{
public:
  using R = Sp::R;            ///< Type of real values.
  /** Default constructor. */
  MVect() = default;
  /**
  Constructs and initializes multi-vector from basis element.

  @param be Basis element.
  @param r Optional factor for basis element.
  */
  MVect(const BElem<Sp>& be, R r=1.0) {(*this)[be]=r;}
  /**
  Constructs and initializes multi-vector of certain grade.

  @param nGrade Grade of elements.
  @param lst Initialization list for elements of this grade.
  */
  MVect(nat nGrade, const std::initializer_list<R>& lst) { AddGrade(nGrade,lst);}
  /** Multi-vector representation in co-space. */
  MVect<typename Sp::Cospace> CoMVect() const;
  /**
  Adds components of certain grade. If elements already exist - they are replaced.

  @param nGrade Grade of elements.
  @param lst Initialization list for elements of this grade.
  */
  void AddGrade(nat nGrade, const std::initializer_list<R>& lst);
  /**
  Picks up components of certain grade.

  @param nGrade Grade of elements.
  @return Multi-vector, containing components of certain grade.
  */
  MVect Grade(nat nGrade) const;
  /** Calculates inverse multi-vector. */
  MVect Inv() const;
  /** Calculates hat involution of multi-vector. */
  MVect Hat() const;
  /** Calculates reverse of multi-vector. */
  MVect Rev() const;
  /** Calculates Clifford conjugation of multi-vector. */
  MVect Conj() const;
  /**
  @name Addition and multiplication operators.

  @param r Scaling factor.
  @param mv1 Multi-vector operand.
  */
///@{
  /** Calculates scaled multi-vector. */
  MVect& operator *=(R r);
  /** Calculates scaled multi-vector. */
  MVect& operator /=(R r);
  /** Calculates sum with multi-vector. */
  MVect& operator +=(const MVect& mv1);
  /** Calculates difference with multi-vector. */
  MVect& operator -=(const MVect& mv1);
  /** Calculates interior product with multi-vector. */
  MVect& operator |=(const MVect& mv1)
  {
    *this=(*this)|mv1;
    return (*this);
  }
  /** Calculates exterior product with multi-vector. */
  MVect& operator ^=(const MVect& mv1)
  {
    *this=(*this)^mv1;
    return (*this);
  }
  /** Calculates geometric product with multi-vector. */
  MVect& operator *=(const MVect& mv1)
  {
    *this=(*this)*mv1;
    return (*this);
  }
  /** Calculates division by multi-vector. */
  MVect& operator /=(const MVect& mv1)
  {
    return (*this)*=mv1.Inv();
  }
///@}
  /** Removes zero elements. */
  void Cleanup(R s)
  {
    for (auto it=this->begin(); it != this->end();)
    {
      if (!CmpReal0(it->second,s))
        it=this->erase(it);
      else
        ++it;
    }
  }
  /**
  @name Setters for certain multi-vector components.

  @param val Value to assign.
  */
///@{
  /** Point weight in projective space. */
  void w(real val) requires (Sp::Proj == 1) {(*this)[BElem<Sp>(Sp::vect,Sp::weight)]=val;}
  /** Positive weight in conformal space. */
  void p(real val) requires (Sp::Proj > 1) {(*this)[BElem<Sp>(Sp::vect,Sp::weight)]=val;}
  /** Negative weight in conformal space. */
  void n(real val) requires (Sp::Proj > 1) {(*this)[BElem<Sp>(Sp::vect,Sp::weight+1)]=val;}
  /** Time in Minkovsky space or equivalent to x coordinate. */
  void t(real val) {(*this)[BElem<Sp>(Sp::vect,Sp::time)]=val;}
  /** X coordinate. */
  void x(real val) {(*this)[BElem<Sp>(Sp::vect,Sp::vect0)]=val;}
  /** Y coordinate. */
  void y(real val) {(*this)[BElem<Sp>(Sp::vect,Sp::vect0+1)]=val;}
  /** Z coordinate. */
  void z(real val) {(*this)[BElem<Sp>(Sp::vect,Sp::vect0+2)]=val;}
  /**
  Arbitrary vector component.

  @param idx Component index.
  */
  void v(nat idx, real val) {(*this)[BElem<Sp>(Sp::vect,Sp::vect0+idx)]=val;}
///@}
  /**
  @name Returns value of certain multi-vector component.
  */
///@{
  /** Point weight in projective space. */
  R w() const requires (Sp::Proj == 1) {auto it=this->find(BElem<Sp>(Sp::vect,Sp::weight)); return (it == this->cend()? 0 : it->second);}
  /** Point weight in conformal space. */
  R w() const requires (Sp::Proj > 1) {return p()+n();}
  /** Positive weight in conformal space. */
  R p() const requires (Sp::Proj > 1) {auto it=this->find(BElem<Sp>(Sp::vect,Sp::weight)); return it == this->cend()? 0 : it->second;}
  /** Negative weight in conformal space. */
  R n() const requires (Sp::Proj > 1) {auto it=this->find(BElem<Sp>(Sp::vect,Sp::weight+1)); return it == this->cend()? 0 : it->second;}
  /** Weight at infinity in conformal space. */
  R i() const requires (Sp::Proj > 1) {return 0.5*(n()-p());}
  /** Time in Minkovsky space. */
  R t() const {auto it=this->find(BElem<Sp>(Sp::vect,Sp::time)); return it == this->cend()? 0 : it->second; }
  /** X coordinate. */
  R x() const {auto it=this->find(BElem<Sp>(Sp::vect,Sp::vect0)); return it == this->cend()? 0 : it->second; }
  /** Y coordinate. */
  R y() const {auto it=this->find(BElem<Sp>(Sp::vect,Sp::vect0+1)); return it == this->cend()? 0 : it->second; }
  /** Z coordinate. */
  R z() const {auto it=this->find(BElem<Sp>(Sp::vect,Sp::vect0+2)); return it == this->cend()? 0 : it->second; }
  /**
  Arbitrary vector component.

  @param idx Component index.
  */
  R v(nat idx) const {auto it=this->find(BElem<Sp>(Sp::vect,Sp::vect0+idx)); return it == this->cend()? 0 : it->second; }
  /** Redefinition of base class operator to return zero for non-existing basis element. */
  R operator [](const BElem<Sp>& be) const {auto it=this->find(be); return it == this->cend()? 0 : it->second; }
  /** Redefinition of base class operator to return zero for non-existing basis element. */
  R& operator [](const BElem<Sp>& be) {return std::map<BElem<Sp>,R>::operator [](be);}
///@}
};
/*===========================================================================*\
*/ /**
  Multi-vector representation in co-space.

  @tparam Sp Space of multi-vectors.
  @return Representation of this multi-vector in co-space basis.
  @see    MVect<Sp>
*/
template<typename Sp> MVect<typename Sp::Cospace> MVect<Sp>::CoMVect() const
{
  MVect<typename Sp::Cospace> mv;
  for (auto it=this->cbegin(); it != this->cend(); ++it)
  {
    int nProd=1;
    nat nBit=0;
    nat nBits=it->first.Bits();
    for (nat nMask=1; nMask<=nBits; ++nBit,nMask <<= 1)
    {
      if (nMask & nBits)
        nProd *= Sp::Sign(nBit);
    }
    mv[BElem<typename Sp::Cospace>(nBits)]=nProd*it->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Adds components of certain grade. If elements already exist - they are replaced.

  @tparam Sp Space of multi-vectors.
  @param nGrade Grade of elements.
  @param lst Initialization list for elements of this grade.
  @see    MVect<Sp>
*/
template<typename Sp> void MVect<Sp>::AddGrade(nat nGrade, const std::initializer_list<R>& lst)
{
  R r;
  nat i=0;
  nat cnt=Sp::Components(nGrade);
  for (auto it=lst.begin(); (i < cnt) && (it != lst.end()); ++it,++i)
  {
    r=*it;
    BElem<Sp> be(nGrade,i);
    if (r != 0.0)
      (*this)[be]=r;
    else if (this->find(be) != this->end())
      this->erase(be);
  }
}
/*===========================================================================*\
*/ /**
  Picks up components of certain grade.

  @tparam Sp Space of multi-vectors.
  @param nGrade Grade of elements.
*/
template<typename Sp> MVect<Sp> MVect<Sp>::Grade(nat nGrade) const
{
  MVect<Sp> mv;
  for (auto it=this->cbegin(); (it != this->cend()) && (it->first.Grade() <= nGrade); ++it)
  {
    if (it->first.Grade() == nGrade)
      mv[it->first]=it->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates inverse multi-vector.

  @tparam Sp Space of multi-vectors.
  @return Inverse of multi-vector.
*/
template<typename Sp> MVect<Sp> MVect<Sp>::Inv() const
{
  if (this->cbegin() == this->cend())
    throw std::overflow_error("Division by empty multi-vector.");

  auto it=this->cbegin();
  if ((it->first.Grade() == this->crbegin()->first.Grade()) && ((it->first.Grade() < 2) || (it->first.Grade() > Sp::Dim-2)))
  // This is a blade.
  {
    MVect mvr=this->Rev();
    R r=((*this)*mvr)[0];
    if (r == 0.0)
      throw std::overflow_error("Division by zero blade.");
    return mvr/r;
  }
/*
  nat nOddEven=0;
  auto it=this->cbegin();
  for (; it != this->cend(); ++it)
    nOddEven |= it->first.Grade() & 1? 2 : 1;
  if (nOddEven != 3)
  // Should be a versor.
  {
//    std::cout << "V*";
    MVect mvr=Rev();
    R r=((*this)*mvr)[0];
    if (r == 0.0)
      throw std::overflow_error("Division by zero versor.");
    return mvr/r;
  }
*/
  // General solution.
  nat r,c;
  const nat s=1<<Sp::Dim;
  // Create matrix representation of multi-vector.
  tens::Matrix<typename Sp::R> m({s,-(int)s});
  for (r=0; r<s; ++r)
  {
    auto p0=BElem<Sp>(r)*BElem<Sp>(r);
    m[r][r]=(*this)[p0.first];
    for (c=r+1; c<s; c++)
    {
      auto p=BElem<Sp>(r)*BElem<Sp>(c);
      m[r][c]=p0.second*p.second*(*this)[p.first];
      p=p.first*p.first;
      m[c][r]=p.second*m[r][c];
    }
  }
  // Solve equation system to find the first column of inverse matrix.
  tens::Matrix<typename Sp::R> x({s,-1});
  x[0][0]=1.0;
  R rDet=tens::SolveGauss(m,x);
  if (rDet == 0.0)
    throw std::overflow_error("Division by zero multi-vector.");
  // Reconstruct inverse multi-vector from the first column of matrix solution.
  R sc=0.0;
  MVect<Sp> mvi;
  for (r=0; r<s; ++r)
  {
    if (x[r][0] == 0.0)
      continue;
    auto p=BElem<Sp>(r)*BElem<Sp>(r);
    mvi[r]=p.second*x[r][0];
    sc += std::abs(x[r][0]);
  }
  mvi.Cleanup(sc);
  return mvi;
}
/*===========================================================================*\
*/ /**
  Calculates hat involution of multi-vector.

  @tparam Sp Space of multi-vectors.
  @return Hat involution of multi-vector.
*/
template<typename Sp> MVect<Sp> MVect<Sp>::Hat() const
{
  nat k;
  MVect<Sp> mv;
  for (auto it=this->cbegin(); it != this->cend(); ++it)
  {
    k=it->first.Grade();
    mv[it->first]=(k & 1? -1 : 1)*it->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates reverse of multi-vector.

  @tparam Sp Space of multi-vectors.
  @return Reverted multi-vector.
*/
template<typename Sp> MVect<Sp> MVect<Sp>::Rev() const
{
  nat k;
  MVect<Sp> mv;
  for (auto it=this->cbegin(); it != this->cend(); ++it)
  {
    k=it->first.Grade();
    mv[it->first]=(k*(k-1) & 2? -1 : 1)*it->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates Clifford conjugation of multi-vector.

  @tparam Sp Space of multi-vectors.
  @return Clifford conjugation of multi-vector.
*/
template<typename Sp> MVect<Sp> MVect<Sp>::Conj() const
{
  nat k;
  MVect<Sp> mv;
  for (auto it=this->cbegin(); it != this->cend(); ++it)
  {
    k=it->first.Grade();
    mv[it->first]=(k*(k+1) & 2? -1 : 1)*it->second;
  }
  return mv;
}
///@}
/**
  @addtogroup mvectop Addition and multiplication operators.
  @ingroup mvect
  @tparam Sp Space of multi-vectors.
  @param r Multiplication factor.
  @param mv1 Multi-vector operand.
  @see    MVect<Sp>
 */
///@{
/** Calculates scaled multi-vector. */
template<typename Sp> MVect<Sp>& MVect<Sp>::operator *=(R r)
{
  if (r == 0.0)
  {
    this->clear();
    return *this;
  }
  for (auto it=this->begin(); it != this->end(); ++it)
    it->second *= r;
  return *this;
}
/** Calculates scaled multi-vector. */
template<typename Sp> MVect<Sp>& MVect<Sp>::operator /=(R r)
{
  if (r == 0.0)
    throw std::overflow_error("Multi-vector divide by zero.");
  r=1.0/r;
  for (auto it=this->begin(); it != this->end(); ++it)
    it->second *= r;
  return *this;
}
/** Calculates sum with multi-vector. */
template<typename Sp> MVect<Sp>& MVect<Sp>::operator +=(const MVect& mv1)
{
  typename Sp::R s=0.0;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    if (it1->second != 0)
    {
      (*this)[it1->first] += it1->second;
      s += std::abs(it1->second);
    }
  }
  this->Cleanup(s);
  return *this;
}
/** Calculates difference with multi-vector. */
template<typename Sp> MVect<Sp>& MVect<Sp>::operator -=(const MVect& mv1)
{
  typename Sp::R s=0.0;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    if (it1->second != 0)
    {
      (*this)[it1->first] -= it1->second;
      s += std::abs(it1->second);
    }
  }
  this->Cleanup(s);
  return *this;
}
///@}
/**
  @addtogroup belemop1 Operators on of basis elements.
  @ingroup belem
  @tparam Sp Space of multi-vectors.
  @param be Basis element.
  @see    BElem<Sp>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates right complement of basis element.

  @return  Returns right complement of basis element together with its sign.
*/
template<typename Sp> std::pair<BElem<Sp>,int> operator ~(const BElem<Sp>& be)
{
  nat nPerm,nBits2,nMask;
  nat nRes=~be.Bits() & Sp::I;

  nPerm=nBits2=0;
  for (nMask=1; nMask && (nMask<=be.Bits()); nMask<<=1)
  {
    if (nMask & be.Bits())
    // Amount of permutations with this bit.
      nPerm+=nBits2;
    else
    // Count of bits for second operand.
      nBits2++;
  }
  // Each permutation changes the sign.
  return {nRes,nPerm & 1? -1 : 1};
}
/*===========================================================================*\
*/ /**
  Calculates left complement of basis element.

  @return  Returns left complement of basis element together with its sign.
*/
template<typename Sp> std::pair<BElem<Sp>,int> operator !(const BElem<Sp>& be)
{
  nat nPerm,nBits2,nMask;
  nat nRes=~be.Bits() & Sp::I;

  nPerm=nBits2=0;
  for (nMask=1; nMask && (nMask<=nRes); nMask<<=1)
  {
    if (nMask & nRes)
    // Amount of permutations with this bit.
      nPerm+=nBits2;
    else
    // Count of bits for second operand.
      nBits2++;
  }
  // Each permutation changes the sign.
  return {nRes,nPerm & 1? -1 : 1};
}
/*===========================================================================*\
*/ /**
  Outputs basis element to the stream.

  @param be Basis element.
  @param os Output stream.
  @return  Returns output stream.
*/
template<typename Sp> std::ostream& operator <<(std::ostream& os, const BElem<Sp>& be)
{
  nat n=0;
  nat nBits=be.Bits();
  if (!nBits)
    return os << "1";
  if (nBits == Sp::I)
    return os << "I" << Sp::Dim;
  nat nMask=1;
  if constexpr (Sp::Proj == 1)
  {
    if (nBits & 1)
      os << "o";
    nMask=2;
  }
  if constexpr (Sp::Proj > 1)
  {
    if (nBits & 1)
      os << "p";
    if (nBits & 2)
      os << "n";
    nMask=4;
  }
  if (nMask > nBits)
    return os;
  os << "e";
  for (; nMask <= nBits; ++n,nMask <<= 1)
  {
    if (!(nBits & nMask))
      continue;
    os << n;
  }
  return os;
}
///@}
/**
  @addtogroup belemop2 Products of basis elements.
  @ingroup belem
  @tparam Sp Space of multi-vectors.
  @tparam Cosp Space OR co-space of multi-vectors.
  @param be1,be2 Operands.
  @see    BElem<Sp>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates external product of basis elements.

  @return  Returns external product of basis elements together with sign of the product.
*/
template<typename Sp,typename Cosp> inline std::pair<BElem<Sp>,int> operator ^(const BElem<Cosp>& be1, const BElem<Sp>& be2)
{
  if (be1.Bits() & be2.Bits())
  // External product of any element with itself is 0.
    return {Sp::l,0};
  // Otherwise identical to geometric product.
  return be1*be2;
}
/*===========================================================================*\
*/ /**
  Calculates inner product of basis elements.

  @return  Returns inner product of basis elements.
*/
template<typename Sp,typename Cosp> inline std::pair<BElem<Sp>,int> operator |(const BElem<Cosp>& be1, const BElem<Sp>& be2)
{
  if (!(be1.Bits() & be2.Bits()))
  // Inner product of independent elements is 0.
    return {Sp::l,0};

  // Otherwise identical to geometric product.
  return be1*be2;
}
/*===========================================================================*\
*/ /**
  Calculates geometric product of basis elements.

  @return  Returns geometric product of basis elements together with sign of the product.
*/
template<typename Sp> std::pair<BElem<Sp>,int> operator *(const BElem<Sp>& be1, const BElem<Sp>& be2)
{
  nat nRes,nPerm,nBit,nBits2,nMask;
  int nProd=1;

  // We have to check all basis elements.
  nRes=be1.Bits()|be2.Bits();

  nPerm=nBit=nBits2=0;
  for (nMask=1; nMask && (nMask<=nRes); nMask<<=1)
  {
    if (nMask & be1.Bits())
    {
      // In any case permutations should be counted.
      //   Amount of permutations with this bit.
      nPerm+=nBits2;
      if (nMask & be2.Bits())
      // This basis element will not exist in product.
      {
        nRes &= ~nMask;
        // This bit should not be permuted with current, but with following bits.
        nBits2++;
        nProd *= Sp::Sign(nBit);
      }
    }
    else if (nMask & be2.Bits())
    // Count of bits for second operand.
      nBits2++;
    nBit++;
  }
  // Each permutation changes the sign.
  return {nRes,nProd*(nPerm & 1? -1.0 : 1.0)};
}
/*===========================================================================*\
*/ /**
  Calculates geometric product of co-space basis element with space basis element.

  @return  Returns geometric product of basis elements together with sign of the product.
*/
template<typename Sp> std::pair<BElem<Sp>,int> operator *(const BElem<typename Sp::Cospace>& be1, const BElem<Sp>& be2)
{
  nat nRes,nPerm,nBit,nBits2,nMask;
  int nProd=1;

  // We have to check all basis elements.
  nRes=be1.Bits()|be2.Bits();

  nPerm=nBit=nBits2=0;
  for (nMask=1; nMask && (nMask<=nRes); nMask<<=1)
  {
    if (nMask & be1.Bits())
    {
      // In any case permutations should be counted.
      //   Amount of permutations with this bit.
      nPerm+=nBits2;
      if (nMask & be2.Bits())
      // This basis element will not exist in product.
      {
        nRes &= ~nMask;
        // This bit should not be permuted with current, but with following bits.
        nBits2++;
      }
      else
        nProd *= Sp::Cospace::Sign(nBit);
    }
    else if (nMask & be2.Bits())
    // Count of bits for second operand.
      nBits2++;
    nBit++;
  }
  // Each permutation changes the sign.
  return {nRes,nProd*(nPerm & 1? -1.0 : 1.0)};
}
///@}
/**
  @addtogroup mvectop1 Operators on multi-vectors.
  @ingroup mvect
  @tparam Sp Space of multi-vectors.
  @param mv1 Multi-vector.
  @see    MVect<Sp>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates right complement of multi-vector.

  @return  Returns right complement of multi-vector.
*/
template<typename Sp> MVect<Sp> operator ~(const MVect<Sp>& mv1)
{
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    auto p=~it1->first;
    mv[p.first]=p.second*it1->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates left complement of multi-vector.

  @return  Returns left complement of multi-vector.
*/
template<typename Sp> MVect<Sp> operator !(const MVect<Sp>& mv1)
{
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    auto p=!it1->first;
    mv[p.first]=p.second*it1->second;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates negative multi-vector.

  @return  Returns negative multi-vector.
*/
template<typename Sp> MVect<Sp> operator -(const MVect<Sp>& mv1)
{
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
    mv[it1->first]=-it1->second;
  return mv;
}
/*===========================================================================*\
*/ /**
  Outputs multi-vector to the stream.

  @param mv1 Multi-vector.
  @param os Output stream.
  @return  Returns output stream.
*/
template<typename Sp> std::ostream& operator <<(std::ostream& os, const MVect<Sp>& mv1)
{
  for (auto it=mv1.cbegin(); it != mv1.cend(); ++it)
  {
    if ((it != mv1.cbegin()) && (it->second >= 0.0))
      os << "+";
    os << it->second << "*" << it->first;
  }
  return os;
}
///@}
/**
  @addtogroup mvects Scaling of multi-vectors.
  @ingroup mvect
  @tparam Sp Space of multi-vectors.
  @param mv1 Multi-vector.
  @param r Scalar.
  @see    MVect<Sp>
 */
///@{
/** Calculates scaled multi-vector. */
template<typename Sp> MVect<Sp> inline operator *(const MVect<Sp>& mv1, typename Sp::R r)
{
  MVect<Sp> mv(mv1);
  mv *= r;
  return mv;
}
/** Calculates scaled multi-vector. */
template<typename Sp> inline MVect<Sp> operator *(typename Sp::R r, const MVect<Sp>& mv1)
{
  MVect<Sp> mv(mv1);
  mv *= r;
  return mv;
}
/** Calculates divided multi-vector. */
template<typename Sp> inline MVect<Sp> operator /(const MVect<Sp>& mv1, typename Sp::R r)
{
  MVect<Sp> mv(mv1);
  mv /= r;
  return mv;
}
///@}
/**
  @addtogroup mvectop2 Sums and products of multi-vectors.
  @ingroup mvect
  @tparam Sp Space of multi-vectors.
  @tparam Cosp Space OR co-space of multi-vectors.
  @param mv1,mv2 Multi-vector operands.
  @see    MVect<Sp>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates interior product of multi-vectors.

  @return  Returns interior product.
*/
template<typename Sp,typename Cosp> MVect<Sp> operator |(const MVect<Cosp>& mv1, const MVect<Sp>& mv2)
{
  typename Sp::R r;
  typename Sp::R s=0.0;
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    for (auto it2=mv2.cbegin(); it2 != mv2.cend(); ++it2)
    {
      auto p=it1->first|it2->first;
      r=p.second*it1->second*it2->second;
      if (r != 0.0)
      {
        mv[p.first] += r;
        s += std::abs(r);
      }
    }
  }
  mv.Cleanup(s);
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates exterior product of multi-vectors.

  @return  Returns exterior product.
*/
template<typename Sp,typename Cosp> MVect<Sp> operator ^(const MVect<Cosp>& mv1, const MVect<Sp>& mv2)
{
  typename Sp::R r;
  typename Sp::R s=0.0;
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    for (auto it2=mv2.cbegin(); it2 != mv2.cend(); ++it2)
    {
      auto p=it1->first^it2->first;
      r=p.second*it1->second*it2->second;
      if (r != 0.0)
      {
        mv[p.first] += r;
        s += std::abs(r);
      }
    }
  }
  mv.Cleanup(s);
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates geometric product of multi-vectors.

  @return  Returns geometric product.
*/
template<typename Sp,typename Cosp> MVect<Sp> operator *(const MVect<Cosp>& mv1, const MVect<Sp>& mv2)
{
  typename Sp::R r;
  typename Sp::R s=0.0;
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    for (auto it2=mv2.cbegin(); it2 != mv2.cend(); ++it2)
    {
      auto p=it1->first*it2->first;
      r=p.second*it1->second*it2->second;
      if (r != 0.0)
      {
        mv[p.first] += r;
        s += std::abs(r);
      }
    }
  }
  mv.Cleanup(s);
  return mv;
}
/** Calculates sum of multi-vectors. */
template<typename Sp> inline MVect<Sp> operator +(const MVect<Sp>& mv1, const MVect<Sp>& mv2)
{
  MVect<Sp> mv(mv1);
  mv += mv2;
  return mv;
}
/** Calculates difference of multi-vectors. */
template<typename Sp> inline MVect<Sp> operator -(const MVect<Sp>& mv1, const MVect<Sp>& mv2)
{
  MVect<Sp> mv(mv1);
  mv -= mv2;
  return mv;
}
/** Divides arbitrary multi-vectors. */
template<typename Sp,typename Cosp> inline MVect<Sp> operator /(const MVect<Cosp>& mv1, const MVect<Sp>& mv2)
{
  return mv1*mv2.Inv();
}
/** Calculates Grassmann regressive product (meet). */
template<typename Sp> inline MVect<Sp> Meet(const MVect<Sp>& mv1, const MVect<Sp>& mv2)
{
  return ~((~mv1)^(~mv2));
}
/*===========================================================================*\
*/ /**
  Calculates Grassmann regressive product (meet), using alternative dual spaces.

  @param mv1,mv2 Multi-vector operands.
  @param dual Method to calculate dual basis elements.
  @return  Returns regressive product (meet).
*/
template<typename Sp> inline MVect<Sp> Meet(const MVect<Sp>& mv1, const MVect<Sp>& mv2, auto dual)
{
  return dual(dual(mv1)^dual(mv2));
}
///@}
/**
  @addtogroup mvectlin Linear transformations of multi-vectors.
  @ingroup mvect
  @tparam Sp Space of multi-vectors.
  @param m Linear transformation matrix.
  @param mv1,mv2 Multi-vector operands.
  @see    MVect<Sp>
 */
///@{
/*===========================================================================*\
*/ /**
  Calculates linear transformation of basis element.

  @param m Linear transformation matrix.
  @param be Basis element.
  @return  Returns multi-vector of transformed basis element.
*/
template<typename Sp> MVect<Sp> operator *(const tens::Matrix<typename Sp::R>& m,const BElem<Sp>& be)
{
  if ((m.DSize(0) != Sp::Dim) || (m.DSize(1) != Sp::Dim))
    throw std::domain_error("Wrong matrix dimensions in operator *(Matrix,BElem).");
  MVect<Sp> mv(0);
  nat r,nMask;
  nat c=0;
  for (nMask=1; nMask && (nMask<=be.Bits()); ++c,nMask<<=1)
  {
    if (!(nMask & be.Bits()))
      continue;
    MVect<Sp> me;
    for (r=0; r<Sp::Dim; ++r)
    {
      if (m[r][c] != 0.0)
        me[(nat)1<<r]=m[r][c];
    }
    mv=mv^me;
  }
  return mv;
}
/*===========================================================================*\
*/ /**
  Calculates linear transformation of multi-vector.

  @param m Linear transformation matrix.
  @param mv1 Multi-vector operand.
  @return  Returns transformed multi-vector.
*/
template<typename Sp> MVect<Sp> operator *(const tens::Matrix<typename Sp::R>& m,const MVect<Sp>& mv1)
{
  if ((m.DSize(0) != Sp::Dim) || (m.DSize(1) != Sp::Dim))
    throw std::domain_error("Wrong matrix dimensions in operator *(Matrix,MVect).");
  MVect<Sp> mv;
  for (auto it1=mv1.cbegin(); it1 != mv1.cend(); ++it1)
  {
    if (it1->second == 0.0)
      continue;
    MVect<Sp> me=m*it1->first;
    mv=mv+it1->second*me;
  }
  return mv;
}
///@}
/*===========================================================================*\
*/ /**
  @addtogroup duals Dual multi-vectors.
  @ingroup mvect

  @param mv Multi-vector operand.
*/
///@{
/** Grassmann right complement ~mv. */
template<typename Sp> inline MVect<Sp> DualR(const MVect<Sp>& mv) {return ~mv;}
/** Grassmann left complement !mv. */
template<typename Sp> inline MVect<Sp> DualL(const MVect<Sp>& mv) {return !mv;}
/** Right dual mv*I. */
template<typename Sp> inline MVect<Sp> DualIR(const MVect<Sp>& mv) {return mv*MVect<Sp>(Sp::I);}
/** Left dual I*mv. */
template<typename Sp> inline MVect<Sp> DualIL(const MVect<Sp>& mv) {return MVect<Sp>(Sp::I)*mv;}
///@}
};
};
#endif

