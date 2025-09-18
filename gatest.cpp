/*===========================================================================*\
*/ /**
  @file gatest.cpp

  Tests for geometric algebra library.

  @date Created: 18.04.2024
  @date Last Revision:

  @author E.Yanenko
*/
#include <chrono>
#include <iostream>
#include "elib/geom.hpp"

using namespace ey;
using namespace ey::ga;
using namespace ey::tens;

typedef std::chrono::duration<float> float_seconds;
// Range for multi-vector components values.
constexpr int c_nRange=1000;

/*===========================================================================*\
*/ /**
  Create random multi-vector, containing specified grades.

  @param nGrades Bit-field of necessary grades, default all.
  @param bSpace For vectors create only space components.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> MVect<Sp> RandomMVect(nat nGrades = Sp::I, bool bSpace=true)
{
  MVect<Sp> mv;
  // For vectors we will need only space components.
  nat i0=bSpace && (nGrades == (1<<Sp::vect))? Sp::Proj : 0;

  for (nat nGrade=0; nGrade <= Sp::Dim; nGrade++)
  {
    // Skip unnecessary grades.
    if (!(((nat)1<<nGrade) & nGrades))
      continue;

    nat cnt=Sp::Components(nGrade);
    for (nat i=i0; i<cnt; i++)
    {
      int r=std::rand()%(2*c_nRange+1)-c_nRange;
      // Do not set 0 values, non-existing components are anyhow zero.
      if (r)
        mv[BElem<Sp>(nGrade,i)]=r;
    }
  }
  return mv;
}

/*===========================================================================*\
*/ /**
  Create random blade of specified grade.

  @param nGrade Necessary grade.
  @param bSpace For vectors create only space components.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> MVect<Sp> RandomBlade(nat nGrade, bool bSpace=false)
{
  if (nGrade > Sp::Dim)
    nGrade=Sp::Dim;
  if (nGrade < 2)
    return RandomMVect<Sp>((nat)1<<nGrade,bSpace);

  // Random vector.
  MVect<Sp> mv=RandomMVect<Sp>((nat)1<<Sp::vect,bSpace);
  for (nat n=2; n <= nGrade; n++)
    mv ^= RandomMVect<Sp>((nat)1<<Sp::vect,bSpace);
  return mv;
}

/*===========================================================================*\
*/ /**
  Create random point in projection space.

  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> inline MVect<Sp> RandomPoint() requires (Sp::Proj == 1)
{
  MVect<Sp> mv=RandomMVect<Sp>((nat)1<<Sp::vect,true);
  mv.w(1.0);
  return mv;
}

/*===========================================================================*\
*/ /**
  Create random point in conformal space.

  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> MVect<Sp> RandomPoint() requires (Sp::Proj == 2)
{
  MVect<Sp> mv;
  real s2=0.0;
  nat cnt=Sp::Components(Sp::vect);
  for (nat i=Sp::Proj; i<cnt; i++)
  {
    int r=std::rand()%(2*c_nRange+1)-c_nRange;
    // Do not set 0 values, non-existing components are anyhow zero.
    if (r)
    {
      mv[BElem<Sp>(Sp::vect,i)]=r;
      s2+=r*r;
    }
  }
  mv.p(0.5-0.5*s2);
  mv.n(0.5+0.5*s2);
  return mv;
}

/*===========================================================================*\
*/ /**
  Create random matrix.

  @param nDim Amount of matrix rows and columns.
  @tparam R Matrix value type.
*/
template<typename R> Matrix<R> RandomMatrix(int nDim)
{
  Matrix<R> m({nDim,-nDim});
  for (int r=0; r < nDim; r++)
  {
    for (int c=0; c < nDim; c++)
    {
      real v=2.0*(std::rand()%(2*c_nRange+1)-c_nRange)/c_nRange;
      m[r][c]=v;
    }
  }
  return m;
}

/*===========================================================================*\
*/ /**
  Check correctness basis elements.

  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestBasis()
{
  nat nSum=0;
  for (nat nGrade=0; nGrade <= Sp::Dim; nGrade++)
  {
    // Amount of basis elements of a certian grade.
    nat cnt=Sp::Components(nGrade);
    for (nat i=0; i<cnt; i++)
    {
      // Create corresponding basis element.
      BElem<Sp> be(nGrade,i);
      // Check that correct index is calculated.
      nat idx=Sp::Index(be.Bits());
      if (i != idx)
      {
        std::cout << "Grade " << nGrade << ": " << i << " != " << idx << std::endl;
        throw std::runtime_error("Error basis element");
      }
    }
    nSum += cnt;
  }
  // Check that total amount of basis elements is correct.
  if (nSum != ((nat)1<<Sp::Dim))
  {
    std::cout << "Total: " << nSum << " != " << ((nat)1<<Sp::Dim) << std::endl;
    throw std::runtime_error("Error basis count");
  }
  std::cout << "Basis elements: " << nSum << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Test blade products.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestBlades(nat nTests = 2000)
{
  nat n;
  for (n=0; n<nTests; n++)
  {
    // Try all blade combinations.
    for (nat g1=0; g1<=Sp::Dim; g1++)
    {
      MVect<Sp> b1=RandomBlade<Sp>(g1, false);
      for (nat g2=0; g2<=Sp::Dim; g2++)
      {
        MVect<Sp> b2=RandomBlade<Sp>(g2, false);
        MVect<Sp> mvg=b1*b2;
        MVect<Sp> mvi=b1|b2;
        MVect<Sp> mve=b1^b2;
        // Check if the sum of interior and exterior products is equal to geometric product.
        if (mvg != mvi+mve)
        {
          std::cout << "A" << g1 << "=(" << b1 << ") B" << g2 << "=(" << b2 << ")" << std::endl;
          std::cout << "(" << mvi << ")+(" << mve << ") != (" << mvg << ")" << std::endl;
          throw std::runtime_error("Error blade product Ai|Bj + Ai^Bj != Ai*Bj");
        }
      }
    }
  }
  std::cout << "Blade products: " << n*(Sp::Dim+1)*(Sp::Dim+1) << " -> OK" <<std::endl;
}

/*===========================================================================*\
*/ /**
  Test blade inverse.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestBladeInverse(nat nTests = 2000)
{
  // Amount of multi-vector components.
  const nat s=1<<Sp::Dim;
  real rThr=(real)s*c_nRange*c_nRange;

  nat n;
  nat cnt=0;
  real rSum=0;
  real rMax=0;
  for (n=0; n<nTests; n++)
  {
    // Try all blades.
    for (nat g1=0; g1<=Sp::Dim; g1++)
    {
      MVect<Sp> b1=RandomBlade<Sp>(g1, false);

      // Check if blade reverse is correct.
      if (b1*b1.Rev() != b1.Rev()*b1)
      {
        std::cout << "A" << g1 << "=(" << b1 << ")" << std::endl;
        std::cout << "A" << g1 << ".Rev()=(" << b1.Rev() << ")" << std::endl;
        std::cout << "A" << g1 << "*A" << g1 << ".Rev()=(" << b1*b1.Rev() << ")" << std::endl;
        std::cout << "A" << g1 << ".Rev()*A" << g1 << "=(" << b1.Rev()*b1 << ")" << std::endl;
        throw std::runtime_error("Error blade reverse Ai*Ai.Rev() != Ai.Rev()*Ai");
      }
      MVect<Sp> mvInv;
      try
      {
        mvInv=b1.Rev()/(b1*b1.Rev())[0];
      }
      catch (...)
      {
        // Zero multi-vector, skip.
        continue;
      }
      // Check multi-vector is not near 0.
      real rTest=0;
      for (auto it=mvInv.begin(); it != mvInv.end(); ++it)
        rTest = std::max(rTest,std::abs(it->second));
      if (rTest > 0.5)
      {
        // Skip this one.
        n--;
        continue;
      }
      // Check if product is equal 1.
      MVect<Sp> d=mvInv*b1;
      // Calculate rounding error.
      real err=std::abs(d[0]-1.0);
      for (auto it=++d.begin(); it != d.end(); ++it)
        err += std::abs(it->second);
      cnt++;
      rSum += err;
      rMax=std::max(rMax, err);
      if (CmpReal0(err,rThr))
      {
        std::cout << "M=(" << b1 << ")" << std::endl;
        std::cout << "Inv=(" << mvInv << ")" << std::endl;
        std::cout << "M*Inv=(" << d << ")" << std::endl;
        std::cout << "err=" << err << " > thr=" << std::numeric_limits<real>::epsilon()*rThr << std::endl;
        throw std::runtime_error("Error blade inverse");
      }
    }
  }
  std::cout << "Blade inversions: " << cnt <<", error av/max: " << rSum/cnt << "/" << rMax << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Test products.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestProducts(nat nTests = 2000)
{
  nat n;
  for (n=0; n<nTests; n++)
  {
    MVect<Sp> mv1=RandomMVect<Sp>();
    MVect<Sp> mv2=RandomMVect<Sp>();

    MVect<Sp> mvg=mv1*mv2;
    MVect<Sp> mvi=mv1|mv2;
    MVect<Sp> mve=mv1^mv2;

    // Check if the sum of interior and exterior products is equal to geometric product.
    if (mvg != mvi+mve)
    {
      std::cout << "==== Error A|B + A^B != A*B ====" << std::endl;
      std::cout << "A=(" << mv1 << ") B=(" << mv2 << ")" << std::endl;
      std::cout << "(" << mvi << ")+(" << mve << ") != (" << mvg << ")" << std::endl;
      throw std::runtime_error("Error product A|B + A^B != A*B");
    }
  }
  std::cout << "Multi-vector products: " << n << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Test regressive products.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestMeet(nat nTests = 2000)
{
  nat n;
  MVect<Sp> I(Sp::I);
  for (n=0; n<nTests; n++)
  {
    MVect<Sp> mv1=RandomMVect<Sp>();
    MVect<Sp> mv2=RandomMVect<Sp>();
    // In odd dimension left and right complements should be equal.
    if ((Sp::Dim & 1) && (~mv1 != !mv1))
    {
      std::cout << "~A=(" << ~mv1 << ") !A=(" << !mv1 << ")" << std::endl;
      throw std::runtime_error("Error odd dimension ~A != !A");
    }
    // Double left and double right complements should be equal to original wit sign of I*I.
    if (~~mv1 != !!mv1)
    {
      std::cout << "~~A=(" << ~~mv1 << ") !!A=(" << !!mv1 << ")" << std::endl;
      throw std::runtime_error("Error ~~A != !!A");
    }
    for (nat g1=0; g1<=Sp::Dim; g1++)
    {
      MVect<Sp> b1=mv1.Grade(g1);
      if (((Sp::Dim & 1) || !(g1 & 1)) && (~b1 != !b1))
      {
        std::cout << "~A" << g1 << "=(" << ~b1 << ") !A" << g1 << "=(" << !b1 << ")" << std::endl;
        throw std::runtime_error("Error even grade or odd dimension  ~Aj != !Aj");
      }
      for (nat g2=0; g2<=Sp::Dim; g2++)
      {
        MVect<Sp> b2=mv2.Grade(g2);
        if (((Sp::Dim & 1) || !(g2 & 1)) && (~b2 != !b2))
        {
          std::cout << "~A" << g2 << "=(" << ~b2 << ") !A" << g2 << "=(" << !b2 << ")" << std::endl;
          throw std::runtime_error("Error even grade or odd dimension  ~Aj != !Aj");
        }
        // Grassmann's regressive product with right complement ~(~A^~B).
        MVect<Sp> mv3=Meet(b1,b2);
        // Regressive product with left complement !(!A^!B).
        MVect<Sp> mv4=Meet(b1,b2,DualL<Sp>);
        // Regressive product using right dual space ((A*I)^(B*I))*I.
        MVect<Sp> mv5=Meet(b1,b2,DualIR<Sp>);
        // Regressive product using left dual space I*((I*A)^(I*B)).
        MVect<Sp> mv6=Meet(b1,b2,DualIL<Sp>);
        // Check if they are all equal.
        if (mv3 != mv4)
        {
          std::cout << "A" << g1 << "=(" << b1 << ") B" << g2 << "=(" << b2 << ")" << std::endl;
          std::cout << "(" << mv3 << ") != (" << mv4 << ")" << std::endl;
          throw std::runtime_error("Error regressive product ~(~Aj^~Bj) != !(!Aj^!Bj)");
        }
        if (mv5 != mv6)
        {
          std::cout << "A" << g1 << "=(" << b1 << ") B" << g2 << "=(" << b2 << ")" << std::endl;
          std::cout << "(" << mv5 << ") != (" << mv6 << ")" << std::endl;
          throw std::runtime_error("Error regressive product ((Aj*I)^(Bj*I))*I != I*((I*Aj)^(I*Bj))");
        }
        if ((mv3 != mv5) && (mv3 != -mv5))
        {
          std::cout << "A" << g1 << "=(" << b1 << ") B" << g2 << "=(" << b2 << ")" << std::endl;
          std::cout << "(" << mv3 << ") != (" << mv5 << ")" << std::endl;
          std::cout << "I*I = " << "(" << I*I << ")" << std::endl;
          throw std::runtime_error("Error regressive product ~(~Aj^~Bj) != ((Aj*I)^(Bj*I))*I");
        }
      }
    }
  }
  std::cout << "Regressive products: " << n*(Sp::Dim+1)*(Sp::Dim+1) << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Test multi-vector inversion.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestInverse(nat nTests = 2000)
{
  // Amount of multi-vector components.
  const nat s=1<<Sp::Dim;
  real rThr=(real)s*s*c_nRange*c_nRange;

  nat n;
  nat cnt=0;
  real rSum=0;
  real rMax=0;
  for (n=0; n<nTests; n++)
  {
    MVect<Sp> mv=RandomMVect<Sp>();
    MVect<Sp> mvInv;
    try
    {
      mvInv=mv.Inv();
    }
    catch (...)
    {
      // Zero multi-vector, skip.
      continue;
    }
    // Check multi-vector is not near 0.
    real rTest=0;
    for (auto it=mvInv.begin(); it != mvInv.end(); ++it)
      rTest = std::max(rTest,std::abs(it->second));
    if (rTest > 0.5)
    {
      // Skip this one.
      n--;
      continue;
    }
    // Check if product is equal 1.
    MVect<Sp> d=mvInv*mv;
    // Calculate rounding error.
    real err=std::abs(d[0]-1.0);
    for (auto it=++d.begin(); it != d.end(); ++it)
      err += std::abs(it->second);
    cnt++;
    rSum += err;
    rMax=std::max(rMax, err);
    if (CmpReal0(err,rThr))
    {
      std::cout << "M=(" << mv << ")" << std::endl;
      std::cout << "Inv=(" << mvInv << ")" << std::endl;
      std::cout << "M*Inv=(" << d << ")" << std::endl;
      std::cout << "err=" << err << " > thr=" << std::numeric_limits<real>::epsilon()*rThr << std::endl;
      throw std::runtime_error("Error inverse");
    }
  }
  std::cout << "Inversions: " << cnt <<", error av/max: " << rSum/cnt << "/" << rMax << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Test linear transformations.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestLinear(nat nTests = 20000)
{
  // Amount of multi-vector components.
  const nat s=1<<Sp::Dim;
  real rThrD=(real)s*s;
  real rThrL=s;

  nat n;
  real rDet;
  real rSumD=0;
  real rMaxD=0;
  real rSumL=0;
  real rMaxL=0;
  real errD;
  real errL;
  for (n=0; n<nTests; n++)
  {
    Matrix<typename Sp::R> m=RandomMatrix<typename Sp::R>(Sp::Dim);
    // Linear transformation of pseudo-scalar is always pseudo-scalar with factor equal to matrix determinant.
    MVect<Sp> d=m*BElem<Sp>(Sp::I);
    // Note: matrix with 0 columns, we need only determinant.
    Matrix<typename Sp::R> mi({Sp::Dim,0});
    // Find determinant with Gauss method.
    rDet=SolveGauss(m,mi);

    // Compare results.
    errD=std::abs(d[Sp::I]-rDet);
    rSumD += errD;
    rMaxD=std::max(rMaxD, errD);
    if (CmpReal0(errD,(real)s*s))
    {
      std::cout << "Matr*I=(" << d << ")" << std::endl;
      std::cout << "Det=" << rDet << std::endl;
      std::cout << "err=" << errD << " > thr=" << std::numeric_limits<real>::epsilon()*rThrD << std::endl;
      throw std::runtime_error("Error determinant");
    }
    // Create random space vector.
    MVect mv=RandomMVect<Sp>((nat)1<<Sp::vect);
    // Geometric algebra transform.
    d=m*mv;
    // The same vector in matrix form.
    Vector<typename Sp::R> v({Sp::Dim});
    for (nat i=0; i<Sp::Dim; i++)
      v[i]=mv[BElem<Sp>(((nat)1<<Sp::vect)+i)];
    // Matrix multiplication.
    v=m*v;
    // Compare results.
    errL=0;
    for (nat i=0; i<Sp::Dim; i++)
      errL += std::abs(v[i] - mv[BElem<Sp>(((nat)1<<Sp::vect)+i)]);
    rSumL += errL;
    rMaxL=std::max(rMaxL, errL);
    if (CmpReal0(errL,rThrL))
    {
      std::cout << "Matr*V=(" << d << ")" << std::endl;
      std::cout << "err=" << errL << " > thr=" << std::numeric_limits<real>::epsilon()*rThrL << std::endl;
      throw std::runtime_error("Error linear transformation");
    }
  }
  std::cout << "Determinants: " << n << ", error av/max: " << rSumD/n << "/" << rMaxD << " -> OK" << std::endl;
  std::cout << "Linear transformations: " << n << ", error av/max: " << rSumL/n << "/" << rMaxL << " -> OK" << std::endl;
}

/*===========================================================================*\
*/ /**
  Convert multi-vector from e+, e- representation to o, inf.
  Used for printing to be more readable.

  @param mv Multi-vector in e+, e- representation.
  @tparam Sp Space of multi-vectors.
  @return Multi-vector in o, inf representation.
*/
template<typename Sp> MVect<Sp> Convert2OriginInfinity(const MVect<Sp>& mv) requires (Sp::Proj == 2)
{
  MVect<Sp> mvo;
  for (auto it=mv.cbegin(); it != mv.cend(); ++it)
  {
    MVect<Sp> mve;
    nat nb=it->first.Bits();
    nat nr=nb & ~3;
    // p*n = o*i
    if (!(nb & 3) || ((nb & 3) == 3))
      mve[it->first]=it->second;
    else
    {
      mve[1|nr]=it->second;
      mve[2|nr]=(nb & 1? -0.5 : 0.5)*it->second;
    }
    mvo+=mve;
  }
  return mvo;
}

/*===========================================================================*\
*/ /**
  Do not change multi-vector in projection space.

  @param mv Multi-vector.
  @tparam Sp Space of multi-vectors.
  @return The same multi-vector.
*/
template<typename Sp> constexpr MVect<Sp> Convert2OriginInfinity(const MVect<Sp>& mv) requires (Sp::Proj == 1)
{
  return mv;
}

/*===========================================================================*\
*/ /**
  Test geometry.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
  @note For geometry wee need projective space.
*/
template<typename Sp> void TestPlaneIntersections(nat nTests = 20000) requires (Sp::Proj > 0)
{
  nat n;
  real r,err;
  nat cntInt=0;
  nat cntPts=0;
  real rSumDist=0;
  real rMinDist=0;
  real rMaxDist=0;
  real rSumErr=0;
  real rMaxErr=0;
  const nat nDim=Sp::Pos+Sp::Neg;
  for (n=0; n<nTests; n++)
  {
    MVect<Sp> pt[2*nDim];
    for (nat i=0; i<2*nDim; ++i)
      pt[i]=RandomPoint<Sp>();

    // Construct two hyper-planes from random points.
    MVect<Sp> b1=pt[0];
    for (nat j=1; j<nDim; ++j)
      b1 ^= pt[j];
    MVect<Sp> b2=pt[nDim];
    for (nat j=nDim+1; j<2*nDim; ++j)
      b2 ^= pt[j];
    // In conformal space we need infinity.
    if constexpr (Sp::Proj>1)
    {
      b1 ^= Sp::Infinity();
      b2 ^= Sp::Infinity();
    }
    // Calculate intersection.
    MVect<Sp> b0=Meet(b1,b2);
    cntInt++;

    for (nat i=0; i<2*nDim; ++i)
    {
      real s=0.0;
      // Construct hyperplane from intersection and one of the points.
      MVect<Sp> bConstr=b0^pt[i];
      if ((bConstr == MVect<Sp>()) || !CmpReal0((bConstr*bConstr.Rev())[0]))
      // This point is exactly on the intersection.
        err=0.0;
      else
      {
        // Point belongs to this hyperplane.
        MVect<Sp>& bWith=i < nDim? b1 : b2;
        // Relation between original and constructed hyperplane.
        MVect<Sp> d=bWith/bConstr;
        // Original and constructed hyperplanes should be congruent, even if original planes were parallel.
        // In this case intersection is just a direction, not bounded to origin.
        if (!CmpReal0(d[0]))
        {
          std::cout << "P" << i << "=" << Convert2OriginInfinity(pt[i]) << std::endl;
          std::cout << "B1=" << Convert2OriginInfinity(b1) << std::endl;
          std::cout << "B2=" << Convert2OriginInfinity(b2) << std::endl;
          std::cout << "B0=" << Convert2OriginInfinity(b0) << std::endl;
          std::cout << "B0^P"<< i << "=" << Convert2OriginInfinity(bConstr) << std::endl;
          std::cout << "d=" << Convert2OriginInfinity(d) << std::endl;
          throw std::runtime_error("Error blades Bj and B0^Pi different");
        }
        // Calculate error.
        r=d[0];
        err=0.0;
        for (auto it=++d.begin(); it != d.end(); ++it)
          err += std::abs(it->second);
        err /= std::abs(r);

        // Not necessary, but we calculate also distance from this point to intersection
        // (it can be point, plane, hyperplane).
        MVect<Sp> I(Sp::I);
        MVect<Sp> bn;
        // Calculate orthogonal bivector to intersection.
        if constexpr (Sp::Proj == 1)
          // For projective space it is simply dual.
          bn=b0*I;
        else
          // For conformal we first map to projective, take projective dual and flatten the bivector.
          bn=(b0|Sp::Origin())*(I|Sp::Origin())^Sp::Infinity();

        // Bind this bivector to the point, it is now crossing the point and orthogonal to intersection.
        d=pt[i]^bn;
        // Intersect it with intersection. This should be flat point.
        d=Meet(d,b0);
        if constexpr (Sp::Proj == 2)
        // Convert flat point to space point.
          d |= Sp::Origin();

        if (d.w() == 0)
        // No intersection, planes were parallel.
          continue;
        // Normalize.
        d /= d.w();
        // Calculate distance.
        d -= pt[i];
        s=std::sqrt(std::abs((d*d)[0]));
      }
      // Statistic distance.
      rSumDist += s;
      if (!n)
        rMinDist=s;
      else
        rMinDist=std::min(rMinDist,s);
      rMaxDist=std::max(rMaxDist,s);
      // Statistic error.
      rSumErr += err;
      rMaxErr=std::max(rMaxErr, err);
      cntPts++;
    }
  }
  std::cout << "Intersections/points " << cntInt << "/" << cntPts << ", error av/max: " << rSumErr/cntPts << "/" << rMaxErr  << " -> OK" << std::endl;
  std::cout << "Distance intersection-point av/min/max: " << rSumDist/cntPts << "/" << rMinDist << "/" << rMaxDist << std::endl;
}

/*===========================================================================*\
*/ /**
  Perform all tests for certain space.

  @param nTests Tests to perform.
  @tparam Sp Space of multi-vectors.
*/
template<typename Sp> void TestSpace(nat nTests = 2000)
{
  auto t=std::chrono::system_clock::now();
  std::cout << "--- Space " << Sp::Pos << "," << Sp::Neg << "," << Sp::Proj <<std::endl;

  TestBasis<Sp>();
  TestBlades<Sp>();
  TestBladeInverse<Sp>();

  TestProducts<Sp>();
  TestMeet<Sp>();
  TestInverse<Sp>();
  TestLinear<Sp>();

  // For plane intersections wee need projective space.
  if constexpr (Sp::Proj > 0)
  {
    // Set random to compare results in similar spaces.
    srand(1);
    TestPlaneIntersections<Sp>();
  }

  auto dt=std::chrono::system_clock::now()-t;
  std::cout << std::chrono::duration_cast<float_seconds>(dt) << std::endl << std::endl;
}

int main()
{
  auto tm=std::chrono::system_clock::now();
  std::cout << std::uppercase;
  try
  {
    // Euclidean-type.
    TestSpace<Space<2>>();
    TestSpace<Space<3>>();
    TestSpace<Space<4>>();
    TestSpace<Space<5>>();
    TestSpace<Space<6>>();
    // Minkovsky-type.
    TestSpace<Space<1,1>>();
    TestSpace<Space<1,2>>();
    TestSpace<Space<1,3>>();
    TestSpace<Space<1,4>>();
    TestSpace<Space<1,5>>();
    // Euclidean-type projection.
    TestSpace<Space<2,0,1>>();
    TestSpace<Space<3,0,1>>();
    TestSpace<Space<4,0,1>>();
    TestSpace<Space<5,0,1>>();
    // Minkovsky-type projection.
    TestSpace<Space<1,1,1>>();
    TestSpace<Space<1,2,1>>();
    TestSpace<Space<1,3,1>>();
    TestSpace<Space<1,4,1>>();
    // Euclidean-type conformal.
    TestSpace<Space<2,0,2>>();
    TestSpace<Space<3,0,2>>();
    TestSpace<Space<4,0,2>>();
    // Minkovsky-type conformal.
    TestSpace<Space<1,1,2>>();
    TestSpace<Space<1,2,2>>();
    TestSpace<Space<1,3,2>>();
  }
  catch (std::exception &e)
  {
    std::cout << std::endl << "***************** " << e.what() << " *****************" << std::endl << std::endl;
  }
  auto dtm=std::chrono::system_clock::now()-tm;
  std::cout << "Total time: " << std::chrono::duration_cast<float_seconds>(dtm) << std::endl;
  return 0;
}

