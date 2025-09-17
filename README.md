# Geometric Algebra and Tensors (geometric_algebra)

Compact header-only with no dependencies C++20 geometric algebra library with support for Euclidean-type, Minkovsky-type, projection and conformal spaces.

These less than 2 thousand lines of code unveil the power of geometric algebra, which rapidly becomes popularity. 

It is a pure header library. It was tested under gcc 12.1 and Visual Studio 2022. Visual Studio is a bit (ca 8%) slower and, due to absence of long double, produces at least 3 orders of magnitude bigger rounding errors.

Here is a short introduction to the library.

## Basics

### Spaces

- Usual N-dimensional space: Space\<N\>
- Minkovsky-type with negative dimensions, e.g. Einsteins space-time: Space<1,3>
- Projection spaces with 1 additional dimension for geometry calculations: Space<N,M,1>
- Conformal spaces with 2 additional dimensions for geometry calculations also with spheres: Space<N,M,2>
- Each space has corresponding co-space: Space::Cospace
- For each space scalar l (like 1) and pseudo-scalar I are defined 

### Basis elements

Basis elements BElem<Space> use effective bit representation, which allows up to 27 dimensions in 32-bit mode and 59 dimensions in 64-bit mode. Operations on basis elements:

- Grassmann left and right complement: !be, ~be
- Inner product, basis elements belong to the same space, or to space and co-space: be1 | be2
- Outer product, basis elements belong to the same space, or to space and co-space: be1 ^ be2
- Geometric product, basis elements belong to the same space, or to space and co-space: be1 * be2
- Linear transformation: Matrix * be

### Multi-vectors

Multi-vectors MVect<Space> are implemented as map of basis elements to corresponding factors and support following operations:

- Inversions for arbitrary (!) multi-vectors: mv.Inv()
- Hat involution: mv.Hat()
- Reversion: mv.Rev()
- Clifford conjugation: mv.Conj()
- Representation in co-space: mv.CoMVect()
- Left and right complement: mv, ~mv
- Negative: -mv
- Addition and subtraction: mv1 + mv2, mv1 - mv2
- Multiplication and division by factor: mv1 * factor, mv1 / factor
- Interior product, multi-vectors belong to the same space, or to space and co-space: mv1 | mv2
- Exterior product, multi-vectors belong to the same space, or to space and co-space: mv1 ^ mv2
- Geometric product, multi-vectors belong to the same space, or to space and co-space: mv1 * mv2
- Geometric division, multi-vectors belong to the same space, or to space and co-space: mv1 / mv2
- Linear transformation: Matrix * mv 
- Grassmann regressive product: Meet(mv1, mv2)
- Grassmann left and right complement: DualL(mv), DualR(mv)
- Left and right dual spaces: DualIL(mv), DualIR(mv)

### Tensors

Tensors are here to implement linear transformations and multi-vector inversions.

**Note:** tensor interface is not jet stable and may be subject to change.

- **SubTensor**
  - Access to tensor components, mainly to have a possibility to access elements like in multi-dimensional array, like: t[i][j][k]
  - New tensor can be constructed from any sub-tensor: t2(t1[i][j])
- **Tensor** supports arbitrary number of upper and lower indices and following operations:
  - Negative: -t
  - Addition and subtraction: t1 + t2, t1 - t2
  - Multiplication and division by factor: t1 * factor, t1 / factor
  - Tensor product, convolution last index t1 and first index t2, one of them should be upper, another - lower: t1 * t2
- **Vector** is Tensor with 1 index, upper or lower
- **Matrix** is Tensor with 2 indices, lower and upper
  - Unary square matrix: l(dim)
  - Transposed matrix: T(M)
  - Gauss elimination: SolveGauss(a, x)
  
## Applications

Geometric algebra enormously simplifies many non-trivial calculations in arbitrary spaces. Here are only some examples of what is possible.

### Reflections

- Any vector a can be reflected against hyperplane, orthogonal to vector v: -v * a * v
- Any blade of dimension k is reflected similarly, up to the sign: (k & 1? -1 : 1) * v * Bk * v 

### Rotations

Rotation can be performed with two reflections, first against hyperplane, orthogonal to vector v1 and then to vector v2. Rotation angle is then double angle between vectors. Product R = v1 * v2 is called rotor. And we need also R.Rev() = v2 * v1.

- Any multi-vector is rotated by: R * mv * R.Rev()
- Back rotation: R.Rev() * mv * R

### Projective geometry

- Points here contain coordinate origin o (this is additional projective dimension): A = o + v
- Line can be defined with two points or point and vector: L = A ^ B, L = A ^ v
- Plane can be defined with three points, two points and vector or point and two vectors: P = A ^ B ^ C, P = A ^ B ^ v, P = A ^ v1 ^ v2
- Intersection of two lines can be 0, point or line: P = Meet(L1, L2)
- Intersection of two planes can be 0, line or plane: L = Meet(P1, P2)
- Intersection of line with plane can be 0, point or line: P = Meet(L, P)

### Conformal geometry

- Points in conformal geometry get one more additional projection - infinity: A = o + v + 0.5 * v * v * inf
- Internally instead of o and inf equivalent basis p and n is used: o = 0.5 * (p + n), inf = n - p
- Circle can be defined with 3 points: Circ = A ^ B ^ C
- To define line, the third point must be at infinity: L = A ^ B ^ Space::Infinity()
- Sphere can be defined with 4 points: S = A ^ B ^ C ^ D
- To define a plane, the forth point must be at infinity: P = A ^ B ^ C ^ Space::Infinity()
- Intersections of all kinds are calculated with the same Grassmann regressive product: Meet(S, P)

