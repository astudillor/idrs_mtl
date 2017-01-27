// Copyright 2017 Reinaldo Astudillo and Martin van Gijzen
//
// Permission is hereby granted, free of charge, to any person obtaining a 
// copy of this software and associated documentation files (the "Software"), 
// to deal in the Software without restriction, including without limitation 
// the rights to use, copy, modify, merge, publish, distribute, sublicense, 
// and/or sell copies of the Software, and to permit persons to whom the 
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in 
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
// DEALINGS IN THE SOFTWARE.
// 
// IDR(s) for solving systems of linear equations, biorthogonal version.
//
// Reference: 
// Martin B. van Gijzen and Peter Sonneveld, Algorithm 913: An Elegant IDR(s) 
// Variant that Efficiently Exploits Bi-orthogonality Properties.
// ACM Transactions on Mathematical Software, Vol. 38, No. 1, pp. 5:1-5:19, 2011

#ifndef ITL_BIDR_S_INCLUDE
#define ITL_BIDR_S_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/operation/random.hpp>
#include <boost/numeric/mtl/operation/orth.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/matrix/strict_upper.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>

namespace itl {
using namespace mtl;

/// Induced Dimension Reduction (IDR(s)) 
template < class LinearOperator, class HilbertSpaceX, class HilbertSpaceB, 
	   class Preconditioner, class Iteration >
int idr_s(const LinearOperator& A, HilbertSpaceX& x, const HilbertSpaceB& b, 
	     const Preconditioner& M, Iteration& iter, size_t s)
{
    typedef typename mtl::Collection<HilbertSpaceX>::value_type Scalar;
    typedef HilbertSpaceX                                       Vector;
    mtl::vampir_trace<7010> tracer;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");
    if (s < 1) s= 1;

    const Scalar                zero= math::zero(Scalar());
    const Scalar                one= math::one(Scalar());
    const double angle = 0.7; 
    Vector                      v(resource(x)), t(resource(x)), r(b - A * x);
    mat::multi_vector<Vector>   G(Vector(resource(x), zero), s), U(Vector(resource(x), zero), s), P(Vector(resource(x), zero), s);
    dense_vector<Scalar>   f(s), c(s);
    mat::dense2D<Scalar>        Ms(s, s); Ms = one;
    double om = 1.0;
    random(P);
    orth(P);
    // Main iteration loop, build G-spaces:
    while (! iter.finished(r)) {
         // New righ-hand size for small system:
         f = trans(P) * r;
         for (size_t k= 0; k < s; k++) {
            // Solve small system and make v orthogonal to P
            irange range(k, s);
            c[range] = lower_trisolve(Ms[range][range], f[range]);
            // Preconditioning
            v = solve(M, r - G.vector(range)*c[range]);
            // Compute new U(:,k) and G(:,k), G(:,k) is in space G_j
            U.vector(k) = U.vector(range)*c[range] + om*v;
            // Matvec
            G.vector(k) = A*U.vector(k);
            // Bi-Orthogonalise the new basis vectors
            for (size_t i= 0; i < k; i++) {
                Scalar alpha = dot(P.vector(i), G.vector(k))/Ms[i][i];
                G.vector(k) = G.vector(k) - alpha*G.vector(i);
                U.vector(k) = U.vector(k) - alpha*U.vector(i);
            }
            // New column of M = P'*G  (first k-1 entries are zero)
            for (size_t i= k; i < s; i++) {
                Ms[i][k] = dot(P.vector(i), G.vector(k));
            }
            // Make r orthogonal to g_i, i = 1..k 
            Scalar beta = f[k] / Ms[k][k];
            r = r - beta*G.vector(k);
            x = x + beta*U.vector(k);
            if ((++iter).finished(r))
                return iter;
            if (k < s-1) {
               for (size_t i= k+1; i < s; i++) {
                   f[i] = f[i] - beta*Ms[i][k];
               }
            }
         }
         // Now we have sufficient vectors in G_j to compute residual in G_j+1
         // Note: r is already perpendicular to P so v = r
         // Preconditioning
         v = solve(M, r);
         // Matvec
         t = A*v;
         // Computation of a new omega
         double ns = two_norm(r);
         double nt = two_norm(t);
         double ts = dot(t, r);
         double rho = fabs(ts/(nt*ns));
         om=ts/(nt*nt);
         if (rho < angle)
             om = om*angle/rho;
         r = r - om*t;
         x = x + om*v;
         ++iter;
    }
    return iter;
}


} // namespace itl

#endif // ITL_BIDR_S_INCLUDE
