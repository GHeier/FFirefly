/*
 * spline.h
 *
 * simple cubic spline interpolation library without external
 * dependencies
 *
 * ---------------------------------------------------------------------
 * Copyright (C) 2011, 2014, 2016, 2021 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 *
 */


#ifndef TK_SPLINE_H
#define TK_SPLINE_H

#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#ifdef HAVE_SSTREAM
#include <sstream>
#include <string>
#endif // HAVE_SSTREAM

// not ideal but disable unused-function warnings
// (we get them because we have implementations in the header file,
// and this is because we want to be able to quickly separate them
// into a cpp file if necessary)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"

// unnamed namespace only because the implementation is in this
// header file and we don't want to export symbols to the obj files

namespace tk
{

// spline interpolation
class spline
{
public:
    // spline types
    enum spline_type {
        linear = 10,            // linear interpolation
        cspline = 30,           // cubic splines (classical C^2)
        cspline_hermite = 31    // cubic hermite splines (local, only C^1)
    };

    // boundary condition type for the spline end-points
    enum bd_type {
        first_deriv = 1,
        second_deriv = 2,
        not_a_knot = 3
    };

protected:
    std::vector<double> m_x,m_y;            // x,y coordinates of points
    // interpolation parameters
    // f(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
    // where a_i = y_i, or else it won't go through grid points
    std::vector<double> m_b,m_c,m_d;        // spline coefficients
    double m_c0;                            // for left extrapolation
    spline_type m_type;
    bd_type m_left, m_right;
    double  m_left_value, m_right_value;
    bool m_made_monotonic;
    void set_coeffs_from_b();               // calculate c_i, d_i from b_i
    size_t find_closest(double x) const;    // closest idx so that m_x[idx]<=x

public:
    // default constructor: set boundary condition to be zero curvature
    // at both ends, i.e. natural splines
    spline(): m_type(cspline),
        m_left(second_deriv), m_right(second_deriv),
        m_left_value(0.0), m_right_value(0.0), m_made_monotonic(false)
    {
        ;
    }
    spline(const std::vector<double>& X, const std::vector<double>& Y,
           spline_type type = cspline,
           bool make_monotonic = false,
           bd_type left  = second_deriv, double left_value  = 0.0,
           bd_type right = second_deriv, double right_value = 0.0
          ):
        m_type(type),
        m_left(left), m_right(right),
        m_left_value(left_value), m_right_value(right_value),
        m_made_monotonic(false) // false correct here: make_monotonic() sets it
    {
        this->set_points(X,Y,m_type);
        if(make_monotonic) {
            this->make_monotonic();
        }
    }


    // modify boundary conditions: if called it must be before set_points()
    void set_boundary(bd_type left, double left_value,
                      bd_type right, double right_value);

    // set all data points (cubic_spline=false means linear interpolation)
    void set_points(const std::vector<double>& x,
                    const std::vector<double>& y,
                    spline_type type=cspline);

    // adjust coefficients so that the spline becomes piecewise monotonic
    // where possible
    //   this is done by adjusting slopes at grid points by a non-negative
    //   factor and this will break C^2
    //   this can also break boundary conditions if adjustments need to
    //   be made at the boundary points
    // returns false if no adjustments have been made, true otherwise
    bool make_monotonic();

    // evaluates the spline at point x
    double operator() (double x) const;
    double deriv(int order, double x) const;

    // solves for all x so that: spline(x) = y
    std::vector<double> solve(double y, bool ignore_extrapolation=true) const;

    // returns the input data points
    std::vector<double> get_x() const { return m_x; }
    std::vector<double> get_y() const { return m_y; }
    double get_x_min() const { assert(!m_x.empty()); return m_x.front(); }
    double get_x_max() const { assert(!m_x.empty()); return m_x.back(); }

#ifdef HAVE_SSTREAM
    // spline info string, i.e. spline type, boundary conditions etc.
    std::string info() const;
#endif // HAVE_SSTREAM

};



namespace internal
{

// band matrix solver
class band_matrix
{
private:
    std::vector< std::vector<double> > m_upper;  // upper band
    std::vector< std::vector<double> > m_lower;  // lower band
public:
    band_matrix() {};                             // constructor
    band_matrix(int dim, int n_u, int n_l);       // constructor
    ~band_matrix() {};                            // destructor
    void resize(int dim, int n_u, int n_l);      // init with dim,n_u,n_l
    int dim() const;                             // matrix dimension
    int num_upper() const
    {
        return (int)m_upper.size()-1;
    }
    int num_lower() const
    {
        return (int)m_lower.size()-1;
    }
    // access operator
    double & operator () (int i, int j);            // write
    double   operator () (int i, int j) const;      // read
    // we can store an additional diagonal (in m_lower)
    double& saved_diag(int i);
    double  saved_diag(int i) const;
    void lu_decompose();
    std::vector<double> r_solve(const std::vector<double>& b) const;
    std::vector<double> l_solve(const std::vector<double>& b) const;
    std::vector<double> lu_solve(const std::vector<double>& b,
                                 bool is_lu_decomposed=false);

};

double get_eps();

std::vector<double> solve_cubic(double a, double b, double c, double d,
                                int newton_iter=0);

} // namespace internal






namespace internal
{

// band_matrix implementation
// -------------------------

band_matrix::band_matrix(int dim, int n_u, int n_l)
{
    resize(dim, n_u, n_l);
}
void band_matrix::resize(int dim, int n_u, int n_l)
{
    assert(dim>0);
    assert(n_u>=0);
    assert(n_l>=0);
    m_upper.resize(n_u+1);
    m_lower.resize(n_l+1);
    for(size_t i=0; i<m_upper.size(); i++) {
        m_upper[i].resize(dim);
    }
    for(size_t i=0; i<m_lower.size(); i++) {
        m_lower[i].resize(dim);
    }
}
int band_matrix::dim() const
{
    if(m_upper.size()>0) {
        return m_upper[0].size();
    } else {
        return 0;
    }
}


// defines the new operator (), so that we can access the elements
// by A(i,j), index going from i=0,...,dim()-1
double & band_matrix::operator () (int i, int j)
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
    if(k>=0)    return m_upper[k][i];
    else        return m_lower[-k][i];
}
double band_matrix::operator () (int i, int j) const
{
    int k=j-i;       // what band is the entry
    assert( (i>=0) && (i<dim()) && (j>=0) && (j<dim()) );
    assert( (-num_lower()<=k) && (k<=num_upper()) );
    // k=0 -> diagonal, k<0 lower left part, k>0 upper right part
    if(k>=0)    return m_upper[k][i];
    else        return m_lower[-k][i];
}
// second diag (used in LU decomposition), saved in m_lower
double band_matrix::saved_diag(int i) const
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}
double & band_matrix::saved_diag(int i)
{
    assert( (i>=0) && (i<dim()) );
    return m_lower[0][i];
}

// LR-Decomposition of a band matrix
void band_matrix::lu_decompose()
{
    int  i_max,j_max;
    int  j_min;
    double x;

    // preconditioning
    // normalize column i so that a_ii=1
    for(int i=0; i<this->dim(); i++) {
        assert(this->operator()(i,i)!=0.0);
        this->saved_diag(i)=1.0/this->operator()(i,i);
        j_min=std::max(0,i-this->num_lower());
        j_max=std::min(this->dim()-1,i+this->num_upper());
        for(int j=j_min; j<=j_max; j++) {
            this->operator()(i,j) *= this->saved_diag(i);
        }
        this->operator()(i,i)=1.0;          // prevents rounding errors
    }

    // Gauss LR-Decomposition
    for(int k=0; k<this->dim(); k++) {
        i_max=std::min(this->dim()-1,k+this->num_lower());  // num_lower not a mistake!
        for(int i=k+1; i<=i_max; i++) {
            assert(this->operator()(k,k)!=0.0);
            x=-this->operator()(i,k)/this->operator()(k,k);
            this->operator()(i,k)=-x;                         // assembly part of L
            j_max=std::min(this->dim()-1,k+this->num_upper());
            for(int j=k+1; j<=j_max; j++) {
                // assembly part of R
                this->operator()(i,j)=this->operator()(i,j)+x*this->operator()(k,j);
            }
        }
    }
}
// solves Ly=b
std::vector<double> band_matrix::l_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_start;
    double sum;
    for(int i=0; i<this->dim(); i++) {
        sum=0;
        j_start=std::max(0,i-this->num_lower());
        for(int j=j_start; j<i; j++) sum += this->operator()(i,j)*x[j];
        x[i]=(b[i]*this->saved_diag(i)) - sum;
    }
    return x;
}
// solves Rx=y
std::vector<double> band_matrix::r_solve(const std::vector<double>& b) const
{
    assert( this->dim()==(int)b.size() );
    std::vector<double> x(this->dim());
    int j_stop;
    double sum;
    for(int i=this->dim()-1; i>=0; i--) {
        sum=0;
        j_stop=std::min(this->dim()-1,i+this->num_upper());
        for(int j=i+1; j<=j_stop; j++) sum += this->operator()(i,j)*x[j];
        x[i]=( b[i] - sum ) / this->operator()(i,i);
    }
    return x;
}

std::vector<double> band_matrix::lu_solve(const std::vector<double>& b,
        bool is_lu_decomposed)
{
    assert( this->dim()==(int)b.size() );
    std::vector<double>  x,y;
    if(is_lu_decomposed==false) {
        this->lu_decompose();
    }
    y=this->l_solve(b);
    x=this->r_solve(y);
    return x;
}

// machine precision of a double, i.e. the successor of 1 is 1+eps
double get_eps()
{
    //return std::numeric_limits<double>::epsilon();    // __DBL_EPSILON__
    return 2.2204460492503131e-16;                      // 2^-52
}

// solutions for a + b*x = 0
std::vector<double> solve_linear(double a, double b)
{
    std::vector<double> x;      // roots
    if(b==0.0) {
        if(a==0.0) {
            // 0*x = 0
            x.resize(1);
            x[0] = 0.0;   // any x solves it but we need to pick one
            return x;
        } else {
            // 0*x + ... = 0, no solution
            return x;
        }
    } else {
        x.resize(1);
        x[0] = -a/b;
        return x;
    }
}

// solutions for a + b*x + c*x^2 = 0
std::vector<double> solve_quadratic(double a, double b, double c,
                                    int newton_iter=0)
{
    if(c==0.0) {
        return solve_linear(a,b);
    }
    // rescale so that we solve x^2 + 2p x + q = (x+p)^2 + q - p^2 = 0
    double p=0.5*b/c;
    double q=a/c;
    double discr = p*p-q;
    const double eps=0.5*internal::get_eps();
    double discr_err = (6.0*(p*p)+3.0*fabs(q)+fabs(discr))*eps;

    std::vector<double> x;      // roots
    if(fabs(discr)<=discr_err) {
        // discriminant is zero --> one root
        x.resize(1);
        x[0] = -p;
    } else if(discr<0) {
        // no root
    } else {
        // two roots
        x.resize(2);
        x[0] = -p - sqrt(discr);
        x[1] = -p + sqrt(discr);
    }

    // improve solution via newton steps
    for(size_t i=0; i<x.size(); i++) {
        for(int k=0; k<newton_iter; k++) {
            double f  = (c*x[i] + b)*x[i] + a;
            double f1 = 2.0*c*x[i] + b;
            // only adjust if slope is large enough
            if(fabs(f1)>1e-8) {
                x[i] -= f/f1;
            }
        }
    }

    return x;
}

// solutions for the cubic equation: a + b*x +c*x^2 + d*x^3 = 0
// this is a naive implementation of the analytic solution without
// optimisation for speed or numerical accuracy
// newton_iter: number of newton iterations to improve analytical solution
// see also
//   gsl: gsl_poly_solve_cubic() in solve_cubic.c
//   octave: roots.m - via eigenvalues of the Frobenius companion matrix
std::vector<double> solve_cubic(double a, double b, double c, double d,
                                int newton_iter)
{
    if(d==0.0) {
        return solve_quadratic(a,b,c,newton_iter);
    }

    // convert to normalised form: a + bx + cx^2 + x^3 = 0
    if(d!=1.0) {
        a/=d;
        b/=d;
        c/=d;
    }

    // convert to depressed cubic: z^3 - 3pz - 2q = 0
    // via substitution: z = x + c/3
    std::vector<double> z;              // roots of the depressed cubic
    double p = -(1.0/3.0)*b + (1.0/9.0)*(c*c);
    double r = 2.0*(c*c)-9.0*b;
    double q = -0.5*a - (1.0/54.0)*(c*r);
    double discr=p*p*p-q*q;             // discriminant
    // calculating numerical round-off errors with assumptions:
    //  - each operation is precise but each intermediate result x
    //    when stored has max error of x*eps
    //  - only multiplication with a power of 2 introduces no new error
    //  - a,b,c,d and some fractions (e.g. 1/3) have rounding errors eps
    //  - p_err << |p|, q_err << |q|, ... (this is violated in rare cases)
    // would be more elegant to use boost::numeric::interval<double>
    const double eps = internal::get_eps();
    double p_err = eps*((3.0/3.0)*fabs(b)+(4.0/9.0)*(c*c)+fabs(p));
    double r_err = eps*(6.0*(c*c)+18.0*fabs(b)+fabs(r));
    double q_err = 0.5*fabs(a)*eps + (1.0/54.0)*fabs(c)*(r_err+fabs(r)*3.0*eps)
                   + fabs(q)*eps;
    double discr_err = (p*p) * (3.0*p_err + fabs(p)*2.0*eps)
                       + fabs(q) * (2.0*q_err + fabs(q)*eps) + fabs(discr)*eps;

    // depending on the discriminant we get different solutions
    if(fabs(discr)<=discr_err) {
        // discriminant zero: one or two real roots
        if(fabs(p)<=p_err) {
            // p and q are zero: single root
            z.resize(1);
            z[0] = 0.0;             // triple root
        } else {
            z.resize(2);
            z[0] = 2.0*q/p;         // single root
            z[1] = -0.5*z[0];       // double root
        }
    } else if(discr>0) {
        // three real roots: via trigonometric solution
        z.resize(3);
        double ac = (1.0/3.0) * acos( q/(p*sqrt(p)) );
        double sq = 2.0*sqrt(p);
        z[0] = sq * cos(ac);
        z[1] = sq * cos(ac-2.0*M_PI/3.0);
        z[2] = sq * cos(ac-4.0*M_PI/3.0);
    } else if (discr<0.0) {
        // single real root: via Cardano's fromula
        z.resize(1);
        double sgnq = (q >= 0 ? 1 : -1);
        double basis = fabs(q) + sqrt(-discr);
        double C = sgnq * pow(basis, 1.0/3.0); // c++11 has std::cbrt()
        z[0] = C + p/C;
    }
    for(size_t i=0; i<z.size(); i++) {
        // convert depressed cubic roots to original cubic: x = z - c/3
        z[i] -= (1.0/3.0)*c;
        // improve solution via newton steps
        for(int k=0; k<newton_iter; k++) {
            double f  = ((z[i] + c)*z[i] + b)*z[i] + a;
            double f1 = (3.0*z[i] + 2.0*c)*z[i] + b;
            // only adjust if slope is large enough
            if(fabs(f1)>1e-8) {
                z[i] -= f/f1;
            }
        }
    }
    // ensure if a=0 we get exactly x=0 as root
    // TODO: remove this fudge
    if(a==0.0) {
        assert(z.size()>0);     // cubic should always have at least one root
        double xmin=fabs(z[0]);
        size_t imin=0;
        for(size_t i=1; i<z.size(); i++) {
            if(xmin>fabs(z[i])) {
                xmin=fabs(z[i]);
                imin=i;
            }
        }
        z[imin]=0.0;        // replace the smallest absolute value with 0
    }
    std::sort(z.begin(), z.end());
    return z;
}


} // namespace internal


} // namespace tk



#pragma GCC diagnostic pop

#endif /* TK_SPLINE_H */
