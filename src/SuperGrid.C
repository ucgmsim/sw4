//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 
#include "SuperGrid.h"
#include "Require.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;

SuperGrid::SuperGrid()
{
  m_left = false;
  m_right = false;
  m_x0=0.;
  m_x1=1.;
  m_width=0.1;
  m_const_width=0.;
  m_epsL = 1e-4;
  m_tw_omega = 1.2;
}

void SuperGrid::print_parameters() const
{
   printf("SuperGrid parameters left=%i, right=%i, x0=%e, x1=%e, width=%e, transition=%e epsL=%e\n", m_left, m_right, m_x0, m_x1, m_width, m_trans_width,m_epsL);
}

void SuperGrid::define_taper(bool left, float_sw4 leftStart, bool right, float_sw4 rightEnd, float_sw4 width)
{
  m_left = left;
  m_x0 = leftStart;
  m_right = right;
  m_x1 = rightEnd;
  m_width = width;
// always use the full width for the transition, making m_const_width=0
// experimenting ...
  m_trans_width = 0.5*width;
//  m_trans_width = 1.0*width;
//  m_trans_width = transWidth;
  m_const_width = m_width - m_trans_width;
  
// sanity checks
  if (m_left || m_right)
  {
     float_sw4 dlen = m_x1-m_x0;
     CHECK_INPUT(m_width > 0., "The supergrid taper width must be positive, not = " << m_width);
     CHECK_INPUT(m_width < dlen, "The supergrid taper width must be smaller than the domain, not = " << m_width);
     CHECK_INPUT(m_trans_width > 0., "The supergrid taper transition width must be positive, not = " << m_trans_width);
     CHECK_INPUT(m_const_width >= 0., "The supergrid const_width = width - trans_width must be non-negative, not = " << m_const_width);
  }
  
  if (m_left && m_right)
  {
    if (m_x0+m_width > m_x1-m_width)
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left and right must be separated. Here x0+width = " << m_x0+m_width << 
		  " and x1-width = " << m_x1-m_width);
    }
    
  }
  else if( m_left )
  {
    if (m_x0+m_width > m_x1 )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the left must be smaller than the domain. Here x0+width = " << m_x0+m_width << 
		  " and x1 = " << m_x1);
    }
  }    
  else if( m_right )
  {
    if (m_x0 > m_x1-m_width )
    {
      print_parameters();
      CHECK_INPUT(false, "The supergrid taper functions at the right must be smaller than the domain. Here x0 = " << m_x0 << 
		  " and x1-width = " << m_x1-m_width );
    }
  }
}

float_sw4 SuperGrid::dampingCoeff(float_sw4 x) const
{
  float_sw4 phi = stretching(x);
// should be equivalent to PsiAux/phi
//  float_sw4 f=(1-phi)/phi/(1-m_epsL);
  float_sw4 f = PsiDamp(x)/phi;
// replaced PsiAux by PsiDamp, which goes to one faster
  
  return f;
}

float_sw4 SuperGrid::stretching( float_sw4 x ) const
{ // this function satisfies 0 < epsL <= f <= 1
  return 1-(1-m_epsL)*PsiAux(x); // PsiAux(x) = psi(x) in our papers
}

float_sw4 SuperGrid::PsiAux(float_sw4 x) const
{ // PsiAux = psi in our papers
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  float_sw4 f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 <= x <= m_width
    f=Psi0( (m_x0 + m_width - x)/m_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 
    f=Psi0( (x - (m_x1-m_width) )/m_width);
  return f;
}

float_sw4 SuperGrid::PsiDamp(float_sw4 x) const
{ // PsiAux = psi in our papers
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  float_sw4 f=0.;
  if (m_left && x < m_x0+m_width)
// the following makes the damping transition in 0 < const_width <= x <= const_width+trans_width = m_width
// constant damping in 0 <= x <= const_width
    f=Psi0( (m_x0 + m_width - x)/m_trans_width); 
  else if (m_right && x > m_x1-m_width)
// the following makes the damping transition in m_x1-m_width < x < m_x1 - const_width < m_x1
// constant damping in m_x1 - const_width <= x <= m_x1
    f=Psi0( (x - (m_x1-m_width) )/m_trans_width);
  return f;
}

float_sw4 SuperGrid::linTaper(float_sw4 x) const
{ 
// this function is zero for m_x0+m_width <= x <= m_x1-m_width
// and one for x=m_x0 and x=m_x1
  float_sw4 f=0.;
  if (m_left && x < m_x0+m_width)
//  linear taper from 0 to 1
    f= (m_x0 + m_width - x)/m_width; 
  else if (m_right && x > m_x1-m_width)
// linear taper from 0 to 1
    f= (x - (m_x1-m_width) )/m_width;
  return f;
}


// used for damping coefficient
float_sw4 SuperGrid::Psi0(float_sw4 xi) const
{
   float_sw4 f;
   if (xi<=0.)
      f = 0;
   else if (xi>=1.)
      f = 1.0;
   else
//    f=xi*xi*xi*(10 - 15*xi + 6*xi*xi);
//    f = fmin + (1.-fmin)*xi*xi*xi*(10 - 15*xi + 6*xi*xi);
// C4 function
//    f = fmin + (1.-fmin)* xi*xi*xi*xi*xi*( 
//      126 - 420*xi + 540*xi*xi - 315*xi*xi*xi + 70*xi*xi*xi*xi );
// Skewed C4 fcn (p3)
//    f = xi*xi*xi*xi*xi*(-14.0 + 70.0*xi - 90.0*xi*xi + 35.0*xi*xi*xi);
// C5 function (currently the default stretching function)  (p1)
// Evaluated in double regardless of float_sw4: the bracketed polynomial sums
// terms up to ~3465 that nearly cancel as xi->1, which in float32 loses
// enough precision to make the "monotonically increasing" damping profile
// non-monotonic (confirmed: ~4e-4 absolute error, small overshoot above 1).
// This only runs at setup time to build the static per-axis stretching
// arrays, so the extra double arithmetic has no runtime cost.
   {
      double xid = xi;
      double fd =  xid*xid*xid*xid*xid*xid*(
        462-1980*xid+3465*xid*xid-3080*xid*xid*xid+1386*xid*xid*xid*xid-252*xid*xid*xid*xid*xid);
      f = (float_sw4)fd;
   }
// one-sided C5 fcn (p2)
//     f =  xi*xi*xi*xi*xi*xi*(84.0 - 216.0*xi + 189.0*xi*xi - 56.0*xi*xi*xi);
   
   return f;
}

float_sw4 SuperGrid::cornerTaper( float_sw4 x ) const
{ // this function is 1 in the interior and tapers linearly to 1/2 in the SG layers
  const float_sw4 cmin=0.33;
  return 1.0 - (1.0-cmin)*linTaper(x);
}

float_sw4 SuperGrid::tw_stretching( float_sw4 x ) const
{
   return 1 + 0.5*sin(m_tw_omega*x);
}
  
void SuperGrid::set_twilight( float_sw4 omega )
{
   m_tw_omega = omega;
}

void SuperGrid::set_eps( float_sw4 new_eps )
{
   m_epsL = new_eps;
}

//-----------------------------------------------------------------------
// Spectral radius of the 1-D supergrid damping operator, per unit damping
// coefficient (the 'dc' of the supergrid command).
//
// Why it matters. After each predictor-corrector step SW4 subtracts
//    dc * D (u^n - u^{n-1}),   D = diag(str) * B^T diag(a) B
// (addsgd4_ci: B = second difference, a = dcx; addsgd6_ci: B = third
// difference, a = dcx averaged to half points). Freezing coefficients, a mode
// with -dt^2 L -> mu >= 0 and dc*D -> d >= 0 obeys
//    z^2 - (2 - g - d) z + (1 - d) = 0,   g = mu - mu^2/12,
// whose roots stay in the unit disk iff 0 <= d < 2 and 0 <= g <= 4 - 2d.
// The CFL limit (mu <= 12) only guarantees 0 <= g <= 3, so the damping is
// harmless at every stable time step iff d <= 1/2. For d > 1/2 the modes with
// mu in 6 -/+ sqrt(24d-12) grow: the run turns unstable once the time step is
// large enough for mu to reach that band, and for d >= 2 at any time step.
//
// d is not bounded by the constant-coefficient symbol (16 at 4th order, 64 at
// 6th): the product str(i)*a(i+1) ~ phi(i)/phi(i+1) and the stretching phi
// falls like (distance to the boundary)^6 towards epsL, so the neighbour
// ratios, and with them the largest eigenvalue, grow as the layer gets
// narrower in grid points (4th order: ~16.4 at 60 points, 17.7 at 30, 20 at
// 20, 33 at 12, 59 at 6). With the default dc=0.02 a 12-point layer has
// d ~ 0.67 and grows without bound at CFL >= 1.
//
// D is self-adjoint in the 1/str-weighted inner product, so its eigenvalues
// are those of the symmetric banded matrix G = S^(1/2) B^T A B S^(1/2),
// S = diag(str). The largest one is found by bisection: sigma > lambda_max(G)
// iff sigma*I - G has a Cholesky factorisation. Density is taken constant
// (the kernel's rho(i+-1)/rho(i) factors are ~1 for smooth material).
double SuperGrid::damping_spectral_radius( float_sw4 xmin, float_sw4 h, int n, int order ) const
{
   if( !is_active() || n < 1 )
      return 0;
   const int p = order == 6 ? 3 : 2;          // half band width of G
   // difference stencil rows m (B_m u = sum_q w[q] u_{m+q-off})
   const double w4[3] = {1, -2, 1}, w6[4] = {-1, 3, -3, 1};
   const double* w = order == 6 ? w6 : w4;
   const int nw = order == 6 ? 4 : 3, off = 1;  // stencil covers m-1 .. m-1+nw-1
   // grid index i = 1..n; rows m whose stencil touches 1..n: m = 2-nw+off .. n+off
   const int mlo = 2-nw+off, mhi = n+off;
   std::vector<double> a(mhi-mlo+1), sq(n+1);
   for( int m=mlo ; m <= mhi ; m++ )
   {
      if( order == 6 ) // half point m+1/2
         a[m-mlo] = 0.5*( dampingCoeff(xmin+(m-1)*h) + dampingCoeff(xmin+m*h) );
      else
         a[m-mlo] = dampingCoeff(xmin+(m-1)*h);
   }
   for( int i=1 ; i <= n ; i++ )
      sq[i] = sqrt( (double)stretching(xmin+(i-1)*h) );

   // lower band of G: g[(i-1)*(p+1)+(i-j)] = G(i,j), 0 <= i-j <= p
   std::vector<double> g((size_t)n*(p+1), 0.0);
   for( int m=mlo ; m <= mhi ; m++ )
      for( int q1=0 ; q1 < nw ; q1++ )
      {
         int i = m - off + q1;
         if( i < 1 || i > n ) continue;
         for( int q2=0 ; q2 <= q1 ; q2++ )
         {
            int j = m - off + q2;
            if( j < 1 ) continue;
            g[(size_t)(i-1)*(p+1)+(i-j)] += w[q1]*a[m-mlo]*w[q2]*sq[i]*sq[j];
         }
      }
   // Gershgorin bound as the upper end of the bisection interval
   double hi = 0;
   {
      std::vector<double> rs(n+1, 0.0);
      for( int i=1 ; i <= n ; i++ )
         for( int k=0 ; k <= p && i-k >= 1 ; k++ )
         {
            double v = fabs(g[(size_t)(i-1)*(p+1)+k]);
            rs[i] += v;
            if( k > 0 ) rs[i-k] += v;
         }
      for( int i=1 ; i <= n ; i++ )
         hi = rs[i] > hi ? rs[i] : hi;
   }
   if( hi == 0 )
      return 0;
   double lo = 0;
   std::vector<double> L((size_t)n*(p+1));
   for( int it=0 ; it < 60 && hi-lo > 1e-6*hi ; it++ )
   {
      double sigma = 0.5*(lo+hi);
      bool pd = true;
      for( int i=1 ; i <= n && pd ; i++ )
         for( int j=(i-p > 1 ? i-p : 1) ; j <= i ; j++ )
         {
            double s = (i == j ? sigma : 0) - g[(size_t)(i-1)*(p+1)+(i-j)];
            for( int k=(i-p > 1 ? i-p : 1) ; k < j ; k++ )
               s -= L[(size_t)(i-1)*(p+1)+(i-k)]*L[(size_t)(j-1)*(p+1)+(j-k)];
            if( i == j )
            {
               if( s <= 0 ) { pd = false; break; }
               L[(size_t)(i-1)*(p+1)] = sqrt(s);
            }
            else
               L[(size_t)(i-1)*(p+1)+(i-j)] = s/L[(size_t)(j-1)*(p+1)];
         }
      if( pd )
         hi = sigma;
      else
         lo = sigma;
   }
   return hi;
}
