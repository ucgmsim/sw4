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
//-----------------------------------------------------------------------
//  Adds 4th order artificial disssipation for super-grid damping layers
//
//-----------------------------------------------------------------------

#include "EW.h"


//-----------------------------------------------------------------------
// Supergrid damping is a sponge-layer operator, but the kernels below were
// sweeping the whole subdomain: EW::addSuperGridDamping passes
// m_iStart[g]..m_kEnd[g] and the loops ran that full extent minus their stencil
// margin. Every term in each body carries a dcx/dcy/dcz factor, and
// SuperGrid::PsiDamp returns a literal 0. outside the taper -- its own comment
// says "this function is zero for m_x0+m_width <= x <= m_x1-m_width" -- so for a
// point far from every sponge face the whole bracket is exactly 0.0 and
// `up -= birho*0.0` is the identity.
//
// Measured on HPC3, production settings: the supergrid phase costs the SAME
// fraction of runtime at gp=12 and gp=30 (4.35% vs 4.33%). A 2.5x thicker
// sponge for the same price is the signature -- thickness never entered the
// loop bounds.
//
// Fraction of points computing an exact zero, by geometry:
//     nx=235  gp=30   46%        nx=601   gp=30   76%
//     nx=235  gp=12   73%        nx=1001  gp=30   85%
// It gets worse as the grid refines, the sponge being a fixed point count.
//
// This is bit-exact, not an approximation: we drop additions of exactly zero.

// Find the index interval over which dc is zero across the whole stencil reach.
// Returns false if there is none, in which case the caller sweeps everything.
static bool sgd_dead_range( const float_sw4* dc, int first, int last, int reach,
			    int& a, int& b )
{
   if( dc == 0 || last < first )
      return false;
   auto active = [&]( int idx ) -> bool {
      for( int o = -reach ; o <= reach ; o++ )
      {
	 int t = idx + o;
	 if( t >= first && t <= last && dc[t-first] != 0 )
	    return true;
      }
      return false;
   };
   a = first;  while( a <= last  && active(a) ) a++;
   b = last;   while( b >= first && active(b) ) b--;
   if( a > b )
      return false;
// Verify the candidate interval really is dead throughout. The sponge profile
// makes this hold, but checking is O(n) once per call and removes the
// assumption: a dc array with a nonzero patch in the middle would otherwise
// silently lose that contribution.
   for( int idx = a ; idx <= b ; idx++ )
      if( active(idx) )
	 return false;
   return true;
}

// Decompose the box [i0,i1]x[j0,j1]x[k0,k1] minus the dead box
// [ia,ib]x[ja,jb]x[ka,kb] into up to six disjoint rectangular slabs. Each slab
// is contiguous in i, so the innermost loop stays unit-stride and vectorises as
// before. Returns the slab count, or one full-box slab when there is nothing to
// skip.
static int sgd_slabs( int i0,int i1, int j0,int j1, int k0,int k1,
		      bool dead, int ia,int ib, int ja,int jb, int ka,int kb,
		      int sl[][6] )
{
   int n = 0;
   auto add = [&]( int p,int q,int r,int s,int t,int u ) {
      if( p <= q && r <= s && t <= u )
      {
	 sl[n][0]=p; sl[n][1]=q; sl[n][2]=r; sl[n][3]=s; sl[n][4]=t; sl[n][5]=u; n++;
      }
   };
   if( !dead )
   {
      add( k0,k1, j0,j1, i0,i1 );
      return n;
   }
   add( k0, ka-1, j0, j1, i0, i1 );          // below the dead box in k
   add( kb+1, k1, j0, j1, i0, i1 );          // above it in k
   add( ka, kb, j0, ja-1, i0, i1 );          // dead k, below in j
   add( ka, kb, jb+1, j1, i0, i1 );          // dead k, above in j
   add( ka, kb, ja, jb, i0, ia-1 );          // dead k and j, below in i
   add( ka, kb, ja, jb, ib+1, i1 );          // dead k and j, above in i
   return n;
}

//-----------------------------------------------------------------------
void EW::addsgd4_ci( int ifirst, int ilast, int jfirst, int jlast,
		     int kfirst, int klast,
		     float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		     float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		     float_sw4* __restrict__ a_dcz, float_sw4* __restrict__ a_strx, 
		     float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz,
		     float_sw4* __restrict__ a_cox,  float_sw4* __restrict__ a_coy,
		     float_sw4* __restrict__ a_coz,
		     float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]

      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);

// AP: The for c loop could be inside the for i loop. The simd, ivdep pragmas should be outside the inner-most loop
// Restrict the sweep to where the damping coefficients are nonzero; see
// sgd_dead_range above.
      int sl[6][6], ia,ib,ja,jb,ka,kb;
      bool dx = sgd_dead_range( a_dcx, ifirst, ilast, 1, ia, ib );
      bool dy = sgd_dead_range( a_dcy, jfirst, jlast, 1, ja, jb );
      bool dz = sgd_dead_range( a_dcz, kfirst, klast, 1, ka, kb );
      int nsl = sgd_slabs( ifirst+2, ilast-2, jfirst+2, jlast-2,
			   kfirst+2, klast-2, dx && dy && dz,
			   ia, ib, ja, jb, ka, kb, sl );
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
      for( int s=0 ; s < nsl ; s++ )
// nowait: the slabs are disjoint by construction and the components are
// distinct, so no worksharing region here depends on another. Without it
// this kernel would pay up to 18 barriers per call instead of one.
#pragma omp for collapse(2) nowait
      for( int k=sl[s][0]; k <= sl[s][1] ; k++ )
	 for( int j=sl[s][2]; j <= sl[s][3] ; j++ )
	    //#pragma simd
#pragma ivdep
#pragma omp simd
	    for( int i=sl[s][4]; i <= sl[s][5] ; i++ )
	    {
	       float_sw4 birho=beta/rho(i,j,k);
	       {
		  up(c,i,j,k) -= birho*( 
		  // x-differences
		   strx(i)*coy(j)*coz(k)*(
       rho(i+1,j,k)*dcx(i+1)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
      -2*rho(i,j,k)*dcx(i)  *
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
      +rho(i-1,j,k)*dcx(i-1)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
      -rho(i+1,j,k)*dcx(i+1)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
      +2*rho(i,j,k)*dcx(i)  *
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
      -rho(i-1,j,k)*dcx(i-1)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
      stry(j)*cox(i)*coz(k)*(
      +rho(i,j+1,k)*dcy(j+1)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
      -2*rho(i,j,k)*dcy(j)  *
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
      +rho(i,j-1,k)*dcy(j-1)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
      -rho(i,j+1,k)*dcy(j+1)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
      +2*rho(i,j,k)*dcy(j)  *
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
      -rho(i,j-1,k)*dcy(j-1)*
                   (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) +
       strz(k)*cox(i)*coy(j)*(
// z-differences
      +rho(i,j,k+1)*dcz(k+1)* 
                 ( u(c,i,j,k+2) -2*u(c,i,j,k+1)+ u(c,i,j,k  )) 
      -2*rho(i,j,k)*dcz(k)  *
                 ( u(c,i,j,k+1) -2*u(c,i,j,k  )+ u(c,i,j,k-1))
      +rho(i,j,k-1)*dcz(k-1)*
                 ( u(c,i,j,k  ) -2*u(c,i,j,k-1)+ u(c,i,j,k-2)) 
      -rho(i,j,k+1)*dcz(k+1)*
                 (um(c,i,j,k+2)-2*um(c,i,j,k+1)+um(c,i,j,k  )) 
      +2*rho(i,j,k)*dcz(k)  *
                 (um(c,i,j,k+1)-2*um(c,i,j,k  )+um(c,i,j,k-1)) 
      -rho(i,j,k-1)*dcz(k-1)*
                 (um(c,i,j,k  )-2*um(c,i,j,k-1)+um(c,i,j,k-2)) ) 
					 );

	       }
	    }
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6_ci( int ifirst, int ilast, int jfirst, int jlast,
		     int kfirst, int klast,
		     float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		     float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		     float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		     float_sw4* __restrict__ a_dcz, float_sw4* __restrict__ a_strx,
		     float_sw4* __restrict__ a_stry, float_sw4* __restrict__ a_strz,
		     float_sw4* __restrict__ a_cox,  float_sw4* __restrict__ a_coy,
		     float_sw4* __restrict__ a_coz, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
#define strz(k) a_strz[(k-kfirst)]
#define dcz(k) a_dcz[(k-kfirst)]
#define coz(k) a_coz[(k-kfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
// Restrict the sweep to where the damping coefficients are nonzero; see
// sgd_dead_range above.
      int sl[6][6], ia,ib,ja,jb,ka,kb;
      bool dx = sgd_dead_range( a_dcx, ifirst, ilast, 2, ia, ib );
      bool dy = sgd_dead_range( a_dcy, jfirst, jlast, 2, ja, jb );
      bool dz = sgd_dead_range( a_dcz, kfirst, klast, 2, ka, kb );
      int nsl = sgd_slabs( ifirst+3, ilast-3, jfirst+3, jlast-3,
			   kfirst+3, klast-3, dx && dy && dz,
			   ia, ib, ja, jb, ka, kb, sl );
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
      for( int s=0 ; s < nsl ; s++ )
// nowait: the slabs are disjoint by construction and the components are
// distinct, so no worksharing region here depends on another. Without it
// this kernel would pay up to 18 barriers per call instead of one.
#pragma omp for collapse(2) nowait
      for( int k=sl[s][0]; k <= sl[s][1] ; k++ )
	 for( int j=sl[s][2]; j <= sl[s][3] ; j++ )
	    //#pragma simd
#pragma ivdep
#pragma omp simd
	    for( int i=sl[s][4]; i <= sl[s][5] ; i++ )
	    {
	       float_sw4 birho=0.5*beta/rho(i,j,k);
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*coz(k)*(
// x-differences
         (rho(i+2,j,k)*dcx(i+2)+rho(i+1,j,k)*dcx(i+1))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)+rho(i,j,k)*dcx(i))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)+rho(i-1,j,k)*dcx(i-1))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
       - (rho(i-1,j,k)*dcx(i-1)+rho(i-2,j,k)*dcx(i-2))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
                 ) +  stry(j)*cox(i)*coz(k)*(
// y-differences
         (rho(i,j+2,k)*dcy(j+2)+rho(i,j+1,k)*dcy(j+1))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
      -3*(rho(i,j+1,k)*dcy(j+1)+rho(i,j,k)*dcy(j))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
      +3*(rho(i,j,k)*dcy(j)+rho(i,j-1,k)*dcy(j-1))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
       - (rho(i,j-1,k)*dcy(j-1)+rho(i,j-2,k)*dcy(j-2))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
                 ) +  strz(k)*cox(i)*coy(j)*(
// z-differences
         ( rho(i,j,k+2)*dcz(k+2) + rho(i,j,k+1)*dcz(k+1) )*(
         u(c,i,j,k+3)- 3*u(c,i,j,k+2)+ 3*u(c,i,j,k+1)- u(c,i,  j,k) 
      -(um(c,i,j,k+3)-3*um(c,i,j,k+2)+3*um(c,i,j,k+1)-um(c,i,  j,k)) )
      -3*(rho(i,j,k+1)*dcz(k+1)+rho(i,j,k)*dcz(k))*(
         u(c,i,j,k+2) -3*u(c,i,j,k+1)+ 3*u(c,i,  j,k)- u(c,i,j,k-1) 
      -(um(c,i,j,k+2)-3*um(c,i,j,k+1)+3*um(c,i,  j,k)-um(c,i,j,k-1)) )
      +3*(rho(i,j,k)*dcz(k)+rho(i,j,k-1)*dcz(k-1))*(
         u(c,i,j,k+1)- 3*u(c,i,  j,k)+ 3*u(c,i,j,k-1)-u(c,i,j,k-2) 
      -(um(c,i,j,k+1)-3*um(c,i,  j,k)+3*um(c,i,j,k-1)-um(c,i,j,k-2)) )
       - (rho(i,j,k-1)*dcz(k-1)+rho(i,j,k-2)*dcz(k-2))*(
         u(c,i,  j,k) -3*u(c,i,j,k-1)+ 3*u(c,i,j,k-2)- u(c,i,j,k-3)
      -(um(c,i,  j,k)-3*um(c,i,j,k-1)+3*um(c,i,j,k-2)-um(c,i,j,k-3)) )
					     )  );
	       }
	    }
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef strz
#undef dcz
#undef coz
   }
}

//-----------------------------------------------------------------------
void EW::addsgd4c_ci( int ifirst, int ilast, int jfirst, int jlast,
		      int kfirst, int klast,
		      float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		      float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		      float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		      float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		      float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
		      float_sw4* __restrict__ a_coy, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]

      const size_t ni   =     (ilast-ifirst+1);
      const size_t nij  =  ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
// Restrict the sweep to where the damping coefficients are nonzero; see
// sgd_dead_range above.
      int sl[6][6], ia,ib,ja,jb,ka,kb;
      bool dx = sgd_dead_range( a_dcx, ifirst, ilast, 1, ia, ib );
      bool dy = sgd_dead_range( a_dcy, jfirst, jlast, 1, ja, jb );
// No dcz in the curvilinear kernels: the z-direction carries no supergrid
// damping here, so a point is dead whenever x and y are, at every k.
      bool dz = true; ka = kfirst+2; kb = klast-2;
      int nsl = sgd_slabs( ifirst+2, ilast-2, jfirst+2, jlast-2,
			   kfirst+2, klast-2, dx && dy && dz,
			   ia, ib, ja, jb, ka, kb, sl );
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
      for( int s=0 ; s < nsl ; s++ )
// nowait: the slabs are disjoint by construction and the components are
// distinct, so no worksharing region here depends on another. Without it
// this kernel would pay up to 18 barriers per call instead of one.
#pragma omp for collapse(2) nowait
      for( int k=sl[s][0]; k <= sl[s][1] ; k++ )
	 for( int j=sl[s][2]; j <= sl[s][3] ; j++ )
	    //#pragma simd
#pragma ivdep
#pragma omp simd
	    for( int i=sl[s][4]; i <= sl[s][5] ; i++ )
	    {
	       float_sw4 irhoj=beta/(rho(i,j,k)*jac(i,j,k));
	       {
		  up(c,i,j,k) -= irhoj*( 
		  // x-differences
		   strx(i)*coy(j)*(
		   rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   ( u(c,i+2,j,k) -2*u(c,i+1,j,k)+ u(c,i,  j,k))
		   -2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   ( u(c,i+1,j,k) -2*u(c,i,  j,k)+ u(c,i-1,j,k))
		   +rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   ( u(c,i,  j,k) -2*u(c,i-1,j,k)+ u(c,i-2,j,k)) 
		   -rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)*
                   (um(c,i+2,j,k)-2*um(c,i+1,j,k)+um(c,i,  j,k)) 
		   +2*rho(i,j,k)*dcx(i)*jac(i,j,k)*
                   (um(c,i+1,j,k)-2*um(c,i,  j,k)+um(c,i-1,j,k)) 
		   -rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)*
                   (um(c,i,  j,k)-2*um(c,i-1,j,k)+um(c,i-2,j,k)) ) +
// y-differences
		   stry(j)*cox(i)*(
		    +rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   ( u(c,i,j+2,k) -2*u(c,i,j+1,k)+ u(c,i,j,  k)) 
		    -2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   ( u(c,i,j+1,k) -2*u(c,i,j,  k)+ u(c,i,j-1,k))
		    +rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
                   ( u(c,i,j,  k) -2*u(c,i,j-1,k)+ u(c,i,j-2,k)) 
		    -rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)*
                   (um(c,i,j+2,k)-2*um(c,i,j+1,k)+um(c,i,j,  k)) 
		    +2*rho(i,j,k)*dcy(j)*jac(i,j,k)*
                   (um(c,i,j+1,k)-2*um(c,i,j,  k)+um(c,i,j-1,k)) 
		    -rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)*
		    (um(c,i,j,  k)-2*um(c,i,j-1,k)+um(c,i,j-2,k)) ) );
	       }
	    } 
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}

//-----------------------------------------------------------------------
void EW::addsgd6c_ci(  int ifirst, int ilast, int jfirst, int jlast,
		       int kfirst, int klast,
		       float_sw4* __restrict__ a_up, float_sw4* __restrict__ a_u,
		       float_sw4* __restrict__ a_um, float_sw4* __restrict__ a_rho,
		       float_sw4* __restrict__ a_dcx, float_sw4* __restrict__ a_dcy,
		       float_sw4* __restrict__ a_strx, float_sw4* __restrict__ a_stry, 
		       float_sw4* __restrict__ a_jac, float_sw4* __restrict__ a_cox,
		       float_sw4* __restrict__ a_coy, float_sw4 beta )
{
   if( beta != 0 )
   {
#define rho(i,j,k) a_rho[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define up(c,i,j,k) a_up[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define u(c,i,j,k)   a_u[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define um(c,i,j,k) a_um[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)+(c)*npts]
#define jac(i,j,k) a_jac[(i-ifirst)+ni*(j-jfirst)+nij*(k-kfirst)]
#define strx(i) a_strx[(i-ifirst)]
#define dcx(i) a_dcx[(i-ifirst)]
#define cox(i) a_cox[(i-ifirst)]
#define stry(j) a_stry[(j-jfirst)]
#define dcy(j) a_dcy[(j-jfirst)]
#define coy(j) a_coy[(j-jfirst)]
      const size_t ni = ilast-ifirst+1;
      const size_t nij = ni*(jlast-jfirst+1);
      const size_t npts = nij*(klast-kfirst+1);
// Restrict the sweep to where the damping coefficients are nonzero; see
// sgd_dead_range above.
      int sl[6][6], ia,ib,ja,jb,ka,kb;
      bool dx = sgd_dead_range( a_dcx, ifirst, ilast, 2, ia, ib );
      bool dy = sgd_dead_range( a_dcy, jfirst, jlast, 2, ja, jb );
// No dcz in the curvilinear kernels: the z-direction carries no supergrid
// damping here, so a point is dead whenever x and y are, at every k.
      bool dz = true; ka = kfirst+3; kb = klast-3;
      int nsl = sgd_slabs( ifirst+3, ilast-3, jfirst+3, jlast-3,
			   kfirst+3, klast-3, dx && dy && dz,
			   ia, ib, ja, jb, ka, kb, sl );
#pragma omp parallel
      {
      for( int c=0 ; c < 3 ; c++ )
      for( int s=0 ; s < nsl ; s++ )
// nowait: the slabs are disjoint by construction and the components are
// distinct, so no worksharing region here depends on another. Without it
// this kernel would pay up to 18 barriers per call instead of one.
#pragma omp for collapse(2) nowait
      for( int k=sl[s][0]; k <= sl[s][1] ; k++ )
	 for( int j=sl[s][2]; j <= sl[s][3] ; j++ )
	    //#pragma simd
#pragma ivdep
#pragma omp simd
	    for( int i=sl[s][4]; i <= sl[s][5] ; i++ )
	    {
	       float_sw4 birho=0.5*beta/(rho(i,j,k)*jac(i,j,k));
	       {
		 up(c,i,j,k) += birho*( 
       strx(i)*coy(j)*(
// x-differences
      (rho(i+2,j,k)*dcx(i+2)*jac(i+2,j,k)+rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k))*(
         u(c,i+3,j,k) -3*u(c,i+2,j,k)+ 3*u(c,i+1,j,k)- u(c,i, j,k) 
      -(um(c,i+3,j,k)-3*um(c,i+2,j,k)+3*um(c,i+1,j,k)-um(c,i, j,k)) )
      -3*(rho(i+1,j,k)*dcx(i+1)*jac(i+1,j,k)+rho(i,j,k)*dcx(i)*jac(i,j,k))*(
         u(c,i+2,j,k)- 3*u(c,i+1,j,k)+ 3*u(c,i, j,k)- u(c,i-1,j,k)
      -(um(c,i+2,j,k)-3*um(c,i+1,j,k)+3*um(c,i, j,k)-um(c,i-1,j,k)) )
      +3*(rho(i,j,k)*dcx(i)*jac(i,j,k)+rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k))*(
         u(c,i+1,j,k)- 3*u(c,i,  j,k)+3*u(c,i-1,j,k)- u(c,i-2,j,k) 
      -(um(c,i+1,j,k)-3*um(c,i, j,k)+3*um(c,i-1,j,k)-um(c,i-2,j,k)) )
      - (rho(i-1,j,k)*dcx(i-1)*jac(i-1,j,k)+rho(i-2,j,k)*dcx(i-2)*jac(i-2,j,k))*(
         u(c,i, j,k)- 3*u(c,i-1,j,k)+ 3*u(c,i-2,j,k)- u(c,i-3,j,k) 
      -(um(c,i, j,k)-3*um(c,i-1,j,k)+3*um(c,i-2,j,k)-um(c,i-3,j,k)) )
		       ) +  stry(j)*cox(i)*( 
// y-differences
     (rho(i,j+2,k)*dcy(j+2)*jac(i,j+2,k)+rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k))*(
         u(c,i,j+3,k) -3*u(c,i,j+2,k)+ 3*u(c,i,j+1,k)- u(c,i,  j,k)
      -(um(c,i,j+3,k)-3*um(c,i,j+2,k)+3*um(c,i,j+1,k)-um(c,i,  j,k)) )
     -3*(rho(i,j+1,k)*dcy(j+1)*jac(i,j+1,k)+rho(i,j,k)*dcy(j)*jac(i,j,k))*(
         u(c,i,j+2,k) -3*u(c,i,j+1,k)+ 3*u(c,i,  j,k)- u(c,i,j-1,k) 
      -(um(c,i,j+2,k)-3*um(c,i,j+1,k)+3*um(c,i,  j,k)-um(c,i,j-1,k)) )
     +3*(rho(i,j,k)*dcy(j)*jac(i,j,k)+rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k))*(
         u(c,i,j+1,k)- 3*u(c,i, j,k)+ 3*u(c,i,j-1,k)- u(c,i,j-2,k) 
      -(um(c,i,j+1,k)-3*um(c,i, j,k)+3*um(c,i,j-1,k)-um(c,i,j-2,k)) )
     - (rho(i,j-1,k)*dcy(j-1)*jac(i,j-1,k)+rho(i,j-2,k)*dcy(j-2)*jac(i,j-2,k))*(
         u(c,i, j,k)- 3*u(c,i,j-1,k)+  3*u(c,i,j-2,k)- u(c,i,j-3,k) 
      -(um(c,i, j,k)-3*um(c,i,j-1,k)+ 3*um(c,i,j-2,k)-um(c,i,j-3,k)) )
					     )  );
	       }
	    }
      }
#undef rho
#undef up
#undef u
#undef um
#undef strx
#undef dcx
#undef cox
#undef stry
#undef dcy
#undef coy
#undef jac
   }
}
