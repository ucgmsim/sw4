#include <algorithm>

#include "Sarray.h"
#include "EW.h"
#include "CurvilinearInterface2.h"
#include "TestTwilight.h"
#include "TestEcons.h"
#include "F77_FUNC.h"
#include "GridGenerator.h"

extern "C" {
   void F77_FUNC(dgetrf,DGETRF)(int*,  int*, double*, int*, int*, int*);
   void F77_FUNC(dgetri,DGETRI)(int*, double*, int*, int*, double*, int*, int* );
   void F77_FUNC(dgetrs,DGETRS)(char*, int*, int*, double*, int*, int*, double*, int*, int*);
   void F77_FUNC(sgetrf,SGETRF)(int*,  int*, float*, int*, int*, int*);
   void F77_FUNC(sgetri,SGETRI)(int*, float*, int*, int*, float*, int*, int* );
   void F77_FUNC(sgetrs,SGETRS)(char*, int*, int*, float*, int*, int*, float*, int*, int*);
}

// Precision-dispatched LAPACK wrappers: the factorized mass blocks are stored
// in float_sw4, so the LU routines must match that type in both builds.
static inline void getrf_sw4(int* m, int* n, double* a, int* lda, int* ipiv, int* info)
{ F77_FUNC(dgetrf,DGETRF)(m, n, a, lda, ipiv, info); }
static inline void getrf_sw4(int* m, int* n, float* a, int* lda, int* ipiv, int* info)
{ F77_FUNC(sgetrf,SGETRF)(m, n, a, lda, ipiv, info); }
static inline void getri_sw4(int* n, double* a, int* lda, int* ipiv, double* work, int* lwork, int* info)
{ F77_FUNC(dgetri,DGETRI)(n, a, lda, ipiv, work, lwork, info); }
static inline void getri_sw4(int* n, float* a, int* lda, int* ipiv, float* work, int* lwork, int* info)
{ F77_FUNC(sgetri,SGETRI)(n, a, lda, ipiv, work, lwork, info); }

void bndryOpNoGhostc( float_sw4 *acof_no_gp, float_sw4 *ghcof_no_gp, float_sw4 *sbop_no_gp );

void curvilinear4sgwind( int, int, int, int, int, int, int, int, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, int*, float_sw4*, float_sw4*, float_sw4*, float_sw4*,
                         float_sw4*, float_sw4*, float_sw4*, int, char );

//-----------------------------------------------------------------------
CurvilinearInterface2::CurvilinearInterface2( int a_gc, EW* a_ew )
{
   m_reltol= 1e-6;
   m_abstol= 1e-6;
   m_maxit = 30;
   m_lhs_cof_built = false;
   m_gc = a_gc;
   m_gf = a_gc+1;
   m_ew = a_ew;
   m_etest = a_ew->create_energytest();
   m_tw    = a_ew->create_twilight();
   m_psource = a_ew->get_point_source_test();
   m_nghost = 5;
   a_ew->GetStencilCoefficients( m_acof, m_ghcof, m_bop, m_bope, m_sbop );
   bndryOpNoGhostc( m_acof_no_gp, m_ghcof_no_gp, m_sbop_no_gp );
   for( int s=0 ; s < 4 ; s++ )
      m_isbndry[s] = true;
   m_use_attenuation   = a_ew->usingAttenuation();
   m_number_mechanisms = a_ew->getNumberOfMechanisms();
   m_memory_is_allocated = false;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::bnd_zero( Sarray& u, int npts )
{
// Homogeneous Dirichet at boundaries on sides. Do not apply at upper and lower boundaries.
   for( int s=0 ; s < 4 ; s++ )
      if( m_isbndry[s] )
      {
         int kb=u.m_kb, ke=u.m_ke, jb=u.m_jb, je=u.m_je, ib=u.m_ib, ie=u.m_ie;
         if( s == 0 )
            ie = ib+npts-1;
         if( s == 1 )
            ib = ie-npts+1;
         if( s == 2 )
            je = jb+npts-1;
         if( s == 3 )
            jb = je-npts+1;
         for(int c=1 ; c <= u.m_nc ; c++)
         for( int k=kb ; k <= ke ; k++ )
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
                  u(c,i,j,k)=0;
      }
}
 
//-----------------------------------------------------------------------
void CurvilinearInterface2::copy_str( float_sw4* dest, float_sw4* src,
				      int offset, int n, int nsw )
{
 //
 // Copy supergrid stretching function array into an array with different 
 // number of ghost points. The new values are filled in by constant extrapolation.
 // Input:  src    - Old array
 //         offset - Number of added points points in new array at the
 //                  lower end (<0 means fewer points)
 //         n      - Size of new array 
 //         nsw    - Size of old array
 //
 // Output: dest   - New array
 // Note: Giving n as input implicitly determines the number of new points at the upper end.
 //
   if( offset >= 0 )
   {
      for( int i=0; i < nsw ;i++)
         dest[i+offset] = src[i];
      for( int i=0; i< offset ; i++ )
         dest[i] = dest[offset];
      for( int i=nsw; i < n ; i++ )
         dest[i] = dest[nsw-1];
   }
   else
   {
      for( int i=0; i < min(n,nsw-offset) ;i++)
         dest[i] =src[i+offset];
   }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::init_arrays( vector<float_sw4*>& a_strx,
					 vector<float_sw4*>& a_stry,
                                         vector<Sarray>& a_rho, 
                                         vector<Sarray>& a_mu, 
                                         vector<Sarray>& a_lambda )
{
   if( m_memory_is_allocated )
   {
      delete[] m_strx_c, m_stry_c, m_strx_f, m_stry_f;
      delete[] m_mass_block, m_ipiv_block;
   }
   m_memory_is_allocated = true;
   for( int s=0 ; s < 4; s++ )
      m_isbndry[s] = m_ew->getLocalBcType( m_gc, s ) != bProcessor;

   m_ib = m_ew->m_iStartInt[m_gc]-m_nghost;
   m_ie = m_ew->m_iEndInt[m_gc]+m_nghost;
   m_jb = m_ew->m_jStartInt[m_gc]-m_nghost;
   m_je = m_ew->m_jEndInt[m_gc]+m_nghost;
   m_ibf= m_ew->m_iStartInt[m_gf]-m_nghost;
   m_ief= m_ew->m_iEndInt[m_gf]+m_nghost;
   m_jbf= m_ew->m_jStartInt[m_gf]-m_nghost;
   m_jef= m_ew->m_jEndInt[m_gf]+m_nghost;
   m_nkf= m_ew->m_global_nz[m_gf];

   m_kb  = 0;
   m_ke  = 8;
   m_kbf = m_nkf-7;
   m_kef = m_nkf+1;

   m_strx_c = new float_sw4[m_ie-m_ib+1];
   m_stry_c = new float_sw4[m_je-m_jb+1];
   m_strx_f = new float_sw4[m_ief-m_ibf+1];
   m_stry_f = new float_sw4[m_jef-m_jbf+1];

   int ndif = m_nghost-(m_ew->m_iStartInt[m_gc]-m_ew->m_iStart[m_gc]);
   int nsw=m_ew->m_iEnd[m_gc] - m_ew->m_iStart[m_gc]+1;
   copy_str( m_strx_c, a_strx[m_gc], ndif, m_ie-m_ib+1, nsw );
   communicate_array1d( m_strx_c, m_ie-m_ib+1, 0, m_nghost );

   ndif = m_nghost-(m_ew->m_iStartInt[m_gf]-m_ew->m_iStart[m_gf]);
   nsw  = m_ew->m_iEnd[m_gf] - m_ew->m_iStart[m_gf]+1;
   copy_str( m_strx_f, a_strx[m_gf], ndif, m_ief-m_ibf+1, nsw );
   communicate_array1d( m_strx_f, m_ief-m_ibf+1, 0, m_nghost );

   ndif = m_nghost-(m_ew->m_jStartInt[m_gc]-m_ew->m_jStart[m_gc]);
   nsw  = m_ew->m_jEnd[m_gc] - m_ew->m_jStart[m_gc]+1;
   copy_str( m_stry_c, a_stry[m_gc], ndif, m_je-m_jb+1, nsw );
   communicate_array1d( m_stry_c, m_je-m_jb+1, 1, m_nghost );

   ndif = m_nghost-(m_ew->m_jStartInt[m_gf]-m_ew->m_jStart[m_gf]);
   nsw=m_ew->m_jEnd[m_gf] - m_ew->m_jStart[m_gf]+1;
   copy_str( m_stry_f, a_stry[m_gf], ndif, m_jef-m_jbf+1, nsw );
   communicate_array1d( m_stry_f, m_jef-m_jbf+1, 1, m_nghost );

   m_rho_c.define(m_ib,m_ie,m_jb,m_je,1,1);
   m_rho_f.define(m_ibf,m_ief,m_jbf,m_jef,m_nkf,m_nkf);

   m_mu_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_lambda_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_jac_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);

   m_mu_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_lambda_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_jac_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);

   m_x_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_y_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_z_c.define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_met_c.define(4,m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   m_ew->m_gridGenerator->generate_grid_and_met( m_ew, m_gc, m_x_c, m_y_c, m_z_c, m_jac_c, m_met_c, false );   
   m_met_c.insert_intersection(m_ew->mMetric[m_gc]);
   m_jac_c.insert_intersection(m_ew->mJ[m_gc]);

   communicate_array( m_met_c, true );
   communicate_array( m_jac_c, true );

   m_x_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_y_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_z_f.define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_met_f.define(4,m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   m_ew->m_gridGenerator->generate_grid_and_met( m_ew, m_gf, m_x_f, m_y_f, m_z_f, m_jac_f, m_met_f, false );
   m_met_f.insert_intersection(m_ew->mMetric[m_gf]);
   m_jac_f.insert_intersection(m_ew->mJ[m_gf]);

   communicate_array( m_met_f, true );
   communicate_array( m_jac_f, true );

   if( m_tw != 0 )
   {
      m_tw->get_rho(m_rho_c,m_x_c,m_y_c,m_z_c);
      m_tw->get_rho(m_rho_f,m_x_f,m_y_f,m_z_f);
      m_tw->get_mula(m_mu_c,m_lambda_c,m_x_c,m_y_c,m_z_c);
      m_tw->get_mula(m_mu_f,m_lambda_f,m_x_f,m_y_f,m_z_f);
   }
   else 
   {
      //      m_rho_c.insert_intersection(m_ew->mRho[m_gc]);
      //      m_rho_f.insert_intersection(m_ew->mRho[m_gf]);
      //      m_mu_c.insert_intersection(m_ew->mMu[m_gc]);
      //      m_mu_f.insert_intersection(m_ew->mMu[m_gf]);
      //      m_lambda_c.insert_intersection(m_ew->mLambda[m_gc]);
      //      m_lambda_f.insert_intersection(m_ew->mLambda[m_gf]);

      m_rho_c.insert_intersection(a_rho[m_gc]);
      m_rho_f.insert_intersection(a_rho[m_gf]);
      m_mu_c.insert_intersection(a_mu[m_gc]);
      m_mu_f.insert_intersection(a_mu[m_gf]);
      m_lambda_c.insert_intersection(a_lambda[m_gc]);
      m_lambda_f.insert_intersection(a_lambda[m_gf]);

      int extra_ghost = m_nghost - m_ew->getNumberOfGhostPoints();
      if( extra_ghost > 0 )
      {
	 if( m_etest != 0 ) 
	 {
    	    int sides[6]={1,1,1,1,0,0};
            m_etest->get_rhobnd(m_rho_c,extra_ghost,sides);
            m_etest->get_rhobnd(m_rho_f,extra_ghost,sides);
            m_etest->get_mulabnd(m_mu_c,m_lambda_c,extra_ghost,sides);
            m_etest->get_mulabnd(m_mu_f,m_lambda_f,extra_ghost,sides);
	 }
	 else
	 {
	    m_rho_c.extrapolij(extra_ghost);
	    m_rho_f.extrapolij(extra_ghost);
	    m_mu_c.extrapolij(extra_ghost);
	    m_mu_f.extrapolij(extra_ghost);
	    m_lambda_c.extrapolij(extra_ghost);
	    m_lambda_f.extrapolij(extra_ghost);
	 }
      }
   }
   communicate_array( m_rho_c, true );
   communicate_array( m_mu_c, true );
   communicate_array( m_lambda_c, true );

   communicate_array( m_rho_f, true );
   communicate_array( m_mu_f, true );
   communicate_array( m_lambda_f, true );

   if( m_use_attenuation )
      init_arrays_att();

  // Matrix only defined at interior points
   m_Mass_block.define(9,m_ib+m_nghost,m_ie-m_nghost,m_jb+m_nghost,m_je-m_nghost,1,1);
   interface_block( m_Mass_block );

 // Repackage Mass_block into array of fortran order.
   int nimb = (m_Mass_block.m_ie-m_Mass_block.m_ib+1);
   size_t msize = nimb*(m_Mass_block.m_je-m_Mass_block.m_jb+1);
   m_mass_block = new float_sw4[9*msize];
   for( int j=m_jb+m_nghost ; j <= m_je-m_nghost ; j++ )
      for( int i=m_ib+m_nghost ; i <= m_ie-m_nghost ; i++ )
      {
	 size_t ind = (i-(m_ib+m_nghost)) + nimb*(j-(m_jb+m_nghost));
	 for( int c=1 ; c <= 9 ;c++)
	   m_mass_block[c-1+9*ind] = m_Mass_block(c,i,j,1);
      }
   int three    = 3;
   int info     = 0;
   m_ipiv_block = new int[3*msize];
   int lwork=9;
   float_sw4* work=new float_sw4[lwork];
   for( size_t ind=0 ; ind < msize; ind++ )
   {
      getrf_sw4(&three, &three, &m_mass_block[9*ind], &three,
                &m_ipiv_block[3*ind], &info );
      if( info != 0)
      {
	 int j = ind/m_Mass_block.m_ni+m_Mass_block.m_jb;
	 int i = ind + m_Mass_block.m_ib - m_Mass_block.m_ni*(j - m_Mass_block.m_jb);
         std::cerr << "LU Fails at (i,j) equals" << i << "," << j
                   << " info = " << info << " " << m_Mass_block(info+3*(info-1), i, j,1)
                      << "\n";
         for (int l = 1; l <= 3; l++) 
            for (int m = 1; m <= 3; m++)
	      std::cerr << m_Mass_block(m +3*(l-1), i, j,1) << ",";
	    std::cerr << "\n";
      }
      getri_sw4(&three, &m_mass_block[9*ind], &three,
                &m_ipiv_block[3*ind], work, &lwork, &info );
      if( info != 0)
      {
	 int j = ind/m_Mass_block.m_ni+m_Mass_block.m_jb;
	 int i = ind + m_Mass_block.m_ib - m_Mass_block.m_ni*(j - m_Mass_block.m_jb);
         std::cerr << "DGETRI Fails at (i,j) equals" << i << "," << j
                   << " info = " << info << " " << m_Mass_block(info+3*(info-1), i, j,1)
                      << "\n";
         for (int l = 1; l <= 3; l++) 
            for (int m = 1; m <= 3; m++)
	      std::cerr << m_Mass_block(m +3*(l-1), i, j,1) << ",";
	    std::cerr << "\n";
      }
   }
   delete[] work;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::init_arrays_att()
{
   // Attenuation material setup 
   if( m_use_attenuation )
   {
      m_muve_c.resize(m_number_mechanisms);
      m_lambdave_c.resize(m_number_mechanisms);
      for( int a=0 ; a < m_number_mechanisms; a++ )
      {
         m_muve_c[a].define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
         m_lambdave_c[a].define(m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
      }
      m_muve_f.resize(m_number_mechanisms);
      m_lambdave_f.resize(m_number_mechanisms);
      for( int a=0 ; a < m_number_mechanisms; a++ )
      {
         m_muve_f[a].define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
         m_lambdave_f[a].define(m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
      }
      if( m_tw != 0 )
      {
         // Note, twilight uses only one attenuation mechanism
         m_tw->get_mula_att(m_muve_c[0],m_lambdave_c[0],m_x_c,m_y_c,m_z_c);
         m_tw->get_mula_att(m_muve_f[0],m_lambdave_f[0],m_x_f,m_y_f,m_z_f);
      }
      else
      {
         for( int a=0 ; a < m_number_mechanisms; a++ )
         {
            m_muve_c[a].insert_intersection(m_ew->mMuVE[m_gc][a]);
            m_lambdave_c[a].insert_intersection(m_ew->mLambdaVE[m_gc][a]);
            m_muve_f[a].insert_intersection(m_ew->mMuVE[m_gf][a]);
            m_lambdave_f[a].insert_intersection(m_ew->mLambdaVE[m_gf][a]);
         }
         int extra_ghost = m_nghost - m_ew->getNumberOfGhostPoints();
         if( extra_ghost > 0 )
         {
            if( m_etest != 0 ) 
            {
               int sides[6]={1,1,1,1,0,0};
               for( int a=0 ; a < m_number_mechanisms; a++ )
               {
                  m_etest->get_mulabnd(m_muve_c[a],m_lambdave_c[a],extra_ghost,sides);
                  m_etest->get_mulabnd(m_muve_f[a],m_lambdave_f[a],extra_ghost,sides);
               }
            }
            else
            {
               for( int a=0 ; a < m_number_mechanisms; a++ )
               {
                  m_muve_c[a].extrapolij(extra_ghost);
                  m_muve_f[a].extrapolij(extra_ghost);
                  m_lambdave_c[a].extrapolij(extra_ghost);
                  m_lambdave_f[a].extrapolij(extra_ghost);
               }
            }
         }
      }
      for( int a=0 ; a < m_number_mechanisms; a++ )
      {
         communicate_array( m_muve_c[a], true );
         communicate_array( m_lambdave_c[a], true );         
         communicate_array( m_muve_f[a], true );
         communicate_array( m_lambdave_f[a], true );         
      }
   }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::impose_ic( std::vector<Sarray>& a_U, float_sw4 t,
                                       std::vector<Sarray>& a_F,
                                       std::vector<Sarray*>& a_AlphaVE, 
                                       bool injection_only )
{
   bool force_dirichlet = false; //, check_stress_cont=false;
   //   int fg=0;
   //   if( force_dirichlet )
   //      fg = 1;

   Sarray U_f(3,m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
   Sarray U_c(3,m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
   Sarray F_f(3,m_ibf,m_ief,m_jbf,m_jef,m_nkf,m_nkf);
   Sarray F_c(3,m_ib,m_ie,m_jb,m_je,1,1);
   vector<Sarray> Alpha_c, Alpha_f;

//  1. copy   a_U into U_f and U_c
   U_f.insert_intersection(a_U[m_gf]);
   U_c.insert_intersection(a_U[m_gc]);
   //   F_f.insert_intersection(a_F[m_gf]);
   //   F_c.insert_intersection(a_F[m_gc]);

   if( m_use_attenuation )
   {
      Alpha_c.resize(m_number_mechanisms);
      Alpha_f.resize(m_number_mechanisms);
      //      Alpha_c = new Sarray[m_number_mechanisms];
      //      Alpha_f = new Sarray[m_number_mechanisms];
      for( int a=0 ; a < m_number_mechanisms ; a++)
      {
         Alpha_f[a].define(3,m_ibf,m_ief,m_jbf,m_jef,m_kbf,m_kef);
         Alpha_c[a].define(3,m_ib,m_ie,m_jb,m_je,m_kb,m_ke);
         Alpha_f[a].insert_intersection(a_AlphaVE[m_gf][a]);
         Alpha_c[a].insert_intersection(a_AlphaVE[m_gc][a]);
      }
   }
// 2a. Impose dirichlet conditions at ghost points
   int sides[6]={m_isbndry[0],m_isbndry[1],m_isbndry[2],m_isbndry[3],0,0};
   if( m_tw != 0 )
   {
      m_tw->get_ubnd( U_f, m_x_f, m_y_f, m_z_f, t, m_nghost, sides );
      m_tw->get_ubnd( U_c, m_x_c, m_y_c, m_z_c, t, m_nghost, sides );
      if( m_use_attenuation )
      {
         m_tw->get_bnd_att( Alpha_f[0], m_x_f, m_y_f, m_z_f, t, m_nghost, sides );
         m_tw->get_bnd_att( Alpha_c[0], m_x_c, m_y_c, m_z_c, t, m_nghost, sides );
      }
      if( force_dirichlet )
      {
      // Debug
         sides[0]=sides[1]=sides[2]=sides[3]=sides[4]=0;
         sides[5]=1;
         m_tw->get_ubnd( U_f, m_x_f, m_y_f, m_z_f, t, m_nghost+1, sides );
         sides[0]=sides[1]=sides[2]=sides[3] = sides[5]=0;
         sides[4]=1;
         m_tw->get_ubnd( U_c, m_x_c, m_y_c, m_z_c, t, m_nghost+1, sides );
         a_U[m_gc].copy_kplane2(U_c,0); // have computed U_c:s ghost points
         a_U[m_gc].copy_kplane2(U_c,1); // have computed U_c:s ghost points
         a_U[m_gf].copy_kplane2(U_f,m_nkf);    // .. and U_f:s interface points
         a_U[m_gf].copy_kplane2(U_f,m_nkf+1);    // .. and U_f:s ghost point
         return;
      // End debug
      }
   }
   else if( m_etest != 0 )
   {
      m_etest->get_ubnd( U_f, m_nghost, sides );
      m_etest->get_ubnd( U_c, m_nghost, sides );
   }
   else if( m_psource != 0 )
   {
      m_psource->ubnd( U_f, m_x_f, m_y_f, m_z_f, t, m_ew->mGridSize[m_gf], m_nghost, sides );
      m_psource->ubnd( U_c, m_x_c, m_y_c, m_z_c, t, m_ew->mGridSize[m_gc], m_nghost, sides );
   }
   else
   {
      bnd_zero( U_c, m_nghost );
      bnd_zero( U_f, m_nghost );
      if( m_use_attenuation )
         for( int a=0 ; a < m_number_mechanisms ; a++)
         {
            bnd_zero( Alpha_c[a], m_nghost );
            bnd_zero( Alpha_f[a], m_nghost );           
         }
   }

// 3. Inject U_f := U_c on interface
   communicate_array( U_c, true );
   injection( U_f, U_c );
   communicate_array( U_f, true );

   if( m_use_attenuation )
      for( int a=0 ; a < m_number_mechanisms ; a++)
      {
         communicate_array(Alpha_c[a],true);
         injection(Alpha_f[a],Alpha_c[a]);
         communicate_array( Alpha_f[a], true );
      }

   if( injection_only )
   {
      a_U[m_gf].copy_kplane2(U_f,m_nkf);
      if( m_use_attenuation )
      {
         for( int a=0 ; a < m_number_mechanisms ; a++ )
            a_AlphaVE[m_gf][a].copy_kplane2(Alpha_f[a],m_nkf);
      }
      return;
   }
// 4. Solve equation for stress continuity, formulated as lhs*x+rhs=0, where x are uc's ghost points at k=0.

   // 4.a-b The ghost-plane unknown (xd), fixed forcing term (rhsd), lhs
   // (lhsd) and residual (resd) are held in double regardless of float_sw4,
   // because the residual floor here is 1 ULP of float_sw4 -- unattainable
   // at the requested reltol/abstol in single precision (see
   // convergence_plan.md). Coefficients (metric, jacobian, mu, lambda, rho,
   // mass-block inverse) are read at their native float_sw4 precision and
   // promoted; rhsd is assembled by interface_rhs_d, which likewise
   // promotes each interior float_sw4 result before combining it, so rhsd
   // carries no additional float_sw4 rounding beyond what the underlying
   // (float_sw4-stored) interior solution already has.
   DPlane xd(3,m_ib,m_ie,m_jb,m_je), rhsd(3,m_ib,m_ie,m_jb,m_je);
   DPlane lhsd(3,m_ib,m_ie,m_jb,m_je), resd(3,m_ib,m_ie,m_jb,m_je);
   interface_rhs_d( rhsd, U_c, U_f, F_c, F_f, Alpha_c, Alpha_f );
#pragma omp parallel for collapse(2)
   for( int c=1 ; c <= 3 ;c++)
      for( int j=m_jb ; j <= m_je ; j++ )
         for( int i=m_ib ; i <= m_ie ; i++ )
            xd(c,i,j)   = U_c(c,i,j,0);
   interface_lhs_d( lhsd, xd );

   // Initial residual
   double maxresloc=0;
#pragma omp parallel for collapse(2) reduction(max:maxresloc)
   for( int c=1 ; c <= 3 ;c++)
     for( int j=m_jb+5 ; j <= m_je-5 ; j++ )
       for( int i=m_ib+5 ; i <= m_ie-5 ; i++ )
	 {
	   resd(c,i,j) = lhsd(c,i,j)+rhsd(c,i,j);
	   if( abs(resd(c,i,j)) > maxresloc )
	     maxresloc = abs(resd(c,i,j));
	 }
   double maxres=maxresloc;
   MPI_Allreduce( &maxresloc, &maxres, 1, MPI_DOUBLE, MPI_MAX, m_ew->m_cartesian_communicator );

   // 4.c Jacobi iteration
   float_sw4 scalef=(m_ew->m_global_nx[m_gc]-1)*(m_ew->m_global_ny[m_gc]-1); //scale residual to be size O(1).
   //   int maxit = 30;
   //   float_sw4 reltol=1e-6, abstol=1e-6;
   int iter = 0;
   int nimb = m_Mass_block.m_ie-m_Mass_block.m_ib+1;
   // Block Jacobi, lhs*x+rhs=0 and lhs=M+N --> M*xp+N*x+rhs=0 --> M*(xp-x)+lhs*x+rhs=0
   //        --> xp-x=-inv(M)*(lhs*x+rhs) --> xp = x - inv(M)*(lhs*x+rhs)
   double maxres0 = maxres;
   double relax = 1.0;
   while( maxres > m_reltol*maxres0 && scalef*maxres > m_abstol && iter <= m_maxit )
   {
      iter++;
#pragma omp parallel for
      for( int j=m_Mass_block.m_jb ; j <= m_Mass_block.m_je ; j++ )
         for( int i=m_Mass_block.m_ib ; i <= m_Mass_block.m_ie ; i++ )
	 {
	    size_t ind=(i-m_Mass_block.m_ib)+nimb*(j-m_Mass_block.m_jb);
            double x1, x2, x3;
            double b1=resd(1,i,j), b2=resd(2,i,j), b3=resd(3,i,j);
            x1 = m_mass_block[9*ind  ]*b1+
                 m_mass_block[9*ind+3]*b2+
                 m_mass_block[9*ind+6]*b3;
            x2 = m_mass_block[9*ind+1]*b1+
                 m_mass_block[9*ind+4]*b2+
                 m_mass_block[9*ind+7]*b3;
            x3 = m_mass_block[9*ind+2]*b1+
                 m_mass_block[9*ind+5]*b2+
                 m_mass_block[9*ind+8]*b3;
	    xd(1,i,j) -= relax*x1;
	    xd(2,i,j) -= relax*x2;
	    xd(3,i,j) -= relax*x3;
	 }

  // 4.d Communicate xd here (only the k=0 ghost plane)
      communicate_plane_d( xd );
      interface_lhs_d( lhsd, xd );

// 4.e. Compute residual and its norm
      maxresloc=0;
#pragma omp parallel for collapse(2) reduction(max:maxresloc)
      for( int c=1 ; c <= 3 ;c++)
	  for( int j=m_jb+5 ; j <= m_je-5 ; j++ )
	     for( int i=m_ib+5 ; i <= m_ie-5 ; i++ )
	     {
	        resd(c,i,j) = lhsd(c,i,j)+rhsd(c,i,j);
	        if( abs(resd(c,i,j)) > maxresloc )
	           maxresloc = abs(resd(c,i,j));
	     }
      MPI_Allreduce( &maxresloc, &maxres, 1, MPI_DOUBLE, MPI_MAX, m_ew->m_cartesian_communicator);
   }
   if( (maxres > m_reltol*maxres0 && scalef*maxres > m_abstol) && m_ew->getRank()==0 )
   {
      std::cout << "WARNING, no convergence in curvilinear interface, res = "
                << maxres << " reltol= " << m_reltol << " initial res = " << maxres0
                << std::endl;
      std::cout << "     scaled res = " << scalef*maxres << " abstol= " << m_abstol
                << std::endl;
   }
   // Write the converged double-precision ghost-plane solution back into
   // U_c at its native storage precision (float_sw4).
   for( int c=1 ; c <= 3 ;c++)
      for( int j=m_jb ; j <= m_je ; j++ )
         for( int i=m_ib ; i <= m_ie ; i++ )
            U_c(c,i,j,0) = (float_sw4)xd(c,i,j);

// 5. Copy U_c and U_f back to a_U, only k=0 for U_c and k=n3f for U_f.
   a_U[m_gc].copy_kplane2(U_c,0);     // have computed U_c:s ghost points
   a_U[m_gf].copy_kplane2(U_f,m_nkf);   // .. and U_f:s interface points
   if( m_use_attenuation )
   {
      for( int a=0 ; a < m_number_mechanisms ; a++ )
         a_AlphaVE[m_gf][a].copy_kplane2(Alpha_f[a],m_nkf);
   }
}


//-----------------------------------------------------------------------
void CurvilinearInterface2::injection(Sarray &u_f, Sarray &u_c )
{
  // Injection at the interface

  const float_sw4 a= 9.0/16;
  const float_sw4 b=-1.0/16;
  //  const int ngh = m_nghost;
  int i1=u_c.m_ib+m_nghost-1, i2=u_c.m_ie-m_nghost+1;
  int j1=u_c.m_jb+m_nghost-1, j2=u_c.m_je-m_nghost+1;
  if( m_isbndry[0] )
     i1++;
  if( m_isbndry[1] )
     i2 -= 2;
  if( m_isbndry[2] )
     j1++;
  if( m_isbndry[3] )
     j2 -= 2;
  
  for (int l = 1; l <= u_c.m_nc; l++) 
    //    for (int j = u_c.m_jb+ngh-1; j <= u_c.m_je-ngh+1; j++)
    //      for (int i = u_c.m_ib+ngh-1; i <= u_c.m_ie-ngh+1; i++) 
    for (int j = j1; j <= j2; j++)
      for (int i = i1; i <= i2; i++) 
      {
        u_f(l, 2 * i - 1, 2 * j - 1, m_nkf) = u_c(l, i, j, 1);
        u_f(l, 2 * i, 2 * j - 1, m_nkf) =
            b * u_c(l, i - 1, j, 1) + a * u_c(l, i, j, 1) +
            a * u_c(l, i + 1, j, 1) + b * u_c(l, i + 2, j, 1);
        u_f(l, 2 * i - 1, 2 * j, m_nkf) =
            b * u_c(l, i, j - 1, 1) + a * u_c(l, i, j, 1) +
            a * u_c(l, i, j + 1, 1) + b * u_c(l, i, j + 2, 1);
        u_f(l, 2 * i, 2 * j, m_nkf) =
            b * ( b * u_c(l, i - 1, j - 1, 1) +
                  a * u_c(l, i,     j - 1, 1) +
                  a * u_c(l, i + 1, j - 1, 1) +
                  b * u_c(l, i + 2, j - 1, 1)) +
            a * ( b * u_c(l, i - 1, j, 1) + 
                  a * u_c(l, i,     j, 1) +
                  a * u_c(l, i + 1, j, 1) +
                  b * u_c(l, i + 2, j, 1)) +
            a * ( b * u_c(l, i - 1, j + 1, 1) +
                  a * u_c(l, i,     j + 1, 1) +
                  a * u_c(l, i + 1, j + 1, 1) +
                  b * u_c(l, i + 2, j + 1, 1)) +
            b * ( b * u_c(l, i - 1, j + 2, 1) +
                  a * u_c(l, i,     j + 2, 1) +
                  a * u_c(l, i + 1, j + 2, 1) +
                  b * u_c(l, i + 2, j + 2, 1));
      }
  if( m_isbndry[1] )
  {
     int i=i2+1; 
     for (int l = 1; l <= u_c.m_nc; l++) 
        for (int j = j1; j <= j2; j++)
        {
           u_f(l, 2 * i - 1, 2 * j - 1, m_nkf) = u_c(l, i, j, 1);
           u_f(l, 2 * i - 1, 2 * j,     m_nkf) =
              b * u_c(l, i, j - 1, 1) + a * u_c(l, i, j, 1) +
              a * u_c(l, i, j + 1, 1) + b * u_c(l, i, j + 2, 1);
        }
  }
  if( m_isbndry[3] )
  {
     int j=j2+1; 
     for (int l = 1; l <= u_c.m_nc; l++) 
        for (int i = i1; i <= i2; i++)
        {
           u_f(l, 2 * i - 1, 2 * j - 1, m_nkf) = u_c(l, i, j, 1);
           u_f(l, 2 * i,     2 * j - 1, m_nkf) =
              b * u_c(l, i - 1, j, 1) + a * u_c(l, i, j, 1) +
              a * u_c(l, i + 1, j, 1) + b * u_c(l, i + 2, j, 1);
        }
  }
  if( m_isbndry[3] && m_isbndry[1] )
  {
     int i=i2+1; 
     int j=j2+1; 
     for (int l = 1; l <= u_c.m_nc; l++)
        u_f(l, 2 * i - 1, 2 * j - 1, m_nkf) = u_c(l, i, j, 1);
  }
}

//-----------------------------------------------------------------------
 void CurvilinearInterface2::interface_block( Sarray& matrix )
{
   const float_sw4 w1=17.0/48;
   matrix_Lu( m_ib, m_jb, matrix, m_met_c, m_jac_c, m_mu_c, m_lambda_c,
	      m_strx_c, m_stry_c, m_ghcof[0] );

   for( int c=1 ; c <= 9; c++ )
      for( int j=matrix.m_jb ; j <= matrix.m_je ; j++ )
         for( int i=matrix.m_ib ; i <= matrix.m_ie ; i++ )
	   matrix(c,i,j,1) /= m_rho_c(i,j,1);

   Sarray alpha(m_ibf,m_ief,m_jbf,m_jef,m_nkf,m_nkf);
   for( int j=alpha.m_jb ; j <= alpha.m_je ; j++ )
     for( int i=alpha.m_ib ; i <= alpha.m_ie ; i++ )
       alpha(i,j,m_nkf) = w1*m_jac_f(i,j,m_nkf)*m_rho_f(i,j,m_nkf)/(m_strx_f[i-m_ibf]*m_stry_f[j-m_jbf]);
   if( !m_tw && !m_psource )
      bnd_zero(alpha,m_nghost);
   restprol2D( matrix, alpha, 1, m_nkf );

 // Add -B(uc) contribution to block matrix
   mat_icstresses_curv( m_ib, m_jb, matrix, 1, m_met_c, m_mu_c, m_lambda_c,
			m_strx_c, m_stry_c, m_sbop );
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::interface_lhs_d( DPlane& lhsd, DPlane& xd )
{
   const double w1=17.0/48;
   lhs_Lu_d( xd, lhsd );

#pragma omp parallel for collapse(2)
   for( int c=1 ; c <= 3; c++ )
      for( int j=lhsd.jb ; j <= lhsd.je ; j++ )
         for( int i=lhsd.ib ; i <= lhsd.ie ; i++ )
	    lhsd(c,i,j) /= m_rho_c(i,j,1);
   if( !m_tw && !m_psource )
      bnd_zero_d(lhsd,m_nghost);

   // Reused rather than heap-allocated per call; DPlane's sizing constructor
   // zero-fills, so the buffers are zeroed here to keep that guarantee for
   // whatever prolongate2D_d / lhs_icstresses_curv_d leave untouched.
   if( m_prol_scratch.nc == 0 )
      m_prol_scratch = DPlane(3,m_ibf,m_ief,m_jbf,m_jef);
   else
      std::fill( m_prol_scratch.data.begin(), m_prol_scratch.data.end(), 0.0 );
   DPlane& prollhsd = m_prol_scratch;
   prolongate2D_d( lhsd, prollhsd );
#pragma omp parallel for collapse(2)
   for( int c=1 ; c <= 3 ;c++)
      for( int j=prollhsd.jb ; j <= prollhsd.je ; j++ )
         for( int i=prollhsd.ib ; i <= prollhsd.ie ; i++ )
	   prollhsd(c,i,j) = w1*m_jac_f(i,j,m_nkf)*m_rho_f(i,j,m_nkf)*prollhsd(c,i,j)/
	     (m_strx_f[i-m_ibf]*m_stry_f[j-m_jbf]);
   if( !m_tw && !m_psource )
      bnd_zero_d(prollhsd,m_nghost);
   restrict2D_d( lhsd, prollhsd );

   if( m_bc_scratch.nc == 0 )
      m_bc_scratch = DPlane(3,lhsd.ib,lhsd.ie,lhsd.jb,lhsd.je);
   else
      std::fill( m_bc_scratch.data.begin(), m_bc_scratch.data.end(), 0.0 );
   DPlane& Bc = m_bc_scratch;
   lhs_icstresses_curv_d( xd, Bc );
#pragma omp parallel for collapse(2)
   for( int c=1 ; c <= 3; c++ )
      for( int j=lhsd.jb ; j <= lhsd.je ; j++ )
         for( int i=lhsd.ib ; i <= lhsd.ie ; i++ )
	   lhsd(c,i,j) -= Bc(c,i,j);
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::interface_rhs_d( DPlane& rhsd, Sarray& uc, Sarray& uf,
                                             Sarray& fc, Sarray& ff,
                             vector<Sarray>& Alpha_c, vector<Sarray>& Alpha_f )
{
   // Assembles the fixed forcing term of lhs(x)+rhs=0 that interface_lhs_d
   // solves against. The interior stencil evaluations (curvilinear4sgwind,
   // compute_icstresses_curv) run at float_sw4 precision -- they read uc/uf,
   // which are themselves stored at float_sw4 precision throughout the rest
   // of the simulation, so recomputing them in double would only duplicate
   // those kernels without reducing error. But once each float_sw4 result
   // (L(uc)/rho, L(uf), B(uf), B(uc)) is produced, it is promoted to double
   // immediately and all further interface-local assembly -- prolongation,
   // restriction, the linear combination, and the accumulation into rhsd --
   // runs in double via the same kernels used by interface_lhs_d
   // (prolongate2D_d, restrict2D_d, bnd_zero_d), instead of rounding to
   // float_sw4 at each intermediate step as the old float-only interface_rhs
   // did. This keeps the fixed rhsd fed to the double Jacobi solve as
   // accurate as the float_sw4 inputs allow.
   Sarray utmp(3,uc.m_ib,uc.m_ie,uc.m_jb,uc.m_je,0,0);
   Sarray rhs(3,m_ib,m_ie,m_jb,m_je,1,1);

   const double w1=17.0/48;
//  1. Set ghost points to zero, and save the old value to restore after, so
//     that the routine does not change uc.
   for( int c=1 ; c <= 3; c++ )
      for( int j=uc.m_jb ; j <= uc.m_je ; j++ )
         for( int i=uc.m_ib ; i <= uc.m_ie ; i++ )
	   {
	       utmp(c,i,j,0) = uc(c,i,j,0);
               uc(c,i,j,0)   = 0;
	   }

   int onesided[6]={0,0,0,0,1,1};
// 2. Compute L(uc)/rhoc (float_sw4, unavoidable -- interior stencil on uc)
   curvilinear4sgwind( m_ib, m_ie, m_jb, m_je, m_kb, m_ke, 1, 1, uc.c_ptr(),
		       m_mu_c.c_ptr(), m_lambda_c.c_ptr(),
                       m_met_c.c_ptr(), m_jac_c.c_ptr(), rhs.c_ptr(),
                       onesided, m_acof, m_bope, m_ghcof, m_acof_no_gp,
                       m_ghcof_no_gp, m_strx_c, m_stry_c, 8, '=');
   if( m_use_attenuation )
      for( int a=0 ; a < m_number_mechanisms ; a++ )
         curvilinear4sgwind( m_ib, m_ie, m_jb, m_je, m_kb, m_ke, 1, 1, Alpha_c[a].c_ptr(),
                             m_muve_c[a].c_ptr(), m_lambdave_c[a].c_ptr(),
                             m_met_c.c_ptr(), m_jac_c.c_ptr(), rhs.c_ptr(),
                             onesided, m_acof_no_gp, m_bope, m_ghcof_no_gp,
                             m_acof_no_gp, m_ghcof_no_gp, m_strx_c, m_stry_c, 8, '-');

   for( int c=1 ; c <= 3; c++ )
      for( int j=rhs.m_jb ; j <= rhs.m_je ; j++ )
         for( int i=rhs.m_ib ; i <= rhs.m_ie ; i++ )
            rhs(c,i,j,1) /= m_rho_c(i,j,1);

// Promote to double: no further rounding to float_sw4 happens from here on.
   for( int c=1 ; c <= 3; c++ )
      for( int j=rhs.m_jb ; j <= rhs.m_je ; j++ )
         for( int i=rhs.m_ib ; i <= rhs.m_ie ; i++ )
            rhsd(c,i,j) = rhs(c,i,j,1);
   if( !m_tw && !m_psource )
      bnd_zero_d(rhsd,m_nghost);

// 3. Compute prolrhsd := p(L(uc)/rhoc)  (double)
   DPlane prolrhsd(3,m_ibf,m_ief,m_jbf,m_jef);
   prolongate2D_d( rhsd, prolrhsd );

// 4. Compute L(uf) (float_sw4, unavoidable -- interior stencil on uf), promote to double
   Sarray Luf(3,m_ibf,m_ief,m_jbf,m_jef,m_nkf,m_nkf);
   curvilinear4sgwind( m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef, m_nkf, m_nkf, uf.c_ptr(),
		       m_mu_f.c_ptr(), m_lambda_f.c_ptr(), m_met_f.c_ptr(), m_jac_f.c_ptr(), Luf.c_ptr(),
                       onesided, m_acof, m_bope, m_ghcof, m_acof_no_gp,
                       m_ghcof_no_gp, m_strx_f, m_stry_f, m_nkf, '=');
   if( m_use_attenuation )
      for( int a=0 ; a < m_number_mechanisms ; a++ )
         curvilinear4sgwind( m_ibf, m_ief, m_jbf, m_jef, m_kbf, m_kef, m_nkf, m_nkf, Alpha_f[a].c_ptr(),
		       m_muve_f[a].c_ptr(), m_lambdave_f[a].c_ptr(), m_met_f.c_ptr(),
                       m_jac_f.c_ptr(), Luf.c_ptr(),
                       onesided, m_acof_no_gp, m_bope, m_ghcof_no_gp, m_acof_no_gp,
                       m_ghcof_no_gp, m_strx_f, m_stry_f, m_nkf, '-');

   DPlane lufd(3,m_ibf,m_ief,m_jbf,m_jef);
   for( int c=1 ; c <= 3 ;c++)
      for( int j=Luf.m_jb ; j <= Luf.m_je ; j++ )
         for( int i=Luf.m_ib ; i <= Luf.m_ie ; i++ )
            lufd(c,i,j) = Luf(c,i,j,m_nkf);

// 5. Compute B(uf) (float_sw4, unavoidable -- interior stencil on uf), promote to double
   Sarray Bf(3,m_ibf,m_ief,m_jbf,m_jef,m_nkf,m_nkf);
   compute_icstresses_curv( uf, Bf, m_nkf, m_met_f, m_mu_f, m_lambda_f,
			    m_strx_f, m_stry_f, m_sbop_no_gp, '=' );
   if( m_use_attenuation )
      for( int a=0 ; a < m_number_mechanisms ; a++ )
         compute_icstresses_curv( Alpha_f[a], Bf, m_nkf, m_met_f, m_muve_f[a], m_lambdave_f[a],
                                  m_strx_f, m_stry_f, m_sbop_no_gp, '-' );

   DPlane bfd(3,m_ibf,m_ief,m_jbf,m_jef);
   for( int c=1 ; c <= 3 ;c++)
      for( int j=Bf.m_jb ; j <= Bf.m_je ; j++ )
         for( int i=Bf.m_ib ; i <= Bf.m_ie ; i++ )
            bfd(c,i,j) = Bf(c,i,j,m_nkf);

// 6. Form term prolrhsd := r(w1*J[gf]*(rhof*p(L(uc)/rhoc)-L(uf))+B(uf))  (double)
   for( int c=1 ; c <= 3 ;c++)
      for( int j=prolrhsd.jb ; j <= prolrhsd.je ; j++ )
         for( int i=prolrhsd.ib ; i <= prolrhsd.ie ; i++ )
            prolrhsd(c,i,j) = w1*(double)m_jac_f(i,j,m_nkf)*( (double)m_rho_f(i,j,m_nkf)*prolrhsd(c,i,j)-
	       lufd(c,i,j))/((double)m_strx_f[i-m_ibf]*(double)m_stry_f[j-m_jbf])+bfd(c,i,j);
   if( !m_tw && !m_psource )
      bnd_zero_d(prolrhsd,m_nghost);
   restrict2D_d( rhsd, prolrhsd );

// 7. Compute B(uc) (float_sw4, unavoidable -- interior stencil on uc), and
//    form rhsd := rhsd - B(uc) = r(w1*J[gf]*(rhof*p(L(uc)/rhoc)-L(uf))+B(uf))-B(uc)  (double)
   Sarray Bc(rhs);
   compute_icstresses_curv( uc, Bc, 1,  m_met_c, m_mu_c, m_lambda_c,
			    m_strx_c, m_stry_c, m_sbop, '=' );
   if( m_use_attenuation )
      for( int a=0 ; a < m_number_mechanisms ; a++ )
         compute_icstresses_curv( Alpha_c[a], Bc, 1, m_met_c, m_muve_c[a], m_lambdave_c[a],
                                  m_strx_c, m_stry_c, m_sbop_no_gp, '-' );

   for( int c=1 ; c <= 3 ;c++)
      for( int j=rhsd.jb ; j <= rhsd.je ; j++ )
         for( int i=rhsd.ib ; i <= rhsd.ie ; i++ )
	   rhsd(c,i,j) -= (double)Bc(c,i,j,1);

// 8. Restore ghost point values to U.
   for( int c=1 ; c <= 3; c++ )
      for( int j=uc.m_jb ; j <= uc.m_je ; j++ )
         for( int i=uc.m_ib ; i <= uc.m_ie ; i++ )
            uc(c,i,j,0) = utmp(c,i,j,0);
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::compute_icstresses_curv( Sarray& a_Up, Sarray& B, int kic,
						     Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
						     float_sw4* a_str_x, float_sw4* a_str_y, 
                                                     float_sw4* sbop, char op )
{
   const float_sw4 a1=2.0/3, a2=-1.0/12;
   const bool upper = (kic == 1);
   const int k=kic;
   const int kl = upper ? 1 :-1;
   const int ifirst = a_Up.m_ib;
   const int jfirst = a_Up.m_jb;
#define str_x(i) a_str_x[(i-ifirst)]   
#define str_y(j) a_str_y[(j-jfirst)]   
   float_sw4 sgn=1;
   if( op == '=' )
   {
      B.set_value(0.0);
      sgn = 1;
   }
   if( op == '-' )
   {
      sgn = -1;
   }


#pragma omp parallel for
   for( int j=B.m_jb+2 ; j <= B.m_je-2 ; j++ )
#pragma omp simd
      for( int i=B.m_ib+2 ; i <= B.m_ie-2 ; i++ )
      {
	 float_sw4 uz, vz, wz;	 
	 uz = vz = wz = 0;
         for( int m=0 ; m <= 5 ; m++ )
         {
            uz += sbop[m]*a_Up(1,i,j,k+kl*(m-1));
            vz += sbop[m]*a_Up(2,i,j,k+kl*(m-1));
            wz += sbop[m]*a_Up(3,i,j,k+kl*(m-1));
         }
         uz *=kl;
         vz *=kl;
         wz *=kl;

         // Normal terms
         float_sw4 m2 = str_x(i)*a_metric(2,i,j,k);
         float_sw4 m3 = str_y(j)*a_metric(3,i,j,k);
         float_sw4 m4 = a_metric(4,i,j,k);
         float_sw4 un = m2*uz+m3*vz+m4*wz;
         float_sw4 mnrm = m2*m2+m3*m3+m4*m4;
         float_sw4 B1, B2, B3;

         B1 = a_mu(i,j,k)*mnrm*uz + (a_mu(i,j,k)+a_lambda(i,j,k))*m2*un;
         B2 = a_mu(i,j,k)*mnrm*vz + (a_mu(i,j,k)+a_lambda(i,j,k))*m3*un;
         B3 = a_mu(i,j,k)*mnrm*wz + (a_mu(i,j,k)+a_lambda(i,j,k))*m4*un;

         // Tangential terms
         // p-derivatives
         float_sw4 up1=str_x(i)*(a2*(a_Up(1,i+2,j,k)-a_Up(1,i-2,j,k))+a1*(a_Up(1,i+1,j,k)-a_Up(1,i-1,j,k)));   
         float_sw4 up2=str_x(i)*(a2*(a_Up(2,i+2,j,k)-a_Up(2,i-2,j,k))+a1*(a_Up(2,i+1,j,k)-a_Up(2,i-1,j,k)));   
         float_sw4 up3=str_x(i)*(a2*(a_Up(3,i+2,j,k)-a_Up(3,i-2,j,k))+a1*(a_Up(3,i+1,j,k)-a_Up(3,i-1,j,k)));
         B1 += a_metric(1,i,j,k)*( (2*a_mu(i,j,k)+a_lambda(i,j,k))*m2*up1 + a_mu(i,j,k)*(m3*up2 + m4*up3));
         B2 += a_metric(1,i,j,k)*( a_lambda(i,j,k)*m3*up1 + a_mu(i,j,k)*m2*up2 );
         B3 += a_metric(1,i,j,k)*( a_lambda(i,j,k)*m4*up1 + a_mu(i,j,k)*m2*up3 );
         
         // q-derivatives
         float_sw4 uq1=str_y(j)*(a2*(a_Up(1,i,j+2,k)-a_Up(1,i,j-2,k))+a1*(a_Up(1,i,j+1,k)-a_Up(1,i,j-1,k)));
         float_sw4 uq2=str_y(j)*(a2*(a_Up(2,i,j+2,k)-a_Up(2,i,j-2,k))+a1*(a_Up(2,i,j+1,k)-a_Up(2,i,j-1,k)));
         float_sw4 uq3=str_y(j)*(a2*(a_Up(3,i,j+2,k)-a_Up(3,i,j-2,k))+a1*(a_Up(3,i,j+1,k)-a_Up(3,i,j-1,k)));
         B1 += a_metric(1,i,j,k)*( a_lambda(i,j,k)*m2*uq2 + a_mu(i,j,k)*m3*uq1);
         B2 += a_metric(1,i,j,k)*( (2*a_mu(i,j,k)+a_lambda(i,j,k))*m3*uq2 +a_mu(i,j,k)*(m2*uq1 + m4*uq3));
         B3 += a_metric(1,i,j,k)*( a_lambda(i,j,k)*m4*uq2 + a_mu(i,j,k)*m3*uq3);

         float_sw4 isgxy = 1.0/(str_x(i)*str_y(j));
         B1 *= isgxy;
         B2 *= isgxy;
         B3 *= isgxy;

         B(1,i,j,k) += sgn*B1;
         B(2,i,j,k) += sgn*B2;
         B(3,i,j,k) += sgn*B3;
      }
#undef str_x
#undef str_y
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::lhs_icstresses_curv_d( DPlane& xd, DPlane& Bc )
{
   // As lhs_icstresses_curv (float version, removed), but accumulates in
   // double and reads the unknown ghost-plane from xd instead of a Sarray.
   // Only ever called with kic==1 (upper==true) from interface_lhs_d.
   const int k=1;
   const int ifirst = m_ib;
   const int jfirst = m_jb;
#define str_x_d(i) m_strx_c[(i-ifirst)]
#define str_y_d(j) m_stry_c[(j-jfirst)]

#pragma omp parallel for
   for( int j=Bc.jb ; j <= Bc.je ; j++ )
#pragma omp simd
      for( int i=Bc.ib ; i <= Bc.ie ; i++ )
      {
	 double uz = m_sbop[0]*xd(1,i,j);
	 double vz = m_sbop[0]*xd(2,i,j);
	 double wz = m_sbop[0]*xd(3,i,j);

         // Normal terms
         double m2 = str_x_d(i)*(double)m_met_c(2,i,j,k);
         double m3 = str_y_d(j)*(double)m_met_c(3,i,j,k);
         double m4 = (double)m_met_c(4,i,j,k);
         double un   = m2*uz + m3*vz + m4*wz;
         double mnrm = m2*m2 + m3*m3 + m4*m4;
         double mu = m_mu_c(i,j,k), la = m_lambda_c(i,j,k);

         Bc(1,i,j) = mu*mnrm*uz + (mu+la)*m2*un;
         Bc(2,i,j) = mu*mnrm*vz + (mu+la)*m3*un;
         Bc(3,i,j) = mu*mnrm*wz + (mu+la)*m4*un;

         double isgxy = 1.0/(str_x_d(i)*str_y_d(j));
         Bc(1,i,j) *= isgxy;
         Bc(2,i,j) *= isgxy;
         Bc(3,i,j) *= isgxy;
      }
#undef str_x_d
#undef str_y_d
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::lhs_Lu_d( DPlane& xd, DPlane& lhsd )
{
   // As lhs_Lu (float version, removed), but accumulates in double and
   // reads the unknown ghost-plane from xd instead of a Sarray.
   const int ifirst = m_ib;
   const int jfirst = m_jb;
#define strx_d(i) m_strx_c[(i-ifirst)]
#define stry_d(j) m_stry_c[(j-jfirst)]
   // The six stress coefficients and ijac are functions of mu, lambda, the
   // metric, the Jacobian and the supergrid stretching only -- none of which
   // change in time -- so they are built once here instead of on every Jacobi
   // iteration of every solve of every timestep. Same expressions in the same
   // order, so the stored values are bit-identical to the ones the inline
   // version produced.
   // Rebuild if a future call site ever passes different bounds; today both
   // call sites are impose_ic's single lhsd, sized (3,m_ib,m_ie,m_jb,m_je).
   if( m_lhs_cof_built && ( m_lhs_cof.ib != lhsd.ib || m_lhs_cof.ie != lhsd.ie ||
                            m_lhs_cof.jb != lhsd.jb || m_lhs_cof.je != lhsd.je ) )
      m_lhs_cof_built = false;
   if( !m_lhs_cof_built )
   {
      m_lhs_cof = DPlane( 7, lhsd.ib, lhsd.ie, lhsd.jb, lhsd.je );
#pragma omp parallel for
      for( int j=lhsd.jb; j <= lhsd.je ;j++ )
	 for( int i=lhsd.ib; i <= lhsd.ie ;i++ )
	 {
	    double mu = m_mu_c(i,j,1), la = m_lambda_c(i,j,1);
	    double met2 = m_met_c(2,i,j,1), met3 = m_met_c(3,i,j,1), met4 = m_met_c(4,i,j,1);
	    double sx = strx_d(i), sy = stry_d(j);
	    m_lhs_cof(1,i,j) = ((2*mu+la)*met2*sx*met2*sx + mu*(met3*sy*met3*sy+met4*met4));
	    m_lhs_cof(2,i,j) = ((2*mu+la)*met3*sy*met3*sy + mu*(met2*sx*met2*sx+met4*met4));
	    m_lhs_cof(3,i,j) = ((2*mu+la)*met4*met4 + mu*(met2*sx*met2*sx+met3*sy*met3*sy));
	    m_lhs_cof(4,i,j) = (mu+la)*met2*met3*sx*sy;
	    m_lhs_cof(5,i,j) = (mu+la)*met2*met4*sx;
	    m_lhs_cof(6,i,j) = (mu+la)*met3*met4*sy;
	    m_lhs_cof(7,i,j) = m_ghcof[0]/(double)m_jac_c(i,j,1);
	 }
      m_lhs_cof_built = true;
   }
#pragma omp parallel for
   for( int j=lhsd.jb; j <= lhsd.je ;j++ )
#pragma omp simd
      for( int i=lhsd.ib; i <= lhsd.ie ;i++ )
      {
         double mucofu2 = m_lhs_cof(1,i,j), mucofv2 = m_lhs_cof(2,i,j);
         double mucofw2 = m_lhs_cof(3,i,j), mucofuv = m_lhs_cof(4,i,j);
         double mucofuw = m_lhs_cof(5,i,j), mucofvw = m_lhs_cof(6,i,j);
         double ijac    = m_lhs_cof(7,i,j);
         lhsd(1,i,j) = (mucofu2*xd(1,i,j) + mucofuv*xd(2,i,j) + mucofuw*xd(3,i,j))*ijac;
	 lhsd(2,i,j) = (mucofuv*xd(1,i,j) + mucofv2*xd(2,i,j) + mucofvw*xd(3,i,j))*ijac;
         lhsd(3,i,j) = (mucofuw*xd(1,i,j) + mucofvw*xd(2,i,j) + mucofw2*xd(3,i,j))*ijac;
      }
#undef strx_d
#undef stry_d
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::prolongate2D_d( DPlane& Uc, DPlane& Uf )
{
   const double i16 = 1.0/16;
   const double i256 = 1.0/256;
   int ib1, ie1, ib2, ie2;
   if( Uf.ib % 2 == 0 )
      ib1 = Uf.ib/2+1;
   else
      ib1 = (Uf.ib+1)/2;
   ib1 = max(Uc.ib,ib1);
   if( Uf.ie % 2 == 0 )
      ie1 = Uf.ie/2;
   else
      ie1 = (Uf.ie+1)/2;
   ie1 = min(Uc.ie,ie1);

   if( Uf.ib % 2 == 0 )
      ib2 = Uf.ib/2;
   else
      ib2 = (Uf.ib+1)/2;
   ib2 = max(Uc.ib+1,ib2);
   if( Uf.ie % 2 == 0 )
      ie2 = Uf.ie/2;
   else
      ie2 = (Uf.ie-1)/2;
   ie2 = min(Uc.ie-2,ie2);

   int jb1, je1, jb2, je2;
   if( Uf.jb % 2 == 0 )
      jb1 = Uf.jb/2+1;
   else
      jb1 = (Uf.jb+1)/2;
   jb1 = max(Uc.jb,jb1);
   if( Uf.je % 2 == 0 )
      je1 = Uf.je/2;
   else
      je1 = (Uf.je+1)/2;
   je1 = min(Uc.je,je1);

   if( Uf.jb % 2 == 0 )
      jb2 = Uf.jb/2;
   else
      jb2 = (Uf.jb+1)/2;
   jb2 = max(Uc.jb+1,jb2);
   if( Uf.je % 2 == 0 )
      je2 = Uf.je/2;
   else
      je2 = (Uf.je-1)/2;
   je2 = min(Uc.je-2,je2);

#pragma omp parallel
   {
   for( int c=1 ; c <= Uf.nc ;c++)
#pragma omp for
      for( int j=jb1 ; j <= je1 ; j++ )
#pragma omp simd
         for( int i=ib1 ; i <= ie1 ; i++ )
            Uf(c,2*i-1,2*j-1) = Uc(c,i,j);
   for( int c=1 ; c <= Uf.nc ;c++)
#pragma omp for
      for( int j=jb2 ; j <= je2 ; j++ )
#pragma omp simd
         for( int i=ib1 ; i <= ie1 ; i++ )
            Uf(c,2*i-1,2*j  ) = i16*(-Uc(c,i,j-1)+9*(Uc(c,i,j)+Uc(c,i,j+1))-Uc(c,i,j+2));
   for( int c=1 ; c <= Uf.nc ;c++)
#pragma omp for
      for( int j=jb1 ; j <= je1 ; j++ )
#pragma omp simd
         for( int i=ib2 ; i <= ie2 ; i++ )
            Uf(c,2*i,  2*j-1) = i16*(-Uc(c,i-1,j)+9*(Uc(c,i,j)+Uc(c,i+1,j))-Uc(c,i+2,j));
   for( int c=1 ; c <= Uf.nc ;c++)
#pragma omp for
      for( int j=jb2 ; j <= je2 ; j++ )
#pragma omp simd
         for( int i=ib2 ; i <= ie2 ; i++ )
            Uf(c,2*i,  2*j  ) = i256*
               ( Uc(c,i-1,j-1)-9*(Uc(c,i,j-1)+Uc(c,i+1,j-1))+Uc(c,i+2,j-1)
           + 9*(-Uc(c,i-1,j)+9*(Uc(c,i,j)+Uc(c,i+1,j))-Uc(c,i+2,j)
                -Uc(c,i-1,j+1)+9*(Uc(c,i,j+1)+Uc(c,i+1,j+1))-Uc(c,i+2,j+1))
                +Uc(c,i-1,j+2)-9*(Uc(c,i,j+2)+Uc(c,i+1,j+2))+Uc(c,i+2,j+2));
}
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::restrict2D_d( DPlane& Uc, DPlane& Uf )
{
   int icb, ice, jcb, jce;
   if( Uf.ib % 2 == 0 )
      icb = Uf.ib/2+2;
   else
      icb = (Uf.ib+1)/2+2;
   icb = max(Uc.ib,icb);
   if( Uf.ie % 2 == 0 )
      ice = Uf.ie/2-1;
   else
      ice = (Uf.ie-1)/2-1;
   ice = min(Uc.ie,ice);

   if( Uf.jb % 2 == 0 )
      jcb = Uf.jb/2+2;
   else
      jcb = (Uf.jb+1)/2+2;
   jcb = max(Uc.jb,jcb);
   if( Uf.je % 2 == 0 )
      jce = Uf.je/2-1;
   else
      jce = (Uf.je-1)/2-1;
   jce = min(Uc.je,jce);

   const double i1024 = 4.0/1024; // Multiply r:=4*r
#pragma omp parallel
   for (int c=1; c <= Uf.nc; c++)
#pragma omp for
      for( int jc= jcb ; jc <= jce ; jc++ )
#pragma omp simd
         for( int ic= icb ; ic <= ice ; ic++ )
         {
            int i=2*ic-1, j=2*jc-1;
            Uc(c,ic,jc)  = i1024*(
                    Uf(c,i-3,j-3)-9*Uf(c,i-3,j-1)-16*Uf(c,i-3,j)-9*Uf(c,i-3,j+1)+Uf(c,i-3,j+3)
               +9*(-Uf(c,i-1,j-3)+9*Uf(c,i-1,j-1)+16*Uf(c,i-1,j)+9*Uf(c,i-1,j+1)-Uf(c,i-1,j+3))
              +16*(-Uf(c,i,  j-3)+9*Uf(c,i,  j-1)+16*Uf(c,i,  j)+9*Uf(c,i,  j+1)-Uf(c,i,  j+3))
               +9*(-Uf(c,i+1,j-3)+9*Uf(c,i+1,j-1)+16*Uf(c,i+1,j)+9*Uf(c,i+1,j+1)-Uf(c,i+1,j+3)) +
                    Uf(c,i+3,j-3)-9*Uf(c,i+3,j-1)-16*Uf(c,i+3,j)-9*Uf(c,i+3,j+1)+Uf(c,i+3,j+3) );
         }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::bnd_zero_d( DPlane& u, int npts )
{
// Homogeneous Dirichlet at boundaries on sides (2D-plane analog of bnd_zero).
   for( int s=0 ; s < 4 ; s++ )
      if( m_isbndry[s] )
      {
         int jb=u.jb, je=u.je, ib=u.ib, ie=u.ie;
         if( s == 0 )
            ie = ib+npts-1;
         if( s == 1 )
            ib = ie-npts+1;
         if( s == 2 )
            je = jb+npts-1;
         if( s == 3 )
            jb = je-npts+1;
         for(int c=1 ; c <= u.nc ; c++)
            for( int j=jb ; j <= je ; j++ )
               for( int i=ib ; i <= ie ; i++ )
                  u(c,i,j)=0;
      }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_plane_d( DPlane& u )
{
//
// General ghost point exchange at processor boundaries, double-precision
// single-plane analog of communicate_array (which operates on float_sw4
// Sarrays with an arbitrary k-range).
//
  const int ng = m_nghost;
  const int ni = (u.ie-u.ib+1);
  const int nj = (u.je-u.jb+1);
  double *sbuf1, *sbuf2, *rbuf1, *rbuf2;

  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1=503, tag2=504;

  size_t npts1 = (size_t)ng*nj;
  size_t npts2 = (size_t)ni*ng;
  size_t nptsmax = max(npts1,npts2);
  double* tmp = new double[4*nptsmax*u.nc];
  sbuf1 = &tmp[0];
  rbuf1 = &tmp[  nptsmax*u.nc];
  sbuf2 = &tmp[2*nptsmax*u.nc];
  rbuf2 = &tmp[3*nptsmax*u.nc];

// i-direction communication
  MPI_Irecv( rbuf1, npts1*u.nc, MPI_DOUBLE, m_ew->m_neighbor[1], tag1,
	     m_ew->m_cartesian_communicator, &req1 );
  MPI_Irecv( rbuf2, npts1*u.nc, MPI_DOUBLE, m_ew->m_neighbor[0], tag2,
	     m_ew->m_cartesian_communicator, &req2 );
  if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb ; j <= u.je; j++ )
        for(int i=u.ib+ng ; i <= u.ib+2*ng-1; i++ )
	{
	   size_t ind = i-(u.ib+ng)+ng*(j-u.jb);
	   sbuf1[ind+npts1*(c-1)]= u(c,i,j);
        }
  MPI_Isend( sbuf1, npts1*u.nc, MPI_DOUBLE, m_ew->m_neighbor[0], tag1,
	     m_ew->m_cartesian_communicator, &req3 );
  if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb ; j <= u.je; j++ )
        for(int i=u.ie-2*ng+1 ; i <= u.ie-ng; i++ )
	{
	   size_t ind = i-(u.ie-2*ng+1)+ng*(j-u.jb);
	   sbuf2[ind+npts1*(c-1)]= u(c,i,j);
        }
  MPI_Isend( sbuf2, npts1*u.nc, MPI_DOUBLE, m_ew->m_neighbor[1], tag2,
	     m_ew->m_cartesian_communicator, &req4);
  MPI_Wait( &req1, &status );
  if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb ; j <= u.je; j++ )
        for(int i=u.ie-ng+1 ; i <= u.ie; i++ )
	{
	   size_t ind = i-(u.ie-ng+1)+ng*(j-u.jb);
	   u(c,i,j) = rbuf1[ind+npts1*(c-1)];
        }
  MPI_Wait( &req2, &status );
  if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb ; j <= u.je; j++ )
        for(int i=u.ib ; i <= u.ib+ng-1; i++ )
	{
	   size_t ind = i-u.ib+ng*(j-u.jb);
	   u(c,i,j) = rbuf2[ind+npts1*(c-1)];
        }

  MPI_Wait( &req3, &status );
  MPI_Wait( &req4, &status );

// j-direction communication
  MPI_Irecv( rbuf1, npts2*u.nc, MPI_DOUBLE, m_ew->m_neighbor[3], tag1,
	     m_ew->m_cartesian_communicator, &req1 );
  MPI_Irecv( rbuf2, npts2*u.nc, MPI_DOUBLE, m_ew->m_neighbor[2], tag2,
	     m_ew->m_cartesian_communicator, &req2 );
  if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb+ng ; j <= u.jb+2*ng-1; j++ )
        for(int i=u.ib ; i <= u.ie; i++ )
	{
	   size_t ind = i-u.ib+ni*(j-(u.jb+ng));
	   sbuf1[ind+npts2*(c-1)]= u(c,i,j);
        }
  MPI_Isend( sbuf1, npts2*u.nc, MPI_DOUBLE, m_ew->m_neighbor[2], tag1,
	     m_ew->m_cartesian_communicator, &req3 );
  if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.je-2*ng+1 ; j <= u.je-ng; j++ )
        for(int i=u.ib ; i <= u.ie; i++ )
	{
	   size_t ind = i-u.ib+ni*(j-(u.je-2*ng+1));
	   sbuf2[ind+npts2*(c-1)]= u(c,i,j);
        }
  MPI_Isend( sbuf2, npts2*u.nc, MPI_DOUBLE, m_ew->m_neighbor[3], tag2,
	     m_ew->m_cartesian_communicator, &req4);
  MPI_Wait( &req1, &status );
  if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.je-ng+1 ; j <= u.je; j++ )
        for(int i=u.ib ; i <= u.ie; i++ )
	{
	   size_t ind = i-u.ib + ni*(j-(u.je-ng+1));
	   u(c,i,j) = rbuf1[ind+npts2*(c-1)];
        }
  MPI_Wait( &req2, &status );
  if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.nc ; c++ )
     for( int j=u.jb ; j <= u.jb+ng-1; j++ )
        for(int i=u.ib ; i <= u.ie; i++ )
	{
	   size_t ind = i-u.ib+ni*(j-u.jb);
	   u(c,i,j) = rbuf2[ind+npts2*(c-1)];
        }

  MPI_Wait( &req3, &status );
  MPI_Wait( &req4, &status );
  delete[] tmp;
}

//-----------------------------------------------------------------------
 void CurvilinearInterface2::mat_icstresses_curv( int ib, int jb, Sarray& a_mat, int kic,
                              Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
                              float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop )
{
   // As compute_icstresses_curv, but evaluates the matrix multiplying the ghost point part 
  //   const float_sw4 a1=2.0/3, a2=-1.0/12;
   const bool upper = (kic == 1);
   const int k=kic;
   // const int kl = upper ? 1 :-1;
#define str_x(i) a_str_x[(i-ib)]   
#define str_y(j) a_str_y[(j-jb)]   

   float_sw4 sb = sbop[0];
   if( !upper )
     sb = -sb;

#pragma omp parallel for
   for( int j=a_mat.m_jb ; j <= a_mat.m_je ; j++ )
#pragma omp simd
      for( int i=a_mat.m_ib ; i <= a_mat.m_ie ; i++ )
      {
         float cof = sb/(str_x(i)*str_y(j));
         // Normal terms
         float_sw4 m2 = str_x(i)*a_metric(2,i,j,k);
         float_sw4 m3 = str_y(j)*a_metric(3,i,j,k);
         float_sw4 m4 = a_metric(4,i,j,k);
         //         float_sw4 un   = m2*uz + m3*vz + m4*wz;
         float_sw4 mnrm = m2*m2 + m3*m3 + m4*m4;

         //         a_lhs(1,i,j,k) = a_mu(i,j,k)*mnrm*uz + (a_mu(i,j,k)+a_lambda(i,j,k))*m2*un;
         a_mat(1,i,j,k) -= cof*(a_mu(i,j,k)*mnrm + (a_mu(i,j,k)+a_lambda(i,j,k))*m2*m2); // dB1/du1
         a_mat(4,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m2*m3; // dB1/du2
         a_mat(7,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m2*m4; // dB1/du3

         //         a_lhs(2,i,j,k) = a_mu(i,j,k)*mnrm*vz + (a_mu(i,j,k)+a_lambda(i,j,k))*m3*un;
         a_mat(2,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m3*m2; // dB2/du1
         a_mat(5,i,j,k) -= cof*(a_mu(i,j,k)*mnrm + (a_mu(i,j,k)+a_lambda(i,j,k))*m3*m3); // dB2/du2
         a_mat(8,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m3*m4; // dB2/du3

         //         a_lhs(3,i,j,k) = a_mu(i,j,k)*mnrm*wz + (a_mu(i,j,k)+a_lambda(i,j,k))*m4*un;
         a_mat(3,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m4*m2; // dB3/du1
         a_mat(6,i,j,k) -=                    cof*(a_mu(i,j,k)+a_lambda(i,j,k))*m4*m3; // dB3/du2
         a_mat(9,i,j,k) -= cof*(a_mu(i,j,k)*mnrm + (a_mu(i,j,k)+a_lambda(i,j,k))*m4*m4); // dB3/du3


         //         float_sw4 isgxy = 1.0/(str_x(i)*str_y(j));
         //         a_lhs(1,i,j,k) *= isgxy;
         //         a_lhs(2,i,j,k) *= isgxy;
         //         a_lhs(3,i,j,k) *= isgxy;
      }
#undef str_x
#undef str_y

}

//-----------------------------------------------------------------------
void CurvilinearInterface2::matrix_Lu( int strib, int strjb, Sarray& a_mat,
				       Sarray& met, Sarray& jac, 
				       Sarray& mu, Sarray& la, 
				       float_sw4* a_str_x, float_sw4* a_str_y,
				       float_sw4 ghcof )
{
#define strx(i) a_str_x[(i-strib)]   
#define stry(j) a_str_y[(j-strjb)]   
#pragma omp parallel for
   for( int j=a_mat.m_jb; j <= a_mat.m_je ;j++ )
#pragma omp simd
      for( int i=a_mat.m_ib; i <= a_mat.m_ie ;i++ )
      {
         float_sw4 ijac = ghcof/jac(i,j,1);
         float_sw4 mucofu2 = ((2*mu(i,j,1)+la(i,j,1))*
				   met(2,i,j,1)*strx(i)*met(2,i,j,1)*strx(i)
				   + mu(i,j,1)*(met(3,i,j,1)*stry(j)*met(3,i,j,1)*stry(j)+
						met(4,i,j,1)*met(4,i,j,1) ));
	 float_sw4 mucofv2 = ((2*mu(i,j,1)+la(i,j,1))*
                                met(3,i,j,1)*stry(j)*met(3,i,j,1)*stry(j)
				   + mu(i,j,1)*( met(2,i,j,1)*strx(i)*met(2,i,j,1)*strx(i)+
						 met(4,i,j,1)*met(4,i,j,1) ) );
	 float_sw4 mucofw2 = ((2*mu(i,j,1)+la(i,j,1))*met(4,i,j,1)*met(4,i,j,1)
                                  + mu(i,j,1)*
				( met(2,i,j,1)*strx(i)*met(2,i,j,1)*strx(i)+
				  met(3,i,j,1)*stry(j)*met(3,i,j,1)*stry(j) ) );
	 float_sw4 mucofuv = (mu(i,j,1)+la(i,j,1))*met(2,i,j,1)*met(3,i,j,1)*strx(i)*stry(j);
	 float_sw4 mucofuw = (mu(i,j,1)+la(i,j,1))*met(2,i,j,1)*met(4,i,j,1)*strx(i);
	 float_sw4 mucofvw = (mu(i,j,1)+la(i,j,1))*met(3,i,j,1)*met(4,i,j,1)*stry(j);
         a_mat(1,i,j,1) = mucofu2*ijac;
         a_mat(4,i,j,1) = mucofuv*ijac;
         a_mat(7,i,j,1) = mucofuw*ijac;

         a_mat(2,i,j,1) = mucofuv*ijac;
         a_mat(5,i,j,1) = mucofv2*ijac;
         a_mat(8,i,j,1) = mucofvw*ijac;

         a_mat(3,i,j,1) = mucofuw*ijac;
         a_mat(6,i,j,1) = mucofvw*ijac;
         a_mat(9,i,j,1) = mucofw2*ijac;
      }
#undef strx
#undef stry
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::restprol2D( Sarray& Uc, Sarray& alpha, int kc, int kf )
{
  //
  // Multiplies the diagonal element of the operator P^T*diag(alpha)*P 
  //
  const float_sw4 a=(9.0/16);
  const float_sw4 b=-1.0/16;
  const float_sw4 a2 = a*a;
  const float_sw4 b2 = b*b;
  const float_sw4 a4 = a2*a2;
  const float_sw4 b4 = b2*b2;
  const float_sw4 a2b2 = a2*b2;
  for( int c=1 ; c <= Uc.m_nc; c++)
    //   for( int j=Uc.m_jb+5; j <= Uc.m_je-5 ;j++ )
    //      for( int i=Uc.m_ib+5; i <= Uc.m_ie-5 ;i++ )
    // Note, this routine is called with Uc = matrix block, which is
    // declared without ghost points, so need to loop over the
    // full i- and j- ranges.
#pragma omp parallel for
   for( int j=Uc.m_jb; j <= Uc.m_je ;j++ )
#pragma omp simd
      for( int i=Uc.m_ib; i <= Uc.m_ie ;i++ )
	{
	  Uc(c,i,j,kc) =( alpha(2*i-1,2*j-1,kf)+
                        a2*(alpha(2*i-1,2*j,  kf) + alpha(2*i-1,2*j-2,kf)+
			    alpha(2*i,  2*j-1,kf) + alpha(2*i-2,2*j-1,kf)  ) +
			b2*(alpha(2*i+2,2*j-1,kf) + alpha(2*i-4,2*j-1,kf)+ 
			    alpha(2*i-1,2*j+2,kf) + alpha(2*i-1,2*j-4,kf)  )+
		        a4*(alpha(2*i  ,2*j,  kf) + alpha(2*i  ,2*j-2,kf)+
			    alpha(2*i-2,2*j  ,kf) + alpha(2*i-2,2*j-2,kf)  )+
			b4*(alpha(2*i+2,2*j+2,kf) + alpha(2*i+2,2*j-4,kf)+ 
			    alpha(2*i-4,2*j+2,kf) + alpha(2*i-4,2*j-4,kf)  )+
                     a2b2*( alpha(2*i,  2*j+2,kf) + alpha(2*i+2,2*j,  kf)+
			    alpha(2*i+2,2*j-2,kf) + alpha(2*i-2,2*j+2,kf)+
			    alpha(2*i-2,2*j-4,kf) + alpha(2*i-4,2*j-2,kf)+
			    alpha(2*i,  2*j-4,kf) + alpha(2*i-4,2*j  ,kf) ) )*Uc(c,i,j,kc);
	}
}
//-----------------------------------------------------------------------
void CurvilinearInterface2::prolongate2D( Sarray& Uc, Sarray& Uf, int kc, int kf )
{
   const float_sw4 i16 = 1.0/16;
   const float_sw4 i256 = 1.0/256;
   int ib1, ie1, ib2, ie2;
   if( Uf.m_ib % 2 == 0 )
      ib1 = Uf.m_ib/2+1;
   else
      ib1 = (Uf.m_ib+1)/2;
   ib1 = max(Uc.m_ib,ib1);
   if( Uf.m_ie % 2 == 0 )
      ie1 = Uf.m_ie/2;
   else
      ie1 = (Uf.m_ie+1)/2;
   ie1 = min(Uc.m_ie,ie1);

   if( Uf.m_ib % 2 == 0 )
      ib2 = Uf.m_ib/2;
   else
      ib2 = (Uf.m_ib+1)/2;
   ib2 = max(Uc.m_ib+1,ib2);
   if( Uf.m_ie % 2 == 0 )
      ie2 = Uf.m_ie/2;
   else
      ie2 = (Uf.m_ie-1)/2;
   ie2 = min(Uc.m_ie-2,ie2);

   int jb1, je1, jb2, je2;
   if( Uf.m_jb % 2 == 0 )
      jb1 = Uf.m_jb/2+1;
   else
      jb1 = (Uf.m_jb+1)/2;
   jb1 = max(Uc.m_jb,jb1);
   if( Uf.m_je % 2 == 0 )
      je1 = Uf.m_je/2;
   else
      je1 = (Uf.m_je+1)/2;
   je1 = min(Uc.m_je,je1);

   if( Uf.m_jb % 2 == 0 )
      jb2 = Uf.m_jb/2;
   else
      jb2 = (Uf.m_jb+1)/2;
   jb2 = max(Uc.m_jb+1,jb2);
   if( Uf.m_je % 2 == 0 )
      je2 = Uf.m_je/2;
   else
      je2 = (Uf.m_je-1)/2;
   je2 = min(Uc.m_je-2,je2);

#pragma omp parallel
   {
   for( int c=1 ; c <= Uf.m_nc ;c++)
#pragma omp for
      for( int j=jb1 ; j <= je1 ; j++ )
#pragma omp simd
         for( int i=ib1 ; i <= ie1 ; i++ )
            Uf(c,2*i-1,2*j-1,kf) = Uc(c,i,j,kc);
   for( int c=1 ; c <= Uf.m_nc ;c++)
#pragma omp for
      for( int j=jb2 ; j <= je2 ; j++ )
#pragma omp simd
         for( int i=ib1 ; i <= ie1 ; i++ )
            Uf(c,2*i-1,2*j,  kf) = i16*(-Uc(c,i,j-1,kc)+9*(Uc(c,i,j,kc)+Uc(c,i,j+1,kc))-Uc(c,i,j+2,kc));
   for( int c=1 ; c <= Uf.m_nc ;c++)
#pragma omp for
      for( int j=jb1 ; j <= je1 ; j++ )
#pragma omp simd
         for( int i=ib2 ; i <= ie2 ; i++ )
            Uf(c,2*i,  2*j-1,kf) = i16*(-Uc(c,i-1,j,kc)+9*(Uc(c,i,j,kc)+Uc(c,i+1,j,kc))-Uc(c,i+2,j,kc));
   for( int c=1 ; c <= Uf.m_nc ;c++)
#pragma omp for
      for( int j=jb2 ; j <= je2 ; j++ )
#pragma omp simd
         for( int i=ib2 ; i <= ie2 ; i++ )
            Uf(c,2*i,  2*j,  kf) = i256*
               ( Uc(c,i-1,j-1,kc)-9*(Uc(c,i,j-1,kc)+Uc(c,i+1,j-1,kc))+Uc(c,i+2,j-1,kc)
           + 9*(-Uc(c,i-1,j,  kc)+9*(Uc(c,i,j,  kc)+Uc(c,i+1,j,  kc))-Uc(c,i+2,j,  kc)  
                -Uc(c,i-1,j+1,kc)+9*(Uc(c,i,j+1,kc)+Uc(c,i+1,j+1,kc))-Uc(c,i+2,j+1,kc))
                +Uc(c,i-1,j+2,kc)-9*(Uc(c,i,j+2,kc)+Uc(c,i+1,j+2,kc))+Uc(c,i+2,j+2,kc));
}
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::restrict2D( Sarray& Uc, Sarray& Uf, int kc, int kf )
{
   int icb, ice, jcb, jce;
   if( Uf.m_ib % 2 == 0 )
      icb = Uf.m_ib/2+2;
   else
      icb = (Uf.m_ib+1)/2+2;
   icb = max(Uc.m_ib,icb);
   if( Uf.m_ie % 2 == 0 )
      ice = Uf.m_ie/2-1;
   else
      ice = (Uf.m_ie-1)/2-1;
   ice = min(Uc.m_ie,ice);

   if( Uf.m_jb % 2 == 0 )
      jcb = Uf.m_jb/2+2;
   else
      jcb = (Uf.m_jb+1)/2+2;
   jcb = max(Uc.m_jb,jcb);
   if( Uf.m_je % 2 == 0 )
      jce = Uf.m_je/2-1;
   else
      jce = (Uf.m_je-1)/2-1;
   jce = min(Uc.m_je,jce);


   const float_sw4 i1024 =4.0/1024; // Multiply r:=4*r 
#pragma omp parallel
   for (int c=1; c <= Uf.m_nc; c++)
#pragma omp for
      for( int jc= jcb ; jc <= jce ; jc++ )
#pragma omp simd
         for( int ic= icb ; ic <= ice ; ic++ )
         {
            int i=2*ic-1, j=2*jc-1;
            Uc(c,ic,jc,kc)  = i1024*( 
                    Uf(c,i-3,j-3,kf)-9*Uf(c,i-3,j-1,kf)-16*Uf(c,i-3,j,kf)-9*Uf(c,i-3,j+1,kf)+Uf(c,i-3,j+3,kf)
               +9*(-Uf(c,i-1,j-3,kf)+9*Uf(c,i-1,j-1,kf)+16*Uf(c,i-1,j,kf)+9*Uf(c,i-1,j+1,kf)-Uf(c,i-1,j+3,kf))
              +16*(-Uf(c,i,  j-3,kf)+9*Uf(c,i,  j-1,kf)+16*Uf(c,i,  j,kf)+9*Uf(c,i,  j+1,kf)-Uf(c,i,  j+3,kf))
               +9*(-Uf(c,i+1,j-3,kf)+9*Uf(c,i+1,j-1,kf)+16*Uf(c,i+1,j,kf)+9*Uf(c,i+1,j+1,kf)-Uf(c,i+1,j+3,kf)) +
                    Uf(c,i+3,j-3,kf)-9*Uf(c,i+3,j-1,kf)-16*Uf(c,i+3,j,kf)-9*Uf(c,i+3,j+1,kf)+Uf(c,i+3,j+3,kf) );
         }
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_array( Sarray& u, bool allkplanes, int kplane )
{
// 
// General ghost point exchange at processor boundaries.
//
// Excplicit copy to buffers, not using fancy MPI-datatypes or sendrecv.
//
  int kb = u.m_kb;
  int ke = u.m_ke;
  if( !allkplanes )
    ke = kb = kplane;
  const int ng = m_nghost;
  const int ni = (u.m_ie-u.m_ib+1);
  const int nj = (u.m_je-u.m_jb+1);
  const int nk = ke-kb+1;
  float_sw4 *sbuf1, *sbuf2, *rbuf1, *rbuf2;

  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1=203, tag2=204;

  size_t npts1 = ng*nj*nk;
  size_t npts2 = ni*ng*nk;
  size_t nptsmax = max(npts1,npts2);
  float_sw4* tmp = new float_sw4[4*nptsmax*u.m_nc];
  sbuf1 = &tmp[0];
  rbuf1 = &tmp[  nptsmax*u.m_nc];
  sbuf2 = &tmp[2*nptsmax*u.m_nc];
  rbuf2 = &tmp[3*nptsmax*u.m_nc];

// i-direction communication  
  MPI_Irecv( rbuf1, npts1*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag1,
	     m_ew->m_cartesian_communicator, &req1 );
  MPI_Irecv( rbuf2, npts1*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag2,
	     m_ew->m_cartesian_communicator, &req2 );
  if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb ; j <= u.m_je; j++ )
           for(int i=u.m_ib+ng ; i <= u.m_ib+2*ng-1; i++ )
	   {
	      size_t ind = i-(u.m_ib+ng)+ng*(j-u.m_jb)+ng*nj*(k-kb);
	      sbuf1[ind+npts1*(c-1)]= u(c,i,j,k);
           }
  MPI_Isend( sbuf1, npts1*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[0], tag1,
	     m_ew->m_cartesian_communicator, &req3 );
  if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb ; j <= u.m_je; j++ )
           for(int i=u.m_ie-2*ng+1 ; i <= u.m_ie-ng; i++ )
	   {
	      size_t ind = i-(u.m_ie-2*ng+1)+ng*(j-u.m_jb)+ng*nj*(k-kb);
	      sbuf2[ind+npts1*(c-1)]= u(c,i,j,k);
           }
  MPI_Isend( sbuf2, npts1*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[1], tag2,
	     m_ew->m_cartesian_communicator, &req4);
  MPI_Wait( &req1, &status );
  if( m_ew->m_neighbor[1] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb ; j <= u.m_je; j++ )
           for(int i=u.m_ie-ng+1 ; i <= u.m_ie; i++ )
	   {
	      size_t ind = i-(u.m_ie-ng+1)+ng*(j-u.m_jb)+ng*nj*(k-kb);
	      u(c,i,j,k) = rbuf1[ind+npts1*(c-1)];
           }  
  MPI_Wait( &req2, &status );
  if( m_ew->m_neighbor[0] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb ; j <= u.m_je; j++ )
           for(int i=u.m_ib ; i <= u.m_ib+ng-1; i++ )
	   {
	      size_t ind = i-u.m_ib+ng*(j-u.m_jb)+ng*nj*(k-kb);
	      u(c,i,j,k) = rbuf2[ind+npts1*(c-1)];
           }

  MPI_Wait( &req3, &status );
  MPI_Wait( &req4, &status );

// j-direction communication  
  MPI_Irecv( rbuf1, npts2*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag1,
	     m_ew->m_cartesian_communicator, &req1 );
  MPI_Irecv( rbuf2, npts2*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag2,
	     m_ew->m_cartesian_communicator, &req2 );
  if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb+ng ; j <= u.m_jb+2*ng-1; j++ )
           for(int i=u.m_ib ; i <= u.m_ie; i++ )
	   {
	      size_t ind = i-u.m_ib+ni*(j-(u.m_jb+ng))+ng*ni*(k-kb);
	      sbuf1[ind+npts2*(c-1)]= u(c,i,j,k);
           }
  MPI_Isend( sbuf1, npts2*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[2], tag1,
	     m_ew->m_cartesian_communicator, &req3 );
  if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_je-2*ng+1 ; j <= u.m_je-ng; j++ )
           for(int i=u.m_ib ; i <= u.m_ie; i++ )
	   {
	      size_t ind = i-u.m_ib+ni*(j-(u.m_je-2*ng+1))+ng*ni*(k-kb);
	      sbuf2[ind+npts2*(c-1)]= u(c,i,j,k);
           }
  MPI_Isend( sbuf2, npts2*u.m_nc, m_ew->m_mpifloat, m_ew->m_neighbor[3], tag2,
	     m_ew->m_cartesian_communicator, &req4);
  MPI_Wait( &req1, &status );
  if( m_ew->m_neighbor[3] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_je-ng+1 ; j <= u.m_je; j++ )
           for(int i=u.m_ib ; i <= u.m_ie; i++ )
	   {
	      size_t ind = i-u.m_ib + ni*(j-(u.m_je-ng+1))+ng*ni*(k-kb);
	      u(c,i,j,k) = rbuf1[ind+npts2*(c-1)];
           }  
  MPI_Wait( &req2, &status );
  if( m_ew->m_neighbor[2] != MPI_PROC_NULL )
  for( int c=1 ; c <= u.m_nc ; c++ )
     for( int k=kb ; k <= ke; k++ )
        for( int j=u.m_jb ; j <= u.m_jb+ng-1; j++ )
           for(int i=u.m_ib ; i <= u.m_ie; i++ )
	   {
	      size_t ind = i-u.m_ib+ni*(j-u.m_jb)+ng*ni*(k-kb);
	      u(c,i,j,k) = rbuf2[ind+npts2*(c-1)];
           }
    
  MPI_Wait( &req3, &status );
  MPI_Wait( &req4, &status );
  delete[] tmp;
}

//-----------------------------------------------------------------------
void CurvilinearInterface2::communicate_array1d( float_sw4* u, int n, int dir, int ngh )
{
  // 
  // Communicate one dimensional array in i- or j-direction
  // Input:  u - The array
  //         n - Number of elements in array
  //         dir- Direction, 0 is i, 1 is j
  //         ngh- Number of points to communicate (ghost points)
  //
  // The first and last ngh points of u are updated.
  //
  MPI_Request req1, req2, req3, req4;
  MPI_Status status;
  int tag1=302, tag2=303;
  int no=2*dir;

  MPI_Irecv( &u[n-ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[1+no], tag1,
	     m_ew->m_cartesian_communicator, &req1 );

  MPI_Irecv( &u[0], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[no], tag2,
	     m_ew->m_cartesian_communicator, &req2 );

  MPI_Isend( &u[ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[no], tag1,
	     m_ew->m_cartesian_communicator, &req3 );

  MPI_Isend( &u[n-2*ngh], ngh, m_ew->m_mpifloat, m_ew->m_neighbor[1+no], tag2,
	     m_ew->m_cartesian_communicator, &req4);

  MPI_Wait( &req1, &status );
  MPI_Wait( &req2, &status );
}
