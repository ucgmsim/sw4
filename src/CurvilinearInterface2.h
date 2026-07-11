#ifndef SW4_CURVILINEARINTERFACE2
#define SW4_CURVILINEARINTERFACE2

#include <vector>

class Sarray;
class EW;

class CurvilinearInterface2
{
   // Fixed-precision 2D (c,i,j) plane buffer used to solve the interface
   // ghost-point equation in double precision even when float_sw4==float,
   // since the block-Jacobi residual is otherwise unattainable at the
   // requested tolerances in single precision (residual floor = 1 ULP).
   struct DPlane
   {
      int ib, ie, jb, je, nc;
      std::vector<double> data;
      DPlane():ib(0),ie(-1),jb(0),je(-1),nc(0) {}
      DPlane( int nc_in, int ib_in, int ie_in, int jb_in, int je_in ) :
         ib(ib_in), ie(ie_in), jb(jb_in), je(je_in), nc(nc_in),
         data( (size_t)nc_in*(ie_in-ib_in+1)*(je_in-jb_in+1), 0.0 )
      {}
      inline double& operator()( int c, int i, int j )
      { return data[ (c-1) + nc*((size_t)(i-ib) + (size_t)(ie-ib+1)*(j-jb)) ]; }
   };

   EW* m_ew;
   TestTwilight* m_tw;
   TestEcons* m_etest;
   TestPointSource* m_psource;

   int m_nghost, m_ib, m_ie, m_jb, m_je, m_ibf, m_ief, m_jbf, m_jef, m_nkf;
   int m_kb, m_ke, m_kbf, m_kef;
   int m_gc, m_gf;
   bool m_isbndry[4]; // side is physical boundary

   float_sw4 m_reltol, m_abstol;
   int m_maxit;

   float_sw4 *m_strx_c, *m_strx_f, *m_stry_c, *m_stry_f;
   Sarray m_Mass_block, m_rho_c, m_rho_f, m_mu_c, m_mu_f, m_lambda_c, m_lambda_f;
   Sarray m_met_c, m_met_f, m_jac_c, m_jac_f, m_x_c, m_x_f, m_y_c, m_y_f, m_z_c, m_z_f;
   bool m_use_attenuation;
   int m_number_mechanisms;
   std::vector<Sarray> m_muve_f, m_lambdave_f, m_muve_c, m_lambdave_c;

   float_sw4* m_mass_block;
   int* m_ipiv_block;
   bool m_memory_is_allocated;

   float_sw4 m_acof[384], m_bope[48], m_ghcof[6], m_acof_no_gp[384], m_ghcof_no_gp[6];
   float_sw4 m_sbop[6], m_sbop_no_gp[6], m_bop[24];

   void interface_block( Sarray& matrix );
   void compute_icstresses_curv( Sarray& a_Up, Sarray& B, int kic,
				 Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
				 float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop, char op );

   void mat_icstresses_curv( int ib, int jb, Sarray& a_mat, int kic,
			     Sarray& a_metric, Sarray& a_mu, Sarray& a_lambda,
			     float_sw4* a_str_x, float_sw4* a_str_y, float_sw4* sbop );

   void matrix_Lu( int ib, int jb, Sarray& a_mat, Sarray& met, Sarray& jac, 
		   Sarray& mu, Sarray& la, 
		   float_sw4* a_str_x, float_sw4* a_str_y, float_sw4 ghcof );

   void injection(Sarray &u_f, Sarray &u_c );
   void restrict2D( Sarray& Uc, Sarray& Uf, int kc, int kf );
   void restprol2D( Sarray& Uc, Sarray& alpha, int kc, int kf );
   void bnd_zero( Sarray& u, int npts );
   void copy_str( float_sw4* dest, float_sw4* src, int offset, int n, int nsw );
   void communicate_array1d( float_sw4* u, int n, int dir, int ngh );
   void communicate_array( Sarray& u, bool allkplanes=true, int kplane=0 );
   void init_arrays_att();

   // Double-precision variants of the interface ghost-point solve kernels
   // (see DPlane above). The block-Jacobi solve for U_c's k=0 ghost plane
   // (interface_lhs_d) and the fixed forcing term it solves against
   // (interface_rhs_d) both run in double; only the interior stencil
   // evaluations they call (curvilinear4sgwind, compute_icstresses_curv,
   // still float_sw4 -- see interface_rhs_d) touch data that is genuinely
   // stored at float_sw4 precision elsewhere in the simulation. Everything
   // else in the timestepping loop stays at float_sw4 precision.
   void interface_lhs_d( DPlane& lhsd, DPlane& xd );
   void interface_rhs_d( DPlane& rhsd, Sarray& uc, Sarray& uf, Sarray& fc, Sarray& ff,
                         std::vector<Sarray>& Alpha_c, std::vector<Sarray>& Alpha_f );
   void lhs_Lu_d( DPlane& xd, DPlane& lhsd );
   void lhs_icstresses_curv_d( DPlane& xd, DPlane& Bc );
   void prolongate2D_d( DPlane& Uc, DPlane& Uf );
   void restrict2D_d( DPlane& Uc, DPlane& Uf );
   void bnd_zero_d( DPlane& u, int npts );
   void communicate_plane_d( DPlane& u );
public:
   CurvilinearInterface2( int a_gc, EW* a_ew );
   CurvilinearInterface2() {}
   void init_arrays( std::vector<float_sw4*>& a_strx, std::vector<float_sw4*>& a_stry,
                     std::vector<Sarray>& a_rho, std::vector<Sarray>& a_mu,
                     std::vector<Sarray>& a_lambda );
   //   void test1( EW* a_ew, int gc, std::vector<Sarray>& a_U );
   //   void test2( EW* a_ew, int gc, std::vector<Sarray>& a_U );
   
   void impose_ic( std::vector<Sarray>& a_U, float_sw4 t,
                   std::vector<Sarray>& a_F, 
                   std::vector<Sarray*>& a_AlphaVE,
                   bool injection_only=false );

   void prolongate2D( Sarray& Uc, Sarray& Uf, int kc, int kf );
};

#endif
