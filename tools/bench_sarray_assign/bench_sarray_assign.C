// Standalone benchmark of the Fortran-order -> C-order copy in Sarray::assign(const double*, int)
// (src/Sarray.C), which EW::get_exact_lamb2 calls every time step in the Lamb test cases.
//
//   bench_sarray_assign [n ...]      cubes of n^3 points with 3 components (default: 32 64 128 192 256 305)
//
// The destination is float (a USE_DOUBLE=OFF build), first touched the way Sarray::set_to_zero
// does it; the source is double, allocated and zeroed by one thread, as in get_exact_lamb2.
// Every variant is checked against the original loop, element by element.
#include <omp.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

typedef float float_sw4;

struct Arr {
   int ni, nj, nk, nc;
   size_t npts;
   float_sw4* data;
};

// the loop in src/Sarray.C (class array in C order, input array in Fortran order)
static void assign_orig( Arr& a, const double* ar )
{
   const int m_ni = a.ni, m_nj = a.nj, m_nk = a.nk, m_nc = a.nc;
   float_sw4* m_data = a.data;
#pragma omp parallel for
   for( int i=0 ; i <m_ni ; i++ )
      for( int j=0 ; j <m_nj ; j++ )
         for( int k=0 ; k <m_nk ; k++ )
            for( int c=0 ; c < m_nc ; c++ )
               m_data[i+m_ni*j+m_ni*m_nj*k+m_ni*m_nj*m_nk*c] = ar[c+m_nc*i+m_nc*m_ni*j+m_nc*m_ni*m_nj*k];
}

// (k,j) split over the threads as in set_to_zero; i innermost: contiguous writes, reads stride nc
static void assign_kjci( Arr& a, const double* ar )
{
   const int m_ni = a.ni, m_nj = a.nj, m_nk = a.nk, m_nc = a.nc;
   const size_t offc = (size_t)m_ni*m_nj*m_nk;
   float_sw4* m_data = a.data;
#pragma omp parallel for collapse(2)
   for( int k=0 ; k <m_nk ; k++ )
      for( int j=0 ; j <m_nj ; j++ )
      {
         const size_t line = (size_t)m_ni*j + (size_t)m_ni*m_nj*k;
         for( int c=0 ; c < m_nc ; c++ )
            for( int i=0 ; i <m_ni ; i++ )
               m_data[line+offc*c+i] = ar[c+m_nc*(line+i)];
      }
}

// as kjci, but the i loop kept scalar (no gathers)
static void assign_kjci_novec( Arr& a, const double* ar )
{
   const int m_ni = a.ni, m_nj = a.nj, m_nk = a.nk, m_nc = a.nc;
   const size_t offc = (size_t)m_ni*m_nj*m_nk;
   float_sw4* m_data = a.data;
#pragma omp parallel for collapse(2)
   for( int k=0 ; k <m_nk ; k++ )
      for( int j=0 ; j <m_nj ; j++ )
      {
         const size_t line = (size_t)m_ni*j + (size_t)m_ni*m_nj*k;
         for( int c=0 ; c < m_nc ; c++ )
#pragma novector
            for( int i=0 ; i <m_ni ; i++ )
               m_data[line+offc*c+i] = ar[c+m_nc*(line+i)];
      }
}

// (k,j) split over the threads; c innermost: contiguous reads, nc write streams
static void assign_kjic( Arr& a, const double* ar )
{
   const int m_ni = a.ni, m_nj = a.nj, m_nk = a.nk, m_nc = a.nc;
   const size_t offc = (size_t)m_ni*m_nj*m_nk;
   float_sw4* m_data = a.data;
#pragma omp parallel for collapse(2)
   for( int k=0 ; k <m_nk ; k++ )
      for( int j=0 ; j <m_nj ; j++ )
      {
         const size_t line = (size_t)m_ni*j + (size_t)m_ni*m_nj*k;
         for( int i=0 ; i <m_ni ; i++ )
            for( int c=0 ; c < m_nc ; c++ )
               m_data[line+offc*c+i] = ar[c+m_nc*(line+i)];
      }
}

typedef void (*assign_fn)( Arr&, const double* );

int main( int argc, char** argv )
{
   std::vector<int> sizes;
   for( int a=1 ; a < argc ; a++ )
      sizes.push_back(atoi(argv[a]));
   if( sizes.empty() )
      sizes = {32, 64, 128, 192, 256, 305};

   const char* names[] = {"orig", "kjci", "kjci-novec", "kjic"};
   assign_fn fns[] = {assign_orig, assign_kjci, assign_kjci_novec, assign_kjic};
   const int nvar = 4;
   const char* rank = getenv("SLURM_PROCID");

   printf("rank %s, %d OpenMP threads\n", rank ? rank : "-", omp_get_max_threads());
   printf("%6s %10s %-11s %10s %10s %9s %s\n", "n", "MB(src)", "variant", "min (ms)", "mean (ms)", "GB/s", "check");
   for( int n : sizes )
   {
      Arr a;
      a.ni = a.nj = a.nk = n;
      a.nc = 3;
      a.npts = (size_t)n*n*n;
      const size_t nval = a.npts*a.nc;
      a.data = new float_sw4[nval];
      // first touch as Sarray::set_to_zero
      for( int c=0 ; c < a.nc ; c++ )
#pragma omp parallel for collapse(2)
         for( int k=0 ; k < n ; k++ )
            for( int j=0 ; j < n ; j++ )
               for( int i=0 ; i < n ; i++ )
                  a.data[i+(size_t)n*j+(size_t)n*n*k+a.npts*c] = 0;
      double* ar = new double[nval];
      for( size_t i=0 ; i < nval ; i++ )   // one thread, as in get_exact_lamb2
         ar[i] = 0;
      for( size_t i=0 ; i < nval ; i++ )
         ar[i] = (double)(i % 100003) + 0.25;

      std::vector<float_sw4> ref(nval);
      assign_orig(a, ar);
      memcpy(ref.data(), a.data, nval*sizeof(float_sw4));

      const int reps = std::max(3, (int)(2e9/(nval*12.0)));   // ~2 GB moved per variant
      const double bytes = nval*(sizeof(double)+sizeof(float_sw4));
      for( int v=0 ; v < nvar ; v++ )
      {
         memset(a.data, 0, nval*sizeof(float_sw4));
         fns[v](a, ar);   // warm-up and check
         bool ok = memcmp(ref.data(), a.data, nval*sizeof(float_sw4)) == 0;
         double tmin = 1e30, tsum = 0;
         for( int r=0 ; r < reps ; r++ )
         {
            double t0 = omp_get_wtime();
            fns[v](a, ar);
            double t = omp_get_wtime() - t0;
            tmin = std::min(tmin, t);
            tsum += t;
         }
         printf("%6d %10.1f %-11s %10.3f %10.3f %9.2f %s\n", n, nval*8/1e6, names[v],
                1e3*tmin, 1e3*tsum/reps, bytes/tmin/1e9, ok ? "ok" : "MISMATCH");
      }
      fflush(stdout);
      delete[] ar;
      delete[] a.data;
   }
   return 0;
}
