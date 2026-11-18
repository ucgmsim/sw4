//-*-c++-*-
// This file is part of SW4, distributed under the GNU General Public License
// version 2; see LICENCE.txt.
#ifndef SW4_SOURCETIMETABLE_H
#define SW4_SOURCETIMETABLE_H

#include <algorithm>
#include <cstddef>
#include <vector>

#include "sw4.h"

// A source time function (kind 0, for EW::Force) and its second derivative
// (kind 1, for EW::Force_tt) at the times EW::solve evaluates them.
//
// A prefiltered discrete (SRF/rupture) time function is a quintic spline with a
// knot every SRF time step (often 0.005 s) over the series padded by the
// filter precursor on both sides, 6 coefficients per knot. solve() only ever
// evaluates it at its own time steps, which are far coarser, so for a plain
// forward run Source::tabulate_time_function stores the values the spline
// gives there and frees the spline. The values are computed by the same call
// GridPointSource makes, so forces are bitwise unchanged.
struct SourceTimeTable
{
   // Sorted, distinct evaluation times, shared by every table (owned by EW).
   const std::vector<float_sw4>* times[2] = {nullptr, nullptr};
   // values[k][m] is the value at times[k][first[k] + m]; outside that range
   // the spline is clamped and gives before[k] or after[k].
   std::size_t first[2] = {0, 0};
   std::vector<float_sw4> values[2];
   float_sw4 before[2] = {0, 0}, after[2] = {0, 0};

   float_sw4 at( int k, std::size_t n ) const
   {
      if( n < first[k] )
         return before[k];
      const std::size_t m = n - first[k];
      return m < values[k].size() ? values[k][m] : after[k];
   }
};

// Index of t in the sorted times, or -1 if t is not one of them.
inline long tabulated_time_index( const std::vector<float_sw4>& times, float_sw4 t )
{
   const auto it = std::lower_bound( times.begin(), times.end(), t );
   return ( it != times.end() && *it == t ) ? static_cast<long>( it - times.begin() ) : -1;
}

#endif
