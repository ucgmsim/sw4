// OMPT tool (OpenMP 5.0 tools interface) that charges OpenMP time to the
// running sw4timer region (see sw4timer.h).
//
// The OpenMP runtime calls ompt_start_tool() at start-up. For every top-level
// parallel region the master thread records its wall time and team size and
// charges them to the innermost running timer region. Every thread in the team
// records its busy time: from the start of its implicit task until it starts
// waiting at a barrier, and again after each wait ends.
//
// Busy time is measured, rather than waiting time, because LLVM's runtime
// reports the end of a worker's final implicit barrier only when the worker is
// released into the next parallel region. A worker's measured wait would
// therefore include the serial code between parallel regions. Its wait
// *begin* is reported promptly, so busy time is exact; waiting is then team
// time minus busy time.
//
// Nested parallel regions are attributed to their enclosing top-level region
// and not tracked separately. Compiled only with SW4_HAVE_OMPT (CMake detects
// omp-tools.h); a runtime without OMPT (GCC's libgomp) never calls
// ompt_start_tool and the report omits the OpenMP table.
#include "sw4timer.h"

#if !defined(SW4_NO_TIMERS) && defined(SW4_HAVE_OMPT)

#include <omp-tools.h>

#include <chrono>

namespace
{
inline double now()
{
   return std::chrono::duration<double>( std::chrono::steady_clock::now().time_since_epoch() ).count();
}

struct ThreadState
{
   sw4timer::OmpThread* data = nullptr;
   int depth = 0;      // implicit-task nesting; only depth 1 is recorded
   int region = -1;    // timer region of the current top-level team
   bool busy = false;  // a busy segment is open
   double busy_begin = 0;
// master of a top-level parallel region
   double par_begin = 0;
   unsigned team_size = 1;
};
thread_local ThreadState ts;

void add_work( double seconds )
{
   if( ts.region < 0 )
      return;
   if( ts.data == nullptr )
      ts.data = sw4timer::omp_register_thread();
   std::vector<double>& w = ts.data->work;
   if( static_cast<int>( w.size() ) <= ts.region )
      w.resize( ts.region + 1, 0.0 );
   w[ts.region] += seconds;
}

void close_busy()
{
   if( ts.busy )
   {
      add_work( now() - ts.busy_begin );
      ts.busy = false;
   }
}

void on_parallel_begin( ompt_data_t*, const ompt_frame_t*, ompt_data_t* parallel_data,
                        unsigned int, int, const void* )
{
   if( ts.depth > 0 )
   {
      parallel_data->value = 0; // nested: part of the enclosing region
      return;
   }
   parallel_data->value = static_cast<uint64_t>( sw4timer::omp_current_region() ) + 1;
   ts.team_size = 1;
   ts.par_begin = now();
}

void on_parallel_end( ompt_data_t* parallel_data, ompt_data_t*, int, const void* )
{
   if( parallel_data == nullptr || parallel_data->value == 0 )
      return;
   sw4timer::omp_parallel_done( static_cast<int>( parallel_data->value - 1 ), now() - ts.par_begin,
                                ts.team_size );
}

void on_implicit_task( ompt_scope_endpoint_t endpoint, ompt_data_t* parallel_data, ompt_data_t*,
                       unsigned int actual_parallelism, unsigned int index, int flags )
{
   if( flags & ompt_task_initial )
      return;
   if( endpoint == ompt_scope_begin )
   {
      if( ++ts.depth != 1 )
         return;
      ts.region = ( parallel_data && parallel_data->value ) ? static_cast<int>( parallel_data->value - 1 ) : -1;
      if( index == 0 )
         ts.team_size = actual_parallelism;
      ts.busy = true;
      ts.busy_begin = now();
   }
   else if( endpoint == ompt_scope_end )
   {
      if( ts.depth == 1 )
      {
         close_busy();
         ts.region = -1;
      }
      if( ts.depth > 0 )
         ts.depth--;
   }
}

void on_sync_region_wait( ompt_sync_region_t, ompt_scope_endpoint_t endpoint, ompt_data_t*,
                          ompt_data_t*, const void* )
{
   if( ts.depth != 1 )
      return;
   if( endpoint == ompt_scope_begin )
      close_busy();
   else if( endpoint == ompt_scope_end && ts.region >= 0 )
   {
// After the final implicit barrier this opens a segment of a few microseconds,
// closed by the implicit-task end that follows immediately.
      ts.busy = true;
      ts.busy_begin = now();
   }
}

int initialize( ompt_function_lookup_t lookup, int, ompt_data_t* )
{
   ompt_set_callback_t set_callback = reinterpret_cast<ompt_set_callback_t>( lookup( "ompt_set_callback" ) );
   if( set_callback == nullptr )
      return 0;
   auto ok = [&]( ompt_callbacks_t event, ompt_callback_t cb )
   {
      const ompt_set_result_t res = set_callback( event, cb );
      return res != ompt_set_error && res != ompt_set_never;
   };
   const bool all =
      ok( ompt_callback_parallel_begin, reinterpret_cast<ompt_callback_t>( &on_parallel_begin ) ) &&
      ok( ompt_callback_parallel_end, reinterpret_cast<ompt_callback_t>( &on_parallel_end ) ) &&
      ok( ompt_callback_implicit_task, reinterpret_cast<ompt_callback_t>( &on_implicit_task ) ) &&
      ok( ompt_callback_sync_region_wait, reinterpret_cast<ompt_callback_t>( &on_sync_region_wait ) );
   if( !all )
      return 0;
   sw4timer::omp_set_active();
   return 1;
}

void finalize( ompt_data_t* )
{
}
} // namespace

extern "C" __attribute__( ( visibility( "default" ) ) ) ompt_start_tool_result_t*
ompt_start_tool( unsigned int, const char* )
{
   static ompt_start_tool_result_t result = { &initialize, &finalize, { 0 } };
   return &result;
}

#endif
