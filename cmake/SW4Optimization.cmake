# ---------------------------------------------------------------------------
# SW4 release-mode optimisation flags
# ---------------------------------------------------------------------------
#
# Everything selected here is intended to be NUMERICALLY INERT: it may change
# how fast the code runs and which instructions are emitted, but not which
# floating-point values come out. That rules out the whole -ffast-math family,
# and not merely on general principle:
#
#   -ffast-math implies -ffinite-math-only, which lets the compiler fold
#   std::isnan(x) to false. SW4 has ~12 isnan() guards -- material.C aborts on
#   "non-finite (NaN/Inf) material value(s)", ProjectMtrl.C counts NaNs in
#   rho/mu/lambda, lbfgs.C tests the line-search step. Fast-math does not make
#   those checks fail loudly; it deletes them. A run that should have aborted
#   instead produces confident garbage.
#
# Also excluded for the same reason: -Ofast, -funsafe-math-optimizations,
# -fassociative-math, -freciprocal-math, -fno-signed-zeros, -ffinite-math-only,
# and Intel's -fp-model=fast (which is icpx's DEFAULT -- see below).
#
# The one deliberate exception is FMA contraction. Both GCC and icpx default to
# contracting a*b+c into a single fused instruction, which rounds once instead
# of twice. That changes results relative to a non-FMA build -- generally for
# the better, but it is not bit-reproducible across architectures. Set
# SW4_STRICT_FP=ON to turn contraction off and get reproducibility instead.
# ---------------------------------------------------------------------------

set(SW4_TARGET "auto" CACHE STRING
    "Deployment target. auto = sniff the environment; native = this machine.")
set_property(CACHE SW4_TARGET PROPERTY STRINGS
    auto native generic
    mn5-gpp cascade
    hpc3-genoa hpc3-milan hpc3-portable
    frontera stampede3 stampede3-spr vista)

option(SW4_STRICT_FP
       "Disable FMA contraction for bit-reproducible results across machines." OFF)
option(SW4_LTO
       "Enable link-time / interprocedural optimisation. Slower builds." OFF)
option(SW4_OPT_REPORT
       "Emit vectorisation/optimisation reports at compile time." OFF)

# --- refuse to co-operate with an unsafe user-supplied flag -----------------
foreach(_bad_flag "-ffast-math" "-Ofast" "-funsafe-math-optimizations"
                  "-fassociative-math" "-ffinite-math-only" "-fp-model=fast"
                  "-fp-model fast")
  if("${CMAKE_CXX_FLAGS} ${SW4_ARCH_FLAGS} ${SW4_EXTRA_RELEASE_FLAGS}" MATCHES "${_bad_flag}")
    message(FATAL_ERROR
      "${_bad_flag} was requested, but it implies finite-math-only, which folds "
      "std::isnan() to false and silently disables SW4's NaN guards in "
      "material.C / ProjectMtrl.C / lbfgs.C. Remove it. If you want maximum "
      "speed, the flags this module selects are already as aggressive as is "
      "safe for this code.")
  endif()
endforeach()

# ---------------------------------------------------------------------------
# 1. Resolve SW4_TARGET
# ---------------------------------------------------------------------------
# Auto-detection is best-effort: it reads site-set environment variables, which
# are present on login nodes but are NOT a substitute for knowing your machine.
# Login-node and compute-node CPUs differ on several of these systems, so
# -march=native is specifically NOT the auto fallback. Verify with
# `cmake -LH | grep SW4_RESOLVED` before trusting a production build.
if(SW4_TARGET STREQUAL "auto")
  set(_t "generic")
  if(DEFINED ENV{BSC_MACHINE})
    if("$ENV{BSC_MACHINE}" MATCHES "mn5")
      set(_t "mn5-gpp")
    endif()
  endif()
  if(DEFINED ENV{TACC_SYSTEM})
    if("$ENV{TACC_SYSTEM}" STREQUAL "frontera")
      set(_t "frontera")
    elseif("$ENV{TACC_SYSTEM}" MATCHES "stampede3")
      set(_t "stampede3")
    elseif("$ENV{TACC_SYSTEM}" MATCHES "vista")
      set(_t "vista")
    endif()
  endif()
  # REANNZ HPC3 and ESNZ Cascade expose no documented site variable, so fall
  # back to the Slurm cluster name and then to the module-system name.
  foreach(_v SLURM_CLUSTER_NAME LMOD_SYSTEM_NAME CLUSTER_NAME)
    if(DEFINED ENV{${_v}})
      string(TOLOWER "$ENV{${_v}}" _cn)
      if(_cn MATCHES "hpc3|reannz|mahuika")
        set(_t "hpc3-portable")   # heterogeneous site: default to the safe build
      elseif(_cn MATCHES "cascade")
        set(_t "cascade")
      endif()
    endif()
  endforeach()
  set(SW4_RESOLVED_TARGET "${_t}" CACHE INTERNAL "" FORCE)
else()
  set(SW4_RESOLVED_TARGET "${SW4_TARGET}" CACHE INTERNAL "" FORCE)
endif()

# ---------------------------------------------------------------------------
# 2. Per-target architecture flags
# ---------------------------------------------------------------------------
# _gnu_arch  : GCC / Clang -march-family flags
# _intel_arch: Intel oneAPI (icpx) -x-family flags
# _zmm       : ON if 512-bit vectors are a win on this part.
#
# On Skylake-SP and Cascade Lake, 512-bit ops trigger license-based frequency
# throttling, so GCC's 256-bit default is usually faster -- _zmm stays OFF
# there. Sapphire Rapids throttles far less, and Zen 4 implements AVX-512 on a
# 256-bit datapath so there is no downclock to avoid; both want 512.
set(_gnu_arch "")
set(_intel_arch "")
set(_zmm OFF)
set(_target_note "")

if(SW4_RESOLVED_TARGET STREQUAL "mn5-gpp" OR
   SW4_RESOLVED_TARGET STREQUAL "stampede3-spr")
  # BSC MareNostrum 5 GPP: 2x Xeon Platinum 8480+, 112 c/node, 8ch DDR5.
  # TACC Stampede3 SPR: 2x Xeon CPU Max 9480 (HBM), 112 c/node.
  # Golden Cove: genuine 2x512-bit FMA units -- the best AVX-512 target here.
  set(_gnu_arch "-march=sapphirerapids" "-mtune=sapphirerapids")
  set(_intel_arch "-xSAPPHIRERAPIDS")
  set(_zmm ON)
  set(_target_note "Sapphire Rapids, 2x512-bit FMA, AVX-512 at full width")

elseif(SW4_RESOLVED_TARGET STREQUAL "cascade" OR
       SW4_RESOLVED_TARGET STREQUAL "hpc3-genoa")
  # ESNZ Cascade: HPE Cray XD2000, AMD 4th Gen EPYC (Genoa/Zen 4).
  # REANNZ HPC3 genoa partition: 2x EPYC 9634, 168 c/node.
  # Zen 4 has the AVX-512 ISA on a 256-bit datapath (double-pumped): half the
  # per-core vector throughput of SPR, but 12 memory channels per socket
  # instead of 8, which suits a bandwidth-bound stencil.
  set(_gnu_arch "-march=znver4" "-mtune=znver4")
  set(_intel_arch "-march=core-avx2")   # icpx -x targets are Intel-only
  set(_zmm ON)
  set(_target_note "Zen 4, AVX-512 on a 256-bit datapath, 12 mem ch/socket")

elseif(SW4_RESOLVED_TARGET STREQUAL "hpc3-milan")
  # REANNZ HPC3 milan partition: 2x EPYC 7713, 128 c/node. Zen 3: AVX2 ONLY.
  set(_gnu_arch "-march=znver3" "-mtune=znver3")
  set(_intel_arch "-march=core-avx2")
  set(_target_note "Zen 3, AVX2 only -- no AVX-512 on this partition")

elseif(SW4_RESOLVED_TARGET STREQUAL "hpc3-portable")
  # HPC3 is heterogeneous: the genoa partition is Zen 4 (AVX-512) and the milan
  # partition is Zen 3 (AVX2 only). A -march=znver4 binary SIGILLs on milan.
  # x86-64-v3 is the common ISA floor; -mtune=znver4 still schedules for the
  # newer part without emitting instructions milan cannot execute.
  set(_gnu_arch "-march=x86-64-v3" "-mtune=znver4")
  set(_intel_arch "-march=core-avx2")
  set(_target_note "x86-64-v3 baseline: one binary that runs on genoa AND milan")

elseif(SW4_RESOLVED_TARGET STREQUAL "frontera")
  # TACC Frontera: 8368 CLX nodes, 2x Xeon Platinum 8280, 56 c/node.
  set(_gnu_arch "-march=cascadelake" "-mtune=cascadelake")
  set(_intel_arch "-xCASCADELAKE")
  set(_target_note "Cascade Lake; 256-bit preferred to avoid AVX-512 downclock")

elseif(SW4_RESOLVED_TARGET STREQUAL "stampede3")
  # Stampede3 mixes SKX (Platinum 8160, 48c), ICX (8380, 80c) and SPR-HBM
  # (Max 9480, 112c). All have AVX-512 but differ in generation; skylake-avx512
  # is the common floor. Use SW4_TARGET=stampede3-spr for an SPR-only build.
  set(_gnu_arch "-march=skylake-avx512" "-mtune=sapphirerapids")
  set(_intel_arch "-xCOMMON-AVX512")
  set(_target_note "SKX/ICX/SPR common AVX-512 floor; see stampede3-spr")

elseif(SW4_RESOLVED_TARGET STREQUAL "vista")
  # TACC Vista GG: 256 nodes, NVIDIA Grace Superchip, 2x72 = 144 cores,
  # Neoverse V2, ~850 GiB/s memory bandwidth. This is ARM: there is no -march,
  # no AVX-512, and icpx does not target it. Highest memory bandwidth of any
  # target here, which is exactly what SW4's stencil wants.
  set(_gnu_arch "-mcpu=neoverse-v2")
  set(_target_note "ARM Neoverse V2 (Grace), SVE2, ~850 GiB/s -- aarch64, not x86")

elseif(SW4_RESOLVED_TARGET STREQUAL "native")
  set(_gnu_arch "-march=native" "-mtune=native")
  set(_intel_arch "-xHost")
  set(_target_note "this build host -- do NOT use if login and compute CPUs differ")

else()  # generic
  set(_target_note "no architecture flags; portable but leaves SIMD on the table")
endif()

# ---------------------------------------------------------------------------
# 3. Numerically inert optimisation flags, per compiler
# ---------------------------------------------------------------------------
set(_cxx_opt "")
set(_fort_opt "")

if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  # icpx defaults to -fp-model=fast, i.e. Intel's out-of-the-box behaviour is
  # already the thing this module exists to prevent. Pin it explicitly.
  list(APPEND _cxx_opt -fp-model=precise -fno-math-errno)
  # The ported kernels are enormous (rhs4th3fort_ci is thousands of
  # instructions); icpx gives up on optimising functions past an internal size
  # limit unless told not to.
  list(APPEND _cxx_opt -qoverride-limits -qopt-prefetch=4)
  if(_zmm)
    list(APPEND _cxx_opt -qopt-zmm-usage=high)
  endif()
  if(SW4_OPT_REPORT)
    list(APPEND _cxx_opt -qopt-report=3)
  endif()
  set(_fort_opt -fp-model=precise)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND _cxx_opt
       -funroll-loops           # was already the project default
       -fno-math-errno          # lets sqrt/pow vectorise; only errno changes
       -fno-trapping-math       # asserts no FP traps; values unaffected
       -fipa-pta                # stronger interprocedural alias analysis
       -fno-semantic-interposition
       -falign-functions=32
       -falign-loops=32)
  if(_zmm)
    list(APPEND _cxx_opt -mprefer-vector-width=512)
  endif()
  if(SW4_OPT_REPORT)
    list(APPEND _cxx_opt -fopt-info-vec-optimized)
  endif()
  set(_fort_opt -funroll-loops -fno-math-errno)

elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  list(APPEND _cxx_opt -funroll-loops -fno-math-errno -fno-trapping-math)
  if(_zmm)
    list(APPEND _cxx_opt -mprefer-vector-width=512)
  endif()
  if(SW4_OPT_REPORT)
    list(APPEND _cxx_opt -Rpass=loop-vectorize)
  endif()

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
  # The sensible CPU compiler on Vista alongside GCC.
  list(APPEND _cxx_opt -Mvect -Munroll -Kieee)   # -Kieee = strict IEEE, no fast-math
  set(_fort_opt -Kieee)
endif()

# --- FP contraction --------------------------------------------------------
if(SW4_STRICT_FP)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    list(APPEND _cxx_opt -fp-model=strict -ffp-contract=off)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    list(APPEND _cxx_opt -ffp-contract=off)
    list(APPEND _fort_opt -ffp-contract=off)
  endif()
endif()

# ---------------------------------------------------------------------------
# 4. Publish
# ---------------------------------------------------------------------------
# SW4_ARCH_FLAGS stays an honoured manual override: if the user set it, we do
# not second-guess them.
if(SW4_ARCH_FLAGS STREQUAL "")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" AND _intel_arch)
    set(SW4_RESOLVED_ARCH_FLAGS ${_intel_arch})
  else()
    set(SW4_RESOLVED_ARCH_FLAGS ${_gnu_arch})
  endif()
else()
  set(SW4_RESOLVED_ARCH_FLAGS ${SW4_ARCH_FLAGS})
  set(_target_note "${_target_note} [overridden by SW4_ARCH_FLAGS]")
endif()

set(SW4_RESOLVED_CXX_OPT  ${_cxx_opt}  CACHE INTERNAL "" FORCE)
set(SW4_RESOLVED_FORT_OPT ${_fort_opt} CACHE INTERNAL "" FORCE)

if(SW4_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg)
  if(_ipo_ok)
    # include() does not open a new scope, so this lands in the caller's.
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
  else()
    message(WARNING "SW4_LTO requested but unsupported: ${_ipo_msg}")
  endif()
endif()

message(STATUS "SW4 target:      ${SW4_RESOLVED_TARGET}  (${_target_note})")
message(STATUS "SW4 arch flags:  ${SW4_RESOLVED_ARCH_FLAGS}")
message(STATUS "SW4 opt flags:   ${SW4_RESOLVED_CXX_OPT}")
message(STATUS "SW4 strict FP:   ${SW4_STRICT_FP}   LTO: ${SW4_LTO}")

# ---------------------------------------------------------------------------
# 5. Heavily-nested kernels: raise the complete-unroll budget
# ---------------------------------------------------------------------------
# The SBP boundary-closure bodies in the curvilinear kernels contain
# `for(q=1..8){ for(m=1..8){ ... } }`. GCC's complete-unroll pass runs AFTER the
# vectoriser, and the early pass that could flatten these is gated by
# --param=max-completely-peeled-insns, default 200 -- far below these bodies. So
# GCC reports "loop nest containing two or more consecutive inner loops cannot
# be vectorized" and the closures come out scalar while the interior loops right
# next to them vectorise fine. innerloop-ani-sgstr-vcc.C comes out entirely
# scalar: the vectoriser never even reaches its 2nd and 3rd loops.
#
# This is the GCC analogue of a documented upstream problem: Johansen et al.,
# CSE 19(5) 2017 record the Intel compiler "silently aborting its vectorisation
# pass" on this same code, fixed by overriding optimisation limits. We already
# pass -qoverride-limits to icpx above; this is the missing GCC counterpart.
#
# Measured with -march=znver4 -mprefer-vector-width=512 (200 -> 4000):
#   curvilinear4sgc.C          zmm  5613 -> 19927   text 120K -> 336K   7.3s -> 34.3s
#   curvilinear4sgwind.C       zmm  5641 -> 19980   text 123K -> 341K   7.4s -> 34.4s
#   ilanisocurvc.C             zmm  3859 -> 15395   text 116K -> 326K  11.2s -> 40.7s
#   innerloop-ani-sgstr-vcc.C  zmm     0 ->  5898   text  21K -> 122K   2.0s -> 14.4s
#
# Applied per-source rather than globally: the compile-time and code-size cost
# is real (~5x on these four files) and there is no reason to pay it on the
# other ~200 translation units. 2000 was measured to be insufficient; 4000 is
# enough on both sapphirerapids and znver4.
#
# NOT bit-exact under the default -ffp-contract=fast: more unrolling exposes
# more FMA-contraction opportunities, which shifts results by ~1e-9 relative.
# That is the one class of deviation this project already treats as
# architecture-dependent and gates behind SW4_STRICT_FP=ON, under which this
# option IS bit-identical.
# rhs4th3fortc.C has the same blocker, from a different construct: multiple
# CONSECUTIVE for(q=1..8) loops in the one-sided closure bodies rather than a
# q/m nest (the m-loops there were already hand-unrolled). Same diagnostic, same
# fix. Verified -march=sapphirerapids -mprefer-vector-width=512, 200 -> 4000:
#   "cannot be vectorized" diagnostics  24 -> 0   (all four closure loops)
#   zmm  9069 -> 38336      packed FP  4942 -> 17206
#   text 351K -> 820K       compile    20.2s -> 73.6s
# This matters more than it looks. EW.C:1023 turns a closure on for every
# bRefInterface, so in a mesh-refinement run each Cartesian level has closures
# on both k-sides, and refinement levels are deliberately thin -- the repo's own
# performance/MR_OMP/rfile/berkeley-mr-cart.in gives nk of 33, 49 and 30. At
# those geometries the closures are 68-81% of the kernel call, not a fringe.
# Note the other Cartesian kernels (rhs4th3windc.C, rhs4th3windc2.C,
# rhs4th3point.C) do NOT have the blocker and are deliberately not listed.
set(SW4_UNROLL_HEAVY_SOURCES
    src/curvilinear4sgc.C
    src/curvilinear4sgwind.C
    src/ilanisocurvc.C
    src/innerloop-ani-sgstr-vcc.C
    src/rhs4th3fortc.C)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(SW4_UNROLL_HEAVY_FLAGS --param=max-completely-peeled-insns=4000)
else()
  set(SW4_UNROLL_HEAVY_FLAGS "")   # icpx already gets -qoverride-limits
endif()

if(SW4_UNROLL_HEAVY_FLAGS AND CMAKE_BUILD_TYPE MATCHES "Release|RelWithDebInfo")
  set_source_files_properties(${SW4_UNROLL_HEAVY_SOURCES}
      PROPERTIES COMPILE_OPTIONS "${SW4_UNROLL_HEAVY_FLAGS}")
  message(STATUS "SW4 unroll:      ${SW4_UNROLL_HEAVY_FLAGS} on 4 heavily-nested kernels")
endif()
