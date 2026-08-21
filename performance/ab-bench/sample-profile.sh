#!/usr/bin/env bash
#
# Sampling profiler for SW4 that degrades gracefully through whatever the site
# actually permits, and says which method it used.
#
#   ./sample-profile.sh --probe                      # just report capabilities
#   ./sample-profile.sh --bin BIN --input CASE.in [--ranks N] [--threads N]
#
# Why not just perf: unprivileged perf needs kernel.perf_event_paranoid <= 2,
# and HPC sites commonly set 3 (or omit perf from compute images entirely).
# ptrace-based stack sampling needs kernel.yama.ptrace_scope <= 1 AND the
# sampler to be an ancestor of the target when it is 1. So we probe, pick the
# best available, and record the choice next to the numbers.
#
set -uo pipefail

BIN=""; INPUT=""; RANKS=1; THREADS="${OMP_NUM_THREADS:-4}"
HZ=99; OUT="prof-$$"; PROBE=0; METHOD="auto"; DURATION=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bin) BIN="$2"; shift 2;;
    --input) INPUT="$2"; shift 2;;
    --ranks) RANKS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --hz) HZ="$2"; shift 2;;
    --outdir) OUT="$2"; shift 2;;
    --method) METHOD="$2"; shift 2;;   # perf | pmp | phases
    --probe) PROBE=1; shift;;
    -h|--help) sed -n '3,16p' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done

# ------------------------------------------------------------------ probe ----
PARANOID="$(cat /proc/sys/kernel/perf_event_paranoid 2>/dev/null || echo unknown)"
PTRACE="$(cat /proc/sys/kernel/yama/ptrace_scope 2>/dev/null || echo 0)"
have() { command -v "$1" >/dev/null 2>&1; }

probe_report() {
  echo "=== sampling capability probe on $(hostname) ==="
  printf "  %-26s %s\n" "perf_event_paranoid" "$PARANOID  (<=2 allows unprivileged user-space sampling)"
  printf "  %-26s %s\n" "yama/ptrace_scope"   "$PTRACE  (0 or 1 allows attaching to your own children)"
  for t in perf eu-stack gdb pstack vtune amduprof hpcrun scorep tau_exec; do
    printf "  %-26s %s\n" "$t" "$(have $t && command -v $t || echo '-')"
  done
  echo "  module-available profilers:"
  (module avail 2>&1 | grep -oiE '\b(vtune|advisor|uprof|hpctoolkit|scorep|tau|papi|arm-forge|map|likwid)[^ ]*' \
     | sort -u | sed 's/^/    /') 2>/dev/null || echo "    (no module command)"
  echo
  echo "  chosen method: $(choose_method)"
  echo "  note: SW4's own per-phase timer (developer reporttiming=1) always"
  echo "        works and needs no privileges; it is enabled regardless."
}

choose_method() {
  if [ "$METHOD" != "auto" ]; then echo "$METHOD"; return; fi
  if have perf && [ "$PARANOID" != unknown ] && [ "$PARANOID" -le 2 ] 2>/dev/null; then
    echo perf
  elif have eu-stack || have gdb; then
    echo pmp
  else
    echo phases
  fi
}

[ "$PROBE" = 1 ] && { probe_report; exit 0; }
[ -z "$BIN" ] || [ -z "$INPUT" ] && { echo "need --bin and --input (or --probe)" >&2; exit 2; }
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
probe_report | tee "$OUT/capability.txt"
M="$(choose_method)"

LAUNCH=(); if [ "$RANKS" -gt 1 ]; then
  if have srun && [ -n "${SLURM_JOB_ID:-}" ]; then LAUNCH=(srun --ntasks="$RANKS" --cpus-per-task="$THREADS")
  else LAUNCH=(mpirun -np "$RANKS" --oversubscribe); [ "$(id -u)" = 0 ] && LAUNCH+=(--allow-run-as-root); fi
fi
export OMP_NUM_THREADS="$THREADS" OMP_PROC_BIND=close OMP_PLACES=cores

case "$M" in
# ---------------------------------------------------------------- perf -------
perf)
  echo "sampling with perf at ${HZ}Hz"
  # Only rank 0 is recorded: N concurrent perf sessions on one node distort each
  # other, and every rank runs the same code path.
  if [ "$RANKS" -gt 1 ]; then
    cat > "$OUT/wrap.sh" <<WEOF
#!/bin/bash
r="\${SLURM_PROCID:-\${OMPI_COMM_WORLD_RANK:-0}}"
if [ "\$r" = "0" ]; then exec perf record -F $HZ -g --call-graph=dwarf -o "$OUT/perf.data" -- "\$@";
else exec "\$@"; fi
WEOF
    chmod +x "$OUT/wrap.sh"
    "${LAUNCH[@]}" "$OUT/wrap.sh" "$BIN" "$INPUT" > "$OUT/run.log" 2>&1
  else
    perf record -F "$HZ" -g --call-graph=dwarf -o "$OUT/perf.data" -- \
      "$BIN" "$INPUT" > "$OUT/run.log" 2>&1
  fi
  if [ -s "$OUT/perf.data" ]; then
    perf report -i "$OUT/perf.data" --stdio --sort=symbol --percent-limit=0.3 \
      -g none --no-children > "$OUT/flat.txt" 2>/dev/null
    perf report -i "$OUT/perf.data" --stdio -g folded --percent-limit=1 \
      > "$OUT/callgraph.txt" 2>/dev/null
    # Demangled C++ template signatures run to hundreds of columns, which makes
    # the raw report unreadable. Keep percentage + a trimmed symbol.
awk '
      /^ *[0-9]+\.[0-9]+%/{
        pct=$1; sub(/^ *[0-9.]+% *(\[[^]]*\] *)?/,""); name=$0
        # the OP template parameter prints as <(char)61>; render it readably
        # BEFORE stripping the argument list, since it contains parentheses
        gsub(/<\(char\)61>/,"<'='>",name)
        gsub(/<\(char\)43>/,"<'+'>",name)
        gsub(/<\(char\)45>/,"<'-'>",name)
        sub(/\(int.*/,"",name); sub(/\(void.*/,"",name); sub(/\(\).*/,"",name)
        sub(/^void /,"",name); sub(/ *\[clone[^]]*\]/,"",name)
        gsub(/std::(vector|allocator|__cxx11::)/,"",name)
        gsub(/ +$/,"",name)
        # Unresolved addresses are almost always the OpenMP runtime spinning at
        # a barrier or MPI internals without symbols. Individually they are
        # noise; summed they are worth one line.
        if (name ~ /^0x[0-9a-f]+/) { unres += pct+0; next }
        if (length(name)>72) name=substr(name,1,69) "..."
        printf "  %-8s %s\n", pct, name
      }
      END { if (unres > 0) printf "  %-8s %s\n", sprintf("%.2f%%", unres),
                   "(unresolved addresses, summed - likely libgomp barrier / MPI internals)" }
    ' "$OUT/flat.txt" > "$OUT/top.txt"
    echo "  -> $OUT/top.txt (digest), flat.txt, callgraph.txt"
  else
    echo "  perf produced no data; falling back to phase timings only" >&2
    M=phases
  fi
  ;;
# --------------------------------------------------- poor man's profiler -----
pmp)
  # Periodically dump the target's user-space stacks. No perf_event access
  # needed, only ptrace of our own child -- which is what ptrace_scope=1 allows.
  DUMP=$(have eu-stack && echo eu-stack || echo gdb)
  echo "sampling with the ptrace poor-man's profiler ($DUMP), ~${HZ} samples"
  "${LAUNCH[@]}" "$BIN" "$INPUT" > "$OUT/run.log" 2>&1 &
  APP=$!
  sleep 3   # let setup and the build-up of grids finish
  n=0
  while kill -0 "$APP" 2>/dev/null && [ "$n" -lt "$HZ" ]; do
    # profile the leaf sw4 processes, not the launcher
    for pid in $(pgrep -P "$APP" -f "$(basename "$BIN")" 2>/dev/null || echo "$APP"); do
      if [ "$DUMP" = eu-stack ]; then
        eu-stack -p "$pid" 2>/dev/null
      else
        gdb -p "$pid" -batch -ex "thread apply all bt" 2>/dev/null
      fi
    done >> "$OUT/stacks.txt"
    n=$((n+1)); sleep 0.1
  done
  wait "$APP" 2>/dev/null
  # Flat profile: count leaf-most SW4 frames. Crude, but it answers "which
  # kernel dominates", which is the question.
  grep -oE '\b(rhs4th3fort[a-z0-9_]*|curvilinear4sg[a-z0-9_]*|addsgd[0-9a-z_]*|bcfort[a-z0-9_]*|pack_x|unpack_x|[a-z_]+JacobiOptD|predfort_ci|corrfort_ci|dpdmtfort_ci|satt_ci|MPI_[A-Za-z_]+|memcpy)\b' \
    "$OUT/stacks.txt" 2>/dev/null | sort | uniq -c | sort -rn > "$OUT/flat.txt"
  echo "  -> $OUT/flat.txt ($(wc -l < "$OUT/stacks.txt" 2>/dev/null || echo 0) stack lines)"
  ;;
esac

# ---------------------------------------------------------------- phases -----
# Always available and always worth having: SW4's own per-phase breakdown.
if [ "$M" = phases ] || [ ! -s "$OUT/flat.txt" ]; then
  echo "collecting SW4's own phase breakdown only (no external sampler usable)"
  grep -q "reporttiming" "$INPUT" || echo "developer reporttiming=1" >> "$INPUT"
  "${LAUNCH[@]}" "$BIN" "$INPUT" > "$OUT/run.log" 2>&1
fi
grep -A1 'Total *Div-stress' "$OUT/run.log" | tail -2 > "$OUT/phases.txt" 2>/dev/null
echo
echo "method used: $M"
echo "results in:  $OUT"
[ -s "$OUT/top.txt" ]    && { echo "--- top frames (self time) ---"; head -18 "$OUT/top.txt"; }
[ -s "$OUT/top.txt" ] || { [ -s "$OUT/flat.txt" ] && { echo "--- top frames ---"; head -18 "$OUT/flat.txt"; }; }
[ -s "$OUT/phases.txt" ] && { echo "--- phases ---";      cat "$OUT/phases.txt"; }
