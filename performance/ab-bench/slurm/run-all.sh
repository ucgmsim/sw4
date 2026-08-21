#!/usr/bin/env bash
#
# Archive prior results, submit the full A/B matrix, then watch it.
#
#   ./performance/ab-bench/slurm/run-all.sh -A nesi00213
#
# Everything is submitted at once. An earlier version chained jobs with
# --dependency=afterany to stop two of our own runs sharing a node, but that is
# the wrong trade on a busy cluster: other users' jobs contend regardless, so
# serialising ours controls almost nothing while multiplying wall-clock. Only
# --exclusive would actually isolate a run, and a whole free node is rare enough
# that those jobs effectively never schedule.
#
# The mitigation is to measure the noise instead of trying to remove it: the
# summary now reports each phase's own run-to-run spread, so a contended result
# announces itself rather than masquerading as signal.
#
set -uo pipefail

ACCOUNT="${ACCOUNT:-nesi00213}"
P_GENOA="${P_GENOA:-genoa}"
P_MILAN="${P_MILAN:-milan}"
SIZE="${SIZE:-L}"
REPS="${REPS:-5}"
STEPS="${STEPS:-200}"
MEM="${MEM:-96G}"
WALL="${WALL:-3:00:00}"
MODULES_ENV="${MODULES:-}"
DO_ARCHIVE=1; DO_SUBMIT=1; DO_WATCH=1; DRY=0; INTERVAL="${INTERVAL:-90}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -A|--account) ACCOUNT="$2"; shift 2;;
    --size) SIZE="$2"; shift 2;;
    --reps) REPS="$2"; shift 2;;
    --steps) STEPS="$2"; shift 2;;
    --mem) MEM="$2"; shift 2;;
    --wall) WALL="$2"; shift 2;;
    --no-archive) DO_ARCHIVE=0; shift;;
    --no-submit)  DO_SUBMIT=0; shift;;
    --no-watch)   DO_WATCH=0; shift;;
    --watch-only) DO_ARCHIVE=0; DO_SUBMIT=0; shift;;
    --dry-run) DRY=1; shift;;
    --interval) INTERVAL="$2"; shift 2;;
    -h|--help)
      sed -n '3,13p' "$0" | sed 's/^# \{0,1\}//'
      echo "Options: -A ACCOUNT --size S|M|L|XL --reps N --steps N --mem X --wall H:MM:SS"
      echo "         --no-archive --no-submit --no-watch --watch-only --dry-run --interval SEC"
      exit 0;;
    *) echo "unknown option: $1" >&2; exit 2;;
  esac
done

REPO="$(git rev-parse --show-toplevel 2>/dev/null)"
[ -z "$REPO" ] && { echo "not in a git repo (did .git come across in the rsync?)"; exit 1; }
cd "$REPO" || exit 1
AB="performance/ab-bench/slurm/hpc3-ab.sl"
COMM="performance/ab-bench/slurm/hpc3-comm.sl"
for f in "$AB" "$COMM"; do [ -x "$f" ] || chmod +x "$f" 2>/dev/null; done
STATE="$REPO/.ab-jobs"

# ---------------------------------------------------------------- archive ----
if [ "$DO_ARCHIVE" = 1 ]; then
  STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
  ARC="$REPO/ab-archive/$STAMP"
  shopt -s nullglob
  MOVED=0
  for x in "$REPO"/ab-comm-* "$REPO"/ab-genoa-* "$REPO"/ab-milan-* "$REPO"/ab-local-* \
           "$REPO"/ab-results-* "$REPO"/sw4-ab-*.out "$REPO"/sw4-comm-*.out "$STATE"; do
    [ -e "$x" ] || continue
    mkdir -p "$ARC"
    # detach any worktrees the old run registered, else git keeps a stale ref
    for wt in "$x"/wt-base "$x"/wt-head; do
      [ -d "$wt" ] && git worktree remove --force "$wt" >/dev/null 2>&1
    done
    mv "$x" "$ARC"/ && MOVED=$((MOVED+1))
  done
  shopt -u nullglob
  git worktree prune >/dev/null 2>&1
  if [ "$MOVED" -gt 0 ]; then echo "archived $MOVED item(s) -> ab-archive/$STAMP"
  else echo "nothing to archive"; fi
fi

# ----------------------------------------------------------------- submit ----
# name | partition | precision | ranks | cpus-per-task | script | extra
# name | partition | walltime | nodes | ntasks | cpus-per-task | mem | script | env
#
# Production settings throughout: the prod cases carry real supergrid absorbing
# boundaries (gp=30) and attenuation nmech=3, which is what production runs use
# and what every twilight case silently does not.
#
# SIZE=M is 800k points/rank -> nx~234 at 16 ranks. The previous round used
# SIZE=L (nx=372) on a four-case list and four of six jobs hit the time limit
# two reps into curvi-mr, having built in 14-89s. This is ~4x less work with 3x
# the walltime, and the prod cases have no twilight sin/cos so they are much
# cheaper per point.
PROD_ENV="CASES=prod,prod-mr SIZE=M STEPS=150 REPS=4"
MATRIX=(
  # the headline: production settings, production layout, both precisions
  "P-gen-sp|$P_GENOA|6:00:00|1|16|4|96G|$AB|$PROD_ENV PRECISION=single RANKS=16"
  "P-gen-dp|$P_GENOA|6:00:00|1|16|4|96G|$AB|$PROD_ENV PRECISION=double RANKS=16"
  # AVX2 contrast -- do the gains survive without AVX-512 on production settings?
  "P-mil-sp|$P_MILAN|6:00:00|1|16|4|96G|$AB|$PROD_ENV PRECISION=single RANKS=16"
  # Layout pair. NX is pinned so both solve the SAME problem: per-rank sizing
  # would give the 2-rank job a 1/8-size grid and make the absolute times
  # meaningless, which is exactly what happened in the first attempt.
  "P-lay16|$P_GENOA|6:00:00|1|16|4|96G|$AB|CASES=prod STEPS=150 REPS=4 PRECISION=single RANKS=16 NX=235 TIMEVAL=0.2983"
  "P-lay2|$P_GENOA|6:00:00|1|2|32|96G|$AB|CASES=prod STEPS=150 REPS=4 PRECISION=single RANKS=2 NX=235 TIMEVAL=0.2983"
  # single-precision correctness. STRICT_FP removes FMA re-contraction and one
  # thread removes reduction ordering, so this MUST come back bit-exact.
  # Anything else is a real defect and outranks all remaining optimisation work.
  "P-strict|$P_GENOA|1:00:00|1|4|1|32G|$AB|CASES=prod,prod-mr SIZE=S STEPS=100 REPS=2 PRECISION=single STRICT=1 RANKS=4 THREADS_OVERRIDE=1 FORCE=1"
  # falsifiable prediction: SG cost was measured INDEPENDENT of gp, because the
  # sweep covers the whole grid either way. gp=12 and gp=30 should cost the
  # same. If they differ, the premise behind the windowing work is wrong.
  "P-sg12|$P_GENOA|3:00:00|1|16|4|96G|$AB|CASES=prod SIZE=M STEPS=150 REPS=4 PRECISION=single RANKS=16 SW4_BENCH_SGGP=12"
  "P-sg30|$P_GENOA|3:00:00|1|16|4|96G|$AB|CASES=prod SIZE=M STEPS=150 REPS=4 PRECISION=single RANKS=16 SW4_BENCH_SGGP=30"
  # halo exchange at HEAD across nodes. BASE=e6ccbba pins this to the revert
  # point so it measures the hand-packing, not the reverted non-blocking attempt
  # the script default would have compared.
  "P-comm|$P_GENOA|4:00:00|4|16|4|96G|$COMM|BASE=e6ccbba HEAD_REF=HEAD SIZE=M STEPS=150 REPS=4 PRECISION=single RANKS_PER_NODE=16"
)

if [ "$DO_SUBMIT" = 1 ]; then
  : > "$STATE"
  echo
  printf "%-10s %-7s %-6s %-9s %-9s %s\n" JOB PART LAYOUT WALL SIZE JOBID
  for row in "${MATRIX[@]}"; do
    IFS='|' read -r name part wall nodes ntasks cpt mem script envs <<< "$row"
    if [ "$script" = "$COMM" ]; then
      layout="${nodes}n x${ntasks}x${cpt}"
      geom=(--nodes="$nodes" --ntasks-per-node="$ntasks" --cpus-per-task="$cpt")
    else
      layout="${ntasks}x${cpt}"
      geom=(--nodes="$nodes" --ntasks="$ntasks" --cpus-per-task="$cpt")
    fi
    sz=$(echo "$envs" | grep -o 'SIZE=[A-Z]*' | head -1)
    if [ "$DRY" = 1 ]; then
      echo "DRY: env $envs ${MODULES_ENV:+MODULES='$MODULES_ENV' }sbatch -A $ACCOUNT -p $part -J sw4-$name -t $wall ${geom[*]} --mem=$mem --hint=nomultithread $script"
      continue
    fi
    # shellcheck disable=SC2086
    jid=$(env $envs ${MODULES_ENV:+MODULES="$MODULES_ENV"} \
          sbatch --parsable -A "$ACCOUNT" -p "$part" -J "sw4-$name" -t "$wall" \
              "${geom[@]}" --mem="$mem" --hint=nomultithread "$script" 2>&1) \
      || { echo "  SUBMIT FAILED $name: $jid"; continue; }
    echo "$jid|sw4-$name|$part|$layout|$sz" >> "$STATE"
    printf "%-10s %-7s %-6s %-9s %-9s %s\n" "$name" "$part" "$layout" "$wall" "${sz:-}" "$jid"
  done
  [ "$DRY" = 1 ] && exit 0
  echo
  echo "all submitted at once -- no --dependency chaining. Other users' jobs"
  echo "contend regardless, so serialising ours would cost wall-clock without"
  echo "buying isolation; the per-phase spread in the summary reports the noise"
  echo "instead. state: $STATE"
fi

# ------------------------------------------------------------------ watch ----
[ "$DO_WATCH" = 0 ] && exit 0
[ -f "$STATE" ] || { echo "no $STATE; nothing to watch"; exit 0; }

dash() {
  clear 2>/dev/null || true
  echo "=== sw4 ab-bench $(date '+%F %T') ============================================"
  printf "%-14s %-7s %-11s %-8s %-10s %s\n" JOB PART LAYOUT SIZE STATE RESULT
  echo "---------------------------------------------------------------------------------"
  local pending=0
  while IFS='|' read -r jid name part prec layout; do
    [ -z "$jid" ] && continue
    st=$(squeue -h -j "$jid" -o '%T' 2>/dev/null | head -1)
    if [ -z "$st" ]; then
      st=$(sacct -n -j "$jid" -o State%-12 2>/dev/null | head -1 | tr -d ' ')
      [ -z "$st" ] && st="UNKNOWN"
    else
      pending=$((pending+1))
    fi
    out=$(ls -t "$REPO"/sw4-*-"$name"-"$jid".out "$REPO"/sw4-*"$jid".out 2>/dev/null | head -1)
    res="-"
    if [ -n "$out" ] && [ -f "$out" ]; then
      res=$(python3 - "$out" <<'PYX'
import re, sys
t = open(sys.argv[1], errors="replace").read()
if "REFUSING TO RUN" in t: print("refused: bad allocation"); raise SystemExit
if "FAILED:" in t:
    print("build/module failure - see .out"); raise SystemExit
geo = re.search(r'geometry\s+(\d+) cores -> (\d+) ranks x (\d+) threads', t)
ds  = re.findall(r'^(\S+)\s+\S+\s+divstress\s+\S+\s+\S+\s+([0-9.]+)', t, re.M)
nf  = re.findall(r'noise floor: \+/-([0-9.]+)%', t)
bad = ("DIFFERS BEYOND SPREAD" in t) or ("WARNING:" in t and "skew" in t)
short = "is under" in t
bits=[]
if geo: bits.append(f"{geo.group(2)}x{geo.group(3)}t")
if ds:  bits.append(" ".join(f"{c}:{v}x" for c,v in ds[:4]))
if nf:  bits.append(f"noise+/-{max(float(x) for x in nf):.1f}%")
if short: bits.append("TOO-SHORT")
if bad:   bits.append("CORRECTNESS!")
print("  ".join(bits) if bits else "running...")
PYX
) || res="(parse error)"
    fi
    printf "%-14s %-7s %-11s %-8s %-10s %s\n" "${name#sw4-}" "$part" "$layout" "$prec" "$st" "$res"
  done < "$STATE"
  echo "---------------------------------------------------------------------------------"
  echo "$pending job(s) still queued/running.  refresh ${INTERVAL}s.  Ctrl-C to stop watching."
  return $pending
}

while :; do
  dash; left=$?
  [ "$left" -eq 0 ] && { echo; echo "all jobs finished."; break; }
  sleep "$INTERVAL"
done
