#!/usr/bin/env bash
#
# Archive prior results, submit the full A/B matrix, then watch it.
#
#   ./performance/ab-bench/slurm/run-all.sh -A nesi00213
#
# Jobs are chained with --dependency=afterany within each partition, so no two
# benchmark jobs ever share a node. That matters: an earlier round had two of
# our own jobs co-tenant on c006 and another overlapping the 4-node comm job,
# and for the memory-bound Cartesian cases a neighbour with varying load is
# exactly the quantity being measured. Serialising costs wall-clock and buys
# numbers worth reading.
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
MATRIX=(
  "prod-genoa-sp|$P_GENOA|single|16|4|$AB|"
  "prod-genoa-dp|$P_GENOA|double|16|4|$AB|"
  "omp32-genoa-sp|$P_GENOA|single|2|32|$AB|"
  "comm-genoa-sp|$P_GENOA|single|16|4|$COMM|--nodes=4"
  "prod-milan-sp|$P_MILAN|single|16|4|$AB|"
  "prod-milan-dp|$P_MILAN|double|16|4|$AB|"
)

if [ "$DO_SUBMIT" = 1 ]; then
  : > "$STATE"
  declare -A LASTDEP=()
  echo
  printf "%-16s %-7s %-7s %-9s %s\n" JOB PART PREC LAYOUT JOBID
  for row in "${MATRIX[@]}"; do
    IFS='|' read -r name part prec ranks cpt script extra <<< "$row"
    dep="${LASTDEP[$part]:-}"
    # env vars have to be exported into sbatch's environment, and RANKS vs
    # RANKS_PER_NODE differs per script, so build the call explicitly
    rankvar_name=RANKS; [ "$script" = "$COMM" ] && rankvar_name=RANKS_PER_NODE
    nodes="--nodes=1"; tasks="--ntasks=$ranks"
    [ -n "$extra" ] && { nodes="$extra"; tasks="--ntasks-per-node=$ranks"; }
    depflag=(); [ -n "$dep" ] && depflag=(--dependency=afterany:"$dep")
    if [ "$DRY" = 1 ]; then
      echo "DRY: PRECISION=$prec $rankvar_name=$ranks SIZE=$SIZE sbatch -A $ACCOUNT -p $part -J sw4-$name $nodes $tasks --cpus-per-task=$cpt --mem=$MEM ${depflag[*]} $script"
      continue
    fi
    jid=$(env PRECISION="$prec" "$rankvar_name=$ranks" SIZE="$SIZE" REPS="$REPS" \
              STEPS="$STEPS" ${MODULES_ENV:+MODULES="$MODULES_ENV"} \
          sbatch --parsable -A "$ACCOUNT" -p "$part" -J "sw4-$name" -t "$WALL" \
              $nodes $tasks --cpus-per-task="$cpt" --mem="$MEM" --hint=nomultithread \
              "${depflag[@]}" "$script" 2>&1) || { echo "  SUBMIT FAILED $name: $jid"; continue; }
    LASTDEP[$part]="$jid"
    echo "$jid|sw4-$name|$part|$prec|${ranks}x${cpt}" >> "$STATE"
    printf "%-16s %-7s %-7s %-9s %s\n" "$name" "$part" "$prec" "${ranks}x${cpt}" "$jid"
  done
  [ "$DRY" = 1 ] && exit 0
  echo
  echo "chained with --dependency=afterany within each partition, so no two"
  echo "benchmark jobs share a node. state: $STATE"
fi

# ------------------------------------------------------------------ watch ----
[ "$DO_WATCH" = 0 ] && exit 0
[ -f "$STATE" ] || { echo "no $STATE; nothing to watch"; exit 0; }

dash() {
  clear 2>/dev/null || true
  echo "=== sw4 ab-bench $(date '+%F %T') ============================================"
  printf "%-16s %-7s %-6s %-9s %-10s %s\n" JOB PART PREC LAYOUT STATE RESULT
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
    printf "%-16s %-7s %-6s %-9s %-10s %s\n" "${name#sw4-}" "$part" "$prec" "$layout" "$st" "$res"
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
