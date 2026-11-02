#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
ROOT=$PWD

# Usage: bash run.sh
# Fast workflow check: MODE=pilot SUITE=smoke bash run.sh
# Latest q=7 studies only: SUITE=core bash run.sh
# Resume the exact same sources/settings: RUN_DIR=/absolute/old/run bash run.sh
MODE=${MODE:-research}
SUITE=${SUITE:-paper}
case "$MODE" in
  research) default_seeds=$(seq -s ' ' 1001 1020); default_inputs=16; default_rounds=3 ;;
  pilot) default_seeds='1001 1002'; default_inputs=4; default_rounds=2 ;;
  *) printf 'Unknown MODE: %s\n' "$MODE" >&2; exit 2 ;;
esac
PLAN=${PLAN:-configs/$SUITE.tsv}
[[ -f $PLAN ]] || { printf 'Plan not found: %s\n' "$PLAN" >&2; exit 2; }
SEEDS=${SEEDS:-$default_seeds}
INPUTS=${INPUTS:-$default_inputs}
ROUNDS=${ROUNDS:-$default_rounds}
BATCH=${BATCH:-4}
RANKS=${RANKS:-historical}
ABLATIONS=${ABLATIONS:-1}
GABIDULIN_CONTROLS=${GABIDULIN_CONTROLS:-1}
CAPTURED=${CAPTURED:-1}
ERASURES=${ERASURES:-1}
DUMP_INPUTS=${DUMP_INPUTS:-1}
DENSE_CASES=${DENSE_CASES:-1}
MAX_M=${MAX_M:-512}
PERIOD_CAP=${PERIOD_CAP:-4096}
MAX_SETS=${MAX_SETS:-3}
INFO_ATTEMPTS=${INFO_ATTEMPTS:-7}
MAX_TRIAL_RADIUS=${MAX_TRIAL_RADIUS:-3}
CANDIDATE_POLICY=${CANDIDATE_POLICY:-cyclic}
STAGE_TIMEOUT_SECONDS=${STAGE_TIMEOUT_SECONDS:-1800}
CONSTRUCTION_TIMEOUT_SECONDS=${CONSTRUCTION_TIMEOUT_SECONDS:-7200}
CONSTRUCTION_SEED=${CONSTRUCTION_SEED:-20261027}
CONSTRUCTION_ATTEMPTS=${CONSTRUCTION_ATTEMPTS:-16}
STRICT_M=${STRICT_M:-1}
BOOTSTRAP=${BOOTSTRAP:-2000}
CASE_ORDER_SEED=${CASE_ORDER_SEED:-8422}
read -r -a seed_array <<< "$SEEDS"
(( ${#seed_array[@]} > 0 )) || { printf '%s\n' 'SEEDS must not be empty.' >&2; exit 2; }
for number in "$INPUTS" "$ROUNDS" "$BATCH" "${seed_array[@]}"; do
  [[ $number =~ ^[0-9]+$ ]] || { printf 'Invalid integer: %s\n' "$number" >&2; exit 2; }
done
for number in "$INPUTS" "$ROUNDS" "$BATCH"; do (( number > 0 )) || exit 2; done
unique_seed_count=$(printf '%s\n' "${seed_array[@]}" | sort -u | wc -l)
(( unique_seed_count == ${#seed_array[@]} )) || { printf '%s\n' 'SEEDS contains duplicates.' >&2; exit 2; }

bash build.sh
binary=$ROOT/bin/hhy_bench
if [[ -z ${RUN_DIR:-} ]]; then RUN_DIR=$ROOT/results/$(date -u +%Y%m%dT%H%M%SZ)_$$; fi
mkdir -p -- "$RUN_DIR"
RUN_DIR=$(cd -- "$RUN_DIR" && pwd)
mkdir -p -- "$RUN_DIR/codes" "$RUN_DIR/workers" "$RUN_DIR/preparation"

# A resume must never mix binaries, code inputs, plans or timing settings.
fingerprint_file=$RUN_DIR/settings.new
{
  printf '%s\n' "MODE=$MODE" "SUITE=$SUITE" "SEEDS=$SEEDS" "INPUTS=$INPUTS" "ROUNDS=$ROUNDS" "BATCH=$BATCH" \
    "RANKS=$RANKS" "ABLATIONS=$ABLATIONS" "GABIDULIN_CONTROLS=$GABIDULIN_CONTROLS" "CAPTURED=$CAPTURED" \
    "ERASURES=$ERASURES" "DUMP_INPUTS=$DUMP_INPUTS" "DENSE_CASES=$DENSE_CASES" "MAX_M=$MAX_M" "PERIOD_CAP=$PERIOD_CAP" \
    "MAX_SETS=$MAX_SETS" "INFO_ATTEMPTS=$INFO_ATTEMPTS" "MAX_TRIAL_RADIUS=$MAX_TRIAL_RADIUS" \
    "CANDIDATE_POLICY=$CANDIDATE_POLICY" "CONSTRUCTION_SEED=$CONSTRUCTION_SEED" "CONSTRUCTION_ATTEMPTS=$CONSTRUCTION_ATTEMPTS" "STRICT_M=$STRICT_M" "BOOTSTRAP=$BOOTSTRAP" \
    "CASE_ORDER_SEED=$CASE_ORDER_SEED" "CXX=${CXX:-g++}" "CXXFLAGS=${CXXFLAGS:-default}" "NTL_PREFIX=${NTL_PREFIX:-system}" "CPU=${CPU:-auto}"
  sha256sum "$PLAN" src/*.cpp include/*.hpp Makefile build.sh run.sh data/*.code bin/hhy_bench
} > "$fingerprint_file"
if [[ -f $RUN_DIR/settings.txt ]]; then
  if ! cmp -s "$fingerprint_file" "$RUN_DIR/settings.txt"; then
    printf '%s\n' 'Sources or settings changed. Use a new RUN_DIR; existing results were not changed.' >&2
    exit 2
  fi
  rm -- "$fingerprint_file"
else mv -- "$fingerprint_file" "$RUN_DIR/settings.txt"; fi

# Pin serial workers to one allowed CPU when taskset is available. This avoids
# competing benchmark workers; it does not lock the processor frequency.
if [[ -z ${CPU:-} ]]; then
  allowed=$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)
  CPU=${allowed%%[-,]*}
fi
pin=()
if command -v taskset >/dev/null 2>&1; then pin=(taskset -c "$CPU"); fi
{
  date -u --iso-8601=seconds
  uname -a
  "${CXX:-g++}" --version
  printf 'CPU affinity: %s\n' "$CPU"
  if command -v lscpu >/dev/null 2>&1; then lscpu; fi
  if command -v ldd >/dev/null 2>&1; then ldd "$binary"; fi
  printf '%s\n' 'Frequency was not programmatically locked.'
} > "$RUN_DIR/environment.txt"
printf 'Results: %s\n' "$RUN_DIR"

status=$RUN_DIR/status.tsv
[[ -f $status ]] || printf 'case\tphase\tseed\tstatus\texit_code\n' > "$status"
child_pid=''
cleanup() { if [[ -n $child_pid ]]; then kill -TERM "$child_pid" 2>/dev/null || true; fi; }
trap cleanup EXIT
trap 'cleanup; exit 130' INT
trap 'cleanup; exit 143' TERM

# The preparation has an aggregate limit. Timed workers use a progress watchdog
# so slow but progressing all-rank sweeps are not killed by a short total limit.
run_worker() {
  local directory=$1
  shift
  mkdir -p -- "$directory"
  "$@" > "$directory/stdout.log" 2> "$directory/stderr.log" &
  child_pid=$!
  local last_change=$SECONDS previous='' current='' watchdog=0 result=0
  while kill -0 "$child_pid" 2>/dev/null; do
    sleep 2
    if [[ -f $directory/progress.txt ]]; then
      current=$(stat -c '%Y:%s' "$directory/progress.txt")
      if [[ $current != "$previous" ]]; then previous=$current; last_change=$SECONDS; fi
    fi
    if (( STAGE_TIMEOUT_SECONDS > 0 && SECONDS-last_change > STAGE_TIMEOUT_SECONDS )); then
      printf 'Watchdog timeout: %s\n' "$directory" >&2
      kill -TERM "$child_pid" 2>/dev/null || true
      sleep 2
      kill -KILL "$child_pid" 2>/dev/null || true
      watchdog=1
      break
    fi
  done
  if wait "$child_pid"; then result=0; else result=$?; fi
  child_pid=''
  if (( watchdog )); then return 124; fi
  return "$result"
}

"$binary" --task order --plan "$PLAN" --seed "$CASE_ORDER_SEED" > "$RUN_DIR/execution_plan.tsv"
failures=0
while IFS=$'\t' read -r label family q minimum_M R r b ell s source_M input_file studies historical_ranks; do
  [[ -n $label ]] || continue
  case_ranks=$RANKS
  if [[ $RANKS == historical ]]; then case_ranks=$historical_ranks; fi
  code=$RUN_DIR/codes/$label.code
  prep=$RUN_DIR/preparation/$label
  mkdir -p -- "$prep"
  if [[ ! -f $code ]]; then
    printf 'Prepare %s (%s; sources: %s)\n' "$label" "$family" "$studies"
    if [[ $family == import ]]; then
      cp -- "$ROOT/$input_file" "$code"
    else
      if timeout --signal=TERM --kill-after=10s "${CONSTRUCTION_TIMEOUT_SECONDS}s" \
          "${pin[@]}" "$binary" --task prepare --family "$family" --q "$q" --M "$minimum_M" --R "$R" \
          --r "$r" --b "$b" --ell "$ell" --s "$s" --source-M "$source_M" --max-M "$MAX_M" \
          --period-cap "$PERIOD_CAP" --construction-seed "$CONSTRUCTION_SEED" --construction-attempts "$CONSTRUCTION_ATTEMPTS" \
          --strict-M "$STRICT_M" --label "$label" --output "$code" \
          > "$prep/construction.log" 2>&1; then
        printf '%s\tprepare\t-\tpassed\t0\n' "$label" >> "$status"
      else
        result=$?
        printf '%s\tprepare\t-\tfailed_or_resource_limit\t%s\n' "$label" "$result" >> "$status"
        printf 'Preparation failed: %s; inspect %s\n' "$label" "$prep/construction.log" >&2
        ((failures+=1))
        continue
      fi
    fi
  fi
  if [[ ! -f $prep/VERIFIED ]]; then
    if timeout --signal=TERM --kill-after=10s "${CONSTRUCTION_TIMEOUT_SECONDS}s" \
        "${pin[@]}" "$binary" --task verify --code "$code" > "$prep/verification.log" 2>&1; then
      printf '%s\n' 'Runtime verification passed.' > "$prep/VERIFIED"
    else
      result=$?
      printf '%s\tverify\t-\tfailed\t%s\n' "$label" "$result" >> "$status"
      printf 'Verification failed: %s; inspect %s\n' "$label" "$prep/verification.log" >&2
      ((failures+=1))
      continue
    fi
  fi
  for seed in "${seed_array[@]}"; do
    worker=$RUN_DIR/workers/$label/seed_$seed
    if [[ -f $worker/SUCCESS ]]; then continue; fi
    # An interrupted worker is restarted with the same deterministic inputs.
    # Its old CSV is not appended. Completed workers are skipped.
    printf 'Benchmark %s seed=%s\n' "$label" "$seed"
    if run_worker "$worker" "${pin[@]}" "$binary" --task bench --code "$code" --output "$worker" \
        --seed "$seed" --inputs "$INPUTS" --rounds "$ROUNDS" --batch "$BATCH" --ranks "$case_ranks" \
        --ablations "$ABLATIONS" --gabidulin-controls "$GABIDULIN_CONTROLS" --captured "$CAPTURED" \
        --erasures "$ERASURES" --dump-inputs "$DUMP_INPUTS" --max-sets "$MAX_SETS" --info-attempts "$INFO_ATTEMPTS" \
        --max-trial-radius "$MAX_TRIAL_RADIUS" --candidate-policy "$CANDIDATE_POLICY" --dense-cases "$DENSE_CASES"; then
      printf '%s\tbenchmark\t%s\tpassed\t0\n' "$label" "$seed" >> "$status"
    else
      result=$?
      printf '%s\tbenchmark\t%s\tfailed_or_timeout\t%s\n' "$label" "$seed" "$result" >> "$status"
      printf 'Worker failed: %s; inspect %s\n' "$label" "$worker/stderr.log" >&2
      ((failures+=1))
    fi
  done
done < "$RUN_DIR/execution_plan.tsv"

if "$binary" --task aggregate --root "$RUN_DIR" --expected-seeds "${#seed_array[@]}" --bootstrap "$BOOTSTRAP"; then
  printf 'Open: %s/report/report.html\n' "$RUN_DIR"
else
  printf '%s\n' 'No complete model data, or aggregation failed. Inspect status.tsv and worker logs.' >&2
  ((failures+=1))
fi
printf 'Failed preparation/verification/worker steps: %s\n' "$failures"
if (( failures > 0 )); then exit 1; fi
