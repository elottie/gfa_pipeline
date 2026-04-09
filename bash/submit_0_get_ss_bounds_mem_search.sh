#!/bin/bash
set -euo pipefail

# bash mem_search.sh 8

traits="${1:?usage: mem_search.sh TRAITS}"
low="${2:-400}"     # MB that fails (set appropriately)
high="${3:-600}"    # MB that succeeds (set appropriately)
delta="${4:-10}"    # precision in MB

mkdir -p logs

is_oom() {
  local jobid="$1"
  # Wait for accounting to be available, then check job state.
  # States vary by site; OUT_OF_MEMORY is common for cgroup OOM kills.
  sacct -j "$jobid" --format=State,ExitCode%20 -n -P | head -n 1
}

while (( high - low > delta )); do
  mid=$(( low + (high - low)/2 ))
  echo "Trying traits=$traits mem=${mid}M (bounds: low=$low high=$high)"

  #jobid=$(sbatch --parsable --mem="${mid}M" submit_1_gather_snps.sh)
  jobid=$(sbatch --parsable --mem="${mid}M" submit_0_get_ss_bounds.sh)
  echo "Submitted $jobid"

  # --- WAIT HERE: replace your squeue loop with this ---
  while true; do
    rec=$(sacct -j "$jobid" -X --format=State,ExitCode -n -P | head -n 1)
    state="${rec%%|*}"

    # still pending/running OR sacct hasn't recorded it yet
    if [[ -z "$state" || "$state" == "PENDING" || "$state" == "RUNNING" || \
          "$state" == "CONFIGURING" || "$state" == "COMPLETING" ]]; then
      sleep 5
      continue
    fi

    break
  done
  # --- END WAIT BLOCK ---

  # read final state
  rec=$(sacct -j "$jobid" --format=State,ExitCode,ElapsedRaw -n -P | head -n 1)

  state="${rec%%|*}"
  rest="${rec#*|}"
  exitcode="${rest%%|*}"
  elapsed="${rec##*|}"   # seconds

  echo "Job $jobid state=$state exit=$exitcode elapsed_s=$elapsed"

  if [[ "$state" == "COMPLETED" ]]; then
    high="$mid"
  elif [[ "$state" == "OUT_OF_MEMORY" ]] || grep -qiE "oom|out of memory|killed process" "logs/rtraits_${jobid}.err" "logs/rtraits_${jobid}.out" 2>/dev/null; then
    low="$mid"
  else
    echo "Non-OOM failure; investigate logs. Not updating bounds."
    exit 2
  fi
done

echo "Estimated minimum mem for traits=$traits: ${high}M (precision ${delta}M)"
