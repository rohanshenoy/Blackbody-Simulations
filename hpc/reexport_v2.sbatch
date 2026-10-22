#!/bin/bash -l
# Step 9b-9c (hpc/README.md): the schema-2 re-export of the three solved 2025 R2 jobs, with no solve, and its
# acceptance (HF-050). One AEDT session per project, no solver licence.
#
# Usage (HPC):  mkdir -p /home/rshenoy/BBRSim/outputs   # Slurm needs the --output directory to exist
#               sbatch hpc/reexport_v2.sbatch
#
# Sources (O = /home/rshenoy/BBRSim/outputs) and the fresh jobs written beside them by run_hfss_frequency.py
# --from-solved, which copies each source's project and .aedtresults into the fresh job and never writes the source:
#   crack1     $O/InfParallelPlate_crack1Rohan_500GHz/job_tight1  -> job_tight1_v2      (step 5d)
#   crack2     $O/InfParallelPlate_crack2_500GHz/job_tight1       -> job_tight1_v2      (step 5d)
#   round gap  $O/RoundGap_r50um_2000GHz/job_reference1           -> job_reference1_v2  (step 8)
# Order: crack1 first, with no reference: it is the crack1 500 GHz reference job, which measures the exporter's
# column order and keeps it as job_tight1_v2/scratch/column_order_reference.json (HF-041). If crack1 fails, the job
# stops there, because the others need that file. crack2 and the round gap then run with
# --column-order-reference <that file>; a failure of one does not stop the other.
# Acceptance: bbsim.reexport_check on each pair that re-exported, with --poynting the source's step-8b record:
#   the cracks:    <job_tight1>/diagnose_poynting/poynting.json      (step 8c, the rerun after the .ffd reader fix
#                                                                    9fccac4; RESULTS 8c: it reproduces run 1, kept
#                                                                    as diagnose_poynting_run1/)
#   the round gap: <job_reference1>/diagnose_poynting/poynting.json  (step 8b, its only run)
# reexport_check runs in this job, from this checkout, because it checks the embedded column-order reference against
# this checkout's exporter files (bbsim/extract.py, readers.py, schema.py), as the runner did.
#
# Preflight, before any AEDT session (exit 2, "nothing started"): each source job and its poynting record exist, the
# record was made on that source (its manifest_sha256) and covers every key and both polarizations of it
# (bbsim.reexport_check --preflight), and no fresh job directory exists yet. A failed attempt's job_*_v2 directory
# is removed by hand before a rerun; this script never deletes anything.
# Space: each re-export copies its source's .aedtresults (about 2.5 GB) onto /home NFS, about 7.5 GB for the three;
# check `quota -s` or `df -h /home/rshenoy` first.
# Resources are first guesses: about 10 minutes per project, the copies included; one AEDT seat, no solve.
#
# Expect: three PREFLIGHT: PASS lines; three re-export logs ending with the schema-2 sidecar written, each followed
# by a manifest summary with qualified true for both polarizations, solve_time_s null (no solve) and the source as
# reused_solution; crack1's two column-order cases passed and discriminating, P2/P printed (the legacy tables gave
# 0.713 / 0.725; printed, not pinned); ff_over_flux near 0.179 (crack1) and 0.348 (crack2) at normal incidence;
# REEXPORT: PASS three times, each table with its largest deviations; the last line "== RESULT: PASS" and exit 0.
# Anything else ends "== RESULT: FAIL <what failed>" with exit 1 (exit 2 when nothing was started).
#SBATCH --job-name=bbsim_reexport_v2
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/reexport_v2-%j.out

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
# bb_env.sh is sourced before strict mode, because module and the env scripts are not proven under set -u.
if ! source "$BBSIM_ENV"; then
  echo "ERROR: sourcing $BBSIM_ENV failed" >&2
  echo "== RESULT: FAIL environment: sourcing $BBSIM_ENV failed; nothing started"
  exit 1
fi
# Slurm exports the submitting shell: a stale PYAEDT_DESKTOP_PORT or PYAEDT_PROCESS_ID would make PyAEDT attach to
# another desktop instead of starting its own; bbsim.session.aedt_session also removes them.
unset "${!PYAEDT_@}"
# No set -e: every step's exit status is recorded, and the last line decides the job's.
set -uo pipefail
CHECKOUT=${BBSIM_CHECKOUT:-/home/rshenoy/BBRSim/Blackbody-Simulations}
if ! cd "$CHECKOUT"; then
  echo "== RESULT: FAIL checkout $CHECKOUT not found; nothing started"
  exit 1
fi
PYTHON=${BBSIM_PYTHON:-python}
O=${BBSIM_OUTPUTS:-/home/rshenoy/BBRSim/outputs}
LABELS="crack1 crack2 roundgap"

source_of() {
  case $1 in
    crack1) echo "$O/InfParallelPlate_crack1Rohan_500GHz/job_tight1" ;;
    crack2) echo "$O/InfParallelPlate_crack2_500GHz/job_tight1" ;;
    roundgap) echo "$O/RoundGap_r50um_2000GHz/job_reference1" ;;
  esac
}
job_id_of() { if [ "$1" = roundgap ]; then echo reference1_v2; else echo tight1_v2; fi; }
fresh_of() { echo "$(dirname "$(source_of "$1")")/job_$(job_id_of "$1")"; }
poynting_of() { echo "$(source_of "$1")/diagnose_poynting/poynting.json"; }

# The fresh job's manifest and sidecar, as B1 must show them; informational, the acceptance judges.
summarize() {
  "$PYTHON" - "$1" <<'PY' || echo "   (no summary: the fresh job's manifest or sidecar is unreadable)"
import json
import sys
from pathlib import Path

job = Path(sys.argv[1])
m = json.loads((job / "manifest.json").read_text())
solver = m.get("solver") or {}
print(f"   manifest: qualified {json.dumps(solver.get('qualified'))}; solve_time_s "
      f"{json.dumps((m.get('job') or {}).get('solve_time_s'))} (expect null: no solve); reused_solution "
      f"{solver.get('reused_solution')}")
for e, q in sorted((solver.get("convergence_qualification") or {}).items()):
    print(f"   Ephi={e}: qualified {json.dumps(q.get('qualified'))}, achieved {q.get('achieved')} against target "
          f"{q.get('target')} (passes {q.get('passes_completed')} of {q.get('max_passes')})")
for sidecar in sorted(job.glob("*.dataset.json")):
    record = json.loads(sidecar.read_text())
    print(f"   sidecar {sidecar.name}: schema_version {record.get('schema_version')}")
    for rec in (record.get("far_field") or {}).get("ff_over_flux") or []:
        if rec.get("theta_deg") == 180.0:
            pol = rec.get("polarizations") or {}
            print(f"   ff_over_flux at (phi, theta) = ({rec.get('phi_deg'):g}, 180) deg: Ephi=0 "
                  f"{json.dumps(pol.get('0'))}, Ephi=1 {json.dumps(pol.get('1'))}")
PY
}

echo "host=$(hostname) job=${SLURM_JOB_ID:-none} checkout=$CHECKOUT commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"

# --- preflight: nothing starts unless every source passes -----------------------------------------------------------
PREFLIGHT_FAILED=""
for L in $LABELS; do
  S=$(source_of "$L"); P=$(poynting_of "$L"); F=$(fresh_of "$L")
  echo "== preflight $L: source $S"
  echo "   poynting $P"
  "$PYTHON" -m bbsim.reexport_check --preflight --source "$S" --poynting "$P"
  rc=$?
  [ "$rc" -eq 0 ] || PREFLIGHT_FAILED="$PREFLIGHT_FAILED${PREFLIGHT_FAILED:+; }$L poynting record (exit $rc)"
  if [ -e "$F" ]; then
    echo "   $F exists from an earlier attempt: remove it by hand before a rerun (this script never deletes)"
    PREFLIGHT_FAILED="$PREFLIGHT_FAILED${PREFLIGHT_FAILED:+; }$L fresh job $F exists"
  fi
done
if [ -n "$PREFLIGHT_FAILED" ]; then
  echo "== RESULT: FAIL preflight: $PREFLIGHT_FAILED; nothing started"
  exit 2
fi

# --- re-export: crack1 first; it writes the column-order reference the others embed -----------------------------------
FAILED=""
DONE=""
REFERENCE=$(fresh_of crack1)/scratch/column_order_reference.json
reexport() {
  local label=$1 reference=${2:-}
  echo "== re-export $label: $(source_of "$label") -> $(fresh_of "$label")"
  if [ -n "$reference" ]; then
    "$PYTHON" run_hfss_frequency.py --from-solved "$(source_of "$label")" --job-id "$(job_id_of "$label")" \
      --column-order-reference "$reference"
  else
    "$PYTHON" run_hfss_frequency.py --from-solved "$(source_of "$label")" --job-id "$(job_id_of "$label")"
  fi
}

reexport crack1
rc=$?
if [ "$rc" -ne 0 ] || [ ! -f "$REFERENCE" ]; then
  why="exit $rc"
  [ "$rc" -eq 0 ] && why="no $REFERENCE"
  echo "== RESULT: FAIL crack1 re-export ($why); crack2 and the round gap not started: they need crack1's column-order reference"
  exit 1
fi
DONE="crack1"
summarize "$(fresh_of crack1)"
echo "   column-order reference $REFERENCE:"
"$PYTHON" - "$REFERENCE" <<'PY' || echo "   (the reference is unreadable)"
import json
import sys

reference = json.load(open(sys.argv[1]))
for c in reference.get("cases") or []:
    print(f"   Ephi={c.get('ephi')} at (phi, theta) = ({c.get('phi_deg'):g}, {c.get('theta_deg'):g}) deg: "
          f"P2/P = {c.get('energy_fraction_2'):.4f}, (P1+P2)/P = {c.get('energy_fraction_12'):.4f}, "
          f"discriminating {json.dumps(c.get('discriminating'))}, passed {json.dumps(c.get('passed'))}")
print(f"   reference passed {json.dumps(reference.get('passed'))} (legacy tables: P2/P 0.713 / 0.725; printed, not pinned)")
PY

for L in crack2 roundgap; do
  reexport "$L" "$REFERENCE"
  rc=$?
  if [ "$rc" -eq 0 ]; then
    DONE="$DONE $L"
    summarize "$(fresh_of "$L")"
  else
    FAILED="$FAILED${FAILED:+; }$L re-export (exit $rc)"
  fi
done

# --- acceptance: each pair that re-exported -----------------------------------------------------------------------------
for L in $DONE; do
  echo "== reexport_check $L"
  "$PYTHON" -m bbsim.reexport_check --source "$(source_of "$L")" --fresh "$(fresh_of "$L")" --poynting "$(poynting_of "$L")"
  rc=$?
  [ "$rc" -eq 0 ] || FAILED="$FAILED${FAILED:+; }$L reexport_check (exit $rc)"
done

if [ -n "$FAILED" ]; then
  echo "== RESULT: FAIL $FAILED"
  exit 1
fi
echo "== RESULT: PASS"
exit 0
