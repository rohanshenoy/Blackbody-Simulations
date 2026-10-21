#!/bin/bash -l
# Step 10a (hpc/README.md): one of the first 35-key solves (OR-036), batch C1. A licensed solve: an AEDT seat and an
# HFSS solver increment for its whole run. Run only on Rohan's separate go.
#
# Usage (HPC):  mkdir -p /home/rshenoy/BBRSim/outputs   # Slurm needs the --output directory to exist
#               sbatch hpc/solve_grid35.sbatch configs/<name>_grid35.toml [extra run_hfss_frequency.py args]
# It runs run_hfss_frequency.py --config <config> --job-id grid35_1 [--cores $SLURM_CPUS_PER_TASK] <extra args>, so
# the job is <output root>/<dataset>_<label>/job_grid35_1. --job-id, --config and --from-solved are fixed here and
# refused among the extra args. Order (README 10a): crack1 first with no reference: it is the crack1 500 GHz
# reference job and writes job_grid35_1/scratch/column_order_reference.json (HF-041); then crack2 and the round gap,
# each with --column-order-reference <that file>, each submitted --dependency=afterany on the previous.
#
# Preflight, before any AEDT session (exit 2, "nothing started"): the config loads with the extra args as the runner
# reads them; it declares the 35 keys; the runner's own pre-AEDT check passes (schema 2, both polarizations, a
# theta = 90 deg key, and a column-order reference on every job but crack1 at 500 GHz); a --column-order-reference
# file exists; the prepared project exists; the job directory does not exist yet. This script never deletes. A failed
# attempt's job_grid35_1 whose logs/run.log shows the solve finished holds a licensed solve that a later re-extraction
# can reuse (HF-052): do not delete it; move it aside (mv job_grid35_1 job_grid35_1_failed), keep it and report the
# failure. A rerun solves again, so submit one only on purpose.
#
# After the solve it prints, from the job's manifest and sidecar: per polarization the convergence qualification
# (qualified, achieved delta against the target, passes), the declared key count, the sidecar's transmittance.cell_90
# and, for the reference job, its column-order cases. The last line is "== RESULT: PASS" with exit 0 only if the runner
# exited 0, the manifest records solver.qualified true, the schema-2 sidecar declares 35 keys, its cell_90 equals
# 0.036633054690581604 to 1e-12 relative (theta 90 to 180 deg by 15 deg, OR-036), and the reference job wrote a passed
# column-order reference; otherwise "== RESULT: FAIL <reasons>" with exit 1. It never exits 0 on a failure. When the
# job's logs/run.log shows the solve finished, the FAIL line ends by saying the solved project is kept.
# Resources are first guesses from the 15-key runs times 35/15 in time; memory does not grow with the key count:
# crack1 about 20 min and 2.8 GB, crack2 about 34 min and 4.8 GB (the largest), the round gap about 10 min and 2.5 GB.
#SBATCH --job-name=bbsim_grid35
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/solve_grid35-%j.out
CONFIG=${1:-}
shift || true

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
if [ -z "$CONFIG" ]; then
  echo "usage: sbatch hpc/solve_grid35.sbatch configs/<name>_grid35.toml [extra run_hfss_frequency.py args]"
  echo "== RESULT: FAIL nothing started: no config given"
  exit 2
fi
ARGS=(--config "$CONFIG" --job-id grid35_1)
[ -n "${SLURM_CPUS_PER_TASK:-}" ] && ARGS+=(--cores "$SLURM_CPUS_PER_TASK")
ARGS+=("$@")
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} checkout=$CHECKOUT commit=$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "== preflight: ${ARGS[*]}"

# --- preflight: python reads the config as the runner does; the job directory, or the refusal, comes back in a file --
PRE=$(mktemp "${TMPDIR:-/tmp}/solve_grid35.XXXXXX")
"$PYTHON" - "$PRE" "${ARGS[@]}" <<'PY'
import sys
from pathlib import Path

out = Path(sys.argv[1])


def refuse(reason):
    print(f"preflight: {reason}")
    out.write_text(reason + "\n")
    sys.exit(2)


try:
    from bbsim.config import load_config
    from bbsim.naming import frequency_label
    from bbsim.run_frequency import check_before_aedt, is_reference_job, load_reference, overrides_from_args, parse_args
    from bbsim.sampling import incident_angle_values
except Exception as exc:  # noqa: BLE001 - the reason goes to the verdict line
    refuse(f"cannot import the runner: {type(exc).__name__}: {exc}")

argv = sys.argv[2:]
try:
    args = parse_args(argv)
except SystemExit:
    refuse(f"the runner refuses the arguments {argv}")
config = argv[1]
if args.job_id != "grid35_1" or args.config != Path(config) or args.from_solved is not None:
    refuse("--job-id, --config and --from-solved are fixed by this script (job grid35_1, config " + config
           + f"); the extra arguments change them: job {args.job_id!r}, config {str(args.config)!r}")
try:
    cfg = load_config(args.config, overrides_from_args(args))
except Exception as exc:  # noqa: BLE001
    refuse(f"the config {config} does not load: {type(exc).__name__}: {exc}")
exc_ = cfg.excitation
n_phi = len(incident_angle_values(exc_.phi_lower_deg, exc_.phi_upper_deg, exc_.phi_step_deg))
n_theta = len(incident_angle_values(exc_.theta_lower_deg, exc_.theta_upper_deg, exc_.theta_step_deg))
if n_phi * n_theta != 35:
    refuse(f"the config {config} declares {n_phi} x {n_theta} = {n_phi * n_theta} keys, not the 35 of OR-036")
reference = None
if args.column_order_reference is not None:
    if not Path(args.column_order_reference).is_file():
        refuse(f"the column-order reference {args.column_order_reference} does not exist (crack1's 35-key job writes "
               "it as job_grid35_1/scratch/column_order_reference.json)")
    try:
        reference = load_reference(args.column_order_reference)
    except Exception as exc:  # noqa: BLE001
        refuse(f"the column-order reference {args.column_order_reference} is unreadable: {exc}")
try:
    check_before_aedt(cfg, reference)
except Exception as exc:  # noqa: BLE001
    refuse(f"the runner's pre-AEDT check: {exc}")
if not cfg.project.path.is_file():
    refuse(f"the prepared project {cfg.project.path} does not exist")
job = Path(cfg.output.root) / f"{cfg.project.dataset_id}_{frequency_label(cfg.solver.frequency_ghz)}" / "job_grid35_1"
if job.exists():
    refuse(f"the job directory {job} exists from an earlier attempt: move it aside (mv {job} {job}_failed) and keep "
           "it if its logs/run.log shows the solve finished; a rerun solves again")
writes_reference = reference is None and is_reference_job(cfg)
print(f"preflight: ok; {n_phi} x {n_theta} = 35 keys, schema {cfg.output.schema}, job {job}"
      + ("; the reference job: it writes scratch/column_order_reference.json" if writes_reference
         else f"; column-order reference {args.column_order_reference}"))
out.write_text(f"{job}\n{'yes' if writes_reference else 'no'}\n")
PY
rc=$?
if [ "$rc" -ne 0 ]; then
  why=$(head -n 1 "$PRE")
  rm -f "$PRE"
  echo "== RESULT: FAIL nothing started: ${why:-the preflight exited $rc}"
  exit 2
fi
JOB=$(sed -n 1p "$PRE")
WRITES_REFERENCE=$(sed -n 2p "$PRE")
rm -f "$PRE"

# --- the solve -----------------------------------------------------------------------------------------------------
echo "== solve: run_hfss_frequency.py ${ARGS[*]}"
"$PYTHON" run_hfss_frequency.py "${ARGS[@]}"
rc=$?
echo "== run_hfss_frequency.py exited $rc"

# --- summary and verdict: python prints what the job holds and writes the failed criteria to a file ------------------
REASONS=$(mktemp "${TMPDIR:-/tmp}/solve_grid35.XXXXXX")
echo "== summary: $JOB"
"$PYTHON" - "$JOB" "$WRITES_REFERENCE" "$REASONS" <<'PY'
import json
import math
import sys
from pathlib import Path

job, writes_reference, reasons_file = Path(sys.argv[1]), sys.argv[2] == "yes", Path(sys.argv[3])
CELL_90 = 0.036633054690581604     # cell_90 of theta 90 to 180 deg by 15 deg (OR-036), full precision
reasons = []


def load(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        reasons.append(f"{path.name} unreadable in {path.parent} ({type(exc).__name__})")
        return None


def number(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


manifest = None
if (job / "manifest.json").is_file():
    manifest = load(job / "manifest.json")
else:
    reasons.append(f"no manifest.json in {job}")
if isinstance(manifest, dict):
    solver = manifest.get("solver") or {}
    print(f"   manifest: qualified {json.dumps(solver.get('qualified'))}")
    for e, q in sorted((solver.get("convergence_qualification") or {}).items()):
        q = q or {}
        print(f"   Ephi={e}: qualified {json.dumps(q.get('qualified'))}, achieved {q.get('achieved')} against target "
              f"{q.get('target')} (passes {q.get('passes_completed')} of {q.get('max_passes')})")
    if solver.get("qualified") is not True:
        reasons.append(f"the manifest's solver.qualified is {json.dumps(solver.get('qualified'))}, not true")

sidecars = sorted(job.glob("*.dataset.json")) if job.is_dir() else []
if len(sidecars) != 1:
    reasons.append(f"{len(sidecars)} dataset sidecars in {job}, not 1")
else:
    record = load(sidecars[0])
    if isinstance(record, dict):
        version = record.get("schema_version")
        exc = record.get("excitation") or {}
        phis, thetas = exc.get("incident_phi_deg"), exc.get("incident_theta_deg")
        n = len(phis) * len(thetas) if isinstance(phis, list) and isinstance(thetas, list) else None
        cell = (record.get("transmittance") or {}).get("cell_90")
        print(f"   sidecar {sidecars[0].name}: schema_version {version}; {n} keys, phi {phis} x theta {thetas} deg")
        print(f"   cell_90 {cell!r} (expect {CELL_90!r} to 1e-12)")
        if not (isinstance(version, str) and version.split(".")[0] == "2"):
            reasons.append(f"the sidecar declares schema_version {version!r}, not 2.x")
        if n != 35:
            reasons.append(f"the sidecar declares {n} keys, not 35")
        if not (number(cell) and math.isclose(cell, CELL_90, rel_tol=1e-12, abs_tol=0.0)):
            reasons.append(f"the sidecar's transmittance.cell_90 is {cell!r}, not {CELL_90!r}")

if writes_reference:
    path = job / "scratch" / "column_order_reference.json"
    if not path.is_file():
        reasons.append(f"the reference job wrote no {path}: crack2 and the round gap need it")
    else:
        reference = load(path)
        if isinstance(reference, dict):
            print(f"   column-order reference {path}:")
            for c in reference.get("cases") or []:
                frac = c.get("energy_fraction_2")
                frac = f"{frac:.4f}" if number(frac) else repr(frac)
                print(f"   Ephi={c.get('ephi')} at (phi, theta) = ({c.get('phi_deg')}, {c.get('theta_deg')}) deg: "
                      f"P2/P = {frac}, discriminating {json.dumps(c.get('discriminating'))}, "
                      f"passed {json.dumps(c.get('passed'))}")
            print(f"   reference passed {json.dumps(reference.get('passed'))} (legacy tables: P2/P 0.713 / 0.725; "
                  "printed, not pinned)")
            if reference.get("passed") is not True:
                reasons.append(f"the column-order reference did not pass ({path})")

reasons_file.write_text("".join(r + "\n" for r in reasons))
sys.exit(1 if reasons else 0)
PY
src=$?
FAILED=""
[ "$rc" -eq 0 ] || FAILED="runner exited $rc"
if [ "$src" -ne 0 ]; then
  found=""
  while IFS= read -r line; do
    [ -n "$line" ] && FAILED="$FAILED${FAILED:+; }$line" && found=1
  done < "$REASONS"
  [ -n "$found" ] || FAILED="$FAILED${FAILED:+; }the summary exited $src without a reason"
fi
rm -f "$REASONS"
if [ -n "$FAILED" ]; then
  # A licensed solve that finished is kept whatever failed after it (HF-052): never delete the job.
  if grep -q "solve finished in" "$JOB/logs/run.log" 2>/dev/null; then
    FAILED="$FAILED; the solve finished ($JOB/logs/run.log): the solved project in $JOB/project is kept, do not delete"
    FAILED="$FAILED the job: move it aside (mv $JOB ${JOB}_failed) and report"
  fi
  echo "== RESULT: FAIL $FAILED"
  exit 1
fi
echo "== RESULT: PASS"
exit 0
