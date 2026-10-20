#!/bin/bash -l
# Step 9d (hpc/README.md): assemble a data tree for BBRsim from finished jobs, as steps 5e and 8a laid the trees out:
# <tree_dir>/waveguides/ holds, per job, its sidecar <ID>_<label>.dataset.json and its two <ID>_<label>_Ephi={0,1}
# directories with waveguide.csv, far_field.csv and manifest.json (7 files per dataset), and <tree_dir>/SHA256SUMS
# lists every one of them by its path relative to <tree_dir>. No AEDT and no licence.
#
# Usage (HPC):  sbatch hpc/assemble_tree.sbatch <tree_dir> <expected_files> <job_dir>...
#   T=/resnick/groups/golwala/rshenoy/bbsim/trees; O=/home/rshenoy/BBRSim/outputs
#   sbatch hpc/assemble_tree.sbatch $T/cracks_500GHz_2025r2_v2 14 \
#     $O/InfParallelPlate_crack1Rohan_500GHz/job_tight1_v2 $O/InfParallelPlate_crack2_500GHz/job_tight1_v2
#
# Gate, before anything is copied (OR-027: assembly takes qualified runs only; the per-run override comes with the
# job array): each job's manifest.json must record solver.qualified true, and the job prints one line per
# polarization with its qualification, the achieved delta against the target (a --from-solved job records the
# re-qualification at its source's target, CDX-025); each job holds exactly one sidecar, named by its dataset_id and
# frequency_label, and both Ephi directories with exactly the three files; no dataset appears twice; every sidecar
# declares the same schema major (one schema per tree); the jobs hold exactly <expected_files> files. An existing
# non-empty <tree_dir> is refused: nothing is overwritten or deleted, ever. Then the files are copied,
# SHA256SUMS is written and checked with sha256sum -c, and its line count must equal <expected_files>.
# The last line is "== RESULT: PASS ..." with exit 0, or "== RESULT: FAIL <reason>" with exit 1; argument errors
# exit 2. A copy that fails halfway leaves a partial tree, which the FAIL line names; remove it by hand.
#SBATCH --job-name=bbsim_assemble_tree
#SBATCH --account=golwala
#SBATCH --partition=expansion
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rshenoy/BBRSim/outputs/assemble_tree-%j.out

BBSIM_ENV=${BBSIM_ENV:-/home/rshenoy/BBRSim/bb_env.sh}
# bb_env.sh (for python) is sourced before strict mode, because module and the env scripts are not proven under set -u.
if ! source "$BBSIM_ENV"; then
  echo "ERROR: sourcing $BBSIM_ENV failed" >&2
  echo "== RESULT: FAIL environment: sourcing $BBSIM_ENV failed; nothing copied"
  exit 1
fi
unset "${!PYAEDT_@}"
# No set -e: every step's exit status is checked, and the last line decides the job's. noclobber: > never overwrites.
set -uo pipefail
set -C
PYTHON=${BBSIM_PYTHON:-python}

usage_error() {
  echo "usage: sbatch hpc/assemble_tree.sbatch <tree_dir> <expected_files> <job_dir>..."
  echo "== RESULT: FAIL usage: $1"
  exit 2
}
[ $# -ge 3 ] || usage_error "need a tree directory, the expected file count and at least one job directory"
TREE=$1
EXPECTED=$2
shift 2
case $EXPECTED in
  '' | *[!0-9]* | 0) usage_error "expected_files must be a positive integer, got '$EXPECTED'" ;;
esac
for J in "$@"; do
  [ -d "$J" ] || usage_error "job directory $J does not exist"
done
if [ -e "$TREE" ] && { [ ! -d "$TREE" ] || [ -n "$(ls -A "$TREE")" ]; }; then
  echo "== RESULT: FAIL tree directory $TREE exists and is not empty; nothing is overwritten"
  exit 1
fi
if ! command -v sha256sum >/dev/null; then
  echo "== RESULT: FAIL sha256sum not found; nothing copied"
  exit 1
fi
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} tree=$TREE expected_files=$EXPECTED jobs=$#"

# --- gate, then copy (python reads the JSON; the reason for a refusal comes back through a file) ----------------------
REASON=$(mktemp "${TMPDIR:-/tmp}/assemble_tree.XXXXXX")
"$PYTHON" - "$REASON" "$TREE" "$EXPECTED" "$@" <<'PY'
import json
import shutil
import sys
from pathlib import Path

reason_file, tree, expected = Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3])
jobs = [Path(j) for j in sys.argv[4:]]
EPHI_FILES = ("waveguide.csv", "far_field.csv", "manifest.json")


def stop(problems):
    for p in problems:
        print(f"   refused: {p}")
    more = f" (and {len(problems) - 1} more above)" if len(problems) > 1 else ""
    reason_file.write_text(problems[0] + more)
    sys.exit(1)


def read(path):
    data = json.loads(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("not a JSON object")
    return data


problems, plan, seen, majors = [], [], {}, {}
for job in jobs:
    print(f"-- {job}")
    try:
        manifest = read(job / "manifest.json")
    except (OSError, ValueError) as exc:
        problems.append(f"{job}: manifest.json unreadable ({exc})")
        continue
    solver = manifest.get("solver") if isinstance(manifest.get("solver"), dict) else {}
    quals = solver.get("convergence_qualification")
    if isinstance(quals, dict) and quals:
        for e in sorted(quals):
            q = quals[e] if isinstance(quals[e], dict) else {}
            print(f"   Ephi={e}: qualified {json.dumps(q.get('qualified'))}, achieved {q.get('achieved')!r} against "
                  f"target {q.get('target')!r} (passes {q.get('passes_completed')} of {q.get('max_passes')})")
    else:
        print("   no solver.convergence_qualification recorded")
    if solver.get("qualified") is not True:
        problems.append(f"{job}: solver.qualified is {json.dumps(solver.get('qualified'))}, not true; assembly takes "
                        "qualified runs only (OR-027)")
    sidecars = sorted(job.glob("*.dataset.json"))
    if len(sidecars) != 1:
        problems.append(f"{job}: {len(sidecars)} sidecars (*.dataset.json), not 1")
        continue
    try:
        record = read(sidecars[0])
    except (OSError, ValueError) as exc:
        problems.append(f"{sidecars[0]}: unreadable ({exc})")
        continue
    dataset_id, label, version = record.get("dataset_id"), record.get("frequency_label"), record.get("schema_version")
    stem = f"{dataset_id}_{label}"
    if not (isinstance(dataset_id, str) and isinstance(label, str) and sidecars[0].name == f"{stem}.dataset.json"):
        problems.append(f"{sidecars[0]}: its dataset_id {dataset_id!r} and frequency_label {label!r} do not name it")
        continue
    if not (isinstance(version, str) and version.split(".")[0].isdigit()):
        problems.append(f"{sidecars[0]}: schema_version {version!r} is not a version string")
        continue
    print(f"   dataset {stem}, schema_version {version}")
    if stem in seen:
        problems.append(f"dataset {stem} is in both {seen[stem]} and {job}; a tree holds each dataset once")
        continue
    seen[stem] = job
    majors.setdefault(version.split(".")[0], []).append(stem)
    plan.append((sidecars[0], Path("waveguides") / sidecars[0].name))
    for ephi in (0, 1):
        d = job / f"{stem}_Ephi={ephi}"
        present = sorted(p.name for p in d.iterdir()) if d.is_dir() else []
        missing = [n for n in EPHI_FILES if n not in present]
        extra = [n for n in present if n not in EPHI_FILES]
        if missing or extra:
            problems.append(f"{d}: missing {missing}, unexpected {extra}; a tree directory holds exactly "
                            f"{list(EPHI_FILES)}")
        plan += [(d / n, Path("waveguides") / d.name / n) for n in EPHI_FILES]
if len(majors) > 1:
    problems.append(f"the sidecars declare schema majors {sorted(majors)} ({majors}); one schema per tree")
if not problems and len(plan) != expected:
    problems.append(f"the jobs hold {len(plan)} files for the tree, not the expected {expected}")
if problems:
    stop(problems)

done = 0
try:
    for src, rel in plan:
        dest = tree / rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            raise FileExistsError(f"{dest} exists; nothing is overwritten")
        shutil.copy2(src, dest)
        done += 1
except OSError as exc:
    print(f"   copy failed after {done} of {len(plan)} files: {exc}")
    reason_file.write_text(f"copy failed after {done} of {len(plan)} files ({exc}); remove the partial tree {tree} "
                           "by hand")
    sys.exit(1)
print(f"copied {done} files from {len(jobs)} jobs into {tree / 'waveguides'}")
PY
rc=$?
why=$(cat "$REASON")
rm -f "$REASON"
if [ "$rc" -ne 0 ]; then
  echo "== RESULT: FAIL ${why:-the gate stopped with exit $rc}"
  exit 1
fi

# --- checksums -----------------------------------------------------------------------------------------------------
if ! cd "$TREE"; then
  echo "== RESULT: FAIL cannot enter $TREE after copying"
  exit 1
fi
if ! find waveguides -type f | LC_ALL=C sort | xargs sha256sum > SHA256SUMS; then
  echo "== RESULT: FAIL writing $TREE/SHA256SUMS failed"
  exit 1
fi
if ! sha256sum -c SHA256SUMS; then
  echo "== RESULT: FAIL sha256sum -c SHA256SUMS failed in $TREE (lines above)"
  exit 1
fi
n=$(wc -l < SHA256SUMS | tr -d ' ')
if [ "$n" -ne "$EXPECTED" ]; then
  echo "== RESULT: FAIL SHA256SUMS lists $n files in $TREE, not the expected $EXPECTED"
  exit 1
fi
echo "SHA256SUMS: $n files, all OK; $(du -sh . | awk '{print $1}') in $TREE"
echo "== RESULT: PASS $TREE: $n files from $# jobs"
exit 0
