"""CLI: compare a candidate dataset directory with a reference one."""
from __future__ import annotations

import argparse
from pathlib import Path

from bbsim.compare import Tolerances, compare_datasets, compare_polarization_pair
from bbsim.manifest import write_json


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare HFSS gap exports (waveguide.csv, far_field.csv) against a reference.")
    p.add_argument("candidate", type=Path, help="Directory with the new waveguide.csv and far_field.csv")
    p.add_argument("reference", type=Path, help="Directory with the reference CSVs")
    p.add_argument("--json", type=Path, default=None, help="Write the full comparison record here")
    p.add_argument("--t-atol", type=float, default=Tolerances.t_atol, help="Absolute transmission tolerance")
    p.add_argument("--t-rtol", type=float, default=Tolerances.t_rtol, help="Relative transmission tolerance")
    p.add_argument("--corr-warn", type=float, default=Tolerances.corr_warn, help="Warn if |E| correlation is below this")
    p.add_argument("--pair", nargs=2, type=Path, metavar=("CANDIDATE_EPHI1", "REFERENCE_EPHI1"), default=None,
                   help="Also compare the polarization pair (relative phase, CDX-005); the positional "
                        "directories are then the Ephi=0 ones")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    tol = Tolerances(t_atol=args.t_atol, t_rtol=args.t_rtol, corr_warn=args.corr_warn)
    result = compare_datasets(args.candidate, args.reference, tol)
    print(result.format_table())
    record, passed = result.to_dict(), result.passed
    if args.pair:
        pair = compare_polarization_pair((args.candidate, args.pair[0]), (args.reference, args.pair[1]), tol)
        print()
        print(pair.format_table())
        record["pair"] = pair.to_dict()
        passed = passed and pair.passed
        record["passed"] = passed
    if args.json:
        write_json(args.json, record)
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
