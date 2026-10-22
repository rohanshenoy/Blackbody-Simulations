"""CLI: build a new HFSS project holding one bare base design from a geometry spec, then reopen and verify it.

The solid is created in the canonical pose (propagation along +z, entrance at z = 0) so the runner's
pose rule and BBRsim's frame hold by construction. The input spec is never modified; an existing
output is overwritten only with --force, and never while AEDT's lock file <output>.lock exists.

Exit codes: 0 verified; 2 invalid spec, existing output without --force, or a lock file (no AEDT started);
3 verification differences or a failed AEDT session. Exits 0 and 3 write <output>.build.json.
"""
from __future__ import annotations

import argparse
import dataclasses
import logging
import sys
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bbsim.config import FaceSelector
from bbsim.geometry import GeometryError, check_face_area, select_face_on_bounding_plane
from bbsim.geometry_spec import GeometrySpec, GeometrySpecError, load_geometry_spec, plan_text
from bbsim.inventory import inventory_design, inventory_project
from bbsim.manifest import sha256_file, write_json

log = logging.getLogger("bbsim.build_project")
TOL_MM = 1e-6


class BuildError(RuntimeError):
    """AEDT refused a build step."""


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build an HFSS project with one bare base design from a geometry spec.")
    p.add_argument("--spec", type=Path, required=True, help="Geometry spec TOML (configs/geometries/)")
    p.add_argument("--force", action="store_true", help="Overwrite an existing output project")
    p.add_argument("--dry-run", action="store_true", help="Validate the spec and print the plan; start no AEDT")
    p.add_argument("--aedt-version", default=None, help="Override the spec's aedt_version")
    return p.parse_args(argv)


def verify_build(design: Mapping[str, Any], spec: GeometrySpec) -> list[str]:
    """Differences between an inventoried design and what the spec asks for; empty when the build is right."""
    diffs: list[str] = []
    if design.get("name") != spec.project.design:
        diffs.append(f"design is named {design.get('name')!r}, expected {spec.project.design!r}")
    if design.get("units") != "mm":
        diffs.append(f"model units are {design.get('units')!r}, expected 'mm'")
    objects = list(design.get("objects", []))
    names = [o["name"] for o in objects]
    if names != [spec.project.object]:
        diffs.append(f"objects are {names}, expected exactly [{spec.project.object!r}]")
    for key in ("boundaries", "excitations", "setups"):
        if design.get(key):
            diffs.append(f"base design must be bare but carries {key}: {design[key]}")
    obj = next((o for o in objects if o["name"] == spec.project.object), None)
    if obj is None:
        return diffs
    if str(obj["material"]).lower() != spec.solid.material.lower():
        diffs.append(f"material is {obj['material']!r}, expected {spec.solid.material!r}")
    expected = spec.solid.expected_bounding_box_mm
    if max(abs(float(a) - b) for a, b in zip(obj["bounding_box_mm"], expected)) > TOL_MM:
        diffs.append(f"bounding box {obj['bounding_box_mm']} mm, expected {expected} mm")
    if len(obj["faces"]) != spec.solid.face_count:
        diffs.append(f"{len(obj['faces'])} faces, expected {spec.solid.face_count}")
    for side in ("min", "max"):
        try:
            face = select_face_on_bounding_plane(obj["faces"], expected, FaceSelector("z", side), tol_mm=TOL_MM)
            # Pose: axis x = y = 0. The runner takes the plane-wave and exit-frame origins from these centres.
            if max(abs(c) for c in face.center_mm[:2]) > TOL_MM:
                diffs.append(f"end face z {side}: face {face.id} centre {list(face.center_mm)} mm is off the "
                             "axis x = y = 0")
            check_face_area(face, spec.solid.end_face_area_mm2)
        except GeometryError as exc:
            diffs.append(f"end face z {side}: {exc}")
    return diffs


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    try:
        spec = load_geometry_spec(args.spec)
    except (GeometrySpecError, OSError) as exc:
        print(f"invalid geometry spec {args.spec}: {exc}", file=sys.stderr)
        return 2
    if args.aedt_version:
        spec = dataclasses.replace(spec, project=dataclasses.replace(spec.project, aedt_version=args.aedt_version))
    output = spec.project.output.resolve()
    print(plan_text(spec))
    if args.dry_run:
        print("\nDry run: no project written.")
        return 0
    if output.exists() and not args.force:
        print(f"output exists (use --force to overwrite): {output}", file=sys.stderr)
        return 2
    # AEDT's lock, <project>.aedt.lock, means another AEDT session has the project open or an earlier run died;
    # only the user can tell which, so it is reported, never removed. A leftover <stem>.aedtresults needs
    # nothing: the base design has no setups and run_frequency copies only the .aedt file.
    lock = Path(f"{output}.lock")
    if lock.exists():
        print(f"AEDT lock file exists: {lock}\n  Close the AEDT session that has the project open or, if none does, "
              f"delete the lock left by an earlier run; --force does not remove it.", file=sys.stderr)
        return 2
    output.parent.mkdir(parents=True, exist_ok=True)

    from bbsim.session import aedt_session, aedt_versions  # lazy PyAEDT import

    solid = spec.solid
    versions: dict[str, str] = {}
    built = reopened = error = None
    build_diffs: list[str] = []
    stage = "build"
    try:
        with aedt_session(None, design=spec.project.design, version=spec.project.aedt_version) as hfss:
            versions = aedt_versions(hfss)
            hfss.modeler.model_units = "mm"
            created = hfss.modeler.create_cylinder(orientation="Z", origin=[0.0, 0.0, 0.0], radius=solid.radius_mm,
                                                   height=solid.length_mm, name=spec.project.object,
                                                   material=solid.material)
            if not created:
                raise BuildError("create_cylinder returned False")
            built = inventory_design(hfss)
            build_diffs = verify_build(built, spec)
            hfss.save_project(file_name=str(output), overwrite=True)

        stage = "reopen"
        with aedt_session(output, design=None, version=spec.project.aedt_version) as hfss:
            reopened = inventory_project(hfss)
    except Exception as exc:  # noqa: BLE001 - recorded below and reported as exit 3; the session is already released
        log.exception("%s session failed", stage)
        error = f"{stage} session failed: {type(exc).__name__}: {exc}"

    diffs = [f"after build: {d}" for d in build_diffs]
    if error:
        diffs.append(error)
    else:
        designs = sorted(reopened["designs"])
        if designs != [spec.project.design]:
            diffs.append(f"reopened project designs {designs} != [{spec.project.design!r}]")
        if spec.project.design in reopened["designs"]:
            diffs += [f"after reopen: {d}" for d in verify_build(reopened["designs"][spec.project.design], spec)]
        write_json(output.with_suffix(".inventory.json"), reopened)
    # Written on failure too, so a record of an earlier build never sits beside a project it does not describe.
    write_json(output.with_suffix(".build.json"), {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "spec_file": str(args.spec.resolve()), "spec": dataclasses.asdict(spec), "versions": versions,
        "built_inventory": built, "reopened_inventory": reopened,
        "output_project": str(output), "output_sha256": sha256_file(output) if output.is_file() else None,
        "error": error, "verification_differences": diffs,
    })
    if diffs:
        print("\nBUILD VERIFICATION FAILED:\n  " + "\n  ".join(diffs))
        return 3
    print(f"\nBuild verified: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
