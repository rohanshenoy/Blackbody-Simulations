"""CLI: produce a cleaned copy of the reference project with renamed base designs.

The input project is copied to a temporary file before AEDT opens it, so the
reference is never touched (not even by a lock file). The copy is saved-as to
the output path, which performs the 2023 R2 -> 2025 R2 conversion. Designs not
in the mapping are deleted from the output only. The output is reopened and
inventoried, and retained designs are diffed against the source.
"""
from __future__ import annotations

import argparse
import logging
import shutil
import sys
import tempfile
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bbsim.inventory import compare_retained_designs, inventory_project
from bbsim.manifest import sha256_file, write_json
from bbsim.mapping import DesignMapping, load_mapping, plan_deletions, validate_mapping

log = logging.getLogger("bbsim.prepare_project")


class PreparationError(RuntimeError):
    """A cleanup step did not take effect."""


def parse_args(argv: list[str] | None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create a cleaned AEDT project with renamed base designs.")
    p.add_argument("--input", type=Path, required=True, help="Reference .aedt (never modified)")
    p.add_argument("--output", type=Path, required=True, help="New .aedt to create")
    p.add_argument("--mapping", type=Path, required=True, help="JSON design mapping (configs/design_mapping.json)")
    p.add_argument("--dry-run", action="store_true", help="Inventory and validate only; write no project")
    p.add_argument("--force", action="store_true", help="Overwrite an existing output project")
    p.add_argument("--aedt-version", default="2025.2")
    return p.parse_args(argv)


def validate_paths(input_path: Path, output_path: Path, force: bool) -> None:
    if not input_path.is_file():
        sys.exit(f"input project not found: {input_path}")
    if input_path.resolve() == output_path.resolve():
        sys.exit("input and output must be different files")
    if output_path.exists() and not force:
        sys.exit(f"output exists (use --force to overwrite): {output_path}")


def apply_cleanup(hfss: Any, mapping: DesignMapping, deletions: list[str]) -> list[dict]:
    """Delete unmapped designs and rename retained designs/objects in the open (output) project.

    Every rename is verified, because PyAEDT's Object3d.name setter only warns when it refuses a name.
    """
    actions: list[dict] = []
    hfss.set_active_design(mapping.renames[0].source_design)  # never delete the active design
    for name in deletions:
        hfss.delete_design(name)
        actions.append({"delete_design": name})
    for r in mapping.renames:
        hfss.set_active_design(r.source_design)
        hfss.modeler[r.source_object].name = r.new_object
        if r.new_object not in hfss.modeler.object_names or r.source_object in hfss.modeler.object_names:
            raise PreparationError(f"design {r.source_design!r}: object {r.source_object!r} was not renamed to {r.new_object!r}")
        actions.append({"rename_object": [r.source_design, r.source_object, r.new_object]})
        hfss.rename_design(r.new_design)
        if r.new_design not in hfss.design_list:
            raise PreparationError(f"design {r.source_design!r} was not renamed to {r.new_design!r}")
        actions.append({"rename_design": [r.source_design, r.new_design]})
    hfss.save_project()
    return actions


def _report(source_inv: dict, mapping: DesignMapping, deletions: list[str], errors: list[str]) -> str:
    lines = [f"Source project: {source_inv['project_file']}", f"Designs found: {len(source_inv['designs'])}", ""]
    for name, d in source_inv["designs"].items():
        objs = ", ".join(f"{o['name']}[{o['material']}]" for o in d["objects"])
        lines.append(f"  {name:28} units={d['units']:3} objects: {objs}  setups={d['setups']} boundaries={len(d['boundaries'])}")
    lines += ["", "Renames:"]
    for r in mapping.renames:
        lines.append(f"  {r.source_design}/{r.source_object} -> {r.new_design}/{r.new_object}  (dataset_id {r.dataset_id}, {r.gap_um} um)")
    lines += ["", f"Deletions ({len(deletions)}): {deletions}"]
    if errors:
        lines += ["", "ERRORS:"] + [f"  {e}" for e in errors]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s")
    input_path, output_path = args.input.resolve(), args.output.resolve()
    validate_paths(input_path, output_path, args.force)
    mapping = load_mapping(args.mapping)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    source_sha = sha256_file(input_path)

    from bbsim.session import aedt_session, aedt_versions  # lazy PyAEDT import

    with tempfile.TemporaryDirectory(dir=output_path.parent, prefix=".prepare_tmp_") as tmp:
        work = Path(tmp) / input_path.name
        shutil.copy2(input_path, work)
        with aedt_session(work, design=None, version=args.aedt_version) as hfss:
            versions = aedt_versions(hfss)
            source_inv = inventory_project(hfss)
            write_json(output_path.with_suffix(".inventory.source.json"), source_inv)
            errors = validate_mapping(mapping, source_inv)
            deletions = plan_deletions(mapping, source_inv)
            print(_report(source_inv, mapping, deletions, errors))
            if errors:
                return 2
            if args.dry_run:
                print("\nDry run: no project written.")
                return 0

            hfss.save_project(file_name=str(output_path), overwrite=True)  # 2023 R2 -> 2025 R2 on this copy only
            actions = [{"save_as": str(output_path)}] + apply_cleanup(hfss, mapping, deletions)

    with aedt_session(output_path, design=None, version=args.aedt_version) as hfss:
        prepared_inv = inventory_project(hfss)
    write_json(output_path.with_suffix(".inventory.json"), prepared_inv)

    diffs = compare_retained_designs(source_inv, prepared_inv, mapping.design_map, mapping.object_map)
    expected, actual = set(mapping.design_map.values()), set(prepared_inv["designs"])
    if expected != actual:
        diffs.append(f"prepared project designs {sorted(actual)} != expected {sorted(expected)}")
    for d in prepared_inv["designs"].values():
        if d["boundaries"] or d["excitations"] or d["setups"]:
            diffs.append(f"{d['name']} carries boundaries/excitations/setups; base designs must be bare")

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_project": str(input_path), "source_sha256": source_sha,
        "source_sha256_after": sha256_file(input_path),
        "output_project": str(output_path), "output_sha256": sha256_file(output_path),
        "versions": versions, "renames": [asdict(r) for r in mapping.renames],
        "deleted_designs": deletions, "actions": actions, "verification_differences": diffs,
    }
    if record["source_sha256_after"] != source_sha:
        diffs.append("source project hash changed during preparation")
    write_json(output_path.with_suffix(".mapping.json"), record)
    if diffs:
        print("\nVERIFICATION FAILED:\n  " + "\n  ".join(diffs))
        return 3
    print(f"\nPrepared project verified: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
